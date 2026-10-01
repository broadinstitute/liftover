"""Tests for the liftover service's client-IP parsing and per-IP rate limiting."""

import json
import os
import tempfile
import unittest
from unittest import mock

# The module reads its rate-limit settings from the environment at import time, so point the
# state file somewhere disposable before importing it. Individual tests override the limits by
# patching the module constants directly.
os.environ.setdefault("RATE_LIMIT_STATE_PATH", os.path.join(tempfile.gettempdir(),
                                                            "liftover_rate_limit_import.json"))

import server


class FakeRequest:
    """Minimal stand-in for a Flask request, carrying only the WSGI environ get_user_ip reads."""

    def __init__(self, x_forwarded_for=None):
        self.environ = {} if x_forwarded_for is None else {"HTTP_X_FORWARDED_FOR": x_forwarded_for}


class GetUserIpTests(unittest.TestCase):

    def test_returns_none_when_header_absent(self):
        self.assertIsNone(server.get_user_ip(FakeRequest()))

    def test_returns_none_when_header_empty(self):
        self.assertIsNone(server.get_user_ip(FakeRequest("")))

    def test_single_entry(self):
        self.assertEqual(server.get_user_ip(FakeRequest("1.2.3.4")), "1.2.3.4")

    def test_strips_whitespace(self):
        self.assertEqual(server.get_user_ip(FakeRequest("  1.2.3.4  ")), "1.2.3.4")

    def test_returns_last_entry_which_is_the_one_the_load_balancer_appended(self):
        self.assertEqual(server.get_user_ip(FakeRequest("10.0.0.1, 1.2.3.4")), "1.2.3.4")

    def test_client_cannot_spoof_by_prepending_entries(self):
        # A caller controls everything before the final comma, so a forged prefix must not win:
        # otherwise a fresh prefix per request would mean a fresh rate-limit budget per request.
        spoofed = "9.9.9.9, 8.8.8.8, 7.7.7.7, 1.2.3.4"
        self.assertEqual(server.get_user_ip(FakeRequest(spoofed)), "1.2.3.4")

    def test_ipv6_entry(self):
        self.assertEqual(server.get_user_ip(FakeRequest("10.0.0.1, 2001:db8::1")), "2001:db8::1")

    def test_returns_none_when_last_entry_is_blank(self):
        self.assertIsNone(server.get_user_ip(FakeRequest("1.2.3.4,   ")))


class ExceedsRateLimitTests(unittest.TestCase):

    def setUp(self):
        self.state_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.state_dir.cleanup)
        self.state_path = os.path.join(self.state_dir.name, "state.json")

        limit = mock.patch.object(server, "RATE_LIMIT_MAX_REQUESTS", 3)
        window = mock.patch.object(server, "RATE_LIMIT_WINDOW_SECONDS", 600)
        limit.start()
        window.start()
        self.addCleanup(limit.stop)
        self.addCleanup(window.stop)

    def call(self, ip, now=1000.0):
        return server.exceeds_rate_limit(ip, now=now, state_path=self.state_path)

    def read_state(self):
        with open(self.state_path) as f:
            return json.load(f)

    def test_first_request_is_allowed_and_creates_the_state_file(self):
        self.assertFalse(self.call("1.2.3.4"))
        self.assertEqual(self.read_state(), {"1.2.3.4": [1000.0, 1]})

    def test_allows_exactly_the_limit_then_rejects(self):
        self.assertEqual([self.call("1.2.3.4") for _ in range(5)],
                         [False, False, False, True, True])

    def test_rejected_requests_do_not_increment_the_counter(self):
        for _ in range(10):
            self.call("1.2.3.4")
        # The counter stops at the limit rather than climbing for as long as the client retries.
        self.assertEqual(self.read_state()["1.2.3.4"][1], 3)

    def test_separate_ips_have_separate_budgets(self):
        for _ in range(3):
            self.call("1.2.3.4")
        self.assertTrue(self.call("1.2.3.4"))
        self.assertFalse(self.call("5.6.7.8"))

    def test_window_expiry_resets_the_budget(self):
        for _ in range(3):
            self.call("1.2.3.4", now=1000.0)
        self.assertTrue(self.call("1.2.3.4", now=1000.0))
        self.assertTrue(self.call("1.2.3.4", now=1599.0))
        self.assertFalse(self.call("1.2.3.4", now=1600.0))

    def test_window_start_does_not_slide_while_the_window_is_open(self):
        self.call("1.2.3.4", now=1000.0)
        self.call("1.2.3.4", now=1500.0)
        self.assertEqual(self.read_state()["1.2.3.4"][0], 1000.0)

    def test_expired_entries_are_dropped_from_the_file(self):
        self.call("1.2.3.4", now=1000.0)
        self.call("5.6.7.8", now=1700.0)
        self.assertEqual(list(self.read_state()), ["5.6.7.8"])

    def spend_budget_in_forked_child(self, ip, attempts, release_fd=None, now=1000.0):
        """Fork a child that calls exceeds_rate_limit, and return its pid.

        The child exits with the number of requests it was allowed, so the parent can add up
        what the workers got between them. Forking is what makes these tests meaningful: it is
        how gunicorn creates its workers, and a child's copy of module state is invisible to
        the parent, so only state that genuinely lives in the file crosses back.
        """
        pid = os.fork()
        if pid != 0:
            return pid
        allowed = 0
        try:
            if release_fd is not None:
                os.read(release_fd, 1)  # block until the parent starts everyone at once
            for _ in range(attempts):
                if not server.exceeds_rate_limit(ip, now=now, state_path=self.state_path):
                    allowed += 1
        except BaseException:
            os._exit(255)
        os._exit(allowed)

    def wait_for_allowed_count(self, pid):
        """Reap a child from spend_budget_in_forked_child and return how many it was allowed."""
        status = os.waitpid(pid, 0)[1]
        self.assertEqual(status & 0xFF, 0, "forked worker was killed by a signal")
        self.assertNotEqual(status >> 8, 255, "forked worker raised")
        return status >> 8

    def test_forked_workers_share_one_budget_rather_than_getting_one_each(self):
        # The regression this guards against is a per-process counter: with --preload, module
        # state is copied into each of the 7 gunicorn workers, so an in-memory counter would
        # give this IP a separate budget per worker. Here two children spend 2 of the limit of
        # 3, which is only visible to the parent if the counters really live in the file.
        for _ in range(2):
            self.assertEqual(self.wait_for_allowed_count(
                self.spend_budget_in_forked_child("1.2.3.4", attempts=1)), 1)
        self.assertFalse(self.call("1.2.3.4"))
        self.assertTrue(self.call("1.2.3.4"))

    def test_concurrent_forked_workers_do_not_exceed_the_limit_between_them(self):
        # Started together on purpose: without flock around the read-modify-write, workers read
        # the same count and overwrite each other, and more than the limit gets through.
        release_read_fd, release_write_fd = os.pipe()
        self.addCleanup(os.close, release_read_fd)
        self.addCleanup(os.close, release_write_fd)
        workers = 4
        pids = [self.spend_budget_in_forked_child("1.2.3.4", attempts=5,
                                                  release_fd=release_read_fd)
                for _ in range(workers)]
        os.write(release_write_fd, b"." * workers)
        self.assertEqual(sum(self.wait_for_allowed_count(pid) for pid in pids), 3)

    def test_recovers_from_a_corrupt_state_file(self):
        with open(self.state_path, "w") as f:
            f.write("{not json at all")
        self.assertFalse(self.call("1.2.3.4"))
        self.assertEqual(self.read_state(), {"1.2.3.4": [1000.0, 1]})

    def test_recovers_from_a_state_file_holding_the_wrong_json_type(self):
        with open(self.state_path, "w") as f:
            json.dump([1, 2, 3], f)
        self.assertFalse(self.call("1.2.3.4"))

    def test_malformed_entries_are_dropped_instead_of_raising(self):
        with open(self.state_path, "w") as f:
            json.dump({"1.2.3.4": "not-a-window", "5.6.7.8": [1000.0, 99], "9.9.9.9": [1000.0]}, f)
        self.assertFalse(self.call("1.2.3.4"))
        self.assertTrue(self.call("5.6.7.8"))
        self.assertNotIn("9.9.9.9", self.read_state())

    def test_empty_state_file_is_treated_as_no_state(self):
        open(self.state_path, "w").close()
        self.assertFalse(self.call("1.2.3.4"))

    def test_tracked_ips_are_capped(self):
        with mock.patch.object(server, "RATE_LIMIT_MAX_TRACKED_IPS", 5):
            for i in range(20):
                # Staggered window starts so the eviction has a defined tie-break order.
                self.call(f"10.0.0.{i}", now=1000.0 + i)
        self.assertLessEqual(len(self.read_state()), 5)
        # Among the equal counts, the most recently started windows are the ones kept.
        self.assertIn("10.0.0.19", self.read_state())

    def test_eviction_keeps_the_ip_that_is_at_its_limit(self):
        # An IP part-way through a flood has the oldest window_start of anyone still tracked,
        # because its window opened when the flood began and does not slide. Evicting on
        # window_start alone would therefore drop exactly the counter that is doing the work
        # and hand the flooder a fresh budget.
        with mock.patch.object(server, "RATE_LIMIT_MAX_TRACKED_IPS", 5):
            for _ in range(3):
                self.call("flooder", now=1000.0)
            self.assertTrue(self.call("flooder", now=1000.0))
            for i in range(5):
                self.call(f"10.0.0.{i}", now=1001.0 + i)
            self.assertIn("flooder", self.read_state())
            self.assertTrue(self.call("flooder", now=1006.0))


class RateLimitEndpointTests(unittest.TestCase):

    def setUp(self):
        self.state_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.state_dir.cleanup)
        state_path = os.path.join(self.state_dir.name, "state.json")

        for name, value in (("RATE_LIMIT_MAX_REQUESTS", 2),
                            ("RATE_LIMIT_WINDOW_SECONDS", 600),
                            ("RATE_LIMIT_STATE_PATH", state_path),
                            ("DISABLE_RATE_LIMIT", False),
                            ("BLOCKED_IPS", frozenset())):
            patcher = mock.patch.object(server, name, value)
            patcher.start()
            self.addCleanup(patcher.stop)

        server.app.config["TESTING"] = True
        self.client = server.app.test_client()

    def get(self, path="/liftover/", ip="1.2.3.4", **params):
        """Issue a request carrying an X-Forwarded-For whose last entry is `ip`."""
        headers = {} if ip is None else {"X-Forwarded-For": f"10.0.0.1, {ip}"}
        return self.client.get(path, query_string=params, headers=headers)

    def test_request_under_the_limit_reaches_the_view(self):
        # No "hg" param, so the view answers with its own validation error rather than running
        # liftOver. Reaching that error is what proves the limiter let the request through.
        response = self.get()
        self.assertEqual(response.status_code, 200)
        self.assertIn("hg", json.loads(response.data)["error"])

    def test_request_over_the_limit_gets_429_with_a_parseable_error_body(self):
        for _ in range(2):
            self.get()
        response = self.get()
        self.assertEqual(response.status_code, 429)
        self.assertEqual(json.loads(response.data)["error"], server.RATE_LIMIT_ERROR_MESSAGE)
        self.assertEqual(response.headers["Retry-After"], "600")

    def test_limit_is_enforced_per_ip(self):
        for _ in range(3):
            self.get(ip="1.2.3.4")
        self.assertEqual(self.get(ip="1.2.3.4").status_code, 429)
        self.assertEqual(self.get(ip="5.6.7.8").status_code, 200)

    def test_spoofed_forwarded_for_prefix_does_not_buy_a_new_budget(self):
        for _ in range(3):
            self.client.get("/liftover/", headers={"X-Forwarded-For": "1.1.1.1, 1.2.3.4"})
        response = self.client.get("/liftover/",
                                   headers={"X-Forwarded-For": "2.2.2.2, 3.3.3.3, 1.2.3.4"})
        self.assertEqual(response.status_code, 429)

    def test_normalize_endpoint_is_limited_too(self):
        for _ in range(2):
            self.get(path="/normalize/")
        self.assertEqual(self.get(path="/normalize/").status_code, 429)

    def test_the_two_endpoints_share_one_budget_per_ip(self):
        self.get(path="/liftover/")
        self.get(path="/normalize/")
        self.assertEqual(self.get(path="/liftover/").status_code, 429)

    def test_cors_preflight_does_not_consume_the_budget(self):
        for _ in range(10):
            self.client.options("/liftover/", headers={"X-Forwarded-For": "10.0.0.1, 1.2.3.4"})
        self.assertEqual(self.get(ip="1.2.3.4").status_code, 200)

    def test_catch_all_route_is_not_limited(self):
        for _ in range(10):
            response = self.get(path="/")
        self.assertEqual(response.status_code, 200)

    def test_blocked_ip_is_rejected_on_its_first_request(self):
        with mock.patch.object(server, "BLOCKED_IPS", frozenset({"1.2.3.4"})):
            self.assertEqual(self.get(ip="1.2.3.4").status_code, 429)
            self.assertEqual(self.get(ip="5.6.7.8").status_code, 200)

    def test_blocked_ip_is_rejected_even_when_rate_limiting_is_disabled(self):
        with mock.patch.object(server, "BLOCKED_IPS", frozenset({"1.2.3.4"})), \
             mock.patch.object(server, "DISABLE_RATE_LIMIT", True):
            self.assertEqual(self.get(ip="1.2.3.4").status_code, 429)

    def test_disable_rate_limit_flag_turns_the_limit_off(self):
        with mock.patch.object(server, "DISABLE_RATE_LIMIT", True):
            for _ in range(10):
                response = self.get()
            self.assertEqual(response.status_code, 200)

    def test_request_without_a_forwarded_for_header_is_not_limited(self):
        # Only reachable off Cloud Run, where there is no verified IP to key a limit on.
        for _ in range(10):
            response = self.get(ip=None)
        self.assertEqual(response.status_code, 200)

    def test_limiter_fails_open_when_the_state_file_cannot_be_written(self):
        with mock.patch.object(server, "RATE_LIMIT_STATE_PATH",
                               "/nonexistent-directory/state.json"):
            for _ in range(10):
                response = self.get()
            self.assertEqual(response.status_code, 200)

    def test_cors_headers_are_present_on_a_rate_limited_response(self):
        # The UI reads the error body cross-origin, so the 429 has to carry CORS headers like
        # every other response does.
        for _ in range(2):
            self.get()
        response = self.get()
        self.assertEqual(response.status_code, 429)
        self.assertEqual(response.headers.get("Access-Control-Allow-Origin"), "*")


def write_indexed_fasta(directory, contigs, bases_per_line=5):
    """Write a small FASTA and its .fai index, wrapping each contig's bases at bases_per_line.

    Args:
        directory (str): where to write the files
        contigs (list): (name, sequence) pairs
        bases_per_line (int): line width, kept small so reads cross line boundaries

    Returns:
        str: the FASTA path
    """
    fasta_path = os.path.join(directory, "reference.fa")
    index_lines = []
    with open(fasta_path, "wb") as fasta_file:
        for name, sequence in contigs:
            fasta_file.write(f">{name}\n".encode("ascii"))
            index_lines.append(f"{name}\t{len(sequence)}\t{fasta_file.tell()}\t{bases_per_line}\t{bases_per_line + 1}\n")
            for i in range(0, len(sequence), bases_per_line):
                fasta_file.write(f"{sequence[i:i + bases_per_line]}\n".encode("ascii"))
    with open(f"{fasta_path}.fai", "wt", encoding="UTF-8") as index_file:
        index_file.writelines(index_lines)
    return fasta_path


class ReferenceAlleleCheckTests(unittest.TestCase):

    CHR8_SEQUENCE = "ACGTACGTAAccggttNNNNGATTACA"

    def setUp(self):
        self.fasta_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.fasta_dir.cleanup)
        self.fasta_path = write_indexed_fasta(self.fasta_dir.name, [
            ("chr8", self.CHR8_SEQUENCE),
            ("chrM", "GATCACAGGT"),
        ])
        # The index cache is keyed by path and temp paths are unique, but clear it anyway so a
        # test never sees an index left over from another one.
        server._FASTA_INDEX_BY_FASTA_PATH.clear()
        patcher = mock.patch.dict(server.FASTA_PATHS, {"hg38": self.fasta_path})
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_reads_single_bases_at_every_position(self):
        for pos, base in enumerate(self.CHR8_SEQUENCE.upper(), start=1):
            self.assertEqual(server.fetch_reference_sequence(self.fasta_path, "chr8", pos, 1), base)

    def test_reads_across_line_boundaries_and_upper_cases_soft_masked_bases(self):
        self.assertEqual(server.fetch_reference_sequence(self.fasta_path, "chr8", 4, 9), "TACGTAACC")

    def test_accepts_chrom_with_or_without_chr_prefix(self):
        self.assertEqual(server.fetch_reference_sequence(self.fasta_path, "8", 1, 4), "ACGT")
        self.assertEqual(server.fetch_reference_sequence(self.fasta_path, "CHR8", 1, 4), "ACGT")

    def test_mt_and_m_both_find_the_mitochondrial_contig(self):
        self.assertEqual(server.fetch_reference_sequence(self.fasta_path, "MT", 1, 3), "GAT")
        self.assertEqual(server.fetch_reference_sequence(self.fasta_path, "chrM", 1, 3), "GAT")

    def test_returns_none_outside_the_contig_or_for_an_unknown_contig(self):
        self.assertIsNone(server.fetch_reference_sequence(self.fasta_path, "chr8", 0, 1))
        self.assertIsNone(server.fetch_reference_sequence(self.fasta_path, "chr8", len(self.CHR8_SEQUENCE), 2))
        self.assertIsNone(server.fetch_reference_sequence(self.fasta_path, "chr9", 1, 1))

    def test_matching_ref_returns_none(self):
        self.assertIsNone(server.get_reference_allele_when_input_ref_differs("hg38", "8", 5, "acg"))

    def test_mismatched_ref_returns_the_reference_allele(self):
        self.assertEqual(server.get_reference_allele_when_input_ref_differs("hg38", "8", 5, "TTT"), "ACG")

    def test_n_masked_reference_is_not_reported_as_a_mismatch(self):
        self.assertIsNone(server.get_reference_allele_when_input_ref_differs("hg38", "8", 17, "A"))

    def test_missing_fasta_fails_open(self):
        with mock.patch.dict(server.FASTA_PATHS, {"hg38": os.path.join(self.fasta_dir.name, "missing.fa")}):
            self.assertIsNone(server.get_reference_allele_when_input_ref_differs("hg38", "8", 5, "T"))

    def test_hg19_mitochondrial_variants_are_not_checked(self):
        with mock.patch.dict(server.FASTA_PATHS, {"hg19": self.fasta_path}):
            self.assertIsNone(server.get_reference_allele_when_input_ref_differs("hg19", "MT", 1, "T"))
            self.assertIsNone(server.get_reference_allele_when_input_ref_differs("hg19", "chrM", 1, "T"))
        self.assertEqual(server.get_reference_allele_when_input_ref_differs("hg38", "chrM", 1, "T"), "G")

    def test_non_integer_position_fails_open(self):
        self.assertIsNone(server.get_reference_allele_when_input_ref_differs("hg38", "8", "5x", "T"))


class LiftoverRefMismatchEndpointTests(unittest.TestCase):

    def setUp(self):
        self.fasta_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.fasta_dir.cleanup)
        fasta_path = write_indexed_fasta(self.fasta_dir.name, [("chr8", "ACGTACGTAACCGGTT")])
        server._FASTA_INDEX_BY_FASTA_PATH.clear()
        for patcher in (
                mock.patch.dict(server.FASTA_PATHS, {"hg38": fasta_path}),
                mock.patch.object(server, "DISABLE_RATE_LIMIT", True),
                # The liftover and normalization tools need bcftools, chain files and full
                # references, none of which exist here. Record what they were asked to lift over.
                mock.patch.object(server, "run_variant_liftover_tool", side_effect=self.fake_liftover),
                mock.patch.object(server, "run_bcftools_norm", side_effect=lambda *args, **kwargs: {})):
            patcher.start()
            self.addCleanup(patcher.stop)
        self.lifted_over = []
        server.app.config["TESTING"] = True
        self.client = server.app.test_client()

    def fake_liftover(self, hg, chrom, pos, ref, alt, verbose=False):
        self.lifted_over.append((chrom, pos, ref, alt))
        return {"output_ref": ref, "output_alt": alt, "liftover_tool": "fake"}

    def liftover_variant(self, pos, ref, alt):
        response = self.client.get("/liftover/", query_string={
            "hg": "hg38-to-hg19", "format": "variant", "chrom": "8", "pos": pos, "ref": ref, "alt": alt})
        return json.loads(response.data)

    def test_matching_ref_is_lifted_over_as_entered_without_warnings(self):
        result = self.liftover_variant(5, "A", "G")
        self.assertEqual(self.lifted_over, [("8", "5", "A", "G")])
        self.assertNotIn("warnings", result)

    def test_mismatched_ref_is_replaced_and_reported_as_a_ref_mismatch_warning(self):
        result = self.liftover_variant(5, "T", "G")
        self.assertEqual(self.lifted_over, [("8", "5", "A", "G")])
        self.assertEqual(result["ref"], "A")
        [warning] = result["warnings"]
        self.assertEqual(warning["code"], "REF_MISMATCH")
        self.assertEqual((warning["input_ref"], warning["reference_ref"]), ("T", "A"))
        self.assertIn("REF T does not match the hg38 reference genome", warning["message"])

    def test_wrong_anchor_base_of_an_insertion_is_corrected_in_ref_and_alt(self):
        result = self.liftover_variant(5, "T", "TGG")
        self.assertEqual(self.lifted_over, [("8", "5", "A", "AGG")])
        self.assertEqual((result["ref"], result["alt"]), ("A", "AGG"))
        self.assertIn("so A>AGG was lifted over instead", result["warnings"][0]["message"])

    def test_wrong_bases_of_a_deletion_are_replaced_by_the_reference_ones(self):
        self.liftover_variant(5, "TT", "T")
        self.assertEqual(self.lifted_over, [("8", "5", "AC", "A")])

    def test_multi_allelic_insertions_all_get_the_corrected_anchor(self):
        self.liftover_variant(5, "T", "TG,TGG")
        self.assertEqual(self.lifted_over, [("8", "5", "A", "AG,AGG")])

    def test_ref_equal_to_alt_is_an_error_not_a_liftover(self):
        result = self.liftover_variant(5, "A", "A")
        self.assertEqual(self.lifted_over, [])
        self.assertIn("are the same", result["error"])

    def test_alt_equal_to_the_reference_allele_is_an_error_not_a_liftover(self):
        result = self.liftover_variant(5, "T", "A")
        self.assertEqual(self.lifted_over, [])
        self.assertIn("no variant to lift over", result["error"])


if __name__ == "__main__":
    unittest.main()
