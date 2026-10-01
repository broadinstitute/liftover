from datetime import datetime
import fcntl
import json
import os
import subprocess
import traceback
import tempfile
import time

# flask imports
from flask import Flask, request, Response, send_from_directory
from flask_cors import CORS
from flask_talisman import Talisman


app = Flask(__name__)

CORS(app)

# On Cloud Run, disable Werkzeug's debug PIN / interactive traceback; keep it on for local dev.
DEBUG = not os.environ.get('RUNNING_ON_GOOGLE_CLOUD_RUN')

# force_https=False because Cloud Run's load balancer terminates TLS and forwards plain HTTP to
# the container with X-Forwarded-Proto: https — an app-level redirect would loop with the LB.
Talisman(app, force_https=False)

FASTA_PATHS = {
    "hg19": "/ref/hg19.fa",
    "hg38": "/ref/hg38.fa",
    "t2t": "/ref/chm13v2.0.fa",
}

UCSC_LIFTOVER_TOOL = "UCSC liftover tool"
BCFTOOLS_LIFTOVER_TOOL = "bcftools liftover plugin"

LIFTOVER_EXAMPLE = f"/liftover/?hg=hg19-to-hg38&format=interval&chrom=chr8&start=140300615&end=140300620"

CHAIN_FILE_PATHS = {
    "hg19-to-hg38": "/hg19ToHg38.over.chain.gz",
    "hg38-to-hg19": "/hg38ToHg19.over.chain.gz",
    "hg38-to-t2t": "/hg38ToHs1.over.chain.gz", # replaced hg38-chm13v2.over.chain.gz based on advice from Giulio Genovese
    "t2t-to-hg38": "/hs1ToHg38.over.chain.gz", # replaced chm13v2-hg38.over.chain.gz based on advice from Giulio Genovese

    #"hg38-to-t2t": "/grch38-chm13v2.chain.gz",  # chain files from https://github.com/marbl/CHM13?tab=readme-ov-file#liftover-resources 
    #"t2t-to-hg38": "/chm13v2-grch38.chain.gz",  # 
}

LIFTOVER_REFERENCE_PATHS = {
    "hg19-to-hg38": (FASTA_PATHS["hg19"], FASTA_PATHS["hg38"]),
    "hg38-to-hg19": (FASTA_PATHS["hg38"], FASTA_PATHS["hg19"]),
    "hg38-to-t2t": (FASTA_PATHS["hg38"], FASTA_PATHS["t2t"]),
    "t2t-to-hg38": (FASTA_PATHS["t2t"], FASTA_PATHS["hg38"]),
}


def _env_flag(name, default=False):
    """Parse a boolean environment variable.

    Treats only 1/true/yes/on (case-insensitive) as true and 0/false/no/off/"" as false, so
    NAME=0 reads as false instead of the way bool(os.environ.get(name)) would make any
    non-empty string truthy. Returns `default` when the variable is unset.
    """
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in ("1", "true", "yes", "on")


def _env_int(name, default):
    """Parse an integer environment variable, falling back to `default` if unset or unparsable."""
    value = os.environ.get(name)
    if value is None or not value.strip():
        return default
    try:
        return int(value)
    except ValueError:
        print(f"WARNING: {name}={value!r} is not an integer; using {default}", flush=True)
        return default


# Per-IP rate limiting. One client once sent 968,609 /liftover/ requests over 20 hours
# (~800/minute). Nothing here throttled it: the service simply ran out of capacity (at most 2
# instances x 7 concurrent requests), and Cloud Run then rejected the overflow at the front
# door, which rejects every caller equally: other users saw an 81% error rate for those
# 20 hours, against a 5% baseline.
#
# The default limit is sized from 21 days of production traffic. Excluding security scanners,
# the busiest legitimate IP made 317 requests in 10 minutes, the 99.9th percentile was 74 and
# the 99th was 19, so 300 per 10 minutes throttles roughly 1 IP in 3,760. Against that flood it
# allows 300 x 120 windows = 36,000 requests per instance over the same 20 hours, or 72,000
# service-wide across the two instances (see the per-instance note below), i.e. about 7% of the
# 968,609 that got through.
RATE_LIMIT_MAX_REQUESTS = _env_int("RATE_LIMIT_MAX_REQUESTS", 300)
RATE_LIMIT_WINDOW_SECONDS = _env_int("RATE_LIMIT_WINDOW_SECONDS", 600)
DISABLE_RATE_LIMIT = _env_flag("DISABLE_RATE_LIMIT")

# Comma-separated list of IPs rejected outright, using the same env var name as the
# SpliceAI-lookup scoring services so both are operated the same way.
BLOCKED_IPS = frozenset(ip.strip() for ip in os.environ.get("BLOCKED_IPS", "").split(",") if ip.strip())

# The counters live in a file rather than in module state because this service runs gunicorn
# with --preload and 7 workers: module state is copied into each worker at fork, so a
# per-process counter would give one IP 7 independent budgets on every instance. One flock'd
# file gives every worker on an instance a single shared budget. Instances still do not share
# state, so with --max-instances 2 the effective service-wide ceiling is twice the limit below.
RATE_LIMIT_STATE_PATH = os.environ.get("RATE_LIMIT_STATE_PATH", "/tmp/liftover_rate_limit_state.json")

# Cap on how many IPs are tracked at once, so the state file cannot grow without bound during a
# burst from many addresses. Past the cap the lowest counts are dropped first, so the callers
# closest to their limit keep their counters and the ones handed a fresh budget are those that
# had barely used the old one.
RATE_LIMIT_MAX_TRACKED_IPS = 10000

# Endpoints worth protecting: both shell out to liftOver / bcftools. The catch-all route
# returns a constant string, so it is left unlimited.
RATE_LIMITED_ENDPOINTS = frozenset({"run_liftover", "normalize_variant"})

RATE_LIMIT_ERROR_MESSAGE = (
    "Rate limit exceeded. This server only supports interactive use. To convert large numbers "
    "of variants or intervals, please run the UCSC liftOver tool or the bcftools liftover "
    "plugin locally. Contact us at https://github.com/broadinstitute/liftover/issues if you "
    "have any questions."
)

print(f"Rate limit: {RATE_LIMIT_MAX_REQUESTS} requests per {RATE_LIMIT_WINDOW_SECONDS}s per IP"
      f"{' (DISABLED by DISABLE_RATE_LIMIT)' if DISABLE_RATE_LIMIT else ''}"
      f"{f', {len(BLOCKED_IPS)} blocked IPs' if BLOCKED_IPS else ''}", flush=True)


def error_response(error_message):
    print(f"ERROR: {error_message}")
    return Response(json.dumps({"error": str(error_message)}), status=200, mimetype='application/json')


def reverse_complement(seq):
    # split on "," so multi-allelic alts (e.g. "G,T" from rsIDs with collocated variants) are
    # reverse-complemented per-allele, preserving allele order.
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return ",".join(allele.translate(table)[::-1] for allele in seq.split(","))


# FASTA index (.fai) contents, keyed by FASTA path. Each index is read once per worker process
# and never changes while the service runs.
_FASTA_INDEX_BY_FASTA_PATH = {}


def read_fasta_index(fasta_path):
    """Read the samtools-style .fai index that sits next to a FASTA file.

    Args:
        fasta_path (str): path of the FASTA file. Its index must be at fasta_path + ".fai".

    Returns:
        dict: contig name -> (contig length, byte offset of its first base, bases per line,
            bytes per line), one entry per line of the index.
    """
    if fasta_path not in _FASTA_INDEX_BY_FASTA_PATH:
        fasta_index = {}
        with open(f"{fasta_path}.fai", "rt", encoding="UTF-8") as index_file:
            for line in index_file:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 5:
                    fasta_index[fields[0]] = tuple(int(value) for value in fields[1:5])
        _FASTA_INDEX_BY_FASTA_PATH[fasta_path] = fasta_index
    return _FASTA_INDEX_BY_FASTA_PATH[fasta_path]


def find_fasta_contig_name(fasta_index, chrom):
    """Return the name a FASTA index uses for a chromosome, or None when it has no such contig.

    The references disagree on naming: hg19.fa has "8" and "MT", while hg38.fa and chm13v2.0.fa
    have "chr8" and "chrM". The page and the API accept either spelling for any build.

    Args:
        fasta_index (dict): the output of read_fasta_index
        chrom (str): chromosome name, with or without a "chr" prefix
    """
    chrom_without_prefix = chrom[3:] if chrom.lower().startswith("chr") else chrom
    chrom_without_prefix = chrom_without_prefix.upper()
    if chrom_without_prefix in ("M", "MT"):
        candidate_names = ["MT", "chrM", "M", "chrMT"]
    else:
        candidate_names = [chrom_without_prefix, f"chr{chrom_without_prefix}"]
    for candidate_name in candidate_names:
        if candidate_name in fasta_index:
            return candidate_name
    return None


def fetch_reference_sequence(fasta_path, chrom, pos, length):
    """Read a stretch of the reference genome, using the FASTA's .fai index to seek straight to it.

    Args:
        fasta_path (str): path of an indexed FASTA file
        chrom (str): chromosome name, with or without a "chr" prefix
        pos (int): 1-based position of the first base to read
        length (int): number of bases to read

    Returns:
        str: the bases, upper-cased, or None when the contig isn't in the FASTA or the interval
            runs outside it
    """
    fasta_index = read_fasta_index(fasta_path)
    contig_name = find_fasta_contig_name(fasta_index, chrom)
    if contig_name is None:
        return None
    contig_length, contig_byte_offset, bases_per_line, bytes_per_line = fasta_index[contig_name]
    if pos < 1 or length < 1 or pos + length - 1 > contig_length:
        return None

    # 0-based index of the first base and of the last base, turned into byte offsets by skipping
    # the newline(s) at the end of every full line before them.
    first_base_index = pos - 1
    last_base_index = pos + length - 2
    first_byte = (contig_byte_offset + (first_base_index // bases_per_line) * bytes_per_line
                  + first_base_index % bases_per_line)
    last_byte = (contig_byte_offset + (last_base_index // bases_per_line) * bytes_per_line
                 + last_base_index % bases_per_line)
    with open(fasta_path, "rb") as fasta_file:
        fasta_file.seek(first_byte)
        raw_bytes = fasta_file.read(last_byte - first_byte + 1)
    sequence = raw_bytes.decode("ascii").replace("\n", "").replace("\r", "").upper()
    return sequence if len(sequence) == length else None


def get_reference_allele_when_input_ref_differs(genome, chrom, pos, ref):
    """Check a variant's REF allele against the reference genome it was entered on.

    Deliberately fails open, like check_ref_allele in SpliceAI-lookup's server.py: a missing FASTA
    or index, an unknown contig, a position outside the contig, and a reference that isn't plain
    ACGT (both hg19 and hg38 hard-mask the chrY pseudoautosomal regions to N) all return None, so a
    position the reference can't speak to is lifted over as entered rather than rejected.

    Args:
        genome (str): "hg19", "hg38" or "t2t"
        chrom (str): chromosome name, with or without a "chr" prefix
        pos (int or str): 1-based position of the variant
        ref (str): the REF allele as entered

    Returns:
        str: the reference genome's allele at the same position and length as ref, when it differs
            from ref. None when ref matches, or when there is nothing to compare it to.
    """
    # hg19's mitochondrial contig can't be checked: hg19.fa's MT is the 16,569 bp rCRS, while the
    # hg19ToHg38 chain lifts UCSC hg19's 16,571 bp chrM (NC_001807), whose coordinates drift from
    # rCRS after base 149. Comparing the two would report mismatches for correctly entered variants.
    if genome == "hg19" and chrom.upper().replace("CHR", "") in ("M", "MT"):
        return None
    try:
        reference_allele = fetch_reference_sequence(FASTA_PATHS[genome], chrom, int(pos), len(ref))
    except Exception as e:
        print(f"WARNING: unable to read the {genome} reference at {chrom}:{pos}: {type(e).__name__}: {e}",
              flush=True)
        return None
    if reference_allele is None or any(base not in "ACGT" for base in reference_allele):
        return None
    if reference_allele == ref.upper():
        return None
    return reference_allele


def get_user_ip(request):
    """Return the client IP that Cloud Run's load balancer verified.

    On Cloud Run the X-Forwarded-For header is "<client-supplied>..., <verified-client>", and
    only that final entry is appended by the load balancer. Returning the whole header, or its
    first entry, would let a caller prepend anything to claim a fresh rate-limit budget on every
    request, or to get an innocent IP throttled and logged as the source of a flood. Returns
    None when the header is absent, which off Cloud Run means local development.
    """
    xff = request.environ.get("HTTP_X_FORWARDED_FOR", "")
    if not xff:
        return None
    return xff.rsplit(",", 1)[-1].strip() or None


def _load_rate_limit_state(state_file):
    """Read the JSON counter file, returning {} when it is empty, corrupt, or malformed.

    Entries that are not a [window_start, count] pair of numbers are dropped rather than
    allowed to raise later, since a single bad entry would otherwise disable the limiter for
    as long as the file survives.
    """
    state_file.seek(0)
    raw = state_file.read()
    if not raw.strip():
        return {}
    try:
        state = json.loads(raw)
    except ValueError:
        print("WARNING: rate limit state file was not valid JSON; starting over", flush=True)
        return {}
    if not isinstance(state, dict):
        return {}
    return {
        ip: window for ip, window in state.items()
        if isinstance(window, list) and len(window) == 2
        and all(isinstance(value, (int, float)) for value in window)
    }


def exceeds_rate_limit(user_ip, now=None, state_path=None):
    """Count one request from `user_ip` and return True if it should be rejected.

    Uses a fixed window: the first request from an IP opens a window of
    RATE_LIMIT_WINDOW_SECONDS, and once RATE_LIMIT_MAX_REQUESTS requests land inside it every
    further request is rejected until that window expires. A client straddling a window
    boundary can therefore briefly burst to twice the limit, which is not worth extra state to
    prevent: the goal is to stop one caller from consuming the whole service, not to meter
    exactly.

    Rejected requests are not counted, so an IP's counter stops at the limit instead of climbing
    for as long as the client keeps retrying. All workers on this instance share the one state
    file, serialized with flock.
    """
    if now is None:
        now = time.time()
    if state_path is None:
        state_path = RATE_LIMIT_STATE_PATH

    fd = os.open(state_path, os.O_RDWR | os.O_CREAT, 0o600)
    with os.fdopen(fd, "r+", encoding="UTF-8") as state_file:
        fcntl.flock(state_file, fcntl.LOCK_EX)
        try:
            # Expired windows are dropped on every pass, both so those IPs start fresh and so
            # the file does not accumulate an entry for every IP ever seen.
            state = {
                ip: window for ip, window in _load_rate_limit_state(state_file).items()
                if now - window[0] < RATE_LIMIT_WINDOW_SECONDS
            }

            window_start, count = state.get(user_ip, (now, 0))
            rejected = count >= RATE_LIMIT_MAX_REQUESTS
            if not rejected:
                state[user_ip] = [window_start, count + 1]

            if len(state) > RATE_LIMIT_MAX_TRACKED_IPS:
                # Keep the highest counts, breaking ties toward the most recently opened
                # window. Evicting an entry hands that IP a fresh budget, so the entries that
                # must survive are the ones actually withholding one; dropping a count-of-1
                # entry costs a single request of slack. Sorting on window_start alone would
                # do the opposite, since a sustained flooder's window_start is fixed at the
                # start of its flood and therefore ages into the first-evicted bucket.
                state = dict(sorted(state.items(), key=lambda item: (item[1][1], item[1][0]),
                                    reverse=True)[:RATE_LIMIT_MAX_TRACKED_IPS])

            state_file.seek(0)
            state_file.truncate()
            json.dump(state, state_file)
            state_file.flush()
        finally:
            fcntl.flock(state_file, fcntl.LOCK_UN)

    return rejected


def rate_limit_response():
    """429 whose JSON body matches the {"error": ...} shape the web UI already parses.

    Unlike error_response, which answers 200 for validation errors and is left that way for
    backward compatibility, this deliberately uses a real 429: a client told everything is fine
    has no reason to slow down, and Cloud Run already returns 429 when this service is
    saturated, so callers that handle it already exist. Retry-After is the full window length,
    which over-estimates the wait for a caller whose window is nearly over but never
    under-estimates it.
    """
    return Response(
        json.dumps({"error": RATE_LIMIT_ERROR_MESSAGE}),
        status=429,
        mimetype='application/json',
        headers=[("Retry-After", str(RATE_LIMIT_WINDOW_SECONDS))],
    )


@app.before_request
def enforce_rate_limit():
    """Reject blocked and over-limit callers before any liftover work happens.

    Routing has already run by the time before_request handlers fire, so request.endpoint names
    the view that would handle this request and only the two endpoints that shell out to
    external tools are limited. Returning None lets the request proceed normally.
    """
    if request.endpoint not in RATE_LIMITED_ENDPOINTS:
        return None

    # A CORS preflight does no liftover work and is issued by the browser, not the caller, so
    # charging it against the budget would halve what a cross-origin client gets.
    if request.method == "OPTIONS":
        return None

    user_ip = get_user_ip(request)
    if user_ip is None:
        # No load-balancer-verified client IP, which happens only off Cloud Run.
        return None

    if user_ip in BLOCKED_IPS:
        print(f"RATE LIMIT: rejecting blocked ip {user_ip}", flush=True)
        return rate_limit_response()

    if DISABLE_RATE_LIMIT:
        return None

    try:
        if exceeds_rate_limit(user_ip):
            print(f"RATE LIMIT: {user_ip} exceeded {RATE_LIMIT_MAX_REQUESTS} requests per "
                  f"{RATE_LIMIT_WINDOW_SECONDS}s", flush=True)
            return rate_limit_response()
    except Exception as e:
        # Fail open so a full disk or a permissions problem cannot take the service down, but
        # print loudly: a limiter that is silently off is worse than one that is visibly off.
        print(f"SECURITY: rate-limit check failed (failing open): {type(e).__name__}: {e}",
              flush=True)
        traceback.print_exc()

    return None


def run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=False):
    if hg not in CHAIN_FILE_PATHS or hg not in LIFTOVER_REFERENCE_PATHS:
        raise ValueError(f"Unexpected hg arg value: {hg}")
    chain_file_path = CHAIN_FILE_PATHS[hg]
    source_fasta_path, destination_fasta_path = LIFTOVER_REFERENCE_PATHS[hg]

    with tempfile.NamedTemporaryFile(suffix=".vcf", mode="wt", encoding="UTF-8") as input_file, \
        tempfile.NamedTemporaryFile(suffix=".vcf", mode="rt", encoding="UTF-8") as output_file:

        #  command syntax: liftOver oldFile map.chain newFile unMapped
        chrom = "chr" + chrom.replace("chr", "")

        input_file.write(f"""##fileformat=VCFv4.2
##contig=<ID={chrom},length=100000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
{chrom}	{pos}	.	{ref}	{alt}	60	.	.""")
        input_file.flush()

        command = (
            f"cat {input_file.name} | "
            f"bcftools plugin liftover -- --src-fasta-ref {source_fasta_path} --fasta-ref {destination_fasta_path} --chain {chain_file_path} 2>&1 | "
            f"tail -n 1 > {output_file.name}"
        )

        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT, encoding="UTF-8")

            results = output_file.read().strip()

            # example: chr8	140300616	.	T	G	60	.	.

            result_fields = results.strip().split("\t")
            if verbose:
                print(f"{BCFTOOLS_LIFTOVER_TOOL} {hg} liftover on {chrom}:{pos} {ref}>{alt} returned: {result_fields}", flush=True)

            if len(result_fields) > 4:
                result_fields[1] = int(result_fields[1])

                # bcftools liftover plugin annotates the INFO column with:
                #   FLIP (Flag) — alleles were reverse-complemented (strand flip)
                #   SWAP=N (Integer) — REF/ALT swap; N is which alt became the ref (-1 = new ref from the source-ref allele)
                # See https://github.com/freeseek/score/blob/master/liftover.c
                info_field = result_fields[7] if len(result_fields) > 7 else ""
                info_tags = set(info_field.split(";")) if info_field and info_field != "." else set()
                output_reverse_complemented = "FLIP" in info_tags
                output_ref_alt_swap = None
                for tag in info_tags:
                    if tag.startswith("SWAP="):
                        try:
                            output_ref_alt_swap = int(tag.split("=", 1)[1])
                        except ValueError:
                            pass
                        break

                return {
                    "hg": hg,
                    "chrom": chrom,
                    "start": int(pos) - 1,
                    "end": pos,
                    "output_chrom": result_fields[0],
                    "output_pos": result_fields[1],
                    "output_ref": result_fields[3],
                    "output_alt": result_fields[4],
                    "output_reverse_complemented": output_reverse_complemented,
                    "output_ref_alt_swap": output_ref_alt_swap,
                    "liftover_tool": BCFTOOLS_LIFTOVER_TOOL,
                }

        except Exception as e:
            variant = f"{hg}  {chrom}:{pos} {ref}>{alt}"
            print(f"ERROR in {BCFTOOLS_LIFTOVER_TOOL} for {variant}: {e}")
            print("Falling back on UCSC liftover tool..")
            #traceback.print_exc()
            #raise ValueError(f"liftOver command failed for {variant}: {e}")

        # if bcftools liftover failed, fall back on running UCSC liftover
        chrom = "chr" + chrom.replace("chr", "")
        result = run_UCSC_liftover_tool(hg, chrom, int(pos)-1, pos, verbose=False)
        if result.get("output_strand") == "-":
            result["output_ref"] = reverse_complement(ref)
            result["output_alt"] = reverse_complement(alt)
        else:
            result["output_ref"] = ref
            result["output_alt"] = alt
        return result


def run_UCSC_liftover_tool(hg, chrom, start, end, verbose=False):
    if hg not in CHAIN_FILE_PATHS:
        raise ValueError(f"Unexpected hg arg value: {hg}")
    chain_file_path = CHAIN_FILE_PATHS[hg]

    reason_liftover_failed = ""
    with tempfile.NamedTemporaryFile(suffix=".bed", mode="wt", encoding="UTF-8") as input_file, \
        tempfile.NamedTemporaryFile(suffix=".bed", mode="rt", encoding="UTF-8") as output_file, \
        tempfile.NamedTemporaryFile(suffix=".bed", mode="rt", encoding="UTF-8") as unmapped_output_file:

        #  command syntax: liftOver oldFile map.chain newFile unMapped
        chrom = "chr" + chrom.replace("chr", "")
        input_file.write("\t".join(map(str, [chrom, start, end, ".", "0", "+"])) + "\n")
        input_file.flush()
        command = f"liftOver {input_file.name} {chain_file_path} {output_file.name} {unmapped_output_file.name}"

        try:
            subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT, encoding="UTF-8")
            results = output_file.read()
            if verbose:
                print(f"{UCSC_LIFTOVER_TOOL} {hg} liftover on {chrom}:{start}-{end} returned: {results}", flush=True)

            result_fields = results.strip().split("\t")
            if len(result_fields) > 5:
                result_fields[1] = int(result_fields[1])
                result_fields[2] = int(result_fields[2])

                return {
                    "hg": hg,
                    "chrom": chrom,
                    "pos": int(start) + 1,
                    "start": start,
                    "end": end,
                    "output_chrom": result_fields[0],
                    "output_pos":   int(result_fields[1]) + 1,
                    "output_start": result_fields[1],
                    "output_end":    result_fields[2],
                    "output_strand": result_fields[5],
                    "liftover_tool": UCSC_LIFTOVER_TOOL,
                }
            else:
                reason_liftover_failed = unmapped_output_file.readline().replace("#", "").strip()

        except Exception as e:
            variant = f"{hg}  {chrom}:{start}-{end}"
            print(f"ERROR during liftover for {variant}: {e}")
            traceback.print_exc()
            raise ValueError(f"liftOver command failed for {variant}: {e}")

    if reason_liftover_failed:
        raise ValueError(f"{hg} liftover failed for {chrom}:{start}-{end} {reason_liftover_failed}")
    else:
        raise ValueError(f"{hg} liftover failed for {chrom}:{start}-{end} for unknown reasons")


def run_bcftools_norm(genome_version, chrom, pos, ref, alt, verbose=False):
    if genome_version not in FASTA_PATHS:
        raise ValueError(f"Unexpected genome_version: {genome_version}")

    if len(ref) == len(alt):
        return {
            "normalized_chrom": chrom,
            "normalized_pos": pos,
            "normalized_ref": ref,
            "normalized_alt": alt,
        }

    if genome_version == "hg19":
        chrom = chrom.replace("chr", "")
    else:
        chrom = "chr" + chrom.replace("chr", "")

    with tempfile.NamedTemporaryFile(suffix=".vcf", mode="wt", encoding="UTF-8") as input_file, \
            tempfile.NamedTemporaryFile(suffix=".vcf", mode="rt", encoding="UTF-8") as output_file:

        input_file.write(f"""##fileformat=VCFv4.2        
##contig=<ID={chrom},length=100000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
{chrom}	{pos}	.	{ref}	{alt}	60	.	.""")
        input_file.flush()
        
        fasta_path = FASTA_PATHS[genome_version]

        #results = pysam.bcftools.norm("-f", fasta_path, input_file.name, split_lines=True)
        command = (
            f"cat {input_file.name} | bcftools norm -f {fasta_path} 2>&1 | grep -v total | tail -n 1 > {output_file.name}"
        )

        subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT, encoding="UTF-8")

        last_line_of_output_file = output_file.read().strip()

        # example: chr8	140300616	.	T	G	60	.	.

    result_fields = last_line_of_output_file.strip().split("\t")
    if verbose:
        print(f"bcftools norm {chrom}:{pos} {ref}>{alt} on {genome_version} returned: {result_fields}", flush=True)

    if len(result_fields) < 5:
        raise ValueError(f"bcftools norm failed for {chrom}:{pos} {ref}>{alt}: {last_line_of_output_file}")

    return {
        "normalized_chrom": result_fields[0],
        "normalized_pos": result_fields[1],
        "normalized_ref": result_fields[3],
        "normalized_alt": result_fields[4],
    }


@app.route("/normalize/", methods=['POST', 'GET'])
def normalize_variant():
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" t{os.getpid()}"
    verbose = True

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if not params:
        params.update(request.get_json(force=True, silent=True) or {})

    genome_version = params.get("g")
    if genome_version == "hg37":
        genome_version = "hg19"
    if not genome_version or genome_version not in FASTA_PATHS:
        return error_response(f'"g" param error. It should be set to {" or ".join(FASTA_PATHS)}')

    for key in "chrom", "pos", "ref", "alt":
        if not params.get(key):
            return error_response(f'"{key}" param not specified')

    chrom = params.get("chrom")
    pos = params.get("pos")
    ref = params.get("ref")
    alt = params.get("alt")
    variant_log_string = f"{chrom}:{pos} {ref}>{alt}"

    user_ip = get_user_ip(request)
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" {user_ip} t{os.getpid()}"

    print(f"{logging_prefix}: ======================", flush=True)
    print(f"{logging_prefix}: normalize: {variant_log_string}", flush=True)

    try:
        result = run_bcftools_norm(genome_version, chrom, pos, ref, alt, verbose=True)
    except Exception as e:
        return error_response(e)

    return Response(json.dumps({**params, **result}), mimetype='application/json')


@app.route("/liftover/", methods=['POST', 'GET'])
def run_liftover():
    user_ip = get_user_ip(request)
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" {user_ip} t{os.getpid()}"
    verbose = True

    # check params
    params = {}
    if request.values:
        params.update(request.values)

    if "format" not in params:
        params.update(request.get_json(force=True, silent=True) or {})

    VALID_HG_VALUES = set(CHAIN_FILE_PATHS.keys())
    hg = params.get("hg")
    if not hg or hg not in VALID_HG_VALUES:
        return error_response(f'"hg" param error. It should be set to {" or ".join(VALID_HG_VALUES)}. For example: {LIFTOVER_EXAMPLE}\n')

    VALID_FORMAT_VALUES = ("interval", "variant", "position")
    format = params.get("format", "")
    if not format or format not in VALID_FORMAT_VALUES:
        return error_response(f'"format" param error. It should be set to {" or ".join(VALID_FORMAT_VALUES)}. For example: {LIFTOVER_EXAMPLE}\n')

    chrom = params.get("chrom")
    if not chrom:
        return error_response(f'"chrom" param not specified')

    if format == "interval":
        for key in "start", "end":
            if not params.get(key):
                return error_response(f'"{key}" param not specified')
        start = params.get("start")
        end = params.get("end")
        # liftOver exits 255 on a BED line whose coordinates aren't integers or whose end
        # precedes its start, and that reaches the log as a traceback rather than as something
        # the caller can act on. A digit dropped from the end coordinate is the usual cause.
        try:
            if int(end) < int(start):
                return error_response(f'"end" ({end}) must not be less than "start" ({start})')
        except ValueError:
            return error_response(f'"start" and "end" params must be integers, not "{start}" and "{end}"')
        variant_log_string = f"{start}-{end}"

    elif format == "position":
        try:
            pos = int(params["pos"])
        except Exception as e:
            return error_response(f'"pos" param error: {e}')

        start = pos - 1
        end = pos
        variant_log_string = f"{pos} "
    elif format == "variant":
        for key in "pos", "ref", "alt":
            if not params.get(key):
                return error_response(f'"{key}" param not specified')
        pos = params.get("pos")
        ref = params.get("ref")
        alt = params.get("alt")
        variant_log_string = f"{pos} {ref}>{alt}"

    if verbose:
        print(f"{logging_prefix}: ======================", flush=True)
        print(f"{logging_prefix}: {hg} liftover {format}: {chrom}:{variant_log_string}", flush=True)

    warnings = []
    if format == "variant" and ref.upper() in alt.upper().split(","):
        return error_response(f"REF {ref} and ALT {alt} are the same, so there is no variant to lift over. "
                              f"To lift over the position itself, use format=position.")
    if format == "variant":
        # A REF the reference genome doesn't have makes the bcftools plugin treat the typed REF as
        # a second ALT (8-141310715-T-G on hg38 came back as C>T,G). So the reference allele is
        # lifted over instead, and the mismatch is reported as a warning rather than dropped, the
        # way GeneBe's liftover and variant-relaxed APIs do.
        input_reference_genome = hg.split("-")[0]
        reference_allele = get_reference_allele_when_input_ref_differs(input_reference_genome, chrom, pos, ref)
        if reference_allele:
            # An indel's ALT alleles start with the same anchor base as its REF (8-141310715-T-TA),
            # so the anchor is corrected in them too: GeneBe and Ensembl both read that input as
            # C>CA, where replacing the REF alone would lift over C>TA, a different variant.
            alts = alt.upper().split(",")
            if all(a[0] == ref[0].upper() for a in alts):
                alts = [reference_allele[0] + a[1:] for a in alts]
            corrected_alt = ",".join(alts)
            if reference_allele in alts:
                return error_response(
                    f"REF {ref} does not match the {input_reference_genome} reference genome, which has "
                    f"{reference_allele} at {chrom}:{pos}. That is the ALT allele entered, so there is no "
                    f"variant to lift over.")
            warnings.append({
                "code": "REF_MISMATCH",
                "message": (
                    f"REF {ref} does not match the {input_reference_genome} reference genome, which has "
                    f"{reference_allele} at {chrom}:{pos}, so {reference_allele}>{corrected_alt} was lifted over "
                    f"instead."),
                "input_ref": ref,
                "reference_ref": reference_allele,
            })
            ref = reference_allele
            alt = corrected_alt
            params["ref"] = ref
            params["alt"] = alt

    try:
        if format == "variant":
            result = run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=verbose)
            try:
                normalized_input_variant = run_bcftools_norm(input_reference_genome, chrom, pos, ref, alt, verbose=verbose)
                result.update(normalized_input_variant)
            except Exception as e:
                print(f"WARNING: unable to normalize input variant {chrom}:{pos} {ref}>{alt}: {e}")
        else:
            result = run_UCSC_liftover_tool(hg, chrom, start, end, verbose=verbose)
    except Exception as e:
        return error_response(e)

    if warnings:
        result["warnings"] = warnings
    return Response(json.dumps({**params, **result}), mimetype='application/json')


@app.route('/', strict_slashes=False, defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    return "liftover api"


if __name__ == '__main__':
    app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
