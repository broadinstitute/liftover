from datetime import datetime
import json
import os
import subprocess
import traceback
import tempfile

# flask imports
from flask import Flask, request, Response, send_from_directory
from flask_cors import CORS
from flask_talisman import Talisman


app = Flask(__name__)

CORS(app)

DEBUG = True  #if socket.gethostname() == "spliceai-lookup" else True
if not DEBUG:
    Talisman(app)

HG19_FASTA_PATH = "/ref/hg19.fa"
HG38_FASTA_PATH = "/ref/hg38.fa"
T2T_FASTA_PATH = "/ref/chm13v2.0.fa"

UCSC_LIFTOVER_TOOL = "UCSC liftover tool"
BCFTOOLS_LIFTOVER_TOOL = "bcftools liftover plugin"

LIFTOVER_EXAMPLE = f"/liftover/?hg=hg19-to-hg38&format=interval&chrom=chr8&start=140300615&end=140300620"

CHAIN_FILE_PATHS = {
    "hg19-to-hg38": "/hg19ToHg38.over.chain.gz",
    "hg38-to-hg19": "/hg38ToHg19.over.chain.gz",
    "hg38-to-t2t": "/hg38ToHs1.over.chain.gz", # replaced hg38-chm13v2.over.chain.gz based on advice from Giulio Genovese
    "t2t-to-hg38": "/hs1ToHg38.over.chain.gz", # replaced chm13v2-hg38.over.chain.gz based on advice from Giulio Genovese
}

LIFTOVER_REFERENCE_PATHS = {
    "hg19-to-hg38": (HG19_FASTA_PATH, HG38_FASTA_PATH),
    "hg38-to-hg19": (HG38_FASTA_PATH, HG19_FASTA_PATH),
    "hg38-to-t2t": (HG38_FASTA_PATH, T2T_FASTA_PATH),
    "t2t-to-hg38": (T2T_FASTA_PATH, HG38_FASTA_PATH),
}


def error_response(error_message, source=None):
    response_json = {"error": str(error_message)}
    if source:
        response_json["source"] = source
    return Response(json.dumps(response_json), status=200, mimetype='application/json')


def reverse_complement(seq):
    reverse_complement_map = dict(zip("ACGTN", "TGCAN"))
    return "".join([reverse_complement_map[n] for n in seq[::-1]])


def run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=False):
    if hg not in CHAIN_FILE_PATHS or hg not in LIFTOVER_REFERENCE_PATHS:
        raise ValueError(f"Unexpected hg arg value: {hg}")
    chain_file_path = CHAIN_FILE_PATHS[hg]
    source_fasta_path, destination_fasta_path = LIFTOVER_REFERENCE_PATHS[hg]

    with tempfile.NamedTemporaryFile(suffix=".vcf", mode="wt", encoding="UTF-8") as input_file, \
            tempfile.NamedTemporaryFile(suffix=".vcf", mode="rt", encoding="UTF-8") as output_file:

        #  command syntax: liftOver oldFile map.chain newFile unMapped
        if hg == "hg19-to-hg38":
            chrom = chrom.replace("chr", "")
        else:
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

                return {
                    "hg": hg,
                    "chrom": chrom,
                    "start": int(pos) - 1,
                    "end": pos,
                    "output_chrom": result_fields[0],
                    "output_pos": result_fields[1],
                    "output_ref": result_fields[3],
                    "output_alt": result_fields[4],
                    "liftover_tool": BCFTOOLS_LIFTOVER_TOOL,
                    #"output_strand": "-" if "SWAP=-1" in results else "+",
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
        result["output_ref"] = ref
        result["output_alt"] = alt
        #if result["output_strand"] == "-":
        #    result["output_ref"] = reverse_complement(result["output_ref"])
        #    result["output_alt"] = reverse_complement(result["output_alt"])
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


@app.route("/liftover/", methods=['POST', 'GET'])
def run_liftover():
    logging_prefix = datetime.now().strftime("%m/%d/%Y %H:%M:%S") + f" t{os.getpid()}"
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
        print(f"{logging_prefix}: {request.remote_addr}: ======================", flush=True)
        print(f"{logging_prefix}: {request.remote_addr}: {hg} liftover {format}: {chrom}:{variant_log_string}", flush=True)

    try:
        if format == "variant":
            result = run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=verbose)
        else:
            result = run_UCSC_liftover_tool(hg, chrom, start, end, verbose=verbose)
    except Exception as e:
        return error_response(str(e))

    result.update(params)

    return Response(json.dumps(result), mimetype='application/json')


@app.route('/', strict_slashes=False, defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    return "liftover api"

app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
