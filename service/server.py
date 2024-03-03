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
}

LIFTOVER_REFERENCE_PATHS = {
    "hg19-to-hg38": (FASTA_PATHS["hg19"], FASTA_PATHS["hg38"]),
    "hg38-to-hg19": (FASTA_PATHS["hg38"], FASTA_PATHS["hg19"]),
    "hg38-to-t2t": (FASTA_PATHS["hg38"], FASTA_PATHS["t2t"]),
    "t2t-to-hg38": (FASTA_PATHS["t2t"], FASTA_PATHS["hg38"]),
}


def error_response(error_message):
    print(f"ERROR: {error_message}")
    return Response(json.dumps({"error": str(error_message)}), status=200, mimetype='application/json')


def reverse_complement(seq):
    reverse_complement_map = dict(zip("ACGTN", "TGCAN"))
    return "".join([reverse_complement_map[n] for n in seq[::-1]])


def get_user_ip(request):
    return request.environ.get("HTTP_X_FORWARDED_FOR")


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

    result.update(params)

    return Response(json.dumps(result), mimetype='application/json')


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

    try:
        if format == "variant":
            result = run_variant_liftover_tool(hg, chrom, pos, ref, alt, verbose=verbose)
            try:
                input_reference_genome = hg.split("-")[0]
                normalized_input_variant = run_bcftools_norm(input_reference_genome, chrom, pos, ref, alt, verbose=verbose)
                result.update(normalized_input_variant)
            except Exception as e:
                print(f"WARNING: unable to normalize input variant {chrom}:{pos} {ref}>{alt}: {e}")
        else:
            result = run_UCSC_liftover_tool(hg, chrom, start, end, verbose=verbose)
    except Exception as e:
        return error_response(e)

    result.update(params)

    return Response(json.dumps(result), mimetype='application/json')


@app.route('/', strict_slashes=False, defaults={'path': ''})
@app.route('/<path:path>/')
def catch_all(path):
    return "liftover api"

app.run(debug=DEBUG, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))
