#!/usr/bin/env python3
"""Decode one CNN-Mamba score bundle with AltAnn alone (UniAnn is no longer used).

The bundle holds plus/ and minus/ directories, each with sequence.fna, coding.tsv
and sites.tsv for one oriented sequence. For each orientation, `altann prepare`
writes decoder inputs from the calibrated probabilities without another
temperature, and `altann decode --exact-baseline` takes the global best path from
exact K=1 search and adds local K-best alternatives. An orientation without a
donor whose rescaled score is positive is skipped, as UniAnn's wrapper did.

Outputs: OUTPUT (every AltAnn model in genomic coordinates) with OUTPUT.tsv and
OUTPUT.json, and --ab-initio for INTEGRATE: the global-path (rank-0) models
without kbest_* attributes, so they read as baselines, plus alternatives whose
CDS chain differs from every retained model (merge_altann.py rules).
"""

import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_altann import merge  # noqa: E402


SEGMENT = re.compile(r"^(.*)__(\d+)-(\d+)-(\d+)$")
SITE_HEADER = "chrom\tpos\tstrand\ttype\tmotif\tprob"
# Match `perl -e 'print exp(1)'`, AltAnn's default multiplier.
DEFAULT_MULT = 2.71828182845905
REPORT_HEADER = ("segment\tstrand\treference_gene_index\trank\tgenes_in_path\tgene_in_path\t"
                 "reference_start_1based\treference_end_1based\tcandidate_start_1based\t"
                 "candidate_end_1based\tfixed_interval_score\treference_fixed_interval_score\t"
                 "score_delta\tactual_start_score\tintron_chain\treference_intron_chain\t"
                 "identical_to_reference\tclassification\twindow_start_1based\t"
                 "window_end_1based\tinitial_inter_len\tsplice_site_jaccard\n")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def read_fasta(path):
    identifier, chunks = None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                if identifier is not None or not line[1:].split():
                    raise ValueError(f"Expected exactly one named FASTA record: {path}")
                identifier = line[1:].split()[0]
            elif line.strip():
                if identifier is None:
                    raise ValueError(f"Sequence precedes FASTA header: {path}")
                chunks.append("".join(line.split()).upper())
    sequence = "".join(chunks)
    if not identifier or not sequence:
        raise ValueError(f"Empty FASTA: {path}")
    return identifier, sequence


def coordinates(identifier, length):
    match = SEGMENT.fullmatch(identifier)
    if match:
        chromosome, start, end, total = match.groups()
        start, end, total = int(start), int(end), int(total)
        if end - start != length or not 0 <= start < end <= total:
            raise ValueError(f"Segment header disagrees with sequence: {identifier}")
    else:
        chromosome, start, end, total = identifier, 0, length, length
    return dict(chromosome=chromosome, offset=start, end=end, total=total, length=length,
                segmented=bool(match))


def donor_maximum(path, identifier, length, multiplier):
    """Largest log(p * multiplier) over donors; UniAnn and AltAnn reject a non-positive one."""
    maximum = 0.0
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            if number == 1 and line.rstrip("\r\n") == SITE_HEADER:
                continue
            row = line.split()
            if len(row) not in (6, 7) or row[0] != identifier or row[2] != "+":
                raise ValueError(f"Site table must match the oriented FASTA and use + labels: {path}:{number}")
            if not 1 <= int(row[1]) <= length or row[3] not in {"donor", "acceptor", "start", "stop"}:
                raise ValueError(f"Invalid site position or type: {path}:{number}")
            probability = float(row[-1])
            if not math.isfinite(probability) or not 0 <= probability <= 1:
                raise ValueError(f"Invalid site probability: {path}:{number}")
            if row[3] == "donor":
                maximum = max(maximum, math.log(probability * multiplier + 1e-10))
    return maximum


def call(command, log, label):
    with Path(log).open("w") as handle:
        result = subprocess.run([str(part) for part in command], stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        raise ValueError(f"AltAnn {label} failed ({result.returncode}); see {log}")


def rank_zero_rows(path, destination, baseline):
    """Copy rank-0 (global-path) transcript blocks with decoded seqids.

    With baseline, kbest_* attributes are dropped so INTEGRATE's calibration reads
    the models as baselines (it identifies them by a missing kbest_delta).
    """
    rank, kept = {}, []
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            attrs = dict(item.split("=", 1) for item in fields[8].split(";") if "=" in item)
            if fields[2] in ("transcript", "mRNA"):
                rank[attrs["ID"]] = attrs.get("kbest_rank")
                if attrs.get("kbest_rank") != "0":
                    continue
                if baseline:
                    attrs = {key: value for key, value in attrs.items() if not key.startswith("kbest_")}
            elif all(rank.get(parent) != "0" for parent in attrs.get("Parent", "").split(",")):
                continue
            if baseline:
                fields[0] = unquote(fields[0])
            fields[8] = ";".join(f"{key}={value}" for key, value in attrs.items())
            kept.append("\t".join(fields) + "\n")
    Path(destination).write_text("##gff-version 3\n" + "".join(kept))


def run(args):
    if not 0 <= args.k <= 32767 or args.flank < 0 or args.overlap < 0 or args.margin < 0:
        raise ValueError("K must be between 0 and 32767; flank/overlap/margin must be nonnegative")
    if not math.isfinite(args.mult) or not 1 <= args.mult <= 100:
        raise ValueError("Multiplier must be between 1 and 100")
    source, output = args.input.resolve(), args.output.resolve()
    executable = args.altann_dir.resolve() / "bin" / "altann"
    if not executable.is_file():
        raise ValueError(f"AltAnn executable not found: {executable}")
    work = output.parent / (output.name + ".jobs")
    if work.exists() and any(work.iterdir()):
        raise ValueError(f"AltAnn job directory must be fresh: {work}")
    work.mkdir(parents=True, exist_ok=True)
    # K=0 keeps only the global path; the local search still needs one path per state.
    search_k = max(args.k, 1)
    metadata = dict(schema="chimann_altann_jobs_v2", source_bundle=str(source), decoder="altann",
                    k=args.k, search_k=search_k, flank=args.flank, multiplier=args.mult,
                    metadata_states=args.metadata_states, overlap=args.overlap, margin=args.margin,
                    jobs=[], completed=False)
    if (source / "manifest.json").is_file():
        metadata["prediction_manifest"] = str(source / "manifest.json")
    rows, report_rows, report_header, original = [], [], None, None
    for orientation in ("plus", "minus"):
        fasta = source / orientation / "sequence.fna"
        identifier, sequence = read_fasta(fasta)
        if orientation == "plus":
            original = (identifier, sequence)
        elif identifier != original[0] or sequence != original[1].translate(COMPLEMENT)[::-1]:
            raise ValueError("Minus FASTA must be the reverse complement of the plus FASTA with the same identifier")
        coords = coordinates(identifier, len(sequence))
        directory = work / orientation
        directory.mkdir()
        item = dict(name=orientation, original_strand="+" if orientation == "plus" else "-",
                    sequence_id=identifier, coordinate_mapping=coords)
        maximum = donor_maximum(source / orientation / "sites.tsv", identifier, len(sequence), args.mult)
        item["maximum_donor_log_score"] = maximum
        if maximum <= 0:
            item.update(status="skipped", reason="no_donor_with_positive_rescaled_score")
            metadata["jobs"].append(item)
            continue
        inputs = directory / "inputs"
        prepare = [sys.executable, executable, "prepare", "--fasta", fasta,
                   "--coding-frame", source / orientation / "coding.tsv",
                   "--scores", source / orientation / "sites.tsv", "--mult", repr(args.mult), "--output", inputs]
        call(prepare, directory / "prepare.log", f"prepare {orientation}")
        local_output = directory / "alternatives.gff3"
        decode = [sys.executable, executable, "decode", "--input", inputs, "--exact-baseline",
                  "--k", search_k, "--flank", args.flank, "--threads", args.threads,
                  "--seqid", coords["chromosome"], "--offset", coords["offset"],
                  "--sequence-length", coords["total"], "--short-cds-filter",
                  "--segment-overlap", args.overlap, "--segment-margin", args.margin,
                  "--spill-dir", directory, "--output", local_output]
        if args.metadata_states:
            decode.append("--metadata-states")
        if orientation == "minus":
            decode.append("--reverse")
        call(decode, directory / "altann.log", f"decode {orientation}")
        with local_output.open() as handle:
            rows.extend(line for line in handle if line.strip() and not line.startswith("#"))
        with Path(str(local_output) + ".tsv").open() as handle:
            header = handle.readline()
            if report_header is not None and report_header != header:
                raise ValueError("AltAnn report columns differ between strand jobs")
            report_header = header
            rank = header.rstrip("\n").split("\t").index("rank")
            report_rows.extend(line for line in handle if args.k or line.rstrip("\r\n").split("\t")[rank] == "0")
        item.update(status="complete", commands=[list(map(str, prepare)), list(map(str, decode))],
                    prepare=json.loads((inputs / "prepare.json").read_text()),
                    result=json.loads(Path(str(local_output) + ".json").read_text()))
        metadata["jobs"].append(item)
    output.write_text("##gff-version 3\n" + "".join(rows))
    Path(str(output) + ".tsv").write_text((report_header or REPORT_HEADER) + "".join(report_rows))
    with tempfile.TemporaryDirectory(dir=work) as temporary:
        baseline = Path(temporary) / "baseline.gff3"
        rank_zero_rows(output, baseline, baseline=True)
        if args.k == 0:
            # Only the global path is requested: drop the K=1 search's alternatives everywhere.
            rank_zero_rows(output, Path(temporary) / "global.gff3", baseline=False)
            output.write_text((Path(temporary) / "global.gff3").read_text())
        merged = merge(baseline, output)
    args.ab_initio.write_text("\n".join(merged) + "\n")
    metadata.update(completed=True, output=str(output), ab_initio=str(args.ab_initio.resolve()))
    Path(str(output) + ".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=Path, help="CNN_SCORE bundle with plus/ and minus/")
    parser.add_argument("--altann-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--ab-initio", required=True, type=Path)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("--threads", default="1")
    parser.add_argument("--mult", type=float, default=DEFAULT_MULT)
    parser.add_argument("--metadata-states", action="store_true",
                        help="Exact local K-best (transition metadata in the state)")
    parser.add_argument("--overlap", type=int, default=4_000_000)
    parser.add_argument("--margin", type=int, default=20_000)
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
