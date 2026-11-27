#!/usr/bin/env python3
"""Decode alternatives from the exact oriented UniAnn scores and traceback."""

import argparse
import json
from pathlib import Path
import subprocess
import sys


REPORT_HEADER = ("segment\tstrand\treference_gene_index\trank\tgenes_in_path\tgene_in_path\t"
                 "reference_start_1based\treference_end_1based\tcandidate_start_1based\t"
                 "candidate_end_1based\tfixed_interval_score\treference_fixed_interval_score\t"
                 "score_delta\tactual_start_score\tintron_chain\treference_intron_chain\t"
                 "identical_to_reference\tclassification\twindow_start_1based\t"
                 "window_end_1based\tinitial_inter_len\tsplice_site_jaccard\n")


def run(args):
    if not 0 <= args.k <= 32767 or args.flank < 0 or args.overlap < 0 or args.margin < 0:
        raise ValueError("K must be between 0 and 32767; flank/overlap/margin must be nonnegative")
    source, output = args.input.resolve(), args.output.resolve()
    manifest = json.loads((source / "manifest.json").read_text())
    if manifest.get("schema") != "chimann_uniann_jobs_v1" or not manifest.get("completed"):
        raise ValueError("Expected a complete run_uniann.py job bundle")
    if (args.overlap, args.margin) != (manifest.get("overlap"), manifest.get("margin")):
        raise ValueError("AltAnn overlap/margin must match the UniAnn mapping policy")
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(schema="chimann_altann_jobs_v1", input_manifest=str(source / "manifest.json"),
                    k=args.k, flank=args.flank, overlap=args.overlap, margin=args.margin,
                    jobs=[], completed=False)
    rows, report_rows, report_header = [], [], None
    if args.k:
        executable = args.altann_dir.resolve() / "bin" / "altann"
        if not executable.is_file():
            raise ValueError(f"AltAnn executable not found: {executable}")
        work = output.parent / (output.name + ".jobs")
        if work.exists() and any(work.iterdir()):
            raise ValueError(f"AltAnn job directory must be fresh: {work}")
        work.mkdir(exist_ok=True)
        for item in manifest["jobs"]:
            if item["status"] == "skipped":
                metadata["jobs"].append(dict(name=item["name"], status="skipped", reason=item["reason"]))
                continue
            if item["status"] != "complete" or item["name"] not in ("plus", "minus"):
                raise ValueError("Invalid UniAnn job status/orientation")
            directory = work / item["name"]
            directory.mkdir()
            local_output = directory / "alternatives.gff3"
            command = [sys.executable, str(executable), "decode"]
            # Replay the native DP/BT trace when one was kept; otherwise recompute the path.
            native_trace = manifest.get("trace", True)
            for flag, key in (("--fasta", "fasta"), ("--emissions", "emissions"), ("--gt", "gt"),
                              ("--ag", "ag"), ("--atg", "atg"), ("--stop", "stop"),
                              *((("--log", "log"),) if native_trace else ()),
                              ("--gff", "reference_gff")):
                path = source / item[key]
                if not path.is_file():
                    raise ValueError(f"Missing UniAnn {key}: {path}")
                command.extend([flag, str(path)])
            if not native_trace:
                command.append("--rerun-viterbi")
            coords = item["coordinate_mapping"]
            command.extend(["--k", str(args.k), "--flank", str(args.flank), "--threads", str(args.threads),
                            "--seqid", coords["chromosome"], "--offset", str(coords["offset"]),
                            "--sequence-length", str(coords["total"]), "--short-cds-filter",
                            "--segment-overlap", str(args.overlap), "--segment-margin", str(args.margin),
                            "--output", str(local_output)])
            if item["name"] == "minus":
                command.append("--reverse")
            with (directory / "altann.log").open("w") as handle:
                result = subprocess.run(command, cwd=directory, stdout=handle, stderr=subprocess.STDOUT)
            if result.returncode:
                raise ValueError(f"AltAnn {item['name']} failed ({result.returncode}); see {directory / 'altann.log'}")
            with local_output.open() as handle:
                rows.extend(line for line in handle if line.strip() and not line.startswith("#"))
            with Path(str(local_output) + ".tsv").open() as handle:
                header = handle.readline()
                if report_header is not None and report_header != header:
                    raise ValueError("AltAnn report columns differ between strand jobs")
                report_header = header
                report_rows.extend(handle)
            local_metadata = json.loads(Path(str(local_output) + ".json").read_text())
            metadata["jobs"].append(dict(name=item["name"], status="complete", command=command,
                                         result=local_metadata))
    else:
        metadata["reason"] = "k_zero_disables_alternative_decoding"
    output.write_text("##gff-version 3\n" + "".join(rows))
    Path(str(output) + ".tsv").write_text((report_header or REPORT_HEADER) + "".join(report_rows))
    metadata["completed"] = True
    Path(str(output) + ".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--altann-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("--threads", default="1")
    parser.add_argument("--overlap", type=int, default=4_000_000)
    parser.add_argument("--margin", type=int, default=20_000)
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
