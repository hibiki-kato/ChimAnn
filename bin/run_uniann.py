#!/usr/bin/env python3
"""Run UniAnn on already oriented CNN-Mamba prediction tables."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import quote


SEGMENT = re.compile(r"^(.*)__(\d+)-(\d+)-(\d+)$")
CODING_HEADER = ("chrom\tposition_1based\tstrand\tframe0_probability\t"
                 "frame1_probability\tframe2_probability\tnoncoding_probability\tpredicted_state")
SITE_HEADER = "chrom\tpos\tstrand\ttype\tmotif\tprob"
TRANSCRIPTS = {"transcript", "mRNA"}
# Match `perl -e 'print exp(1)'` in UniAnn and AltAnn's default exactly.
DEFAULT_MULT = 2.71828182845905
SCORES = {"emissions": "out.ps.txt", "gt": "out.gt.txt", "ag": "out.ag.txt",
          "atg": "out.atg.txt", "stop": "out.stop.txt", "log": "out.err"}


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
    return dict(chromosome=chromosome, offset=start, end=end, total=total,
                length=length, segmented=bool(match))


def converter_frame_columns(path):
    """Recognize UniAnn's frame assignment before adapting the local columns."""
    source = "".join(line.split("#", 1)[0] for line in Path(path).read_text().splitlines())
    source = re.sub(r"\s+", "", source)
    formulas = [formula for formula in re.findall(r"\$coding_frame\[([^]]+)\]=", source)
                if "$i" in formula]
    if formulas == ["($n-$i-1)%3"]:
        return formulas[0], (1, 2, 0)
    if formulas == ["($n-$i)%3"]:
        return formulas[0], (0, 1, 2)
    raise ValueError(f"Unrecognized UniAnn coding-frame assignment in {path}: {formulas}")


def prepare_coding(source, destination, identifier, length, frame_columns=(0, 1, 2)):
    """Remove the known header and adapt column order without changing values."""
    count = 0
    with Path(source).open() as handle, Path(destination).open("w") as output:
        for number, line in enumerate(handle, 1):
            if number == 1 and line.rstrip("\r\n") == CODING_HEADER:
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) not in (7, 8) or row[0] != identifier or row[2] != "+":
                raise ValueError(f"Coding table must match the oriented FASTA and use + labels: {source}:{number}")
            count += 1
            if int(row[1]) != count:
                raise ValueError(f"Coding positions must be contiguous and one-based: {source}:{number}")
            probabilities = list(map(float, row[3:7]))
            if any(not math.isfinite(p) or not 0 <= p <= 1 for p in probabilities):
                raise ValueError(f"Invalid coding probability: {source}:{number}")
            if frame_columns != (0, 1, 2):
                row[3:6] = [row[3 + column] for column in frame_columns]
                line = "\t".join(row) + "\n"
            output.write(line if line.endswith("\n") else line + "\n")
    if count != length:
        raise ValueError(f"Coding rows ({count}) differ from FASTA length ({length})")


def donor_maximum(path, identifier, length, multiplier):
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


def stage_runtime(directory, destination, altann_dir=None):
    """Use current source scripts and a decoder binary, never stale bin scripts.

    With altann_dir, AltAnn's decoder stands in for the native binary: same model and
    output, without its whole-sequence DP matrix or DP/BT log.
    """
    directory = Path(directory).resolve()
    destination.mkdir()
    for name in ("uniann.sh", "convert_coding_frame.pl"):
        source = next((p for p in (directory / "scripts" / name, directory / "bin" / name,
                                   directory / name) if p.is_file()), None)
        if source is None:
            raise ValueError(f"Missing UniAnn {name} in {directory}")
        shutil.copy2(source, destination / name)
        (destination / name).chmod(0o755)
    if altann_dir:
        core = Path(altann_dir).resolve() / "build" / "altann-core"
        if not (core.is_file() and os.access(core, os.X_OK)):
            raise ValueError(f"AltAnn decoder is not built: {core}")
        (destination / "uniann").write_text(f'#!/bin/sh\nexec "{core}" "$@" --uniann-gff\n')
        (destination / "uniann").chmod(0o755)
    else:
        native = next((p for p in (directory / "bin" / "uniann", directory / "src" / "uniann",
                                   directory / "uniann") if p.is_file() and os.access(p, os.X_OK)), None)
        if native is None:
            raise ValueError(f"UniAnn native binary is not installed in {directory}")
        (destination / "uniann").symlink_to(native)
    gffread = next((p for p in (directory / "bin" / "gffread", directory / "gffread")
                    if p.is_file() and os.access(p, os.X_OK)), None)
    if gffread:
        (destination / "gffread").symlink_to(gffread)
    elif not shutil.which("gffread"):
        raise ValueError("gffread was not found in UniAnn or PATH")
    return destination / "uniann.sh"


def validate_outputs(directory, sequence, trace="full"):
    """Check emissions and return the final dp row.

    trace is "full" for the native DP/BT dump, "final" for AltAnn's single dp row, or
    "none" when the native dump is switched off; the final scores are then unknown.
    """
    length = len(sequence)
    count = 0
    with (directory / "out.ps.txt").open() as handle:
        for row in handle:
            fields = row.split()
            if len(fields) not in (6, 7) or int(fields[0]) != count:
                raise ValueError("UniAnn emissions must contain one ordered row per base")
            if any(not math.isfinite(float(value)) for value in fields[1:6]):
                raise ValueError("Nonfinite UniAnn emission")
            if len(fields) == 7 and (count >= length or fields[6].upper() != sequence[count]):
                raise ValueError("Emission nucleotide disagrees with FASTA")
            count += 1
    if count != length:
        raise ValueError(f"UniAnn emitted {count} rows for {length} bases")
    count, final_scores = 0, None
    with open(os.devnull if trace == "none" else directory / "out.err") as handle:
        for line in handle:
            row = line.split()
            if len(row) < 2 or row[1] not in ("dp", "bt"):
                continue
            if trace == "final":
                if len(row) != 9 or row[1] != "dp" or int(row[0]) != length - 1:
                    raise ValueError("Expected only the final dp row from the AltAnn decoder")
                final_scores = list(map(int, row[2:]))
                continue
            if len(row) != 9 or int(row[0]) != count // 2 or row[1] != ("dp" if count % 2 == 0 else "bt"):
                raise ValueError("Incomplete or out-of-order UniAnn DP/BT dump")
            values = list(map(int, row[2:]))
            if row[1] == "dp":
                final_scores = values
            elif any(value != -1 if count // 2 == 0 else not 0 <= value < 7 for value in values):
                raise ValueError("Invalid UniAnn traceback predecessor")
            count += 1
    if trace == "final" and final_scores is None:
        raise ValueError("The AltAnn decoder wrote no final dp row")
    if trace == "full" and count != 2 * length:
        raise ValueError(f"UniAnn DP/BT dump has {count} rows; expected {2 * length}")
    for key in ("gt", "ag", "atg", "stop"):
        if not (directory / SCORES[key]).is_file():
            raise ValueError(f"UniAnn did not produce {SCORES[key]}")
    return final_scores


def read_gff(path, identifier, length):
    rows, ids = [], {}
    with Path(path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) != 9 or row[0] != identifier or row[6] != "+":
                raise ValueError(f"UniAnn GFF must match the oriented sequence: {path}")
            if not 1 <= int(row[3]) <= int(row[4]) <= length:
                raise ValueError(f"UniAnn GFF coordinates outside FASTA: {path}")
            if row[2] == "CDS" and row[7] not in ("0", "1", "2"):
                raise ValueError(f"Invalid UniAnn CDS phase: {path}")
            attributes = dict(item.split("=", 1) for item in row[8].split(";") if "=" in item)
            if "ID" in attributes:
                if attributes["ID"] in ids:
                    raise ValueError(f"Duplicate UniAnn GFF ID: {attributes['ID']}")
                ids[attributes["ID"]] = len(rows)
            rows.append((row, attributes))
    for row, attributes in rows:
        for parent in attributes.get("Parent", "").split(","):
            if parent and parent not in ids:
                raise ValueError(f"UniAnn GFF missing Parent: {parent}")
        if row[2] in TRANSCRIPTS and "ID" not in attributes:
            raise ValueError("UniAnn transcript has no ID")
    return rows


def owners(rows):
    by_id = {a["ID"]: (r, a) for r, a in rows if "ID" in a}
    def lookup(row, attributes, visited=()):
        if row[2] in TRANSCRIPTS:
            return {attributes["ID"]}
        result = set()
        for parent in attributes.get("Parent", "").split(","):
            if not parent:
                continue
            if parent in visited:
                raise ValueError("Cyclic GFF Parent relationship")
            result.update(lookup(*by_id[parent], (*visited, parent)))
        return result
    return [lookup(row, attributes) for row, attributes in rows]


def complete_reference(rows, sequence, final_scores):
    """Exclude only the native terminal partial model; retain nonterminal biology."""
    exclusions = []
    owned = owners(rows)
    models = {a["ID"]: r for r, a in rows if r[2] in TRANSCRIPTS}
    # Without final scores, the terminal model's own CDS below decides.
    if not models or (final_scores is not None and final_scores[0] == max(final_scores)):
        return rows, exclusions
    terminal = max(models, key=lambda key: int(models[key][4]))
    cds = sorted((int(r[3]), int(r[4])) for (r, _), parents in zip(rows, owned)
                 if r[2] == "CDS" and terminal in parents)
    if not cds:
        raise ValueError(f"Terminal UniAnn transcript has no CDS: {terminal}")
    coding = "".join(sequence[start - 1:end] for start, end in cds)
    if coding[-3:] in {"TAA", "TAG", "TGA"} and len(coding) % 3 == 0:
        return rows, exclusions
    exclusions.append(dict(transcript=terminal, reason="native_terminal_partial_CDS",
                           detail="Final Viterbi state is coding/intronic; terminal CDS lacks a complete stop or is truncated",
                           oriented_interval_1based=[int(models[terminal][3]), int(models[terminal][4])],
                           final_printed_dp_scores=final_scores, last_codon=coding[-3:], cds_length=len(coding)))
    keep = select_transcripts(rows, set(models) - {terminal})
    return keep, exclusions


def select_transcripts(rows, selected):
    """Keep selected transcript blocks and their gene ancestors."""
    by_id = {a["ID"]: index for index, (_, a) in enumerate(rows) if "ID" in a}
    retained = {index for index, parents in enumerate(owners(rows)) if parents & selected}
    pending = list(retained)
    while pending:
        index = pending.pop()
        for parent in rows[index][1].get("Parent", "").split(","):
            if parent and by_id[parent] not in retained:
                retained.add(by_id[parent])
                pending.append(by_id[parent])
    return [item for index, item in enumerate(rows) if index in retained]


def write_gff(path, rows):
    with Path(path).open("w") as handle:
        handle.write("##gff-version 3\n")
        for row, attributes in rows:
            fields = row.copy()
            fields[8] = ";".join(f"{key}={value}" for key, value in attributes.items()) or "."
            handle.write("\t".join(fields) + "\n")


def map_gff(rows, identifier, orientation, coords, overlap, margin):
    kept = set()
    def transform(row):
        row = row.copy()
        start, end = int(row[3]), int(row[4])
        if orientation == "minus":
            start, end = coords["length"] - end + 1, coords["length"] - start + 1
        row[0], row[3], row[4], row[6] = coords["chromosome"], str(start + coords["offset"]), str(end + coords["offset"]), "-" if orientation == "minus" else "+"
        return row
    for row, attributes in rows:
        if row[2] not in TRANSCRIPTS:
            continue
        mapped = transform(row)
        start, end = int(mapped[3]), int(mapped[4])
        lo = 1 if coords["offset"] == 0 else coords["offset"] + overlap // 2 + 1
        hi = coords["total"] if coords["end"] == coords["total"] else coords["end"] - overlap // 2
        inner_lo = 1 if coords["offset"] == 0 else coords["offset"] + margin + 1
        inner_hi = coords["total"] if coords["end"] == coords["total"] else coords["end"] - margin
        if not coords["segmented"] or lo <= (start + end) // 2 <= hi and start >= inner_lo and end <= inner_hi:
            kept.add(attributes["ID"])
    prefix = quote(f"{identifier}.{orientation}.uniann.", safe="._-")
    result = []
    for row, attributes in select_transcripts(rows, kept):
        attributes = attributes.copy()
        if "ID" in attributes:
            attributes["ID"] = prefix + quote(attributes["ID"], safe="._-")
        if "Parent" in attributes:
            attributes["Parent"] = ",".join(prefix + quote(p, safe="._-") for p in attributes["Parent"].split(","))
        result.append((transform(row), attributes))
    return result


def run(args):
    if not math.isfinite(args.mult) or not 1 <= args.mult <= 100 or args.overlap < 0 or args.margin < 0:
        raise ValueError("Multiplier must be between 1 and 100; overlap/margin must be nonnegative")
    source, outdir, output = args.input.resolve(), args.outdir.resolve(), args.output.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f"UniAnn job directory must be fresh: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    altann_dir = getattr(args, "altann_dir", None)
    wrapper = stage_runtime(args.uniann_dir, outdir / "runtime", altann_dir)
    frame_formula, frame_columns = converter_frame_columns(wrapper.parent / "convert_coding_frame.pl")
    environment = dict(os.environ, PATH=str(wrapper.parent) + os.pathsep + os.environ.get("PATH", ""))
    # A DP/BT log is about 0.5 GB per Mb per strand; AltAnn can recompute the same path.
    trace = "final" if altann_dir else ("none" if getattr(args, "no_trace", False) else "full")
    metadata = dict(schema="chimann_uniann_jobs_v1", decoder="altann" if altann_dir else "uniann",
                    trace=trace == "full",
                    source_bundle=str(source), multiplier=args.mult,
                    overlap=args.overlap, margin=args.margin, jobs=[], completed=False,
                    coding_frame_adapter=dict(input_convention="zero-based codon offsets on each oriented FASTA",
                                              converter_formula=frame_formula,
                                              staged_frame_columns=list(frame_columns),
                                              rotation="left_one" if frame_columns == (1, 2, 0) else "none"))
    if (source / "manifest.json").is_file():
        shutil.copy2(source / "manifest.json", outdir / "prediction.manifest.json")
        metadata["prediction_manifest"] = "prediction.manifest.json"
    combined = []
    original = None
    for orientation in ("plus", "minus"):
        directory = outdir / orientation
        directory.mkdir()
        fasta = directory / "sequence.fna"
        shutil.copy2(source / orientation / "sequence.fna", fasta)
        identifier, sequence = read_fasta(fasta)
        if orientation == "plus":
            original = (identifier, sequence)
        elif identifier != original[0] or sequence != original[1].translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]:
            raise ValueError("Minus FASTA must be the reverse complement of the plus FASTA with the same identifier")
        coords = coordinates(identifier, len(sequence))
        prepare_coding(source / orientation / "coding.tsv", directory / "coding.tsv", identifier,
                       len(sequence), frame_columns)
        shutil.copy2(source / orientation / "sites.tsv", directory / "sites.tsv")
        maximum = donor_maximum(directory / "sites.tsv", identifier, len(sequence), args.mult)
        item = dict(name=orientation, directory=orientation, original_strand="+" if orientation == "plus" else "-",
                    fasta=f"{orientation}/sequence.fna", coordinate_mapping=coords, sequence_id=identifier,
                    reference_gff=f"{orientation}/uniann.complete.gff3", raw_gff=f"{orientation}/sequence.fna.uniann.gff",
                    maximum_donor_log_score=maximum, exclusions=[])
        item.update({key: f"{orientation}/{name}" for key, name in SCORES.items()})
        raw = directory / "sequence.fna.uniann.gff"
        if maximum <= 0:
            item.update(status="skipped", reason="no_donor_with_positive_rescaled_score")
            write_gff(raw, [])
            write_gff(directory / "uniann.complete.gff3", [])
        else:
            # Older wrappers mask pipeline failures with a final success echo.
            # pipefail prevents their && chain from publishing a broken GFF;
            # validation below also requires that expected artifact to exist.
            command = ["bash", "-o", "pipefail", str(wrapper), "-f", "sequence.fna", "-c", "coding.tsv", "-s", "sites.tsv", "-m", str(args.mult)]
            if trace == "none":
                command.append("-n")
            with (directory / "uniann.log").open("w") as handle:
                result = subprocess.run(command, cwd=directory, env=environment, stdout=handle, stderr=subprocess.STDOUT)
            if result.returncode:
                raise ValueError(f"UniAnn {orientation} failed ({result.returncode}); see {directory / 'uniann.log'}")
            final_scores = validate_outputs(directory, sequence, trace)
            rows = read_gff(raw, identifier, len(sequence))
            rows, item["exclusions"] = complete_reference(rows, sequence, final_scores)
            write_gff(directory / "uniann.complete.gff3", rows)
            combined.extend(map_gff(rows, identifier, orientation, coords, args.overlap, args.margin))
            item.update(status="complete", command=command, emission_rows=len(sequence))
            if trace == "full":
                item.update(dp_bt_rows=2 * len(sequence))
        metadata["jobs"].append(item)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_gff(output, combined)
    metadata.update(completed=True, output=str(output))
    (outdir / "manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--uniann-dir", required=True, type=Path)
    parser.add_argument("--altann-dir", type=Path,
                        help="Decode with this AltAnn build instead of the native UniAnn binary")
    parser.add_argument("--no-trace", action="store_true",
                        help="Do not keep the native DP/BT log; AltAnn then recomputes the path")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mult", type=float, default=DEFAULT_MULT)
    parser.add_argument("--overlap", type=int, default=4_000_000)
    parser.add_argument("--margin", type=int, default=20_000)
    args = parser.parse_args(argv)
    try:
        run(args)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
