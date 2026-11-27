#!/usr/bin/env python3
"""Run `psauron -a` on a long sequence in overlapping chunks and stitch the
per-codon probabilities back into one psauron_score.csv row.

Chunking bounds PSAURON's memory use on long sequences.
Its `*_all_prob` columns list one probability
per translated amino acid, *stop codons omitted*; stitching is therefore done
in nucleotide space (codon start -> prob) and re-serialised at the end.
All six frame columns are stitched from each original-sequence chunk.
Reverse probabilities remain in reverse-complement frame order.

  psauron_chunked.py seq.fa [--chunk 4500000] [--overlap 3000] [--out psauron_score.csv]
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# PSAURON probability cells can contain millions of characters.
csv.field_size_limit(sys.maxsize)
DNA = set("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
COLS = ["description", "psauron_is_protein", "in_frame_score", "forward_frame2_score",
        "forward_frame3_score", "reverse_frame1_score", "reverse_frame2_score",
        "reverse_frame3_score", "mean_out_of_frame_score", "in_frame_all_prob",
        "forward_frame2_all_prob", "forward_frame3_all_prob", "reverse_frame1_all_prob",
        "reverse_frame2_all_prob", "reverse_frame3_all_prob"]


def read_one(path):
    cid, chunks = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if cid is not None:
                    raise ValueError("expected a single-sequence FASTA")
                identifiers = line[1:].split()
                if not identifiers:
                    raise ValueError("FASTA header has no sequence identifier")
                cid = identifiers[0]
            elif line.strip():
                if cid is None:
                    raise ValueError("sequence found before FASTA header")
                chunks.append("".join(line.split()))
    sequence = "".join(chunks).upper()
    if cid is None or not sequence:
        raise ValueError("empty FASTA sequence")
    if set(sequence) - DNA:
        raise ValueError("invalid DNA symbol in FASTA")
    return cid, sequence


def run_psauron(seq, workdir, tag):
    fa = Path(workdir) / f"{tag}.fa"
    out = Path(workdir) / f"{tag}.csv"
    fa.write_text(f">{tag}\n{seq}\n")
    subprocess.run(["psauron", "-i", str(fa), "-a", "-c", "-o", str(out)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with out.open() as handle:
        rows = csv.reader(handle)
        header = next((row for row in rows if row and row[0] == "description"), None)
        data = list(rows)
    if header is None or len(data) != 1 or not data[0] or data[0][0] != tag:
        raise ValueError(f"PSAURON did not score chunk {tag}; its first frame may contain no translated non-stop codons")
    row = dict(zip(header, data[0]))
    if any(column not in row for column in COLS[9:15]):
        raise ValueError(f"PSAURON chunk {tag} is missing six-frame probability columns; use -a")
    return [row[column].split(";") if row[column] else [] for column in COLS[9:15]]


def stop_mask(seq_arr, starts):
    """True where the codon starting at each position is a stop."""
    a, b, c = seq_arr[starts], seq_arr[starts + 1], seq_arr[starts + 2]
    t, g = ord("T"), ord("G")
    return (a == t) & (((b == ord("A")) & ((c == ord("A")) | (c == g))) | ((b == g) & (c == ord("A"))))


def stitch(sequence, chunk, overlap, workdir):
    """Return six global-frame arrays, retaining the first overlapping chunk."""
    if chunk < 6 or not 3 <= overlap < chunk:
        raise ValueError("chunk must be at least 6 bases and overlap must be at least 3 bases and smaller than chunk")
    length = len(sequence)
    normalized = sequence.upper().replace("N", "A")
    reverse = normalized.translate(COMPLEMENT)[::-1]
    # PSAURON maps other ambiguous symbols to A when computing its reverse complement.
    reverse = "".join(base if base in "ACGT" else "A" for base in reverse)
    oriented = [np.frombuffer(seq.encode("ascii"), dtype=np.uint8) for seq in (normalized, reverse)]
    probabilities = [np.full(length, np.nan, dtype=np.float32) for _ in range(2)]
    start = 0
    while start < length:
        end = min(start + chunk, length)
        frames = run_psauron(sequence[start:end], workdir, f"c{start}")
        if len(frames) != 6:
            raise ValueError(f"PSAURON chunk c{start} must provide six probability arrays")
        for direction, origin in enumerate((start, length - end)):
            for frame in range(3):
                starts = np.arange(origin + frame, origin + end - start - 2, 3)
                starts = starts[~stop_mask(oriented[direction], starts)]
                try:
                    values = np.asarray(frames[direction * 3 + frame], dtype=np.float32)
                except (TypeError, ValueError) as error:
                    raise ValueError(f"invalid PSAURON probabilities in chunk c{start}") from error
                if not np.all(np.isfinite(values) & (values >= 0) & (values <= 1)):
                    raise ValueError(f"PSAURON probabilities must be finite and between zero and one in chunk c{start}")
                if len(starts) and len(values) != len(starts):
                    raise ValueError(f"PSAURON probability count mismatch in chunk c{start}, frame {direction * 3 + frame}: expected {len(starts)}, found {len(values)}")
                # Native PSAURON uses placeholder protein scores for an empty alternate frame.
                if not len(starts):
                    continue
                unset = np.isnan(probabilities[direction][starts])
                probabilities[direction][starts[unset]] = values[unset]
        if end == length:
            break
        start = end - overlap
        print(f"[psauron_chunked] {end}/{length}", file=sys.stderr)
    columns = []
    for direction in range(2):
        for frame in range(3):
            starts = np.arange(frame, length - 2, 3)
            starts = starts[~stop_mask(oriented[direction], starts)]
            values = probabilities[direction][starts]
            if np.isnan(values).any():
                raise ValueError(f"missing PSAURON probability for frame {direction * 3 + frame}")
            columns.append(";".join(np.char.mod("%g", values)))
    return columns


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--chunk", type=int, default=4_500_000, help="chunk length (rounded to a multiple of 3)")
    ap.add_argument("--overlap", type=int, default=3000, help="overlap between chunks (multiple of 3)")
    ap.add_argument("--out", default="psauron_score.csv")
    a = ap.parse_args(argv)
    chunk = a.chunk - a.chunk % 3
    overlap = a.overlap - a.overlap % 3
    try:
        cid, seq = read_one(a.fasta)
        with tempfile.TemporaryDirectory(prefix="psauron_") as tmp:
            cols = stitch(seq, chunk, overlap, tmp)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        ap.error(str(error))
    with open(a.out, "w", newline="") as fh:
        fh.write(f"psauron_chunked.py {a.fasta} --chunk {chunk} --overlap {overlap}\n")
        fh.write("psauron score: NA\n")
        fh.write("NOTE: all six frames use original-sequence chunks\n")
        w = csv.writer(fh)
        w.writerow(COLS)
        w.writerow([cid, *["NA"] * 8, *cols])


if __name__ == "__main__":
    main()
