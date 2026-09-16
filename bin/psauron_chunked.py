#!/usr/bin/env python3
"""Run `psauron -a` on a long sequence in overlapping chunks and stitch the
per-codon probabilities back into one psauron_score.csv row.

psauron's memory grows linearly with sequence length (~1 GB per Mb), so whole
chromosomes do not fit a GPU. Its `*_all_prob` columns list one probability
per translated amino acid, *stop codons omitted*; stitching is therefore done
in nucleotide space (codon start -> prob) and re-serialised at the end.
Only the three forward-frame columns are filled: the pipeline scores the
minus strand by running on the reverse complement, and UniAnn's
preprocess_psauron_scores.pl reads forward frames only.

  psauron_chunked.py seq.fa [--chunk 4500000] [--overlap 3000] [--out psauron_score.csv]
"""
import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

csv.field_size_limit(sys.maxsize)   # all_prob cells run to millions of characters
STOPS = {"TAA", "TAG", "TGA"}
COLS = ["description", "psauron_is_protein", "in_frame_score", "forward_frame2_score",
        "forward_frame3_score", "reverse_frame1_score", "reverse_frame2_score",
        "reverse_frame3_score", "mean_out_of_frame_score", "in_frame_all_prob",
        "forward_frame2_all_prob", "forward_frame3_all_prob", "reverse_frame1_all_prob",
        "reverse_frame2_all_prob", "reverse_frame3_all_prob"]


def read_one(path):
    cid, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            if cid is not None:
                raise SystemExit("expected a single-sequence FASTA")
            cid = line[1:].split()[0]
        else:
            chunks.append(line.strip())
    return cid, "".join(chunks).upper()


def run_psauron(seq, workdir, tag):
    fa = Path(workdir) / f"{tag}.fa"
    out = Path(workdir) / f"{tag}.csv"
    fa.write_text(f">{tag}\n{seq}\n")
    subprocess.run(["psauron", "-i", str(fa), "-a", "-o", str(out)], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rows = list(csv.reader(open(out)))
    hdr = next(i for i, r in enumerate(rows) if r and r[0] == "description")
    row = dict(zip(rows[hdr], rows[hdr + 1]))
    return [row[c].split(";") if row[c] else [] for c in COLS[9:12]]


STOP_CODES = {"TAA", "TAG", "TGA"}


def stop_mask(seq_arr, starts):
    """True where the codon starting at each position is a stop."""
    a, b, c = seq_arr[starts], seq_arr[starts + 1], seq_arr[starts + 2]
    t, g = ord("T"), ord("G")
    return (a == t) & (((b == ord("A")) & ((c == ord("A")) | (c == g))) | ((b == g) & (c == ord("A"))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta")
    ap.add_argument("--chunk", type=int, default=4_500_000, help="chunk length (rounded to a multiple of 3)")
    ap.add_argument("--overlap", type=int, default=3000, help="overlap between chunks (multiple of 3)")
    ap.add_argument("--out", default="psauron_score.csv")
    a = ap.parse_args()
    chunk = a.chunk - a.chunk % 3
    overlap = a.overlap - a.overlap % 3
    cid, seq = read_one(a.fasta)
    L = len(seq)

    seq_arr = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    prob = np.full(L, np.nan, dtype=np.float32)      # indexed by absolute codon start; first chunk wins
    with tempfile.TemporaryDirectory(prefix="psauron_") as tmp:
        start = 0
        while start < L:
            end = min(start + chunk, L)
            probs = run_psauron(seq[start:end], tmp, f"c{start}")
            for f in range(3):
                starts = np.arange(start + f, end - 2, 3)
                keep = ~stop_mask(seq_arr, starts)
                vals = np.asarray(probs[f][: int(keep.sum())], dtype=np.float32)
                idx = starts[keep][: len(vals)]
                unset = np.isnan(prob[idx])
                prob[idx[unset]] = vals[unset]
            if end == L:
                break
            start = end - overlap
            print(f"[psauron_chunked] {end}/{L}", file=sys.stderr)

    cols = []
    for f in range(3):
        starts = np.arange(f, L - 2, 3)
        starts = starts[~stop_mask(seq_arr, starts)]
        vals = prob[starts]
        vals = vals[~np.isnan(vals)]
        cols.append(";".join(np.char.mod("%g", vals)))
    with open(a.out, "w", newline="") as fh:
        fh.write(f"psauron_chunked.py {a.fasta} --chunk {chunk} --overlap {overlap}\n")
        fh.write("psauron score: NA\n")
        fh.write("NOTE: forward frames only; reverse-frame columns are empty\n")
        w = csv.writer(fh)
        w.writerow(COLS)
        w.writerow([cid, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", *cols, "", "", ""])


if __name__ == "__main__":
    main()
