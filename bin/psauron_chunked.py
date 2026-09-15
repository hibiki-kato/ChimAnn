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
import os
import subprocess
import sys
import tempfile
from pathlib import Path

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


def codon_probs(seq, frame, probs):
    """{codon start (0-based, absolute within seq): prob} for one forward frame."""
    out, it = {}, iter(probs)
    for p in range(frame, len(seq) - 2, 3):
        if seq[p:p + 3] in STOPS:
            continue
        try:
            out[p] = next(it)
        except StopIteration:
            break
    return out


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

    merged = [dict(), dict(), dict()]           # per frame: codon start -> prob
    with tempfile.TemporaryDirectory(prefix="psauron_") as tmp:
        start = 0
        while start < L:
            end = min(start + chunk, L)
            probs = run_psauron(seq[start:end], tmp, f"c{start}")
            for f in range(3):
                # local frame f corresponds to absolute frame (start + f) % 3
                for p, v in codon_probs(seq[start:end], f, probs[f]).items():
                    ap_ = start + p
                    merged[ap_ % 3].setdefault(ap_, v)    # first (upstream) chunk wins in overlaps
            if end == L:
                break
            start = end - overlap
            print(f"[psauron_chunked] {end}/{L}", file=sys.stderr)

    cols = []
    for f in range(3):
        m = merged[f]
        cols.append(";".join(m[p] for p in range(f, L - 2, 3) if seq[p:p + 3] not in STOPS and p in m))
    with open(a.out, "w", newline="") as fh:
        fh.write(f"psauron_chunked.py {a.fasta} --chunk {chunk} --overlap {overlap}\n")
        fh.write("psauron score: NA\n")
        fh.write("NOTE: forward frames only; reverse-frame columns are empty\n")
        w = csv.writer(fh)
        w.writerow(COLS)
        w.writerow([cid, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", *cols, "", "", ""])


if __name__ == "__main__":
    main()
