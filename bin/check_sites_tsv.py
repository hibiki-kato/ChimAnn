#!/usr/bin/env python3
"""Validate a sites.tsv against the sequence it scores (the site-scorer contract).

  check_sites_tsv.py seq.fa sites.tsv

Checks: header, column count, type/strand values, 1-based positions inside the
sequence, motif consistent with the sequence (donor GT / acceptor AG /
start ATG / stop TAA|TAG|TGA on the given strand), prob in (0, 1].
Exit 1 with a message on the first violation; prints a per-type row count otherwise.
"""
import collections
import sys

HEADER = ["chrom", "pos", "strand", "type", "motif", "prob"]
MOTIFS = {"donor": {"GT"}, "acceptor": {"AG"}, "start": {"ATG"}, "stop": {"TAA", "TAG", "TGA"}}
COMP = str.maketrans("ACGTN", "TGCAN")


def die(msg):
    sys.exit(f"check_sites_tsv: {msg}")


fa, tsv = sys.argv[1], sys.argv[2]
seqs, cid = {}, None
for line in open(fa):
    if line.startswith(">"):
        cid = line[1:].split()[0]
        seqs[cid] = []
    else:
        seqs[cid].append(line.strip().upper())
seqs = {c: "".join(p) for c, p in seqs.items()}
counts = collections.Counter()
with open(tsv) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    if hdr != HEADER:
        die(f"header must be {HEADER}, got {hdr}")
    for n, line in enumerate(fh, 2):
        f = line.rstrip("\n").split("\t")
        if len(f) != 6:
            die(f"line {n}: expected 6 columns, got {len(f)}")
        chrom, pos, strand, typ, motif, prob = f
        if chrom not in seqs:
            die(f"line {n}: unknown sequence {chrom}")
        if strand not in "+-" or typ not in MOTIFS:
            die(f"line {n}: bad strand/type {strand} {typ}")
        seq, L = seqs[chrom], len(seqs[chrom])
        p, k = int(pos), len(motif)
        if not 1 <= p <= L:
            die(f"line {n}: pos {p} outside 1..{L}")
        if strand == "+":
            got = seq[p - 1 : p - 1 + k]
        else:                      # pos = L - i (i: 0-based index in the reverse complement)
            i = L - p
            got = seq[max(0, L - i - k) : L - i][::-1].translate(COMP)
        if got != motif or motif not in MOTIFS[typ]:
            die(f"line {n}: motif {motif} for {typ} at {chrom}:{p}{strand} does not match sequence ({got})")
        pr = float(prob)
        if not 0.0 < pr <= 1.0:
            die(f"line {n}: prob {pr} not in (0, 1]")
        counts[(strand, typ)] += 1
if not counts:
    die("no rows")
print("check_sites_tsv: OK " + " ".join(f"{s}{t}={c}" for (s, t), c in sorted(counts.items())), file=sys.stderr)
