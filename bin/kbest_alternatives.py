#!/usr/bin/env python3
"""Turn UniAnn's gene-local k-best GFF into extra ab initio transcripts.

  kbest_alternatives.py <fasta>.kbest.gff --prefix ID_PREFIX [--min-cds 200] > alternatives.gff

Keeps the k >= 2 transcripts (k1 is the global Viterbi path, already in the
main UniAnn output), drops structures identical to one already emitted at the
locus (paths differing only in a neighbouring gene repeat the same transcript),
applies uniann.sh's MIN_CDS rule (1-2 exon transcripts with CDS <= min-cds are
dropped) and prefixes IDs so they are unique across segments and strands.
Output: transcript/exon/CDS rows in the layout of <fasta>.uniann.gff, no header
(meant to be appended to it).
"""
import argparse
import sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("gff")
ap.add_argument("--prefix", required=True)
ap.add_argument("--min-cds", type=int, default=200)
a = ap.parse_args()

tx = {}          # transcript id -> [transcript row, exon rows, cds rows]
order = []
for line in open(a.gff):
    if line.startswith("#") or not line.strip():
        continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] not in ("transcript", "exon", "CDS"):
        continue
    attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
    if p[2] == "transcript":
        tx[attrs["ID"]] = [p, [], []]
        order.append(attrs["ID"])
    else:
        tx[attrs["Parent"]][1 if p[2] == "exon" else 2].append(p)

seen, kept, dup, short, k1 = set(), 0, 0, 0, 0
for tid in order:
    t, exons, cds = tx[tid]
    locus = tid.split(".")[0]
    key = (locus, t[6], tuple((e[3], e[4]) for e in exons), tuple((c[3], c[4]) for c in cds))
    if ".k1." in tid:                      # the global path: record its structure, do not emit
        k1 += 1; seen.add(key); continue
    if key in seen:
        dup += 1; continue
    seen.add(key)
    cds_len = sum(int(c[4]) - int(c[3]) + 1 for c in cds)
    if len(exons) <= 2 and cds_len <= a.min_cds:
        short += 1; continue
    new_id = f"{a.prefix}.{tid}"
    t[5] = "."                             # path scores have 13+ digits and overflow gffread's score buffer
    t[8] = f"ID={new_id}"
    print("\t".join(t))
    for rows in (exons, cds):
        for r in rows:
            r[8] = f"Parent={new_id}"
            print("\t".join(r))
    kept += 1
print(f"kbest_alternatives: kept {kept} alternatives; skipped {k1} k1, {dup} duplicate, {short} short", file=sys.stderr)
