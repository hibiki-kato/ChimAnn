#!/usr/bin/env python3
"""Keep only ab initio transcripts that do not overlap any evidence-based mRNA
on the same strand.

EviAnn's -c path outputs external CDS only at loci where it annotated nothing
(combine_gene_protein_gff.pl: EviAnnE is skipped once a locus has output), but
before that rule fires the external alignments still displace real protein
alignments (containment, intron-chain dedup). Feeding only novel-locus CDS
avoids the displacement while keeping the intended gain.

  novel_cds.py evidence.gff ab_initio.gff > novel.gff
"""
import bisect
import collections
import sys

evidence, query = sys.argv[1], sys.argv[2]
iv = collections.defaultdict(list)
for line in open(evidence):
    p = line.rstrip("\n").split("\t")
    if len(p) > 8 and p[2] in ("mRNA", "transcript"):
        iv[(p[0], p[6])].append((int(p[3]), int(p[4])))
starts, ends = {}, {}
for k, L in iv.items():
    L.sort()
    starts[k] = [a for a, _ in L]
    m, ends[k] = 0, []
    for _, b in L:                      # running max of ends for a sweep
        m = max(m, b); ends[k].append(m)

def overlaps(k, s, e):
    if k not in starts: return False
    i = bisect.bisect_right(starts[k], e)   # candidates start <= e
    return i > 0 and ends[k][i - 1] >= s

keep, drop, tid_keep = 0, 0, set()
for line in open(query):
    p = line.rstrip("\n").split("\t")
    if len(p) > 8 and p[2] == "transcript":
        tid = [x for x in p[8].split(";") if x.startswith("ID=")][0][3:]
        if overlaps((p[0], p[6]), int(p[3]), int(p[4])): drop += 1
        else: keep += 1; tid_keep.add(tid)
for line in open(query):
    p = line.rstrip("\n").split("\t")
    if line.startswith("#") or len(p) < 9:
        sys.stdout.write(line); continue
    attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
    tid = attrs.get("ID") if p[2] == "transcript" else attrs.get("Parent")
    if tid in tid_keep: sys.stdout.write(line)
print(f"novel_cds: kept {keep} transcripts at loci without evidence mRNA, dropped {drop}", file=sys.stderr)
