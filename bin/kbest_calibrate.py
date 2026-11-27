#!/usr/bin/env python3
"""Filter AltAnn alternatives using an evidence-calibrated score cutoff.

  kbest_calibrate.py evidence.gff ab_initio.gff [--fraction 0.5] > ab_initio.filtered.gff

The combined input contains baseline transcripts (the global path) and AltAnn alternatives on both strands.
AltAnn's kbest_delta is its candidate path score minus the reference score over the same local interval.
An alternative is EviAnn-confirmed when its CDS chain matches a transcript at an evidence gene with >= 2 distinct CDS chains.
A baseline must also match a chain at that gene.
The cutoff is fraction * min(delta of confirmed alternatives); alternatives below it are removed.
Baseline transcripts, identified by the absence of kbest_delta, pass through unchanged.
Without confirmed alternatives, all alternatives are dropped.
"""
import argparse
import collections
import sys

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("evidence")
ap.add_argument("uniann")
ap.add_argument("--fraction", type=float, default=0.5, help="cutoff = fraction * min(confirmed delta)")
a = ap.parse_args()


def read(gff):
    """rows, {tid: (seqid, strand, cds chain)}, {tid: attrs}, {tid: parent}"""
    rows, meta, attrs, parent, cds = [], {}, {}, {}, collections.defaultdict(list)
    for line in open(gff):
        p = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(p) < 9:
            rows.append((None, line)); continue
        at = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
        if p[2] in ("mRNA", "transcript"):
            tid = at["ID"]; meta[tid] = (p[0], p[6]); attrs[tid] = at; parent[tid] = at.get("Parent", tid)
        elif p[2] == "CDS":
            for par in at.get("Parent", "").split(","):
                cds[par].append((int(p[3]), int(p[4])))
        rows.append((at.get("ID") if p[2] in ("mRNA", "transcript") else at.get("Parent"), line))
    chain = {tid: meta[tid] + (tuple(sorted(cds[tid])),) for tid in meta if tid in cds}
    return rows, chain, attrs, parent


_, ev_chain, _, ev_parent = read(a.evidence)
chain_gene = {c: ev_parent[t] for t, c in ev_chain.items()}
gene_chains = collections.defaultdict(set)
for c, g in chain_gene.items():
    gene_chains[g].add(c)
alt_genes = {g for g, cs in gene_chains.items() if len(cs) >= 2}

rows, un_chain, un_attrs, _ = read(a.uniann)
# Match UniAnn baselines to EviAnn genes by chromosome, strand, and CDS coordinates.
# Baselines lack kbest_delta; matching does not depend on UniAnn or AltAnn ID formats.
viterbi_genes = {chain_gene[c] for t, c in un_chain.items() if "kbest_delta" not in un_attrs[t] and c in chain_gene}
confirmed = []
for tid, c in un_chain.items():
    if "kbest_delta" not in un_attrs[tid]:
        continue
    g = chain_gene.get(c)
    if g and g in alt_genes and g in viterbi_genes:
        confirmed.append(float(un_attrs[tid]["kbest_delta"]))
if confirmed:
    cutoff = a.fraction * min(confirmed)
else:
    cutoff = float("inf")
keep = set()
n_alt = n_keep = 0
for tid, at in un_attrs.items():
    if "kbest_delta" not in at:
        keep.add(tid); continue
    n_alt += 1
    if float(at["kbest_delta"]) >= cutoff:
        keep.add(tid); n_keep += 1
for tid, line in rows:
    if tid is None or tid in keep:
        sys.stdout.write(line)
print(f"kbest_calibrate: {len(confirmed)} EviAnn-confirmed alternatives, min delta "
      f"{min(confirmed) if confirmed else 'n/a'}, cutoff {cutoff:.1f}; kept {n_keep}/{n_alt} alternatives",
      file=sys.stderr)
