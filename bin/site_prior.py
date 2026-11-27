#!/usr/bin/env python3
"""Raise start/stop probabilities in a sites.tsv at CDS ends of protein alignments.

  site_prior.py segment.fa sites.tsv --gff ALIGN.gff [--gff ...] [--prob 0.9] [--min-identity 90] > sites.tsv

For every model in the GFFs (miniprot alignments: genes with identity=, or EviAnn's
best_unused_proteins), the first CDS base is a start candidate and the last CDS base
(stop codon included or excluded, both tried) a stop candidate. Positions are mapped into
the segment (<seqid>__<start>-<end>-<len>, 0-based half-open), kept only if the sequence
carries the motif on that strand (sites.tsv convention, see check_sites_tsv.py), and their
prob becomes max(prob, --prob); candidates absent from the table are appended.
Mouse chr16, UniAnn alone: CDS F1 18.4 -> 22.5 with EviAnn's best_unused_proteins (2026-09-21).
"""
import argparse
import collections
import re
import sys

SEG = re.compile(r"^(.*)__(\d+)-(\d+)-(\d+)$")
COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    cid, chunks = None, []
    for line in open(path):
        if line.startswith(">"):
            cid = line[1:].split()[0]
        else:
            chunks.append(line.strip().upper())
    return cid, "".join(chunks)


def cds_ends(gff, chrom, min_identity):
    """(genome pos, strand, type) candidates for one sequence."""
    cds, ident = collections.defaultdict(list), {}
    for line in open(gff):
        if not line.startswith(chrom + "\t"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        if f[2] == "CDS":
            cds[re.search(r"Parent=([^;]+)", f[8]).group(1)].append((int(f[3]), int(f[4]), f[6]))
        else:
            m = re.search(r"identity=([\d.]+)", f[8])
            if m:
                ident[re.search(r"ID=([^;]+)", f[8]).group(1)] = float(m.group(1))
    ends = set()
    for parent, c in cds.items():
        if min_identity is not None and parent in ident and ident[parent] < min_identity:
            continue
        c.sort()
        s, e, strand = c[0][0], c[-1][1], c[0][2]
        if strand == "+":
            ends.update({(s, "+", "start"), (e - 2, "+", "stop"), (e + 1, "+", "stop")})
        else:
            ends.update({(e, "-", "start"), (s + 2, "-", "stop"), (s - 1, "-", "stop")})
    return ends


def motif_at(seq, pos, strand, k=3):
    L = len(seq)
    if strand == "+":
        return seq[pos - 1:pos - 1 + k]
    i = L - pos
    return seq[max(0, L - i - k):L - i][::-1].translate(COMP)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("sites")
    ap.add_argument("--gff", action="append", required=True)
    ap.add_argument("--prob", type=float, default=0.9)
    ap.add_argument("--min-identity", type=float, default=None)
    a = ap.parse_args()
    sid, seq = read_fasta(a.fasta)
    m = SEG.match(sid)
    chrom, start0 = (m.group(1), int(m.group(2))) if m else (sid, 0)
    want = {}
    for gff in a.gff:
        for gpos, strand, typ in cds_ends(gff, chrom, a.min_identity):
            p = gpos - start0
            if not 1 <= p <= len(seq):
                continue
            mo = motif_at(seq, p, strand)
            if (typ == "start" and mo == "ATG") or (typ == "stop" and mo in STOPS):
                want[(p, strand, typ)] = mo
    seen, boosted = set(), 0
    out = sys.stdout
    with open(a.sites) as fh:
        out.write(fh.readline())
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (int(f[1]), f[2], f[3])
            if key in want:
                seen.add(key)
                if float(f[5]) < a.prob:
                    f[5] = str(a.prob)
                    boosted += 1
            out.write("\t".join(f) + "\n")
    added = 0
    for key, mo in want.items():
        if key not in seen:
            out.write(f"{sid}\t{key[0]}\t{key[1]}\t{key[2]}\t{mo}\t{a.prob}\n")
            added += 1
    print(f"site_prior: {sid}: {len(want)} alignment-supported sites, {boosted} raised, {added} added", file=sys.stderr)


if __name__ == "__main__":
    main()
