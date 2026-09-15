#!/usr/bin/env python3
"""Minus-strand support for UniAnn (which decodes the + strand only):
run it on the reverse complement, then map results back.

  strand_tools.py revcomp  seq.fa            > seq.rc.fa
  strand_tools.py sites_rc seq.fa sites.tsv  > sites.rc.tsv   (- rows -> + rows in rc coords)
  strand_tools.py flip_gff seq.fa in.gff     > out.gff        (rc coords -> genome coords, strand -)
"""
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def read_one(path):
    cid, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if cid is not None:
                raise SystemExit("expected a single-sequence FASTA")
            cid = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return cid, "".join(seq)


def revcomp(fa):
    cid, seq = read_one(fa)
    rc = seq.translate(COMP)[::-1]
    print(f">{cid}")
    for i in range(0, len(rc), 80):
        print(rc[i:i + 80])


def sites_rc(fa, tsv):
    """sitescore '-' rows carry pos = L - p (p = 0-based index in the rc sequence),
    so the rc 1-based position is L - pos + 1."""
    _, seq = read_one(fa)
    L = len(seq)
    with open(tsv) as fh:
        print(next(fh), end="")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[2] == "-":
                f[1], f[2] = str(L - int(f[1]) + 1), "+"
                print("\t".join(f))


def flip_gff(fa, gff):
    _, seq = read_one(fa)
    L = len(seq)
    for line in open(gff):
        if line.startswith("#") or not line.strip():
            print(line, end="")
            continue
        f = line.rstrip("\n").split("\t")
        s, e = int(f[3]), int(f[4])
        f[3], f[4] = str(L - e + 1), str(L - s + 1)
        f[6] = "-" if f[6] == "+" else "+" if f[6] == "-" else f[6]
        print("\t".join(f))


if __name__ == "__main__":
    cmd, *args = sys.argv[1:]
    {"revcomp": revcomp, "sites_rc": sites_rc, "flip_gff": flip_gff}[cmd](*args)
