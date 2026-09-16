#!/usr/bin/env python3
"""Split a genome into overlapping segments for the per-sequence ab initio track,
and map results back.

  segment_genome.py segment genome.fa --size 25000000 --overlap 4000000 --min-len 100000 --outdir seqs
      one FASTA per segment, id  <seqid>__<start>-<end>-<seqlen>  (0-based half-open, on the + strand)
  segment_genome.py unsegment in.gff > out.gff
      restore seqid and coordinates; keep a transcript only from the segment that
      owns its midpoint (owner interval = segment minus half the overlap on each
      interior side) and only if it lies fully inside the segment with a margin

UniAnn's Viterbi and its PSAURON preprocessing hold whole-sequence arrays
(~0.6 GB per Mb), so mammalian chromosomes must be run in pieces.
"""
import argparse
import re
import sys
from pathlib import Path

SEG = re.compile(r"^(.*)__(\d+)-(\d+)-(\d+)$")


def read_fasta(path):
    cid, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cid is not None:
                    yield cid, "".join(chunks)
                cid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if cid is not None:
        yield cid, "".join(chunks)


def segment(a):
    out = Path(a.outdir)
    out.mkdir(exist_ok=True)
    n = 0
    for cid, seq in read_fasta(a.fasta):
        L = len(seq)
        if L < a.min_len:
            continue
        starts = [0]
        while starts[-1] + a.size < L:
            starts.append(starts[-1] + a.size - a.overlap)
        for s in starts:
            e = min(s + a.size, L)
            name = f"{cid}__{s}-{e}-{L}"
            with open(out / f"{name}.fa", "w") as fh:
                fh.write(f">{name}\n")
                for i in range(s, e, 80):
                    fh.write(seq[i : min(i + 80, e)] + "\n")
            n += 1
    if n == 0:
        sys.exit("segment_genome: no sequence >= --min-len")
    print(f"segment_genome: wrote {n} segments", file=sys.stderr)


def unsegment(a):
    half, margin = a.overlap // 2, a.margin
    keep_tx, rows = set(), []
    for line in open(a.gff):
        if line.startswith("#") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        m = SEG.match(p[0])
        if not m:
            sys.exit(f"unsegment: seqid {p[0]} is not a segment id")
        cid, s, e, L = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
        start, end = int(p[3]) + s, int(p[4]) + s  # genome coordinates, 1-based
        p[0], p[3], p[4] = cid, str(start), str(end)
        attrs = dict(x.split("=", 1) for x in p[8].split(";") if "=" in x)
        tid = attrs.get("ID") if p[2] in ("transcript", "mRNA") else attrs.get("Parent")
        if p[2] in ("transcript", "mRNA"):
            own_lo = 1 if s == 0 else s + half + 1
            own_hi = L if e == L else e - half
            in_lo = 1 if s == 0 else s + margin + 1
            in_hi = L if e == L else e - margin
            mid = (start + end) // 2
            if own_lo <= mid <= own_hi and start >= in_lo and end <= in_hi:
                keep_tx.add(tid)
        rows.append((tid, "\t".join(p)))
    n = 0
    print("##gff-version 3")
    for tid, row in rows:
        if tid in keep_tx:
            print(row)
            n += 1
    print(f"unsegment: kept {len(keep_tx)} transcripts", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("segment")
    s.add_argument("fasta")
    s.add_argument("--size", type=int, default=25_000_000)
    s.add_argument("--overlap", type=int, default=4_000_000)
    s.add_argument("--min-len", type=int, default=100_000)
    s.add_argument("--outdir", default="seqs")
    s.set_defaults(fn=segment)
    u = sp.add_parser("unsegment")
    u.add_argument("gff")
    u.add_argument("--overlap", type=int, default=4_000_000, help="must match the segment step")
    u.add_argument("--margin", type=int, default=20_000, help="drop transcripts this close to an interior segment end")
    u.set_defaults(fn=unsegment)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
