# EviAnn issues observed while integrating ab initio (UniAnn) CDSs

State of upstream `alekseyzimin/eviann`: main = 920e8de (2026-09-18). Numbers below are from the
D. melanogaster EviAnn-paper setup (5 SRA runs, 5 related proteomes, `evaluate.sh` vs
`ref/annotation.fix.gff`) unless stated otherwise; the S. pombe fixture is the one used for PR #3.
Run directories are under `work/chimann/` and `work/eviann_untrusted_cds_rebase/`.

## 1. External CDSs are dropped wherever a protein-only *candidate* exists, and the locus ends up empty

Since 0aef9bd / 2c1e1f9 (`-c` handling rewritten), `eviann.sh` compares the external CDSs with
`gffcompare` against `$GENOME.k.gff` **plus `$GENOME.best_unused_proteins.gff`** and keeps only class
`u|p|o|x|s` (src/eviann.sh, "here we figure out which external CDSs overlap with complete annotations
and then remove them"). `best_unused_proteins` are miniprot-derived candidates, most of which never
reach the output, so an external CDS that shares a junction with such a candidate (class j/k/c/e/i)
is discarded and nothing replaces it.

Evidence (INTEGRATE reruns on the same first pass and the same UniAnn transcripts, Viterbi only,
novel-locus prefilter on, `work/chimann/Dmel_cnn_kbest/eval_ablation_cnn_k1` vs `eval_v3_novel_k1`):

| eviann | Transcript Sn / Pr / F1 | CDS Sn / Pr / F1 | Gene Sn / Pr / F1 | matched ref genes | of which UniAnn-contributed |
| --- | --- | --- | --- | --- | --- |
| 19ddde7 | 51.7 / 59.1 / 55.2 | 68.6 / 74.3 / 71.3 | 77.4 / 91.6 / 83.9 | 11860 | 502 |
| 2c1e1f9 | 51.1 / 58.9 / 54.7 | 67.8 / 74.2 / 70.9 | 76.4 / 91.5 / 83.3 | - | - |
| 920e8de | 51.1 / 58.9 / 54.7 | 67.8 / 74.3 / 70.9 | 76.4 / 91.5 / 83.3 | 11712 | 350 |

Of the 169 reference genes matched on 19ddde7 but not on 920e8de (21 gained the other way):
none is in the EviAnn-only output (all were contributed by the external CDSs); 163 overlap a miniprot
alignment (`palign.fixed.gff`) and 109 overlap a `best_unused_proteins` model; 91 overlap a first-pass
mRNA but only 7 on the same strand; in the 920e8de output 156 of the 169 loci have **no same-strand
mRNA at all** and 13 have an inexact model. Strand split 81 (-) / 88 (+).

A pure-upstream run (eviann 920e8de default options + UniAnn 91477a6 Viterbi, plain `-c`, all UniAnn
transcripts, `work/chimann/Dmel_upstream`): first pass CDS 67.3 / 74.1 / 70.5, gene 75.0 / 92.8 / 83.0 →
with UniAnn CDS 68.3 / 73.8 / 70.9, gene 77.1 / 91.5 / 83.7 (13255 genes, 1093 mRNAs with a UniAnn model as
evidence). The same UniAnn input on 19ddde7 gave gene 77.4 / 91.6 / 83.9 (with the novel-locus prefilter).

Suggested change: test externals against the final annotation only (`k.gff`), or defer the
external-vs-protein-candidate decision to the final selection so that the external survives when the
protein candidate is not output. Related: issue #4.

## 2. Plain `-c` outputs external CDSs with bad splice sites as exact models (PR #3)

S. pombe fixture (`work/shared/eviann/Spom/legacy-inputs`): 200 correct externals (CLEAN), 200 with perturbed
splice sites (ERROR), 50 single-exon (SINGLE), same alignments in every run.

| eviann | CLEAN kept | ERROR kept | SINGLE kept |
| --- | --- | --- | --- |
| 920e8de, `-c` (`run_main_920e8de.py`) | 190/200 | **24/200** | 50/50 |
| 2c1e1f9, `-c` | 191/200 | 24/200 | 50/50 |
| 0aef9bd, `-c` | - | 1/200 | - |
| 2c1e1f9 + PR #3, `-c --untrusted-cds` | - | 0/200 | - |

Externals of class `u|p|o|x|s` are added as `_EXTERNAL` alignments with identity 100 and skip the
splice-site (WAM) filter that every other CDS candidate goes through. Ab initio predictions are
exactly the kind of external input that needs that filter. PR #3 adds it as an opt-in flag; the
default output stays byte-identical to main.

## 3. `-c` either displaces evidence-based models or forfeits externals — the trade-off is unresolved (issue #4)

On 19ddde7, passing *all* UniAnn transcripts to `-c` lost 631 correct CDSs that the EviAnn-only run had
(467 loci vanished): (i) `combine_gene_protein_gff.pl` treats a transcript whose `cmp_ref` is an
`_EXTERNAL` as coding only for class k/=; (ii) real miniprot alignments are marked "contained" in the
`_EXTERNAL` or deduplicated by intron chain in `filter_unused_proteins.pl`; (iii) the final pass
outputs the EviAnnE model only if the locus has no other output (997a6b1) — so both vanish.
0aef9bd..920e8de fixed this by restricting externals to unannotated regions, which is issue 1.
Passing only externals at loci without an evidence mRNA (ChimAnn `bin/novel_cds.py`) was the
workaround for the first behaviour and still gains +0.2 CDS F1 on 920e8de.

## 4. The protein alignment stage cannot fail

```
miniprot ... 2>miniprot.err | convert_miniprot_gff.pl > $GENOME.$PROTEIN.uniq.palign.gff.tmp && \
mv ... && ... touch protein2genome.align.success
```

No `pipefail`: when miniprot is killed (OOM on the 372 MB rodent proteome, and on zebrafish), the
pipeline's exit status is that of `convert_miniprot_gff.pl`, the truncated GFF is moved into place
and `protein2genome.align.success` is written. `eviann.sh` continues, exits 0 and produces an
annotation without protein evidence (or 0 genes, when the truncated file breaks later stages); the
next invocation skips the alignment because the marker exists. Observed on Mmus (2026-09-16) and
Drer. Fork fix: hibiki-kato/eviann f97a9e9 (`set -o pipefail` around miniprot, marker only on
success). ChimAnn works around it by refusing an output with no mRNA.

## 5. Minor

- SNAP cannot be switched off (`AB_INITIO=1`, no flag); `-c` disables it implicitly. Fork: `--no-snap`
  (93d86c1), useful when another ab initio track is supplied.
- `src/` cannot be run from a checkout: `$MYPATH/{stringtie,TransDecoder.Predict,ufasta}` (and the SNAP
  binaries) must exist beside `eviann.sh`. Fork: ec9d216 looks them up on PATH.
- Resolved upstream: gene record overwrite when XLOC/AXLOC share a start (PR #5, fixed differently in 2c1e1f9).
