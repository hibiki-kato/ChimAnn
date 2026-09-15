# ChimAnn stage contracts

`main.nf` chains five stages. Each boundary is a plain file; this page pins
the formats so components can be swapped independently.

```
genome.fa ──> EVIANN ──> <genome>.pseudo_label.gff ──┐
    │                                                ├─> SITESCORE_TRAIN ─> model_dir
    ├─> SPLIT_GENOME ─> <seq>.fa ─┬─> PSAURON ─> psauron_score.csv ─┐
    │                             ├─> SITESCORE_SCORE ─> sites.tsv ───┼─> UNIANN ─> <seq>.fa.uniann.gff
    │                             └─────────────────────────────────┘        │
    └─> INTEGRATE = eviann.sh -c uniann.gff --untrusted-cds ─> chimann.gff <─┘
```

## EVIANN
- Command: `eviann.sh -t N -g genome.fa [-r rnaseq.list] [-p proteins.faa]`
- Output used: `<genome>.pseudo_label.gff` (evidence-based annotation),
  `<genome>.proteins.fasta`.
- `-r` list file must contain absolute paths (EviAnn reads it inside the work dir).

## PSAURON
- Command: `psauron -i seq.fa -a` (the `-a` flag is mandatory: UniAnn reads the
  per-frame per-base probability columns 10–15 of `psauron_score.csv`).
- Input: one sequence per FASTA.

## sitescore (`sitescore/` submodule)
[hibiki-kato/sitescore](https://github.com/hibiki-kato/sitescore): one `SiteModel`
interface (`train` / `load` / `score`), models registered as plug-ins
(`sitescore models`). convmamba is the first plug-in; LLM or other models add a
package + one entry-point line, nothing in ChimAnn changes.

- `sitescore train --model M --genome genome.fa --annotation eviann.gff --out model_dir [--init pretrained_dir]`
  trains (or fine-tunes) on EviAnn's preliminary annotation.
- `sitescore train` ends with a Platt fit on the validation sequences
  (`model_dir/calibration.json`); `sitescore score` emits calibrated probabilities.
- `sitescore score --model-dir model_dir --fasta seq.fa > sites.tsv`

### sites.tsv (consumed by `uniann.sh -s`)
Tab-separated, header row, 1-based coordinates, probabilities in (0, 1]:

```
chrom  pos  strand  type      motif  prob
X      7    +       donor     GT     3.43e-05
```
`type` ∈ {donor, acceptor, start, stop}. `uniann.sh` currently uses `+` strand
rows only and takes the **last** column as the score. Reference file:
`dev/UniAnn/data/dmel/chrX_sites.tsv`.

## UNIANN (per sequence **and strand**)
UniAnn decodes the + strand only. `bin/strand_tools.py` makes the - strand a
second + run: `revcomp` the sequence, `sites_rc` remaps sitescore's `-` rows
(rc pos = L - pos + 1), PSAURON runs on the rc sequence, and `flip_gff` maps
the result back (start' = L - end + 1, end' = L - start + 1, strand -).

- Command: `uniann.sh -f seq.fa -p psauron_score.csv -s sites.tsv -n`
- `uniann.sh` rescales scores by the best donor and aborts unless some donor has
  prob > 1/e; UNIANN emits an empty GFF for such sequences (e.g. mitochondria).
- Output: `seq.fa.uniann.gff`. UniAnn is the `uniann/` submodule
  ([hibiki-kato/UniAnn](https://github.com/hibiki-kato/UniAnn)); run its
  `install.sh` once to build `uniann/bin`. Override with `params.uniann_dir`.

## INTEGRATE (EviAnn, second pass)
- Command: `eviann.sh <same args as EVIANN> -c uniann.gff --untrusted-cds`,
  run in a copy of the first pass's directory so alignment is not repeated
  (`-c` clears `merge.success`, so EviAnn resumes from the merge step).
- UniAnn CDS are *low-trust*: they pass EviAnn's splice-site filtering before
  the final merge. Output: `<genome>.pseudo_label.gff` → `chimann.gff`.
