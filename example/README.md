# ChimAnn example

## Input
| dir | content |
| --- | --- |
| `S_pombe/` | target: `genome.fna` (ASM294v3), `annotation.gff3` (reference, evaluation only), `rna/fastq/` (5 paired RNA-seq runs, not in git), `rnaseq.list` |
| `S_octosporus/` | `protein.faa` (related-species evidence) |
| `S_cryophilus/` | `GCF_000004155.1_SCY4_protein.faa` (alternative evidence) |

## Run
```bash
nextflow run . -profile local_envs -params-file example/params.spom_soct.yaml \
    -work-dir ~/nf_work/chimann_spom -resume
```
Results land in `example/output/spom_soct/` (not in git). Compare with the
reference: `gffcompare -r example/input/S_pombe/annotation.gff3 example/output/spom_soct/chimann.gff`.
