# ChimAnn example

`params.yaml` is a portable template. `params.spom_soct.yaml` uses the existing
S. pombe genome/RNA and S. octosporus protein paths, trains CNN-Mamba on EviAnn
labels with explicit chromosome holdouts, then runs both decoders and final EviAnn integration.
Install upstream EviAnn's local runtime after activating the annotation environment:

```bash
python3 bin/install_eviann.py --source eviann --outdir eviann/bin
```

Initialize the `cnn-mamba` submodule (or set `cnn_mamba_dir` to another clean Bowen checkout) and edit the environment paths before running:

```bash
nextflow run . -profile local_envs --annotation_env /abs/path/annotation_env \
  -params-file example/params.spom_soct.yaml -work-dir /scratch/chimann_spom -resume
```

For an already completed EviAnn run, add `--eviann_resume_dir /abs/path/eviann_run
--eviann_reuse_annotation true`. To reuse a joint-frame model instead of training,
add `--cnn_checkpoint /abs/path/model.pt`; optionally supply its matching
`--cnn_calibration_json`. Model artifacts may remain in another checkout; only source imports use the clean `cnn_mamba_dir`. See the root README for stage contracts.

The current workflow uses unmodified upstream EviAnn: its normal Markov splice filtering and plain `-c` external-CDS merge. Leave `eviann_args` and `integrate_args` empty unless supplying supported upstream options. The fork-only `--no-snap`, `--splice-model` and `--untrusted-cds` flags are unavailable, and the fork's extra external-CDS splice filtering is absent. October 16 validation metrics came from that earlier fork, not this upstream configuration.

The reference `input/S_pombe/annotation.gff3` is for evaluation only. Use `--strict-match -e 1000` for S. pombe; other species keep their existing evaluation settings. Report F1 as `2 * sensitivity * precision / (sensitivity + precision)` from the gffcompare statistics. Transcript uses the native Transcript row, gene uses the native Locus row, and CDS uses the Transcript row after CDS-as-exons conversion:

```bash
gffcompare --strict-match -e 1000 -r example/input/S_pombe/annotation.gff3 example/output/spom_soct/chimann.gff
```

Other species/ablation YAML files record the earlier PSAURON/sitescore workflow;
they need migration to the CNN parameters in `params.yaml` before use.
