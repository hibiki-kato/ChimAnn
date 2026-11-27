<p align="center">
  <img src="docs/assets/logo.svg" width="29%" align="middle" hspace="20" />
  <img src="docs/assets/title.svg" width="50%" align="middle" />
</p>

ChimAnn (Chimeric Annotation software, pronounced Kye-mahn) is experimental genome annotation software combining evidence-based models with UniAnn best paths and AltAnn alternatives.

### Pipeline

```text
EviAnn evidence annotation
    -> CNN-Mamba training on EviAnn pseudo-labels, or supplied checkpoint
    -> coding-frame and site probabilities on both oriented strands
    -> UniAnn emissions, transitions and best-path traceback
    -> AltAnn alternatives using the same processed scores
    -> upstream EviAnn external-CDS merge and final annotation
```

`main.nf` is a Nextflow DSL2 workflow. [Stage contracts](docs/contracts.md) describe score formats and coordinate transformations.

| Stage | Purpose |
| --- | --- |
| `EVIANN` | Assemble RNA evidence, align related-species proteins and produce the initial annotation. Upstream EviAnn uses its normal Markov splice filtering and annotation rules. |
| `CNN_MODEL` | Train CNN-Mamba using only EviAnn pseudo-labels; `cnn_checkpoint` reuses a supplied model. |
| `CNN_SCORE` | Score each segment and its reverse complement with the same model, exporting joint coding-frame and donor/acceptor/start/stop probabilities. |
| `UNIANN` | Convert the probabilities once, run native Viterbi and retain the processed scores and traceback for AltAnn. |
| `ALTANN` | Generate existing K-best alternatives; check reference CDS chains against UniAnn before merging candidates. |
| `INTEGRATE` | Rerun upstream EviAnn with `-c …`, retaining its normal external-CDS integration rules. |

### Installation

Nextflow and a Java runtime launch the workflow. UniAnn needs a C++ compiler; AltAnn also needs CMake.

```bash
git clone https://github.com/hibiki-kato/ChimAnn.git
cd ChimAnn
git submodule update --init uniann eviann altann cnn-mamba
conda env create -f environment.yml
conda activate chimann
python3 bin/install_eviann.py --source eviann --outdir eviann/bin
(cd uniann && ./install.sh)
cmake -S altann -B altann/build -DCMAKE_BUILD_TYPE=Release
cmake --build altann/build --parallel 16
```

The annotation submodules use upstream [EviAnn](https://github.com/alekseyzimin/eviann), upstream [UniAnn](https://github.com/alekseyzimin/UniAnn), and [AltAnn](https://github.com/hibiki-kato/AltAnn). `install_eviann.py` creates a local `eviann/bin` runtime from the unmodified upstream source and resolves its required tools from the annotation environment on `PATH`; `eviann_src` points to that installed runtime. It does not patch upstream files. After updating dependencies with `git submodule update --init --remote uniann eviann altann cnn-mamba`, repeat the installation/build commands before running again. `environment.yml` supplies annotation tools, including SNAP and TransDecoder dependencies. Use `-profile conda` to create the annotation environment, or `-profile local_envs --annotation_env /abs/path/annotation_env` to reuse one.

CNN-Mamba is the `cnn-mamba` submodule, a clean upstream checkout pinned to `809f743`, selected by `cnn_mamba_dir`. Keep experimental edits in other checkouts separate; existing checkpoint and calibration artifacts may still reside there. ChimAnn validates coding tables and strips the header before invoking UniAnn. For converters using legacy `(n-i-1)%3` indexing, it stages frame columns as `[frame1, frame2, frame0]` to obtain the biologically verified `(n-i)%3` assignment. Corrected converters receive unchanged frame order. The manifest records this adapter; probabilities and upstream sources are unchanged.
CNN-Mamba runs through `cnn_python`, an absolute path to Python in an existing GPU environment with CUDA-enabled PyTorch, Mamba, NumPy, pandas and pysam. Follow CNN-Mamba's installation instructions for its compiled CUDA dependencies. The workflow imports its source from `cnn_mamba_dir`; it does not install drivers or rebuild that environment. Verify GPU access in that Python environment before running.

### Usage

All input and environment paths should be absolute. An existing EviAnn run can be reused directly:

```bash
nextflow run . -profile local_envs \
  --annotation_env /abs/path/annotation_env \
  --genome /abs/path/genome.fna \
  --rnaseq /abs/path/original/rnaseq.list \
  --proteins /abs/path/original/protein.faa \
  --eviann_resume_dir /abs/path/completed/eviann_run \
  --eviann_reuse_annotation true \
  --cnn_mamba_dir /abs/path/ChimAnn/cnn-mamba \
  --cnn_python /abs/path/gpu_env/bin/python \
  --cnn_checkpoint /abs/path/joint_frame_checkpoint.pt \
  --outdir /abs/path/results -work-dir /scratch/chimann_work -resume
```

Keep the original genome and protein basenames when resuming: EviAnn locates cached alignments by those names. The original RNA manifest and protein FASTA are still required by its final integration driver. The run directory is copied without BAMs; existing alignments and assembled transcripts are reused. The final pass executes upstream EviAnn's normal merge. Reusing a first-pass annotation retains its original producer's settings; unset `eviann_reuse_annotation` to recompute it with upstream EviAnn.

Training is the default when `cnn_checkpoint` is unset. Replace the checkpoint option above with explicit, disjoint validation and test chromosome lists, leaving at least one chromosome for training. For the three nuclear S. pombe chromosomes, for example:

```bash
  --cnn_validation_chroms NC_003423.3 \
  --cnn_test_chroms NC_003421.2
```

Both lists identify FASTA records, not windows. All training, validation and test targets come from the first-pass EviAnn annotation. Held-out results measure agreement with these pseudo-labels, not accuracy against a reference annotation. The original reference annotation is not a training input.

| Parameter | Meaning | Default |
| --- | --- | --- |
| `genome`, `rnaseq`, `proteins` | Genome FASTA, EviAnn RNA manifest and related-species proteins | Genome required; supply evidence paths for EviAnn |
| `eviann_args` | Additional upstream EviAnn options | empty |
| `eviann_resume_dir` | Completed EviAnn directory used to seed the first pass | unset |
| `eviann_reuse_annotation` | Keep its completed annotation instead of rerunning the first-pass merge | `false` |
| `cnn_mamba_dir` | Clean upstream CNN-Mamba checkout | `cnn-mamba` |
| `cnn_python` | Python executable in the prepared CUDA environment; an absolute path is recommended | `python` |
| `cnn_checkpoint` | Joint four-class checkpoint; skips training | unset |
| `cnn_validation_chroms`, `cnn_test_chroms` | Explicit chromosome holdouts required for training | unset |
| `cnn_window_size`, `cnn_stride` | Model window and overlap stride | `10000`, `5000` |
| `cnn_train_args` | Additional CNN-Mamba training options | empty |
| `cnn_calibration_json` | Matching checkpoint's joint-logit temperature calibration | unset: temperature 1 |
| `cnn_fit_temperature` | Fit a joint temperature on held-out EviAnn pseudo-labels | `false` |
| `min_seq_len` | Shorter sequences retain EviAnn annotation only | `100000` |
| `segment_len`, `segment_overlap`, `segment_margin` | Optional ab initio segment size, overlap and internal edge margin; `segment_len: 0` keeps chromosomes whole | `0`, `4000000`, `20000` |
| `ab_initio_decoder` | Global best path decoder: `uniann` (native binary; about 0.2 GB per Mb with UniAnn PR #4) or `altann` (about 15 MB per Mb, about 3 times slower; needs AltAnn with `--uniann-gff`). Both give identical annotations | `uniann` |
| `uniann_trace` | Keep the native DP/BT log, about 0.5 GB per Mb per strand. Without it AltAnn recomputes the same path | `false` |
| `score_multiplier` | Shared UniAnn site probability multiplier, applied before logarithms; this is not temperature calibration | `e` |
| `altann_kbest` | Paths retained per state; `0` disables alternatives | `10` |
| `altann_flank` | Gene-local flank in bases | `1000` |
| `kbest_delta_fraction` | Optional evidence-confirmed alternative score-drop filter | unset: off |
| `integrate_args` | Additional upstream external-CDS options | empty |
| `integrate_novel_only` | Legacy prefilter limiting candidates to loci without evidence mRNA | `false` |

Temperature calibration is applied once to the CNN joint four-class logits, before softmax. It does not calibrate site heads. UniAnn and AltAnn consume the resulting scores without additional temperature scaling. Supplied calibration must belong to the selected checkpoint; fitting on EviAnn holdouts measures calibration against pseudo-labels. Avoid fitting again when supplying an existing calibration.

AltAnn's K-Viterbi algorithm and candidate-ranking rules are unchanged. `merge_altann.py` checks every retained AltAnn reference CDS chain, strand and phase against UniAnn, then adds distinct complete alternatives. The optional `kbest_delta_fraction` filter concerns candidate score drops; it is separate from neural probability calibration. Without any evidence-confirmed alternatives, that optional filter drops all alternatives. Upstream EviAnn receives all combined predictions by default. Final ORF refinement and external, transcript-only and pseudogene rules can change or drop first-pass CDS. `integrate_novel_only` retains the older novel-locus prefilter.

The scoring adapter accepts inference checkpoints without `optimizer_step`, recording unavailable training progress as missing metadata. This affects provenance only, not logits or scores.

Outputs include the first-pass EviAnn annotation, CNN model/calibration provenance, oriented coding/site tables, UniAnn processed inputs and traceback, AltAnn candidates with score/provenance sidecars, and final `chimann.gff`. Add `-with-report report.html -with-trace trace.txt` for timing and memory. Use `-resume` to retain completed tasks; CNN and decoder source fingerprints invalidate cached scores after code changes.

The old PSAURON and sitescore files remain historical code; they are not active stages or dependencies. Installation initializes only `uniann`, `eviann` and `altann`, leaving unused legacy worktrees untouched. A successful model-scoring or decoder check alone does not establish complete EviAnn integration or annotation accuracy.

The October 16 S. pombe validation used the earlier EviAnn fork, including its CNN splice filter and `--untrusted-cds`; those metrics are historical and do not validate this upstream configuration. Upstream EviAnn lacks `--no-snap`, `--splice-model` and `--untrusted-cds`: the current workflow uses its normal Markov filter and plain `-c` integration, without the fork's extra external-CDS splice-filter step. The October 18 checks reran upstream EviAnn first-pass annotation and final integration on S. pombe using existing alignments and saved full-genome candidates, and verified clean-source CNN GPU inference on 30 kb plus both-strand decoder equivalence on 100 kb. Full-genome CNN scoring through the updated workflow was not regenerated. Shared CDS biology and EviAnn retention/reproducibility limitations remain.
