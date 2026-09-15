<p align="center">
  <img src="docs/assets/logo.svg" width="29%" align="middle" hspace="20" />
  <img src="docs/assets/title.svg" width="50%" align="middle" />
</p>

ChimAnn (Chimeric Annotation software, pronounced Kye-mahn) is a **Experimental** DNA annotation software that combines evidence-based and ab initio methods.

### Software Architecture
<p align="center">
  <img src="docs/assets/pipeline_overall.drawio.svg" width="50%" align="middle" hspace="20" />
</p>

- Alignment is performed using [minimap2](https://github.com/lh3/minimap2), [HISAT2](https://daehwankimlab.github.io/hisat2/), [StringTie](https://github.com/Transcriptome-Analysis-Tools/StringTie), and [miniprot](https://github.com/lh3/miniprot)
- Filtering is performed using [EviAnn](https://github.com/alekseyzimin/EviAnn).
- Ab initio gene prediction is performed using [UniAnn](https://github.com/hibiki-kato/UniAnn), fed by [PSAURON](https://github.com/salzberg-lab/PSAURON) emissions and [sitescore](https://github.com/hibiki-kato/sitescore) site scores.
- Integration is performed by EviAnn: UniAnn's CDS are fed back as low-trust external evidence (`-c … --untrusted-cds`).

### Pipeline
`main.nf` (Nextflow DSL2) runs the stages in order; every boundary is a file whose format is pinned in [docs/contracts.md](docs/contracts.md).

| Stage | Module | Tool |
| --- | --- | --- |
| Alignment + filtering | `modules/eviann.nf` | `eviann.sh` |
| Coding-potential emissions | `modules/psauron.nf` | `psauron -a` (per sequence) |
| Site scores (donor/acceptor/start/stop) | `modules/sitescore.nf` | `sitescore train/score` (model plug-in, convmamba by default) |
| Ab initio prediction | `modules/uniann.nf` | `uniann.sh` per sequence and strand (- via reverse complement, `bin/strand_tools.py`) |
| Integration | `modules/integrate.nf` | `eviann.sh -c uniann.gff --untrusted-cds` (resumes the EviAnn run) |

### Install
Nextflow is the launcher; each stage's tools come from conda environments.
```bash
conda install -c bioconda nextflow                 # launcher, once
git clone --recurse-submodules https://github.com/hibiki-kato/ChimAnn.git
cd ChimAnn
(cd uniann && ./install.sh)                        # builds uniann/bin
```
Two ways to supply the stage environments:

| profile | what it does |
| --- | --- |
| `-profile conda` | Nextflow builds `environment.yml` (EviAnn, gffcompare/gffread, PSAURON) and `sitescore/environment.yml` (PyTorch + mamba-ssm, CUDA) |
| `-profile local_envs` | uses existing envs by path (edit `nextflow.config`); PSAURON and sitescore need `pip install psauron -e sitescore` in the GPU env |

### Usage
```bash
nextflow run ChimAnn -profile local_envs -params-file params.yaml \
    -work-dir /scratch/chimann_work -resume
```
`params.yaml` (see `example/params.yaml`; all paths absolute):

| param | meaning | default |
| --- | --- | --- |
| `genome` | target genome FASTA | required |
| `rnaseq` | EviAnn `-r` list: one sample per line, `R1.fq.gz R2.fq.gz fastq` | none |
| `proteins` | related-species protein FASTA (EviAnn `-p`) | none → EviAnn downloads Swiss-Prot |
| `eviann_args` | extra `eviann.sh` options | `''` |
| `eviann_src` | directory of a newer `eviann.sh` + helper scripts, prepended to PATH (binaries still come from the env). INTEGRATE needs `--untrusted-cds`, which is in the [hibiki-kato/eviann](https://github.com/hibiki-kato/eviann) fork (branch `path-lookup`) but not yet in conda EviAnn 2.0.6 | none |
| `outdir` | results directory | `results` |
| `threads` | CPUs for EviAnn / INTEGRATE | 16 |
| `min_seq_len` | sequences shorter than this skip PSAURON/sitescore/UniAnn (EviAnn still annotates them) | 100000 |
| `psauron_chunk` | nt per psauron call; ~1 GB GPU per Mb | 4500000 |
| `uniann_dir` | UniAnn install (dir with `bin/uniann.sh`) | `uniann/` submodule |
| `uniann_args` | extra `uniann.sh` options | `-n` |
| `sitescore_model` | evaluator plug-in (`sitescore models`) | `convmamba` |
| `sitescore_init` | pretrained `model_dir` to fine-tune from (recommended) | none → train from scratch |
| `sitescore_hparams` | JSON forwarded to the model, e.g. `'{"epochs": 20, "batch_size": 2}'` | model defaults |

Outputs under `outdir/`:

| path | content |
| --- | --- |
| `eviann/<genome>.pseudo_label.gff` | EviAnn evidence-based annotation (first pass) |
| `sitescore/model_dir/` | fine-tuned evaluator (`model.pt`, `train_info.json`) |
| `uniann/<seq>.{plus,minus}.uniann.gff` | UniAnn ab initio predictions per sequence and strand |
| `chimann.gff` | final annotation: EviAnn second pass with UniAnn CDS as low-trust evidence |

Add `-with-report report.html -with-trace trace.txt` for per-task timing and
memory. Rerun with `-resume` after a failure; completed stages are cached.

**GPU**: PSAURON and sitescore share one GPU (`maxForks 1`). PSAURON needs
~1 GB per Mb with `-a`, so sequences are scored in `psauron_chunk` pieces
(default 4.5 Mb) and stitched; on CUDA OOM a task retries once on CPU.
Training with the convmamba plug-in fits in ~6 GB at `batch_size 2`.

### Example
`example/` holds *S. pombe* (genome, RNA-seq, reference GFF) and
*S. octosporus* / *S. cryophilus* proteins; `example/params.spom_soct.yaml`
annotates *S. pombe* with *S. octosporus* proteins.
