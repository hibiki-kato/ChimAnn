# ChimAnn stage contracts

The active workflow has one CNN-Mamba model for coding-frame and site predictions:

```text
genome + RNA/protein evidence -> EVIANN -> evidence.gff + completed run directory
genome + evidence.gff -> CNN_MODEL -> checkpoint + optional joint calibration
                       (cnn_checkpoint skips training)
genome -> SPLIT_GENOME -> segment FASTA
segment + checkpoint -> CNN_SCORE -> plus/minus oriented FASTA, coding.tsv, sites.tsv
oriented inputs -> ALTANN: altann prepare -> decoder inputs (no further temperature)
                          altann decode --exact-baseline -> global best path + K-best candidates
global path (baseline) + distinct AltAnn alternatives -> ab_initio.gff
completed EviAnn run + ab_initio.gff -> INTEGRATE -> chimann.gff
```

## Evidence and training

`EVIANN` runs upstream `eviann.sh -t N -g genome.fna -r rnaseq.list -p protein.faa` with its normal Markov splice filter and annotation rules. RNA manifest entries must contain absolute paths. `eviann_args` and `integrate_args` default to empty. The installed runtime is created by `python3 bin/install_eviann.py --source eviann --outdir eviann/bin`, with required annotation tools available on `PATH`; the upstream source remains unmodified.

`eviann_resume_dir` seeds a copy of completed alignments and transcript assemblies. By default the first-pass merge reruns. With `eviann_reuse_annotation: true`, the nonempty completed `<genome>.pseudo_label.gff` is reused directly; the final integration pass still runs. Reused annotations retain their original provenance; they are not silently regenerated as upstream output. Preserve the genome/protein basenames and original evidence paths so EviAnn resolves its cached files.

`CNN_MODEL` uses only that EviAnn annotation as labels. Explicit `cnn_validation_chroms` and `cnn_test_chroms` are required for training, must be nonempty and disjoint, and must leave training chromosomes. Splits contain whole FASTA records, preventing overlap between training and held-out windows. Model windows default to 10,000 bases with stride 5,000. Both strand orientations, motif labels and local coding-phase targets use CNN-Mamba's preparation rules.

The upstream dataset interface contains `reference_*`/`chess_*` arrays; the ChimAnn adapter fills these from the same EviAnn labels. They are compatibility fields, not independently verified reference truth. Validation, test metrics and optional temperature fitting therefore measure held-out EviAnn pseudo-label consistency.

`cnn_checkpoint` skips training and must contain the joint four-class head as well as the site heads. CNN-Mamba runs with `cnn_python` from the clean `cnn_mamba_dir` checkout, default the `cnn-mamba` submodule. Checkpoints/calibration can be stored elsewhere, including an experimental checkout, without importing that checkout's source. An inference checkpoint may omit `optimizer_step`; the adapter records missing progress without changing model outputs.

## Probabilities and calibration

The joint head's class order is `[non_CDS, frame0, frame1, frame2]`. Frame labels are local codon offsets: offset 0 is the first coding base, followed by 1 and 2. `CNN_SCORE` exports probabilities after `softmax(logits / T)`, with `T=1` unless `cnn_calibration_json` supplies a temperature or `cnn_fit_temperature` requests fitting during model preparation. Calibration belongs to the chosen checkpoint and is applied once to the joint head only. Site heads retain their model probabilities.

`coding.tsv` is tab-separated with the CNN-Mamba header:

```text
chrom	position_1based	strand	frame0_probability	frame1_probability	frame2_probability	noncoding_probability	predicted_state
```

It contains exactly one row per oriented FASTA base, in contiguous 1-based order. The export reorders the internal noncoding-first class probabilities. Chromosome IDs match the oriented FASTA, and both independently oriented exports use `+` labels. ChimAnn's `run_uniann.py` strips the exact header and validates row count, coordinates, strand and finite probabilities before invoking UniAnn. The pinned upstream converter does not handle the header itself.

`sites.tsv` is tab-separated with this header:

```text
chrom	pos	strand	type	motif	prob
```

Site types are donor/acceptor/start/stop, with GT/AG/ATG/TAA,TAG,TGA motifs. In each independently scored oriented FASTA, positions are 1-based motif starts in that sequence and rows passed to the forward UniAnn wrapper use `+`. Minus genomic predictions arise from scoring a reverse-complement FASTA and restoring coordinates later, not from feeding a minus-labelled original-coordinate table into the forward wrapper.

## AltAnn decoder inputs

`altann prepare` converts each orientation independently. Local class `i` at 0-based oriented position `n` scores absolute Viterbi exon frame `(n-i)%3`. No temperature is applied after the model's own calibration; UniAnn 08eaf15's converter applied `sigmoid(logit(p)/3)` again, which is why ChimAnn no longer uses it. The emissions match UniAnn f126cc8, the last upstream converter without a temperature.

`out.ps.txt` contains `position N E0 E1 E2 I nucleotide`: a 0-based position, five log emission columns and the sequence base. The intron column is shared by `I0/I1/I2`; state order is `N E0 E1 E2 I0 I1 I2`. `out.gt.txt`, `out.ag.txt`, `out.atg.txt` and `out.stop.txt` contain sparse 0-based motif-start positions and transformed site scores. The probability multiplier `score_multiplier=e` and donor-based scaling are score transformations, not temperature calibration. An orientation without a donor whose rescaled score is positive is skipped.

The coding emission is `log(e * P(frame_i) * (1-P(non_CDS)) + 1e-10)`, rounded to two decimals, and `log(e * P(non_CDS) + 1e-10)` for noncoding. The extra coding-mass factor is the upstream UniAnn weighting, kept as is.

The global best path comes from exact K=1 metadata-state search over the whole sequence (`--exact-baseline`); local K-best runs use it as their reference. No legacy emission source is substituted when joint-frame predictions are missing.

## Coordinates, segments and alternatives

Both decoders run in oriented 0-based coordinates. GFF uses 1-based inclusive coordinates. For a reverse-complement interval `[a,b]` in 1-based inclusive coordinates on length `L`, genomic segment coordinates are `[L-b+1,L-a+1]` with strand `-`; transcript-order CDS phases are preserved.

Chromosomes are decoded whole by default, as one segment spanning the sequence. A nonzero `segment_len` splits them. Segments use IDs `<seqid>__<start>-<end>-<len>`. After decoding, ChimAnn's coordinate adapters restore chromosome coordinates and retain each transcript from the segment owning its midpoint, with the configured internal edge margin. Sequence identifiers, strands and phases are checked when the strand jobs are merged.

`altann_kbest: 0` keeps only the global path. Otherwise AltAnn retains K paths per state (exact metadata states by default) with gene-local flank and short-CDS filtering. Its GFF3 uses `kbest_rank`, `kbest_locus` and `kbest_delta`; TSV/JSON sidecars hold scores and provenance. `run_altann.py` writes the rank-0 models without `kbest_*` attributes as baselines, and `merge_altann.py` adds distinct complete alternatives only from retained reference loci. Candidate recall is distinct from final annotation precision.

## Final EviAnn integration

`INTEGRATE` copies the first-pass directory without BAMs and calls upstream EviAnn with `-c ab_initio.gff`. Copies are independent: hard links are unsafe because EviAnn truncates files. `-c` clears `merge.success`, so upstream filtering and merging rerun while completed alignments and assemblies are reused. The final annotation must be newly written and nonempty.

Upstream EviAnn has no `--no-snap`, `--splice-model` or `--untrusted-cds` options. Its ordinary external-CDS handling replaces the earlier fork-specific CNN/untrusted-CDS path; in particular, the extra external-CDS splice-filter step is absent. ORF refinement and external, transcript-only and pseudogene rules can change or drop first-pass CDS. `integrate_novel_only: true` retains the optional novel-locus prefilter.

`kbest_delta_fraction`, unset by default, optionally bounds alternative score drops using EviAnn-confirmed alternatives. With no confirmed alternatives, this filter removes all alternatives. It is not neural temperature calibration and does not modify the K-Viterbi search.

The resulting `<genome>.pseudo_label.gff` becomes `chimann.gff`. Report complete end-to-end validation only after this final driver actually runs and its nonempty output is checked.

The October 16 validation used the earlier fork; its metrics do not establish accuracy or complete validation of the current upstream configuration.
