# ChimAnn stage contracts

The active workflow has one CNN-Mamba model for coding-frame and site predictions:

```text
genome + RNA/protein evidence -> EVIANN -> evidence.gff + completed run directory
genome + evidence.gff -> CNN_MODEL -> checkpoint + optional joint calibration
                       (cnn_checkpoint skips training)
genome -> SPLIT_GENOME -> segment FASTA
segment + checkpoint -> CNN_SCORE -> plus/minus oriented FASTA, coding.tsv, sites.tsv
oriented inputs -> UNIANN -> processed scores + native best-path GFF/traceback
processed scores + best path -> ALTANN -> existing K-best candidates [K > 0]
UniAnn + distinct AltAnn alternatives -> ab_initio.gff
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

## UniAnn processed scores and traceback

Each orientation is converted independently. Local offset `i` at 0-based oriented position `n` selects absolute Viterbi exon frame `(n-i)%3`. Upstream UniAnn currently uses `(n-i-1)%3`; ChimAnn detects that formula and stages columns `[frame1, frame2, frame0]`, producing the verified assignment without changing probabilities or patching UniAnn. Corrected `(n-i)%3` converters receive the original order; unknown formulas fail. Header removal, validation and the recorded frame adapter occur in ChimAnn. AltAnn consumes the resulting processed emissions directly.

`out.ps.txt` contains `position N E0 E1 E2 I nucleotide`: a 0-based position, five log emission columns and the sequence base. The intron column is shared by `I0/I1/I2`; state order remains `N E0 E1 E2 I0 I1 I2`. Native stop masks and transition rules remain UniAnn's responsibility. `out.gt.txt`, `out.ag.txt`, `out.atg.txt` and `out.stop.txt` contain sparse 0-based motif-start positions and transformed site scores. UniAnn's probability multiplier `score_multiplier=e` and donor-based scaling are score transformations, not temperature calibration.

The coding converter currently uses `log(e * P(frame_i) * (1-P(non_CDS)) + 1e-10)`, rounded to two decimals, and `log(e * P(non_CDS) + 1e-10)` for noncoding. Bowen exports joint frame probabilities, so the extra coding-mass factor is an upstream UniAnn weighting, not a plain log joint probability. ChimAnn and AltAnn preserve that weighting; its interpretation is for UniAnn to confirm.

The native best path and DP/BT traceback are saved with these exact processed inputs. AltAnn uses them as the reference for its gene-local K-best runs. No legacy emission source is substituted when joint-frame predictions are missing.

## Coordinates, segments and alternatives

Both decoders run in oriented 0-based coordinates. GFF uses 1-based inclusive coordinates. For a reverse-complement interval `[a,b]` in 1-based inclusive coordinates on length `L`, genomic segment coordinates are `[L-b+1,L-a+1]` with strand `-`; transcript-order CDS phases are preserved.

Chromosomes are decoded whole by default, as one segment spanning the sequence. A nonzero `segment_len` splits them. Segments use IDs `<seqid>__<start>-<end>-<len>`. After decoding, ChimAnn's coordinate adapters restore chromosome coordinates and retain each transcript from the segment owning its midpoint, with the configured internal edge margin. Sequence identifiers, strands and phases must agree across UniAnn and AltAnn.

`altann_kbest: 0` disables alternatives. Otherwise AltAnn retains the existing K paths per state and existing gene-local flank and short-CDS filtering. Its GFF3 uses `kbest_rank`, `kbest_locus` and `kbest_delta`; TSV/JSON sidecars hold scores and provenance. `merge_altann.py` requires every AltAnn reference CDS chain to occur in UniAnn on the same strand with identical phases. It keeps UniAnn, removes duplicate references, and adds distinct complete alternatives only from retained reference loci. Candidate recall is distinct from final annotation precision.

## Final EviAnn integration

`INTEGRATE` copies the first-pass directory without BAMs and calls upstream EviAnn with `-c ab_initio.gff`. Copies are independent: hard links are unsafe because EviAnn truncates files. `-c` clears `merge.success`, so upstream filtering and merging rerun while completed alignments and assemblies are reused. The final annotation must be newly written and nonempty.

Upstream EviAnn has no `--no-snap`, `--splice-model` or `--untrusted-cds` options. Its ordinary external-CDS handling replaces the earlier fork-specific CNN/untrusted-CDS path; in particular, the extra external-CDS splice-filter step is absent. ORF refinement and external, transcript-only and pseudogene rules can change or drop first-pass CDS. `integrate_novel_only: true` retains the optional novel-locus prefilter.

`kbest_delta_fraction`, unset by default, optionally bounds alternative score drops using EviAnn-confirmed alternatives. With no confirmed alternatives, this filter removes all alternatives. It is not neural temperature calibration and does not modify the K-Viterbi search.

The resulting `<genome>.pseudo_label.gff` becomes `chimann.gff`. Report complete end-to-end validation only after this final driver actually runs and its nonempty output is checked.

The October 16 validation used the earlier fork; its metrics do not establish accuracy or complete validation of the current upstream configuration.
