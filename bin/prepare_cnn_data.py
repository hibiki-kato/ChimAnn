#!/usr/bin/env python3
"""Prepare chromosome-held-out CNN-Mamba windows from EviAnn pseudo-labels."""

import argparse
from pathlib import Path

import pysam
from cnn_mamba import prepare_droso as bowen


PSEUDO_LABEL_NOTICE = (
    "All labels, including chess_* and reference_* arrays, come from the supplied "
    "EviAnn annotation. Validation/test measure held-out pseudo-label agreement; "
    "they are not gold-reference annotation accuracy."
)


def prepare_dataset(fasta, annotation, outdir, validation_chrom, test_chrom,
                    window_size=10000, stride=5000):
    fasta, annotation, outdir = map(Path, (fasta, annotation, outdir))
    fasta, annotation = fasta.resolve(), annotation.resolve()
    if not fasta.is_file() or not annotation.is_file():
        raise ValueError("--fasta and --annotation must name existing files")
    if not 0 < stride <= window_size:
        raise ValueError("require 0 < --stride <= --window-size")
    with pysam.FastaFile(str(fasta)) as genome:
        chromosomes = list(genome.references)
    validation, test = set(validation_chrom), set(test_chrom)
    if not validation or not test:
        raise ValueError("validation and test chromosome lists must be nonempty")
    if validation & test:
        raise ValueError("validation and test chromosomes must be disjoint")
    unknown = (validation | test) - set(chromosomes)
    if unknown:
        raise ValueError("chromosomes absent from FASTA: " + ", ".join(sorted(unknown)))
    splits = {
        "train": [name for name in chromosomes if name not in validation | test],
        "val": [name for name in chromosomes if name in validation],
        "test": [name for name in chromosomes if name in test],
    }
    if not splits["train"]:
        raise ValueError("at least one chromosome must remain for training")
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("refusing nonempty output directory: " + str(outdir))

    # Reuse Bowen's exact motif, CDS phase/conflict, strand and window rules.
    # His generator constructs dense tracks for only the current chromosome.
    bowen.FASTA = fasta
    bowen.CHROMS = {name: name for name in chromosomes}
    bowen.SPLITS = splits
    bowen.EVIANN_GFF = annotation
    bowen.REFERENCE_GTF = annotation
    meta, exons, cds, _ = bowen.parse_eviann(annotation)
    regions = bowen.coding_regions(meta, exons, cds)
    sites, rejected = bowen.extract_sites(meta, exons, cds)
    indexed = bowen.index_sites(sites)
    phases = bowen.parse_cds_phase_records(annotation)
    manifest = bowen.generate_cell(
        outdir, {"id": "eviann_pseudolabels", "window_bp": window_size,
                 "stride_bp": stride},
        regions, indexed, indexed, phases, phases,
    )
    for split, summary in manifest["summary"].items():
        if not summary["windows"]:
            raise ValueError(
                f"{split} has no eligible windows; check annotations, chromosome "
                "lengths and --window-size; use a fresh output directory to retry"
            )
    manifest.update(
        schema="chimann_eviann_cnn_windows_v1",
        label_notice=PSEUDO_LABEL_NOTICE,
        label_source="EviAnn pseudo-labels only",
        split_policy="disjoint whole FASTA records; no region/window splitting",
        rejected_sites=rejected,
        legacy_reference_arrays="copies of the same EviAnn pseudo-label targets",
    )
    # Both manifests carry the warning, including the upstream-compatible one.
    bowen.atomic_json(outdir / "split_manifest.json", manifest)
    bowen.atomic_json(outdir / "metadata.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--validation-chrom", action="append", required=True)
    parser.add_argument("--test-chrom", action="append", required=True)
    parser.add_argument("--window-size", type=int, default=10000)
    parser.add_argument("--stride", type=int, default=5000)
    args = parser.parse_args()
    try:
        prepare_dataset(args.fasta, args.annotation, args.outdir,
                        args.validation_chrom, args.test_chrom,
                        args.window_size, args.stride)
    except ValueError as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
