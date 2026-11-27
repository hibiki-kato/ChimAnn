"""Exercise actual Bowen helpers on small EviAnn-only chromosome splits."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "prepare_cnn_data.py"
SPEC = importlib.util.spec_from_file_location("prepare_cnn_data", SCRIPT)
prepare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prepare)


class PrepareCnnDataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="chimann-cnn-data-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.fasta = self.directory / "genome.fa"
        self.annotation = self.directory / "eviann.gff3"
        sequence = list("C" * 120)
        sequence[3:23] = "ATGAAAGTCCCCAGCCCTAA"
        reverse_gene = "ATGAAAGTCCCCAGCCCTAA".translate(
            str.maketrans("ACGT", "TGCA"))[::-1]
        sequence[97:117] = reverse_gene
        self.sequence = "".join(sequence)
        self.fasta.write_text("".join(
            f">{chrom}\n{self.sequence}\n" for chrom in ("train", "val", "test")))
        rows = ["##gff-version 3"]
        for chrom in ("train", "val", "test"):
            for strand, intervals, name in (
                    ("+", [(3, 9), (17, 23)], "plus"),
                    ("-", [(97, 103), (111, 117)], "minus")):
                tid = f"{chrom}.{name}"
                rows.append(f"{chrom}\tEviAnn\tmRNA\t{intervals[0][0] + 1}"
                            f"\t{intervals[-1][1]}\t.\t{strand}\t.\tID={tid}")
                for start, end in intervals:
                    for feature, phase in (("exon", "."), ("CDS", "0")):
                        rows.append(f"{chrom}\tEviAnn\t{feature}\t{start + 1}"
                                    f"\t{end}\t.\t{strand}\t{phase}\tParent={tid}")
        self.annotation.write_text("\n".join(rows) + "\n")

    def generate(self, **kwargs):
        options = dict(fasta=self.fasta, annotation=self.annotation,
                       outdir=self.directory / "windows", validation_chrom=["val"],
                       test_chrom=["test"], window_size=60, stride=30)
        options.update(kwargs)
        return prepare.prepare_dataset(**options)

    def test_actual_npz_labels_both_strands_and_pudataset(self):
        from cnn_mamba.train import PUDataset

        manifest = self.generate()
        for split, chrom in (("train", "train"), ("val", "val"), ("test", "test")):
            directory = self.directory / "windows" / split
            files = sorted(directory.glob("*.npz"))
            self.assertEqual(len(files), 2)
            self.assertTrue(all(path.name.startswith(chrom + "_") for path in files))
            for path in files:
                with np.load(path) as data:
                    self.assertEqual(set(data.files), {
                        "sequence", "labels", "start_stop_labels", "chess_labels",
                        "chess_start_stop_labels", "phase_labels", "phase_mask",
                        "reference_phase_labels", "reference_phase_mask",
                        "cds_labels", "reference_cds_labels",
                    })
                    for primary, legacy in (
                            ("labels", "chess_labels"),
                            ("start_stop_labels", "chess_start_stop_labels"),
                            ("phase_labels", "reference_phase_labels"),
                            ("phase_mask", "reference_phase_mask"),
                            ("cds_labels", "reference_cds_labels")):
                        np.testing.assert_array_equal(data[primary], data[legacy])
                    self.assertEqual(data["labels"][9], 1)
                    self.assertEqual(data["labels"][15], 2)
                    self.assertEqual(data["start_stop_labels"][3], 1)
                    self.assertEqual(data["start_stop_labels"][20], 2)
                    positions = np.r_[3:9, 17:23]
                    np.testing.assert_array_equal(data["phase_labels"][positions],
                                                  np.tile([0, 1, 2], 4))
                    np.testing.assert_array_equal(np.flatnonzero(data["phase_mask"]), positions)
                    np.testing.assert_array_equal(np.flatnonzero(data["cds_labels"]), positions)
                    np.testing.assert_array_equal(data["sequence"][3:6], [0, 3, 2])
            dataset = PUDataset(directory, require_phase=True, require_cds=True)
            self.assertEqual(len(dataset), 2)
            self.assertEqual(len(dataset[0]), 13)
            self.assertTrue(all(tuple(array.shape) == (60,) for array in dataset[0]))
        self.assertEqual(manifest["splits"],
                         {"train": ["train"], "val": ["val"], "test": ["test"]})

    def test_manifests_explicitly_identify_all_targets_as_pseudolabels(self):
        self.generate()
        for name in ("metadata.json", "split_manifest.json"):
            metadata = json.loads((self.directory / "windows" / name).read_text())
            self.assertIn("not gold-reference annotation accuracy", metadata["label_notice"])
            self.assertEqual(metadata["sources"]["eviann"], metadata["sources"]["reference"])
            self.assertEqual(metadata["sources"]["eviann"]["path"], str(self.annotation))

    def test_upstream_conflicting_phases_mask_phase_but_retain_cds_union(self):
        with self.annotation.open("a") as handle:
            handle.write("train\tEviAnn\tmRNA\t4\t9\t.\t+\t.\tID=conflict\n"
                         "train\tEviAnn\texon\t4\t9\t.\t+\t.\tParent=conflict\n"
                         "train\tEviAnn\tCDS\t4\t9\t.\t+\t1\tParent=conflict\n")
        self.generate()
        with np.load(self.directory / "windows/train/train_plus_0.npz") as data:
            self.assertFalse(data["phase_mask"][3:9].any())
            self.assertTrue(data["cds_labels"][3:9].all())
            np.testing.assert_array_equal(data["phase_labels"][17:23], [0, 1, 2, 0, 1, 2])

    def test_upstream_noncanonical_splice_is_excluded(self):
        self.fasta.write_text(self.fasta.read_text().replace("GTCCCCAG", "GCCCCCAG"))
        manifest = self.generate()
        self.assertEqual(manifest["rejected_sites"]["splice"], 3)
        with np.load(self.directory / "windows/train/train_plus_0.npz") as data:
            self.assertFalse(data["labels"].any())
            self.assertEqual(data["start_stop_labels"][3], 1)

    def test_chromosome_splits_and_window_parameters_must_be_valid(self):
        for options, message in (
                ({"validation_chrom": ["test"]}, "disjoint"),
                ({"validation_chrom": []}, "nonempty"),
                ({"test_chrom": ["absent"]}, "absent from FASTA"),
                ({"validation_chrom": ["train", "val"]}, "remain for training"),
                ({"window_size": 20}, "stride")):
            with self.subTest(options=options), self.assertRaisesRegex(ValueError, message):
                self.generate(**options)

    def test_empty_window_split_and_nonempty_destination_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "no eligible windows"):
            self.generate(window_size=150, stride=50)
        with self.assertRaisesRegex(ValueError, "nonempty output"):
            self.generate()


if __name__ == "__main__":
    unittest.main()
