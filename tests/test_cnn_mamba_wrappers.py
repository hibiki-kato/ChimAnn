"""Local orchestration and checkpoint-metadata contracts; no GPU inference."""

import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import Mock, patch


ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


runner = load_script("run_cnn_mamba")


class ModelBundleTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="chimann-cnn-wrapper-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.checkpoint = self.root / "published.pt"
        self.checkpoint.write_bytes(b"existing upstream checkpoint bytes\x00")
        self.calibration = self.root / "published-temperature.json"
        self.calibration.write_text('{"temperature": 1.2288177735743935}\n')
        self.client = runner.CNNMamba(self.root / "external-cnn-mamba")
        self.client.call = Mock()
        self.client.module = Mock()

    def model(self, **options):
        defaults = dict(fasta=self.root / "genome.fa", annotation=self.root / "eviann.gff3",
                        outdir=self.root / "model")
        defaults.update(options)
        return self.client.model(**defaults)

    def simulate_training(self):
        def prepare(arguments, log):
            data = Path(arguments[arguments.index("--outdir") + 1])
            data.mkdir()
            (data / "metadata.json").write_text('{"label_source":"EviAnn pseudo-labels"}\n')

        def module(name, arguments, log):
            output = Path(arguments[arguments.index("--output-dir") + 1])
            if name == "train":
                (output / "best_model.pt").write_bytes(b"trained checkpoint")
            elif name == "calibrate_joint_frame":
                output.mkdir()
                (output / "temperature.json").write_text('{"temperature":1.25}\n')
            else:
                self.fail("unexpected upstream module " + name)

        self.client.call.side_effect = prepare
        self.client.module.side_effect = module

    def test_existing_checkpoint_is_copied_without_training_or_default_calibration(self):
        output = self.model(checkpoint=self.checkpoint)
        self.assertEqual((output / "best_model.pt").read_bytes(), self.checkpoint.read_bytes())
        self.assertFalse((output / "temperature.json").exists())
        metadata = json.loads((output / "model.json").read_text())
        self.assertFalse(metadata["trained"])
        self.assertIsNone(metadata["joint_temperature"])
        self.assertEqual(metadata["checkpoint_sha256"], runner.sha256(self.checkpoint))
        self.assertIn("Supplied checkpoint", metadata["label_source"])
        self.assertIn("not inferred", metadata["label_source"])
        self.assertNotIn("EviAnn pseudo-labels", metadata["label_source"])
        self.client.call.assert_not_called()
        self.client.module.assert_not_called()

    def test_existing_calibration_is_copied_verbatim_without_site_temperature(self):
        output = self.model(checkpoint=self.checkpoint, calibration_json=self.calibration)
        self.assertEqual((output / "temperature.json").read_bytes(), self.calibration.read_bytes())
        metadata = json.loads((output / "model.json").read_text())
        self.assertEqual(metadata["joint_temperature"], "temperature.json")
        self.assertIn("site heads unchanged", metadata["temperature_scope"])
        self.client.module.assert_not_called()

    def test_training_passes_disjoint_split_lists_and_required_joint_head(self):
        self.simulate_training()
        output = self.model(validation_chrom=["chrVal1", "chrVal2"],
                            test_chrom=["chrTest"], window_size=120, stride=60,
                            train_args="--epochs 1 --lr 0.0001")
        prepare_args = self.client.call.call_args.args[0]
        self.assertEqual([prepare_args[i + 1] for i, value in enumerate(prepare_args)
                          if value == "--validation-chrom"], ["chrVal1", "chrVal2"])
        self.assertEqual(prepare_args[prepare_args.index("--test-chrom") + 1], "chrTest")
        name, train_args, log = self.client.module.call_args.args
        self.assertEqual(name, "train")
        self.assertEqual(train_args[train_args.index("--mode") + 1], "baseline")
        self.assertGreater(float(train_args[train_args.index("--joint-frame-aux-weight") + 1]), 0)
        self.assertNotIn("--all-splits-as-train", train_args)
        self.assertNotIn("--include-test-in-train", train_args)
        self.assertEqual((output / "training-data.json").read_bytes(),
                         (self.root / "training_data/metadata.json").read_bytes())
        self.assertIn("not gold accuracy", json.loads((output / "model.json").read_text())["label_source"])

    def test_reserved_data_and_head_overrides_cannot_launch_training(self):
        self.simulate_training()
        for index, arguments in enumerate(("--data-root /other", "--output-dir=/other",
                                           "--joint-frame-aux-weight=0",
                                           "--all-splits-as-train", "--include-test-in-train",
                                           "--data-ro /other", "--all-splits-as-tr",
                                           "--joint-frame-aux-weigh=0")):
            with self.subTest(arguments=arguments):
                self.client.call.reset_mock()
                # Preparation is irrelevant to this rejected training request.
                self.client.call.side_effect = None
                with self.assertRaisesRegex(ValueError, "preserve held-out splits"):
                    self.model(outdir=self.root / f"invalid-{index}",
                               validation_chrom=["v"], test_chrom=["t"], train_args=arguments)
                self.client.module.assert_not_called()

    def test_temperature_fit_uses_only_validation_and_preserves_upstream_result(self):
        self.simulate_training()
        output = self.model(validation_chrom=["v"], test_chrom=["t"], fit_temperature=True)
        fit = [call for call in self.client.module.call_args_list
               if call.args[0] == "calibrate_joint_frame"]
        self.assertEqual(len(fit), 1)
        arguments = fit[0].args[1]
        self.assertEqual(arguments[arguments.index("--split") + 1], "val")
        self.assertEqual((output / "temperature.json").read_bytes(),
                         (self.root / "calibration/temperature.json").read_bytes())

    def test_reused_model_cannot_fit_temperature_without_training_data(self):
        with self.assertRaisesRegex(ValueError, "requires training data"):
            self.model(checkpoint=self.checkpoint, fit_temperature=True)
        self.client.call.assert_not_called()
        self.client.module.assert_not_called()

    def test_new_training_cannot_reuse_an_unrelated_checkpoint_temperature(self):
        with self.assertRaisesRegex(ValueError, "calibration-json requires --checkpoint"):
            self.model(calibration_json=self.calibration,
                       validation_chrom=["v"], test_chrom=["t"])
        self.client.call.assert_not_called()
        self.client.module.assert_not_called()
        self.assertFalse((self.root / "model").exists())

    def test_oriented_prediction_calibrates_joint_head_without_recalibrating_site_heads(self):
        fasta = self.root / "genome.fa"
        sequence = "ATGCRYKMBDHVN"
        fasta.write_text(">chr\n" + sequence + "\n")

        def upstream_module(name, arguments, log):
            if name == "score_joint_frame_chromosome":
                output = Path(arguments[arguments.index("--output-dir") + 1])
                output.mkdir()
                (output / "metadata.json").write_text(json.dumps({
                    "temperature": 1.2288177735743935,
                    "temperature_scope": "joint four-class logits only",
                }))

        self.client.module.side_effect = upstream_module
        with patch.object(runner.subprocess, "check_output", side_effect=["abc123\n", ""]):
            result = self.client.predict(fasta, self.checkpoint, self.root / "predictions",
                                         calibration_json=self.calibration)
        manifest = json.loads(result.read_text())
        self.assertEqual(len(manifest["jobs"]), 2)
        for orientation, expected in (("plus", sequence), ("minus", "NBDHVKMRYGCAT")):
            actual = (self.root / "predictions" / orientation / "sequence.fna").read_text()
            self.assertEqual(actual, ">chr\n" + expected + "\n")
        joint = [call for call in self.client.module.call_args_list
                 if call.args[0] == "score_joint_frame_chromosome"]
        self.assertEqual(len(joint), 2)
        self.assertTrue(all("--calibration-json" in call.args[1] for call in joint))
        site = [call.args[0] for call in self.client.call.call_args_list
                if str(call.args[0][0]).endswith("score_cnn_sites.py")]
        self.assertEqual(len(site), 2)
        self.assertTrue(all("--calibration-json" not in arguments
                            and "--temperature" not in arguments for arguments in site))
        self.assertEqual([job["original_strand"] for job in manifest["jobs"]], ["+", "-"])


class SiteMetadataAdapterTests(unittest.TestCase):
    def run_adapter(self, checkpoint, fail=False):
        model = checkpoint["model"]
        original_model = copy.deepcopy(model)
        load = Mock(return_value=checkpoint)
        writes = []
        write = Mock(side_effect=lambda path, payload: writes.append((path, copy.deepcopy(payload))))
        scorer = SimpleNamespace(torch=SimpleNamespace(load=load), atomic_json=write)

        def main():
            saved = scorer.torch.load("checkpoint.pt", map_location="cpu")
            self.assertIs(saved["model"], model)
            self.assertEqual(saved["model"], original_model)
            if fail:
                raise RuntimeError("simulated scoring failure")
            scorer.atomic_json("metadata.json", {"checkpoint_step": saved["optimizer_step"]})
            scorer.atomic_json(".complete.json", {"schema": "complete"})

        scorer.main = main
        package = ModuleType("cnn_mamba")
        package.score_chromosome = scorer
        with patch.dict(sys.modules, {"cnn_mamba": package}):
            adapter = load_script("score_cnn_sites")
        try:
            if fail:
                with self.assertRaisesRegex(RuntimeError, "simulated scoring failure"):
                    adapter.main()
            else:
                adapter.main()
        finally:
            self.assertIs(scorer.torch.load, load)
            self.assertIs(scorer.atomic_json, write)
        self.assertEqual(model, original_model)
        return writes

    def test_public_checkpoint_records_unknown_step_and_known_epoch(self):
        writes = self.run_adapter({"model": {"weight": [1, 2, 3]}, "epoch": 4})
        self.assertEqual(writes[0], ("metadata.json", {"checkpoint_step": None, "checkpoint_epoch": 4}))
        self.assertEqual(writes[1], (".complete.json", {"schema": "complete"}))

    def test_training_checkpoint_preserves_existing_step(self):
        writes = self.run_adapter({"model": {"weight": [4]}, "epoch": 2, "optimizer_step": 123})
        self.assertEqual(writes[0][1], {"checkpoint_step": 123, "checkpoint_epoch": 2})

    def test_unknown_epoch_is_null_and_exception_restores_upstream_functions(self):
        writes = self.run_adapter({"model": {"weight": [7]}})
        self.assertIsNone(writes[0][1]["checkpoint_epoch"])
        self.run_adapter({"model": {"weight": [7]}, "epoch": 4}, fail=True)


if __name__ == "__main__":
    unittest.main()
