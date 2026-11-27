"""Local installation must preserve upstream source bytes and refresh copies."""

import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[1] / "bin/install_eviann.py"
SPEC = importlib.util.spec_from_file_location("install_eviann", SCRIPT)
installer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(installer)


class InstallEviannTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="chimann-install-eviann-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "upstream"
        self.src = self.source / "src"
        self.src.mkdir(parents=True)
        self.runtime = self.source / "bin"
        self.dependencies = self.root / "annotation-env/bin"
        self.dependencies.mkdir(parents=True)
        (self.source / "Makefile.am").write_text(
            "dist_bin_SCRIPTS = src/eviann.sh \\\n"
            " src/helper.pl src/ufasta.pl gffread TransDecoder-v5.tar.gz\n")
        for name in ("eviann.sh", "helper.pl", "ufasta.pl"):
            path = self.src / name
            path.write_text("#!/bin/sh\nprintf '%s\\n' 'upstream " + name + "'\n")
            path.chmod(0o644)
        (self.source / "gffread").write_text("#!/bin/sh\nexit 0\n")
        (self.source / "gffread").chmod(0o755)
        for name in installer.DEPENDENCIES:
            path = self.dependencies / name
            path.write_text("#!/bin/sh\nexit 0\n")
            path.chmod(0o755)
        self.path = patch.dict(os.environ, {"PATH": str(self.dependencies)})
        self.path.start()
        self.addCleanup(self.path.stop)

    def test_exact_scripts_are_executable_copies_and_tools_preserve_real_location(self):
        manifest = installer.install_runtime(self.source, self.runtime)
        for name in ("eviann.sh", "helper.pl", "ufasta.pl"):
            original, installed = self.src / name, self.runtime / name
            self.assertEqual(original.read_bytes(), installed.read_bytes())
            self.assertFalse(installed.is_symlink())
            self.assertEqual(original.stat().st_mode & 0o777, 0o644)
            self.assertEqual(installed.stat().st_mode & 0o777, 0o755)
        self.assertTrue((self.runtime / "gffread").is_symlink())
        self.assertEqual((self.runtime / "gffread").resolve(), self.source / "gffread")
        for name in installer.DEPENDENCIES:
            self.assertTrue((self.runtime / name).is_symlink())
            self.assertEqual((self.runtime / name).resolve(), self.dependencies / name)
        self.assertFalse((self.runtime / "TransDecoder-v5.tar.gz").exists())
        metadata = json.loads(manifest.read_text())
        self.assertEqual(metadata["source"], str(self.source))
        self.assertEqual(metadata["tools"]["eviann.sh"]["kind"], "copy")

    def test_reinstall_refreshes_upstream_copies_and_preserves_unrelated_files(self):
        installer.install_runtime(self.source, self.runtime)
        unrelated = self.runtime / "user-file"
        unrelated.write_text("preserve me\n")
        (self.src / "eviann.sh").write_text("#!/bin/sh\necho new upstream release\n")
        installer.install_runtime(self.source, self.runtime)
        self.assertEqual((self.runtime / "eviann.sh").read_bytes(),
                         (self.src / "eviann.sh").read_bytes())
        self.assertEqual(unrelated.read_text(), "preserve me\n")
        self.assertEqual((self.src / "eviann.sh").stat().st_mode & 0o777, 0o644)

    def test_compiled_ufasta_is_optional_and_symlinked_when_available(self):
        ufasta = self.dependencies / "ufasta"
        ufasta.write_text("#!/bin/sh\nexit 0\n")
        ufasta.chmod(0o755)
        installer.install_runtime(self.source, self.runtime)
        self.assertEqual((self.runtime / "ufasta").resolve(), ufasta)
        self.assertEqual((self.runtime / "ufasta.pl").read_bytes(),
                         (self.src / "ufasta.pl").read_bytes())

    def test_missing_dependency_fails_before_installing_partial_runtime(self):
        (self.dependencies / "TransDecoder.Predict").unlink()
        with self.assertRaisesRegex(ValueError, "TransDecoder.Predict is missing"):
            installer.install_runtime(self.source, self.runtime)
        self.assertFalse(self.runtime.exists())

    def test_source_paths_and_existing_directories_are_never_replaced(self):
        original = (self.src / "eviann.sh").read_bytes()
        with self.assertRaisesRegex(ValueError, "refusing to replace source"):
            installer.install_runtime(self.source, self.src)
        self.assertEqual((self.src / "eviann.sh").read_bytes(), original)
        (self.runtime / "helper.pl").mkdir(parents=True)
        with self.assertRaisesRegex(ValueError, "refusing to replace source or directory"):
            installer.install_runtime(self.source, self.runtime)
        self.assertTrue((self.runtime / "helper.pl").is_dir())


if __name__ == "__main__":
    unittest.main()
