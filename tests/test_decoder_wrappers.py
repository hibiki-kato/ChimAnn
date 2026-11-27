"""Oriented decoder contracts; no network, model inference, or GPU required."""

import argparse
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "bin" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ALTANN = load("run_altann")


class AltAnnWrapperTests(unittest.TestCase):
    """run_altann.py decodes the CNN bundle with AltAnn alone (no UniAnn)."""

    FAKE = r"""#!/usr/bin/env python3
import json
from pathlib import Path
import sys
args = sys.argv[1:]
log = Path(__file__).resolve().parent / 'calls.jsonl'
with log.open('a') as handle:
    handle.write(json.dumps(args) + '\n')
output = Path(args[args.index('--output') + 1])
if args[0] == 'prepare':
    output.mkdir()
    (output / 'prepare.json').write_text(json.dumps({'temperature': None}))
    raise SystemExit(0)
strand = '-' if '--reverse' in args else '+'
tag = 'm' if strand == '-' else 'p'
def block(name, rank, cds, locus='L1'):
    attrs = f'ID={tag}{name};kbest_delta={-rank * 10.0};kbest_rank={rank};kbest_locus={tag}{locus}'
    start, end = cds[0][0], cds[-1][1]
    rows = [f'chr\tUniAnn\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{attrs}']
    length = 0
    for a, b in (cds if strand == '+' else cds[::-1]):
        rows.append(f'chr\tUniAnn\tCDS\t{a}\t{b}\t.\t{strand}\t{(-length) % 3}\tParent={tag}{name}')
        length += b - a + 1
    return rows
rows = block('ref', 0, [(4, 12)]) + block('alt', 1, [(4, 6), (10, 15)]) + block('dup', 2, [(4, 12)])
output.write_text('##gff-version 3\n' + '\n'.join(rows) + '\n')
Path(str(output) + '.tsv').write_text('segment\tstrand\treference_gene_index\trank\n'
                                      'chr\t+\t1\t0\nchr\t+\t1\t1\nchr\t+\t1\t2\n')
Path(str(output) + '.json').write_text(json.dumps({'arguments': args}))
"""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.bundle = self.root / "prediction"
        sequence = "CCCATGAAATAACCC"
        for orientation, dna in (("plus", sequence), ("minus", sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])):
            directory = self.bundle / orientation
            directory.mkdir(parents=True)
            (directory / "sequence.fna").write_text(">chr\n" + dna + "\n")
            (directory / "coding.tsv").write_text("".join(
                f"chr\t{i + 1}\t+\t0.7\t0.1\t0.1\t0.1\tframe0\n" for i in range(len(dna))))
            (directory / "sites.tsv").write_text(ALTANN.SITE_HEADER + "\nchr\t5\t+\tdonor\tGT\t0.9\n")
        self.install = self.root / "altann"
        (self.install / "bin").mkdir(parents=True)
        executable = self.install / "bin" / "altann"
        executable.write_text(self.FAKE)
        executable.chmod(0o755)

    def arguments(self, **changes):
        values = dict(input=self.bundle, altann_dir=self.install, output=self.root / "out" / "x.altann.gff3",
                      ab_initio=self.root / "out" / "x.ab_initio.gff", k=3, flank=1000, threads="2",
                      mult=ALTANN.DEFAULT_MULT, metadata_states=True, overlap=0, margin=0)
        values.update(changes)
        return argparse.Namespace(**values)

    def calls(self):
        return [json.loads(line) for line in (self.install / "bin" / "calls.jsonl").read_text().splitlines()]

    def test_prepares_and_decodes_each_orientation_without_uniann(self):
        metadata = ALTANN.run(self.arguments())
        self.assertTrue(metadata["completed"])
        self.assertEqual([job["status"] for job in metadata["jobs"]], ["complete", "complete"])
        calls = self.calls()
        self.assertEqual([call[0] for call in calls], ["prepare", "decode", "prepare", "decode"])
        for index, (prepare, decode) in enumerate((calls[0:2], calls[2:4])):
            orientation = ("plus", "minus")[index]
            self.assertEqual(prepare[prepare.index("--coding-frame") + 1],
                             str(self.bundle.resolve() / orientation / "coding.tsv"))
            self.assertEqual(float(prepare[prepare.index("--mult") + 1]), ALTANN.DEFAULT_MULT)
            for flag in ("--exact-baseline", "--metadata-states", "--short-cds-filter"):
                self.assertIn(flag, decode)
            self.assertNotIn("--gff", decode)
            self.assertNotIn("--log", decode)
            self.assertEqual(decode[decode.index("--k") + 1], "3")
            self.assertEqual(decode[decode.index("--seqid") + 1], "chr")
            self.assertEqual(decode[decode.index("--sequence-length") + 1], "15")
            self.assertEqual("--reverse" in decode, orientation == "minus")

    def test_ab_initio_keeps_the_global_path_as_baseline_and_adds_distinct_alternatives(self):
        args = self.arguments()
        ALTANN.run(args)
        text = args.ab_initio.read_text()
        transcripts = [line.split("\t")[8] for line in text.splitlines() if "\ttranscript\t" in line]
        self.assertEqual(sorted(t.split(";")[0] for t in transcripts), ["ID=malt", "ID=mref", "ID=palt", "ID=pref"])
        for attributes in transcripts:
            self.assertEqual("kbest_delta" in attributes, "alt" in attributes.split(";")[0])
        self.assertNotIn("dup", text)
        self.assertEqual(args.output.read_text().count("\ttranscript\t"), 6)

    def test_zero_k_keeps_only_the_global_path(self):
        args = self.arguments(k=0, metadata_states=False)
        ALTANN.run(args)
        decodes = [call for call in self.calls() if call[0] == "decode"]
        self.assertTrue(all(call[call.index("--k") + 1] == "1" for call in decodes))
        self.assertTrue(all("--metadata-states" not in call for call in decodes))
        self.assertEqual(args.output.read_text().count("\ttranscript\t"), 2)
        self.assertNotIn("alt", args.ab_initio.read_text())
        report = Path(str(args.output) + ".tsv").read_text().splitlines()
        self.assertEqual([row.split("\t")[3] for row in report[1:]], ["0", "0"])

    def test_low_donor_orientation_is_skipped(self):
        (self.bundle / "minus" / "sites.tsv").write_text(ALTANN.SITE_HEADER + "\nchr\t5\t+\tdonor\tGT\t0.1\n")
        metadata = ALTANN.run(self.arguments())
        self.assertEqual([job["status"] for job in metadata["jobs"]], ["complete", "skipped"])
        self.assertEqual([call[0] for call in self.calls()], ["prepare", "decode"])

    def test_rejects_incorrect_minus_sequence_and_missing_altann(self):
        (self.bundle / "minus" / "sequence.fna").write_text(">chr\nCCCATGAAATAACCC\n")
        with self.assertRaisesRegex(ValueError, "reverse complement"):
            ALTANN.run(self.arguments())
        with self.assertRaisesRegex(ValueError, "AltAnn executable not found"):
            ALTANN.run(self.arguments(altann_dir=self.root / "absent", output=self.root / "o2" / "x.gff3"))

    def test_refuses_stale_job_directory(self):
        args = self.arguments()
        jobs = Path(str(args.output) + ".jobs")
        jobs.mkdir(parents=True)
        (jobs / "old").write_text("x")
        with self.assertRaisesRegex(ValueError, "fresh"):
            ALTANN.run(args)

if __name__ == "__main__":
    unittest.main()
