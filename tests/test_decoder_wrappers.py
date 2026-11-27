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


UNIANN, ALTANN = load("run_uniann"), load("run_altann")


class DecoderWrapperTests(unittest.TestCase):
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
            (directory / "coding.tsv").write_text(UNIANN.CODING_HEADER + "\n" + "".join(
                f"chr\t{i + 1}\t+\t0.7\t0.1\t0.1\t0.1\tframe0\n" for i in range(len(dna))))
            (directory / "sites.tsv").write_text(UNIANN.SITE_HEADER + "\nchr\t5\t+\tdonor\tGT\t0.9\n")
        self.install = self.root / "uniann"
        (self.install / "scripts").mkdir(parents=True)
        (self.install / "bin").mkdir()
        self.executable(self.install / "scripts" / "uniann.sh", '''#!/bin/bash
MYPATH=$(dirname "$0")
while [[ $# -gt 0 ]]; do
 case "$1" in
  -n) NO_TRACE=1; shift; continue;;
  -f) FASTA="$2";; -c) CODING="$2";; -s) SITES="$2";; -m) MULT="$2";;
 esac
 shift 2
done
"$MYPATH/convert_coding_frame.pl" "$FASTA" < "$CODING" > out.ps.txt
touch out.gt.txt out.ag.txt out.atg.txt out.stop.txt
"$MYPATH/uniann" "$FASTA" > "$FASTA.uniann.gff" 2> out.err
echo success
''')
        self.executable(self.install / "scripts" / "convert_coding_frame.pl", '''#!/usr/bin/env perl
open(my $fasta, '<', $ARGV[0]) or die $!;
my $sequence = join('', grep { !/^>/ } <$fasta>);
$sequence =~ s/\\s+//g;
my $n = 0;
while (<STDIN>) {
 for (my $i = 0; $i < 3; $i++) { $coding_frame[($n-$i)%3] = 0; }
 print "$n\\t0\\t0\\t0\\t0\\t0\\t", substr($sequence, $n, 1), "\\n";
 $n++;
}
''')
        self.native = self.install / "bin" / "uniann"
        self.executable(self.native, '''#!/usr/bin/env python3
from pathlib import Path
import sys
sequence = ''.join(Path(sys.argv[1]).read_text().splitlines()[1:])
for n in range(len(sequence)):
 print(str(n) + '\\tdp\\t0\\t-1\\t-2\\t-3\\t-4\\t-5\\t-6', file=sys.stderr)
 print(str(n) + '\\tbt\\t' + '\\t'.join([str(-1 if n == 0 else 0)] * 7), file=sys.stderr)
print('##gff-version 3')
print('chr\\tUniAnn\\ttranscript\\t4\\t12\\t.\\t+\\t.\\tID=t1')
print('chr\\tUniAnn\\texon\\t4\\t12\\t.\\t+\\t.\\tParent=t1')
print('chr\\tUniAnn\\tCDS\\t4\\t12\\t.\\t+\\t0\\tParent=t1')
''')
        self.executable(self.install / "bin" / "gffread", "#!/bin/bash\ncat\n")
        self.args = argparse.Namespace(input=self.bundle, uniann_dir=self.install,
                                       outdir=self.root / "jobs", output=self.root / "best.gff3",
                                       mult=UNIANN.DEFAULT_MULT, overlap=0, margin=0)

    def executable(self, path, source):
        path.write_text(source)
        path.chmod(0o755)

    def test_bundle_keeps_exact_scores_trace_and_local_reference(self):
        metadata = UNIANN.run(self.args)
        self.assertTrue(metadata["completed"])
        self.assertEqual(metadata["coding_frame_adapter"]["rotation"], "none")
        self.assertEqual(metadata["coding_frame_adapter"]["staged_frame_columns"], [0, 1, 2])
        self.assertEqual([job["status"] for job in metadata["jobs"]], ["complete", "complete"])
        for job in metadata["jobs"]:
            self.assertEqual(job["dp_bt_rows"], 30)
            self.assertEqual(job["emission_rows"], 15)
            local = self.args.outdir / job["reference_gff"]
            self.assertIn("\t+\t0\tParent=t1", local.read_text())
            supplied = (self.bundle / job["name"] / "coding.tsv").read_text().splitlines(keepends=True)
            self.assertEqual((self.args.outdir / job["name"] / "coding.tsv").read_text(), "".join(supplied[1:]))
        rows = [line.split("\t") for line in self.args.output.read_text().splitlines() if not line.startswith("#")]
        self.assertEqual({row[6] for row in rows}, {"+", "-"})
        self.assertEqual(len({row[8] for row in rows if row[2] == "transcript"}), 2)

    def test_legacy_converter_rotation_is_recorded_and_source_bundle_is_untouched(self):
        converter = self.install / "scripts" / "convert_coding_frame.pl"
        converter.write_text(converter.read_text().replace("($n-$i)%3", "($n-$i-1)%3"))
        before = {name: (self.bundle / name / "coding.tsv").read_bytes() for name in ("plus", "minus")}
        metadata = UNIANN.run(self.args)
        adapter = metadata["coding_frame_adapter"]
        self.assertEqual(adapter["converter_formula"], "($n-$i-1)%3")
        self.assertEqual(adapter["staged_frame_columns"], [1, 2, 0])
        self.assertEqual(adapter["rotation"], "left_one")
        for name in ("plus", "minus"):
            self.assertEqual((self.bundle / name / "coding.tsv").read_bytes(), before[name])
            row = (self.args.outdir / name / "coding.tsv").read_text().splitlines()[0].split("\t")
            self.assertEqual(row[3:7], ["0.1", "0.1", "0.7", "0.1"])

    def test_masked_native_failure_is_rejected_by_missing_trace(self):
        self.executable(self.native, "#!/bin/bash\nexit 7\n")
        with self.assertRaisesRegex(ValueError, "DP/BT dump"):
            UNIANN.run(self.args)
        self.assertFalse(self.args.output.exists())

    def test_masked_gff_pipeline_failure_cannot_publish_empty_annotation(self):
        wrapper = self.install / "scripts" / "uniann.sh"
        wrapper.write_text(wrapper.read_text().replace(
            '"$MYPATH/uniann" "$FASTA" > "$FASTA.uniann.gff" 2> out.err',
            '"$MYPATH/uniann" "$FASTA" 2> out.err | "$MYPATH/gffread" --tlf | cat > "$FASTA.gff.tmp" && mv "$FASTA.gff.tmp" "$FASTA.uniann.gff"'))
        self.executable(self.install / "bin" / "gffread", "#!/bin/bash\ncat >/dev/null\nexit 19\n")
        with self.assertRaises(FileNotFoundError):
            UNIANN.run(self.args)
        self.assertFalse(self.args.output.exists())

    def test_low_donor_strand_is_explicitly_skipped(self):
        (self.bundle / "minus" / "sites.tsv").write_text("chr\t5\t+\tdonor\tGT\t0.001\n")
        metadata = UNIANN.run(self.args)
        self.assertEqual(metadata["jobs"][1]["status"], "skipped")
        self.assertEqual(metadata["jobs"][1]["reason"], "no_donor_with_positive_rescaled_score")
        self.assertEqual((self.args.outdir / "minus" / "uniann.complete.gff3").read_text(), "##gff-version 3\n")

    def test_rejects_malformed_or_unoriented_tables(self):
        coding = self.bundle / "plus" / "coding.tsv"
        coding.write_text(coding.read_text().replace("\t+\t", "\t-\t"))
        with self.assertRaisesRegex(ValueError, "oriented FASTA"):
            UNIANN.run(self.args)

    def test_refuses_stale_job_directory(self):
        UNIANN.run(self.args)
        with self.assertRaisesRegex(ValueError, "must be fresh"):
            UNIANN.run(self.args)

    def test_rejects_incorrect_minus_sequence(self):
        (self.bundle / "minus" / "sequence.fna").write_text(">chr\n" + "C" * 15 + "\n")
        with self.assertRaisesRegex(ValueError, "reverse complement"):
            UNIANN.run(self.args)

    def test_gene_ancestors_survive_mapping_with_unique_parent_links(self):
        rows = [(["chr", "UniAnn", kind, "4", "12", ".", "+", phase, "."], attrs)
                for kind, phase, attrs in (("gene", ".", {"ID": "g1"}),
                                           ("transcript", ".", {"ID": "t1", "Parent": "g1"}),
                                           ("exon", ".", {"Parent": "t1"}),
                                           ("CDS", "0", {"Parent": "t1"}))]
        mapped = UNIANN.map_gff(rows, "chr", "minus", UNIANN.coordinates("chr", 15), 0, 0)
        self.assertEqual(len(mapped), 4)
        self.assertEqual(mapped[1][1]["Parent"], mapped[0][1]["ID"])
        self.assertEqual(mapped[2][1]["Parent"], mapped[1][1]["ID"])

    def test_terminal_partial_is_explicit_and_nonterminal_model_is_preserved(self):
        def model(identifier, start, end):
            return [(["chr", "UniAnn", kind, str(start), str(end), ".", "+", phase, "."], attributes)
                    for kind, phase, attributes in (("transcript", ".", {"ID": identifier}),
                                                     ("exon", ".", {"Parent": identifier}),
                                                     ("CDS", "0", {"Parent": identifier}))]
        # Deliberately nonterminal invalid coding sequence is retained, not filtered.
        rows = model("earlier", 1, 3) + model("terminal", 7, 14)
        complete, exclusions = UNIANN.complete_reference(rows, "C" * 15, [-5, 10, -5, -5, -5, -5, -5])
        self.assertEqual([item["transcript"] for item in exclusions], ["terminal"])
        self.assertEqual([attrs["ID"] for row, attrs in complete if row[2] == "transcript"], ["earlier"])
        unchanged, exclusions = UNIANN.complete_reference(rows, "C" * 15, [10, 10, -5, -5, -5, -5, -5])
        self.assertEqual(unchanged, rows)
        self.assertEqual(exclusions, [])

    def test_reverse_segment_mapping_ownership_and_phase_match_existing_convention(self):
        coords = UNIANN.coordinates("chr__100-200-300", 100)
        rows = [(["chr__100-200-300", "UniAnn", kind, "21", "60", ".", "+", phase, "."], attrs)
                for kind, phase, attrs in (("transcript", ".", {"ID": "tx"}),
                                           ("exon", ".", {"Parent": "tx"}),
                                           ("CDS", "2", {"Parent": "tx"}))]
        mapped = UNIANN.map_gff(rows, "chr__100-200-300", "minus", coords, 20, 5)
        self.assertEqual(mapped[0][0][3:5], ["141", "180"])
        self.assertEqual(mapped[2][0][6:8], ["-", "2"])
        self.assertEqual(UNIANN.map_gff(rows, "chr__100-200-300", "minus", coords, 100, 5), [])

    def install_altann_decoder(self):
        install = self.root / "altann"
        (install / "build").mkdir(parents=True)
        self.executable(install / "build" / "altann-core", '''#!/usr/bin/env python3
from pathlib import Path
import sys
assert sys.argv[-1] == '--uniann-gff'
sequence = ''.join(Path(sys.argv[1]).read_text().splitlines()[1:])
print(str(len(sequence) - 1) + '\\tdp\\t0\\t-1\\t-2\\t-3\\t-4\\t-5\\t-6', file=sys.stderr)
print('##gff-version 3')
print('chr\\tUniAnn\\ttranscript\\t4\\t12\\t.\\t+\\t.\\tID=t1')
print('chr\\tUniAnn\\texon\\t4\\t12\\t.\\t+\\t.\\tParent=t1')
print('chr\\tUniAnn\\tCDS\\t4\\t12\\t.\\t+\\t0\\tParent=t1')
''')
        return install

    def test_native_decoder_without_trace_needs_no_log_and_altann_recomputes(self):
        self.args.no_trace = True
        metadata = UNIANN.run(self.args)
        self.assertEqual(metadata["decoder"], "uniann")
        self.assertFalse(metadata["trace"])
        self.assertTrue(all("-n" in job["command"] and "dp_bt_rows" not in job for job in metadata["jobs"]))
        (self.args.outdir / "plus" / "out.err").unlink()
        install = self.install_altann_decoder()
        (install / "bin").mkdir()
        self.executable(install / "bin" / "altann", '''#!/usr/bin/env python3
import json
from pathlib import Path
import sys
output = Path(sys.argv[sys.argv.index('--output') + 1])
output.write_text('##gff-version 3\\n')
Path(str(output) + '.tsv').write_text('segment\\tstrand\\n')
Path(str(output) + '.json').write_text(json.dumps({'arguments': sys.argv[1:]}))
''')
        args = argparse.Namespace(input=self.args.outdir, output=self.root / "alternatives.gff3",
                                  altann_dir=install, k=3, flank=1000, threads="2", overlap=0, margin=0)
        for job in ALTANN.run(args)["jobs"]:
            self.assertIn("--rerun-viterbi", job["result"]["arguments"])
            self.assertNotIn("--log", job["result"]["arguments"])

    def test_terminal_partial_is_still_excluded_without_final_scores(self):
        sequence = "ATGAAACCC"
        rows = [(["chr", "UniAnn", "transcript", "1", "9", ".", "+", ".", "ID=t"], {"ID": "t"}),
                (["chr", "UniAnn", "CDS", "1", "9", ".", "+", "0", "Parent=t"], {"Parent": "t"})]
        kept, exclusions = UNIANN.complete_reference(rows, sequence, None)
        self.assertEqual(kept, [])
        self.assertEqual(len(exclusions), 1)
        complete, none = UNIANN.complete_reference(rows, "ATGAAATAA", None)
        self.assertEqual(complete, rows)
        self.assertEqual(none, [])

    def test_altann_decoder_replaces_the_native_binary_and_its_trace(self):
        self.native.unlink()
        self.args.altann_dir = self.install_altann_decoder()
        metadata = UNIANN.run(self.args)
        self.assertEqual(metadata["decoder"], "altann")
        self.assertEqual([job["status"] for job in metadata["jobs"]], ["complete", "complete"])
        self.assertTrue(all("dp_bt_rows" not in job for job in metadata["jobs"]))
        self.assertEqual(len((self.args.outdir / "plus" / "out.err").read_text().splitlines()), 1)
        self.assertIn("ID=", self.args.output.read_text())

    def test_altann_decoder_without_final_row_is_rejected(self):
        self.native.unlink()
        install = self.install_altann_decoder()
        self.executable(install / "build" / "altann-core", "#!/bin/bash\necho '##gff-version 3'\n")
        self.args.altann_dir = install
        with self.assertRaisesRegex(ValueError, "final dp row"):
            UNIANN.run(self.args)

    def test_missing_altann_build_is_reported(self):
        self.args.altann_dir = self.root / "absent"
        with self.assertRaisesRegex(ValueError, "not built"):
            UNIANN.run(self.args)

    def test_altann_stage_recomputes_the_path_without_a_native_trace(self):
        self.native.unlink()
        install = self.install_altann_decoder()
        self.args.altann_dir = install
        UNIANN.run(self.args)
        (install / "bin").mkdir()
        self.executable(install / "bin" / "altann", '''#!/usr/bin/env python3
import json
from pathlib import Path
import sys
output = Path(sys.argv[sys.argv.index('--output') + 1])
output.write_text('##gff-version 3\\n')
Path(str(output) + '.tsv').write_text('segment\\tstrand\\n')
Path(str(output) + '.json').write_text(json.dumps({'arguments': sys.argv[1:]}))
''')
        args = argparse.Namespace(input=self.args.outdir, output=self.root / "alternatives.gff3",
                                  altann_dir=install, k=3, flank=1000, threads="2", overlap=0, margin=0)
        for job in ALTANN.run(args)["jobs"]:
            command = job["result"]["arguments"]
            self.assertIn("--rerun-viterbi", command)
            self.assertNotIn("--log", command)

    def test_segment_size_zero_keeps_sequences_whole(self):
        fasta = self.root / "genome.fa"
        fasta.write_text(">long\n" + "A" * 500 + "\n>short\n" + "C" * 50 + "\n")
        whole, cut = self.root / "whole", self.root / "cut"
        script = ROOT / "bin" / "segment_genome.py"
        subprocess.run(["python3", str(script), "segment", str(fasta), "--size", "0", "--min-len", "100",
                        "--outdir", str(whole)], check=True, capture_output=True)
        self.assertEqual([p.name for p in whole.iterdir()], ["long__0-500-500.fa"])
        subprocess.run(["python3", str(script), "segment", str(fasta), "--size", "300", "--overlap", "100",
                        "--min-len", "100", "--outdir", str(cut)], check=True, capture_output=True)
        self.assertEqual(sorted(p.name for p in cut.iterdir()), ["long__0-300-500.fa", "long__200-500-500.fa"])

    def test_zero_k_needs_no_altann_installation(self):
        UNIANN.run(self.args)
        output = self.root / "alternatives.gff3"
        args = argparse.Namespace(input=self.args.outdir, output=output, altann_dir=self.root / "absent",
                                  k=0, flank=1000, threads="1", overlap=0, margin=0)
        metadata = ALTANN.run(args)
        self.assertTrue(metadata["completed"])
        self.assertEqual(output.read_text(), "##gff-version 3\n")
        self.assertTrue(Path(str(output) + ".tsv").read_text().startswith("segment\tstrand\t"))

    def test_altann_passes_saved_scores_reverse_and_segment_settings(self):
        UNIANN.run(self.args)
        install = self.root / "altann"
        (install / "bin").mkdir(parents=True)
        self.executable(install / "bin" / "altann", '''#!/usr/bin/env python3
import json
from pathlib import Path
import sys
output = Path(sys.argv[sys.argv.index('--output') + 1])
output.write_text('##gff-version 3\\n')
Path(str(output) + '.tsv').write_text('segment\\tstrand\\n')
Path(str(output) + '.json').write_text(json.dumps({'arguments': sys.argv[1:]}))
''')
        output = self.root / "alternatives.gff3"
        args = argparse.Namespace(input=self.args.outdir, output=output, altann_dir=install,
                                  k=3, flank=1000, threads="2", overlap=0, margin=0)
        metadata = ALTANN.run(args)
        self.assertEqual(len(metadata["jobs"]), 2)
        for index, job in enumerate(metadata["jobs"]):
            command = job["result"]["arguments"]
            self.assertIn("--log", command)
            self.assertIn("--emissions", command)
            self.assertIn("--gff", command)
            self.assertIn("--short-cds-filter", command)
            self.assertNotIn("--psauron", command)
            self.assertNotIn("--rerun-viterbi", command)
            self.assertEqual("--reverse" in command, index == 1)
        self.assertTrue(json.loads(Path(str(output) + ".json").read_text())["completed"])


class FrameColumnAdapterTests(unittest.TestCase):
    def test_recognizes_only_known_assignments_ignoring_whitespace_and_comments(self):
        with tempfile.TemporaryDirectory() as directory:
            converter = Path(directory) / "converter.pl"
            for formula, expected in (("( $n - $i - 1 ) % 3", (1, 2, 0)),
                                      ("( $n - $i ) % 3", (0, 1, 2))):
                converter.write_text("# $coding_frame[($n-$i-2)%3] = 0;\n"
                                     + "$coding_frame[" + formula + "] = log($F[$i+3]);\n")
                self.assertEqual(UNIANN.converter_frame_columns(converter)[1], expected)
            for source in ("", "$coding_frame[($n-$i-2)%3] = 0;",
                           "$coding_frame[($n-$i)%3] = 0; $coding_frame[($n-$i-1)%3] = 0;"):
                converter.write_text(source)
                with self.assertRaisesRegex(ValueError, "Unrecognized UniAnn coding-frame"):
                    UNIANN.converter_frame_columns(converter)

    @unittest.skipUnless(shutil.which("perl"), "Perl is required for the real converter contract")
    def test_real_perl_legacy_and_corrected_emissions_match_for_all_local_frames_and_orientations(self):
        # Exercise the actual installed source, changing only its frame-index
        # expression in temporary copies. Calibration and score formulas stay
        # identical; no dependency repository or CNN prediction is changed.
        actual = (ROOT / "uniann" / "scripts" / "convert_coding_frame.pl").read_text()
        expression = r"\(\s*\$n\s*-\s*\$i\s*(?:-\s*1\s*)?\)\s*%\s*3"
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            scripts = {}
            for label, formula in (("legacy", "($n-$i-1)%3"), ("corrected", "($n-$i)%3")):
                source, count = re.subn(expression, lambda _: formula, actual, count=1)
                self.assertEqual(count, 1)
                scripts[label] = directory / (label + ".pl")
                scripts[label].write_text(source)
            forward = "CCCATGAAACCCGGGTAACCC"
            for orientation, sequence in (("plus", forward),
                                          ("minus", forward.translate(str.maketrans("ACGT", "TGCA"))[::-1])):
                fasta = directory / (orientation + ".fna")
                fasta.write_text(">chr\n" + sequence + "\n")
                for phase in range(3):
                    with self.subTest(orientation=orientation, local_frame=phase):
                        source = directory / "coding.tsv"
                        probabilities = ["0.03", "0.03", "0.03", "0.09"]
                        probabilities[phase] = "0.85"
                        source.write_text(UNIANN.CODING_HEADER + "\n" + "".join(
                            "\t".join(["chr", str(n + 1), "+", *probabilities, f"frame{phase}"]) + "\n"
                            for n in range(len(sequence))))
                        original = source.read_bytes()
                        outputs = []
                        for label in ("legacy", "corrected"):
                            staged = directory / (label + ".tsv")
                            _, columns = UNIANN.converter_frame_columns(scripts[label])
                            UNIANN.prepare_coding(source, staged, "chr", len(sequence), columns)
                            result = subprocess.run(["perl", str(scripts[label]), str(fasta)],
                                                    input=staged.read_bytes(), capture_output=True)
                            self.assertEqual(result.returncode, 0, result.stderr.decode())
                            outputs.append(result.stdout)
                        self.assertEqual(outputs[0], outputs[1])
                        self.assertEqual(source.read_bytes(), original)
                        rows = outputs[0].decode().splitlines()
                        self.assertEqual(len(rows), len(sequence))
                        for n, line in enumerate(rows):
                            coding = list(map(float, line.split("\t")[2:5]))
                            expected = (n - phase) % 3
                            if coding[expected] != -1000000:
                                self.assertEqual(coding.index(max(coding)), expected)


if __name__ == "__main__":
    unittest.main()
