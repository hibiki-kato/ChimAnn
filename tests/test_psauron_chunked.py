"""Check six-frame chunk boundaries against whole-sequence codon identities."""
import csv
from contextlib import redirect_stderr
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


MODULE = Path(__file__).resolve().parents[1] / 'bin' / 'psauron_chunked.py'
spec = importlib.util.spec_from_file_location('psauron_chunked', MODULE)
chunked = importlib.util.module_from_spec(spec)
spec.loader.exec_module(chunked)
STOPS = {'TAA', 'TAG', 'TGA'}


def native_sequences(sequence):
    forward = sequence.upper().replace('N', 'A')
    reverse = ''.join({'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(base, 'A')
                      for base in reversed(forward))
    return forward, reverse


def identity(direction, position):
    return .1 + direction * .5 + position / 1000


def native_stub(original, chunk_marker=False):
    """Model native six-frame array layout with identifiable codon scores."""
    def run(sequence, workdir, tag):
        start = int(tag[1:])
        end = start + len(sequence)
        arrays = []
        for direction, oriented in enumerate(native_sequences(sequence)):
            origin = start if direction == 0 else len(original) - end
            for frame in range(3):
                values = [str(identity(direction, origin + position) + (start / 10000 if chunk_marker else 0))
                          for position in range(frame, len(oriented) - 2, 3)
                          if oriented[position:position + 3] not in STOPS]
                arrays.append(values if values else ['.99'] * 12)
        return arrays
    return run


class ChunkedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='psauron-chunks-test-')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.stderr = self.enterContext(redirect_stderr(io.StringIO()))

    def stitch(self, sequence, chunk=21, overlap=6, marker=False):
        with patch.object(chunked, 'run_psauron', side_effect=native_stub(sequence, marker)):
            return chunked.stitch(sequence, chunk, overlap, self.work)

    def assert_whole_identities(self, sequence, columns):
        self.assertEqual(len(columns), 6)
        for direction, oriented in enumerate(native_sequences(sequence)):
            for frame in range(3):
                expected = [identity(direction, position) for position in range(frame, len(oriented) - 2, 3)
                            if oriented[position:position + 3] not in STOPS]
                actual = [float(value) for value in columns[direction * 3 + frame].split(';') if value]
                self.assertEqual(len(actual), len(expected))
                for score, want in zip(actual, expected):
                    self.assertAlmostEqual(score, want, places=6)

    def test_all_length_residues_rephase_reverse_chunk_frames(self):
        for length in (60, 61, 62):
            for size, overlap in ((21, 6), (24, 9), (21, 3)):
                with self.subTest(length=length, chunk=size, overlap=overlap):
                    sequence = ('ATGGACTAAGGGTAGCCCTGACAT' * 3)[:length]
                    self.assert_whole_identities(sequence, self.stitch(sequence, size, overlap))

    def test_native_ambiguity_normalization_keeps_correct_stop_counts(self):
        for length in (60, 61, 62):
            sequence = ('ATNTNATAATRGCCCYTNACGNSWKMBDHV' * 3)[:length]
            with self.subTest(length=length):
                self.assert_whole_identities(sequence, self.stitch(sequence))

    def test_one_chunk_keeps_all_native_frame_probabilities(self):
        sequence = 'ATNTNATAATRGCCCYTNACGNSWKMBDHV' * 2
        self.assert_whole_identities(sequence, self.stitch(sequence, chunk=90, overlap=6))

    def test_first_overlapping_chunk_wins_on_each_strand(self):
        sequence = 'GCC' * 20
        columns = self.stitch(sequence, marker=True)
        self.assertAlmostEqual(float(columns[0].split(';')[5]), identity(0, 15), places=6)
        self.assertAlmostEqual(float(columns[3].split(';')[13]), identity(1, 39), places=6)
        self.assertAlmostEqual(float(columns[0].split(';')[7]), identity(0, 21) + .0015, places=6)
        self.assertAlmostEqual(float(columns[3].split(';')[8]), identity(1, 24) + .0015, places=6)

    def test_empty_alternate_frames_ignore_native_placeholder_scores(self):
        sequence = 'ATAATAATAATAATAATAATAATAAT'
        columns = self.stitch(sequence, chunk=18, overlap=6)
        self.assert_whole_identities(sequence, columns)
        self.assertEqual(columns[1], '')

    def test_probability_count_mismatch_is_rejected_for_both_directions(self):
        sequence = 'GCC' * 20
        for frame in (0, 3):
            for change in ('missing', 'extra'):
                def bad(seq, workdir, tag):
                    arrays = native_stub(sequence)(seq, workdir, tag)
                    arrays[frame] = arrays[frame][:-1] if change == 'missing' else arrays[frame] + ['.5']
                    return arrays
                with self.subTest(frame=frame, change=change), patch.object(chunked, 'run_psauron', side_effect=bad):
                    with self.assertRaisesRegex(ValueError, 'count mismatch'):
                        chunked.stitch(sequence, 21, 6, self.work)

    def test_invalid_probability_is_rejected(self):
        sequence = 'GCC' * 20
        for value in ('nan', 'inf', '-0.1', '1.1', 'broken'):
            def bad(seq, workdir, tag):
                arrays = native_stub(sequence)(seq, workdir, tag)
                arrays[4][0] = value
                return arrays
            with self.subTest(value=value), patch.object(chunked, 'run_psauron', side_effect=bad):
                with self.assertRaises(ValueError):
                    chunked.stitch(sequence, 21, 6, self.work)

    def test_invalid_chunk_bounds_never_start_psauron(self):
        with patch.object(chunked, 'run_psauron') as run:
            for size, overlap in ((0, 3), (3, 3), (6, 0), (6, 2), (6, 6), (6, 9)):
                with self.subTest(chunk=size, overlap=overlap), self.assertRaises(ValueError):
                    chunked.stitch('GCC' * 20, size, overlap, self.work)
            run.assert_not_called()

    def test_cli_writes_four_preamble_lines_and_six_columns(self):
        sequence = 'GCC' * 20 + 'AC'
        fasta, output = self.work / 'sequence.fa', self.work / 'scores.csv'
        fasta.write_text('>original description\n' + sequence + '\n')
        with patch.object(chunked, 'run_psauron', side_effect=native_stub(sequence)):
            chunked.main([str(fasta), '--chunk', '21', '--overlap', '6', '--out', str(output)])
        with output.open() as handle:
            self.assertTrue(handle.readline().startswith('psauron_chunked.py'))
            handle.readline()
            self.assertIn('six frames', handle.readline())
            header, row = list(csv.reader(handle))
        self.assertEqual(header, chunked.COLS)
        self.assertEqual(row[0], 'original')
        self.assertEqual(len(row), 15)
        self.assert_whole_identities(sequence, row[9:])

    def test_cli_failure_preserves_existing_output(self):
        fasta, output = self.work / 'sequence.fa', self.work / 'scores.csv'
        fasta.write_text('>original\nGCCGCC\n')
        output.write_text('existing output\n')
        with patch.object(chunked, 'run_psauron') as run, self.assertRaises(SystemExit):
            chunked.main([str(fasta), '--chunk', '0', '--out', str(output)])
        run.assert_not_called()
        self.assertEqual(output.read_text(), 'existing output\n')

    def test_native_command_requests_all_six_frames_on_cpu(self):
        def write_csv(command, **kwargs):
            self.assertIn('-a', command)
            self.assertIn('-c', command)
            self.assertEqual((self.work / 'chunk.fa').read_text(), '>chunk\nACGNRY\n')
            output = Path(command[command.index('-o') + 1])
            with output.open('w') as handle:
                handle.write('metadata\nmetadata\nmetadata\n')
                writer = csv.writer(handle)
                writer.writerow(chunked.COLS)
                writer.writerow(['chunk', *['NA'] * 8, *[f'.{index + 1}' for index in range(6)]])
        with patch.object(chunked.subprocess, 'run', side_effect=write_csv):
            arrays = chunked.run_psauron('ACGNRY', self.work, 'chunk')
        self.assertEqual(arrays, [[f'.{index + 1}'] for index in range(6)])

    def test_native_skip_has_friendly_error(self):
        def write_csv(command, **kwargs):
            output = Path(command[command.index('-o') + 1])
            with output.open('w') as handle:
                csv.writer(handle).writerow(chunked.COLS)
        with patch.object(chunked.subprocess, 'run', side_effect=write_csv):
            with self.assertRaisesRegex(ValueError, 'did not score chunk'):
                chunked.run_psauron('TAA' * 10, self.work, 'chunk')


if __name__ == '__main__':
    unittest.main()
