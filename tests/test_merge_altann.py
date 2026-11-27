"""Check preservation of UniAnn models and validation of AltAnn alternatives."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from urllib.parse import quote


SCRIPT = Path(__file__).resolve().parents[1] / 'bin/merge_altann.py'


def model(identifier, intervals, strand='+', seqid='chr', rank=None, parent=None):
    attrs = 'ID=' + identifier
    if parent:
        attrs += ';Parent=' + parent
    if rank is not None:
        attrs += (f';kbest_rank={rank};kbest_score=1234.5678901234567'
                  f';kbest_delta=-0.0000000000000123;kbest_locus={seqid}.{strand}.locus1;kbest_origin=' +
                  ('uniann_best' if rank == 0 else 'local_kbest'))
    start, end = min(row[0] for row in intervals), max(row[1] for row in intervals)
    rows = [f'{seqid}\tUniAnn\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t{attrs}']
    for index, (start, end, phase) in enumerate(intervals):
        rows.append(f'{seqid}\tUniAnn\texon\t{start}\t{end}\t.\t{strand}\t.\t'
                    f'ID={identifier}.exon{index};Parent={identifier}')
        rows.append(f'{seqid}\tUniAnn\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t'
                    f'Parent={identifier}')
    return rows


PLUS = [(10, 13, 0), (23, 27, 2)]
MINUS = [(201, 204, 0), (101, 105, 2)]
ALT_PLUS = [(10, 16, 0), (26, 27, 2)]
ALT_MINUS = [(201, 205, 0), (101, 104, 1)]


class MergeAltAnnTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='chimann-merge-altann-')
        self.addCleanup(self.tmp.cleanup)
        self.work = Path(self.tmp.name)
        self.original = [
            'chr\tUniAnn\tgene\t1\t900\t.\t+\t.\tID=gene',
            *model('best.f', PLUS, parent='gene'),
            *model('best.r', MINUS, '-'),
        ]
        self.alt = [*model('reference.f', PLUS, rank=0),
                    *model('candidate.f', ALT_PLUS, rank=1),
                    *model('reference.r', MINUS, '-', rank=0),
                    *model('candidate.r', ALT_MINUS, '-', rank=2)]

    def run_merge(self, original=None, alt=None, check=True):
        original = self.original if original is None else original
        alt = self.alt if alt is None else alt
        paths = [self.work / 'uniann.gff', self.work / 'altann.gff']
        for path, rows in zip(paths, (original, alt)):
            path.write_text('##gff-version 3\n' + '\n'.join(rows) + '\n')
        return subprocess.run([sys.executable, str(SCRIPT), *map(str, paths)],
                              capture_output=True, text=True, check=check)

    def assert_rejected(self, original=None, alt=None, message=None):
        result = self.run_merge(original, alt, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, '')
        if message:
            self.assertIn(message, result.stderr)

    def test_preserves_originals_and_adds_only_alternatives_verbatim(self):
        result = self.run_merge()
        expected = ['##gff-version 3', *self.original,
                    *model('candidate.f', ALT_PLUS, rank=1),
                    *model('candidate.r', ALT_MINUS, '-', rank=2)]
        self.assertEqual(result.stdout.splitlines(), expected)
        self.assertEqual(result.stdout.count('##gff-version 3'), 1)
        self.assertIn('kbest_score=1234.5678901234567', result.stdout)
        self.assertIn('kbest_delta=-0.0000000000000123', result.stdout)
        self.assertIn('\t-\t1\tParent=candidate.r', result.stdout)
        self.assertNotIn('reference.f', result.stdout)

    def test_additional_original_partial_model_is_retained(self):
        original = [*self.original, *model('partial', [(500, 504, 0)])]
        result = self.run_merge(original)
        self.assertEqual(result.stdout.splitlines()[1:1 + len(original)], original)

    def test_altann_sequence_escaping_is_decoded_once_without_changing_attributes(self):
        for seqid in ('chr|x', 'chr%2F', 'chr/x'):
            with self.subTest(seqid=seqid):
                original = model('original%2F', PLUS, seqid=seqid)
                candidate = model('candidate%2F', ALT_PLUS, seqid=seqid, rank=1)
                alt = [*model('reference%2F', PLUS, seqid=seqid, rank=0), *candidate]
                escaped = []
                for row in alt:
                    fields = row.split('\t')
                    fields[0] = quote(seqid, safe='._-:')
                    escaped.append('\t'.join(fields))
                result = self.run_merge(original=original, alt=escaped)
                self.assertEqual(result.stdout.splitlines(), ['##gff-version 3', *original, *candidate])
                self.assertIn('ID=candidate%2F;', result.stdout)
                self.assertIn('Parent=candidate%2F', result.stdout)

    def test_duplicate_cds_chains_are_removed_across_originals_and_alternatives(self):
        alt = [*self.alt, *model('same_as_best', PLUS, rank=3),
               *model('same_as_candidate', ALT_PLUS, rank=4)]
        result = self.run_merge(alt=alt)
        self.assertNotIn('same_as_', result.stdout)
        self.assertIn('ID=candidate.f;', result.stdout)

    def test_deduplication_keeps_distinct_sequences_and_strands(self):
        intervals = [(700, 708, 0)]
        original = [*self.original, *model('other_best', [(40, 48, 0)], seqid='other')]
        alt = [*self.alt, *model('other_reference', [(40, 48, 0)], seqid='other', rank=0),
               *model('chr_plus', intervals, '+', 'chr', rank=3),
               *model('chr_minus', intervals, '-', 'chr', rank=3),
               *model('other_plus', intervals, '+', 'other', rank=3)]
        result = self.run_merge(original=original, alt=alt)
        for identifier in ('chr_plus', 'chr_minus', 'other_plus'):
            self.assertIn('ID=' + identifier + ';', result.stdout)

    def test_reference_requires_matching_chain_and_phase(self):
        mismatch = model('bad_reference', ALT_PLUS, rank=0)
        self.assert_rejected(alt=[*self.alt, *mismatch], message='absent from UniAnn')
        original = [row.replace('\t+\t2\tParent=best.f', '\t+\t1\tParent=best.f')
                    for row in self.original]
        self.assert_rejected(original=original, message='absent from UniAnn')

    def test_transcript_requires_rank_and_consistent_child_ranks(self):
        self.assert_rejected(alt=model('missing_rank', ALT_PLUS), message='kbest_rank')
        for invalid in ('-1', 'one', '1.0'):
            with self.subTest(rank=invalid):
                rows = [row.replace('kbest_rank=1;', f'kbest_rank={invalid};') for row in self.alt]
                self.assert_rejected(alt=rows, message='kbest_rank')
        rows = list(self.alt)
        rows[6] += ';kbest_rank=0'
        self.assert_rejected(alt=rows, message='mixed ranks')

    def test_missing_duplicate_and_colliding_ids_are_rejected(self):
        self.assert_rejected(alt=[*self.alt, self.alt[0]], message='duplicate ID')
        self.assert_rejected(original=[*self.original, self.original[0]], message='duplicate ID')
        rows = [row.replace('candidate.f', 'best.f') for row in self.alt]
        self.assert_rejected(alt=rows, message='collides')
        rows = [row.replace('ID=candidate.f;', '', 1) if '\tmRNA\t' in row else row
                for row in self.alt]
        self.assert_rejected(alt=rows, message='requires an ID')

    def test_missing_parent_and_cross_strand_children_are_rejected(self):
        rows = [row.replace('Parent=candidate.f', 'Parent=missing') for row in self.alt]
        self.assert_rejected(alt=rows, message='missing Parent')
        rows = list(self.alt)
        rows[6] = rows[6].replace('\t+\t', '\t-\t')
        self.assert_rejected(alt=rows, message='different sequence or strand')

    def test_shared_children_and_transcript_parent_blocks_are_rejected(self):
        rows = list(self.alt)
        rows[1] = rows[1].replace('Parent=reference.f', 'Parent=reference.f,candidate.f')
        self.assert_rejected(alt=rows, message='multiple transcript blocks')
        rows = list(self.alt)
        rows[5] += ';Parent=reference.f'
        self.assert_rejected(alt=rows, message='transcript Parent must be a gene')

    def test_empty_cds_and_partial_altann_models_are_rejected(self):
        rows = [row for row in model('no_cds', ALT_PLUS, rank=1) if '\tCDS\t' not in row]
        self.assert_rejected(alt=rows, message='no CDS')
        self.assert_rejected(alt=model('partial', [(500, 504, 0)], rank=1),
                             message='complete coding model')
        self.assert_rejected(alt=model('bad_phase', [(500, 508, 1)], rank=1),
                             message='inconsistent CDS phases')

    def test_header_only_annotations_and_no_alternatives(self):
        result = self.run_merge(original=[], alt=[])
        self.assertEqual(result.stdout, '##gff-version 3\n')
        result = self.run_merge(alt=[])
        self.assertEqual(result.stdout.splitlines(), ['##gff-version 3', *self.original])

    def test_alternatives_without_owned_reference_are_omitted(self):
        result = self.run_merge(alt=model('only_alt', ALT_PLUS, rank=1))
        self.assertEqual(result.stdout.splitlines(), ['##gff-version 3', *self.original])
        rows = [row for row in self.alt if 'reference.r' not in row]
        result = self.run_merge(alt=rows)
        self.assertNotIn('candidate.r', result.stdout)
        self.assertIn('candidate.f', result.stdout)
        result = self.run_merge(original=[], alt=model('only_alt', ALT_PLUS, rank=1))
        self.assertEqual(result.stdout, '##gff-version 3\n')

    def test_alternative_requires_explicit_locus_metadata(self):
        rows = [row.replace(';kbest_locus=chr.+.locus1', '') if 'ID=candidate.f;' in row else row
                for row in self.alt]
        self.assert_rejected(alt=rows, message='requires kbest_locus')

    def test_help_explains_reference_invariant(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--help'],
                                capture_output=True, text=True, check=True)
        self.assertIn('Every AltAnn rank-zero CDS', result.stdout)
        self.assertIn('same genomic coordinates', result.stdout)


if __name__ == '__main__':
    unittest.main()
