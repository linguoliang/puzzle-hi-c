import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from utils.alignment_input import detect_format, normalize_contacts


class TextInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.output = self.root / 'normalized.txt'

    def convert(self, content, filename='reads.txt', **kwargs):
        source = self.root / filename
        if filename.endswith('.gz'):
            with gzip.open(source, 'wt') as handle:
                handle.write(content)
        else:
            source.write_text(content)
        count = normalize_contacts(source, self.output, **kwargs)
        return count, self.output.read_text()

    def test_legacy_juicer_long_and_short_preserve_fields(self):
        short = '0 ctg1 10 7 16 ctg2 20 9'
        count, output = self.convert(short + ' 60 10M A 40 10M A r1 r2\n' + short + '\n')
        self.assertEqual(count, 2)
        self.assertEqual(output, '\t'.join(short.split()) + '\n' + '\t'.join(short.split()) + '\n')

    def test_pairs_standard_gzip_comments_blanks_and_no_final_newline(self):
        count, output = self.convert('#pairs format v1.0\n\nr1 ctg1 10 ctg2 20 + -', 'reads.pairs.gz')
        self.assertEqual(count, 1)
        self.assertEqual(output, '0\tctg1\t10\t0\t16\tctg2\t20\t1\n')

    def test_pairs_reordered_columns_and_chrom_aliases(self):
        count, output = self.convert('#columns: pos2 readID chrom2 strand2 chrom1 strand1 pos1 extra\n'
                                     '20 r1 ctg2 - ctg1 + 10 x\n', 'reads.pairs')
        self.assertEqual(count, 1)
        self.assertEqual(output, '0\tctg1\t10\t0\t16\tctg2\t20\t1\n')

    def test_pairs_filters_unmapped_duplicates_and_nonunique(self):
        header = '#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 pair_type\n'
        rows = [f'r{i} ctg1 10 ctg2 20 + - {kind}\n'
                for i, kind in enumerate(['UU', 'UR', 'RU', 'DD', 'MM', 'MU', 'WW'])]
        rows.append('unmapped ! 0 ctg2 20 - - NU\n')
        count, output = self.convert(header + ''.join(rows), 'reads.pairs')
        self.assertEqual(count, 3)
        self.assertEqual(len(output.splitlines()), 3)

    def test_pairs_mapq_filters_both_ends_and_unknown(self):
        header = '#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 mapq1 mapq2\n'
        rows = [f'r ctg1 10 ctg2 20 + - {a} {b}\n'
                for a, b in [(60, 30), (29, 60), (60, 29), (255, 60)]]
        self.assertEqual(self.convert(header + ''.join(rows), 'reads.pairs', min_mapq=30)[0], 1)

    def test_hicpro_preserves_coordinates_ignores_fragment_names(self):
        count, output = self.convert('r1 ctg1 10 + ctg2 20 - 100 HIC_chr1_1 HIC_chr2_2 60 60\n',
                                     'sample.allValidPairs.gz')
        self.assertEqual(count, 1)
        self.assertEqual(output, '0\tctg1\t10\t0\t16\tctg2\t20\t1\n')

    def test_hicpro_mapq_and_juicer_mapq(self):
        valid = 'r ctg1 10 + ctg2 20 - 100 fragA fragB '
        self.assertEqual(self.convert(valid + '60 60\n' + valid + '10 60\n',
                                      'x.validPairs', min_mapq=30)[0], 1)
        juicer = '0 ctg1 10 0 16 ctg2 20 1 '
        self.assertEqual(self.convert(juicer + '60 10M A 60 10M A r r\n' +
                                     juicer + '60 10M A 10 10M A r r\n', min_mapq=30)[0], 1)

    def test_missing_mapq_is_actionable(self):
        for name, content in [('x.pairs', 'r ctg1 10 ctg2 20 + -\n'),
                              ('x.validPairs', 'r ctg1 10 + ctg2 20 -\n'),
                              ('x.txt', '0 ctg1 10 0 16 ctg2 20 1\n')]:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'min-mapq requires'):
                self.convert(content, name, min_mapq=1)

    def test_invalid_input_does_not_replace_existing_output(self):
        self.output.write_text('keep me\n')
        for record in ['r ctg1 0 ctg2 20 + -', 'r ctg1 xx ctg2 20 + -',
                       'r ctg1 10 ctg2 20 x -', 'too short']:
            with self.subTest(record=record), self.assertRaisesRegex(ValueError, ':2:'):
                self.convert('# comment\n' + record, 'x.pairs')
            self.assertEqual(self.output.read_text(), 'keep me\n')
        self.assertEqual(list(self.root.glob('.contacts-*')), [])

    def test_empty_input_fails(self):
        with self.assertRaisesRegex(ValueError, 'no usable contact'):
            self.convert('# comment\n\n', 'x.pairs')

    def test_duplicate_missing_or_late_column_headers_fail(self):
        for text in ['#columns: readID chr1 pos1 chr2 pos2 strand1 strand1\n',
                     '#columns: readID chr1 pos1\n',
                     'r ctg1 10 ctg2 20 + -\n#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n']:
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, '#columns'):
                self.convert(text, 'x.pairs')

    def test_auto_detection_and_explicit_override(self):
        for name, expected in [('x.BAM', 'bam'), ('x.sam', 'sam'), ('x.pairs.gz', 'pairs'),
                               ('x.allValidPairs', 'validpairs'), ('x.txt', 'juicer')]:
            self.assertEqual(detect_format(name), expected)
        self.assertEqual(self.convert('r ctg1 10 ctg2 20 + -\n', input_format='pairs')[0], 1)
        for name in ['x.matrix', 'x.matrix.gz', 'x.cool', 'x.mcool', 'x.hic', 'x.cram']:
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, 'not supported'):
                detect_format(name)

    def test_matrix_without_extension_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, 'allValidPairs, not an aggregate'):
            self.convert('1 2 100\n')

    def test_weighted_juicer_and_ambiguous_pairs_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'weighted nine-column'):
            self.convert('0 ctg1 10 0 16 ctg2 20 1 10\n')
        with self.assertRaisesRegex(ValueError, '#columns'):
            self.convert('r ctg1 10 ctg2 20 + - DD\n', 'x.pairs')

    def test_same_file_and_invalid_options(self):
        self.output.write_text('unchanged')
        with self.assertRaisesRegex(ValueError, 'different files'):
            normalize_contacts(self.output, self.output)
        self.assertEqual(self.output.read_text(), 'unchanged')
        for opts in [{'input_format': 'matrix'}, {'min_mapq': -1}, {'min_mapq': 256}]:
            with self.subTest(opts=opts), self.assertRaises(ValueError):
                normalize_contacts(self.output, self.root / 'other', **opts)

    def test_missing_pysam_explains_optional_dependency(self):
        with patch.dict('sys.modules', {'pysam': None}):
            with self.assertRaisesRegex(ValueError, 'requirements-alignment.txt'):
                self.convert('', 'reads.bam')


@unittest.skipUnless(importlib.util.find_spec('pysam'), 'optional pysam is not installed')
class AlignmentInputTests(unittest.TestCase):
    setUp = TextInputTests.setUp
    # Binary paths use actual pysam files.
    def write_alignment(self, filename, groups, order='queryname'):
        import pysam
        path = self.root / filename
        header = {'HD': {'VN': '1.6', 'SO': order},
                  'SQ': [{'SN': 'ctg1', 'LN': 1000}, {'SN': 'ctg2', 'LN': 1000}]}
        with pysam.AlignmentFile(str(path), 'wb' if filename.endswith('.bam') else 'w', header=header) as out:
            for name, changes in groups:
                for side, overrides in enumerate(changes):
                    r = pysam.AlignedSegment(out.header)
                    r.query_name = name
                    r.flag = 65 if side == 0 else 145  # Paired, inter-contig, not proper-pair.
                    r.reference_id = side % 2
                    r.reference_start = 9 if side == 0 else 19
                    r.cigarstring = '5S10M2D5M3S'
                    r.query_sequence = 'A' * 23
                    r.mapping_quality = 60
                    r.next_reference_id = 1 - side % 2
                    r.next_reference_start = 19 if side == 0 else 9
                    for key, value in overrides.items():
                        if key == 'SA':
                            r.set_tag('SA', value)
                        else:
                            setattr(r, key, value)
                    out.write(r)
        return path

    def test_bam_and_sam_same_coordinates_once_per_template(self):
        for suffix in ['bam', 'sam']:
            with self.subTest(suffix=suffix):
                source = self.write_alignment('reads.' + suffix, [('read1', [{}, {}])])
                self.assertEqual(normalize_contacts(source, self.output), 1)
                # Reverse uses CIGAR reference span: 19 + 10 + 2 + 5 = 36.
                self.assertEqual(self.output.read_text(), '0\tctg1\t10\t0\t16\tctg2\t36\t1\n')

    def test_bam_filters_pairs_not_individual_mates(self):
        groups = [('a_good', [{}, {}]),
                  ('b_dup', [{}, {'flag': 145 | 1024}]),
                  ('c_qcfail', [{'flag': 65 | 512}, {}]),
                  ('d_unmapped', [{'flag': 65 | 4}, {}]),
                  ('e_low', [{'mapping_quality': 10}, {}]),
                  ('f_orphan', [{}]),
                  ('g_ambiguous', [{}, {}, {'flag': 65}]),
                  ('h_sa', [{'SA': 'ctg2,20,+,23M,60,0;'}, {}]),
                  ('i_supplementary', [{}, {}, {'flag': 65 | 2048}]),
                  ('j_secondary', [{}, {}, {'flag': 65 | 256}]),
                  ('k_unknown_mapq', [{'mapping_quality': 255}, {}]),
                  ('l_same_side', [{}, {'flag': 65}]),
                  ('m_unpaired', [{'flag': 64}, {}]),
                  ('n_mate_unmapped', [{'flag': 65 | 8}, {}])]
        source = self.write_alignment('reads.bam', groups)
        self.assertEqual(normalize_contacts(source, self.output, min_mapq=30), 2)
        self.assertEqual(len(self.output.read_text().splitlines()), 2)

    def test_alignment_requires_name_sorted_header(self):
        for sort in ['coordinate', 'unknown', 'unsorted']:
            source = self.write_alignment('reads.bam', [('read', [{}, {}])], order=sort)
            with self.subTest(sort=sort), self.assertRaisesRegex(ValueError, 'samtools sort -n'):
                normalize_contacts(source, self.output)

    def test_read2_can_appear_before_read1(self):
        source = self.write_alignment('reads.bam', [('read', [{'flag': 145}, {'flag': 65}])])
        self.assertEqual(normalize_contacts(source, self.output), 1)
        self.assertEqual(self.output.read_text(), '0\tctg2\t20\t0\t16\tctg1\t26\t1\n')


if __name__ == '__main__':
    unittest.main()
