"""Small deterministic whole-pipeline checks; no Juicer binary or real dataset needed."""
import gzip
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = ('Bio', 'h5py', 'networkx', 'numpy', 'pandas', 'scipy', 'matplotlib')
HAS_PIPELINE = all(importlib.util.find_spec(module) for module in DEPENDENCIES)


@unittest.skipUnless(HAS_PIPELINE and shutil.which('split') and shutil.which('sort')
                     and (os.cpu_count() or 1) > 1,
                     'pipeline dependencies, GNU split/sort and at least two CPUs are required')
class PipelineInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.root = Path(cls.tmp.name)
        cls.fasta = cls.root / 'contigs.fa'
        cls.fasta.write_text(''.join(f'>ctg{i}\n' + 'ACGT' * 2500 + '\n' for i in range(1, 5)))
        cls.contacts = []
        # Same-bin contacts give nonzero coverage; end links assemble two chromosomes.
        for i in range(1, 5):
            for pos in range(50, 10000, 100):
                cls.contacts.extend([(f'ctg{i}', pos, f'ctg{i}', pos + 10)] * 40)
        for first, second in [(1, 2), (3, 4)]:
            for x in range(0, 600, 20):
                for y in range(0, 600, 20):
                    cls.contacts.append((f'ctg{first}', 10000 - x, f'ctg{second}', 1 + y))
        cls.juicer = cls.root / 'merged_nodups.txt'
        cls.juicer.write_text(''.join(f'0 {a} {x} 0 16 {b} {y} 1\n' for a, x, b, y in cls.contacts))
        cls.pairs = cls.root / 'reads.pairs.gz'
        with gzip.open(cls.pairs, 'wt') as output:
            output.write('#columns: readID chr1 pos1 chr2 pos2 strand1 strand2\n')
            output.writelines(f'r{i:06d} {a} {x} {b} {y} + -\n'
                              for i, (a, x, b, y) in enumerate(cls.contacts))
        cls.validpairs = cls.root / 'sample.allValidPairs'
        cls.validpairs.write_text(''.join(f'r{i:06d} {a} {x} + {b} {y} -\n'
                                         for i, (a, x, b, y) in enumerate(cls.contacts)))
        cls.reference = cls.run_pipeline('juicer', cls.juicer)
        cls.expected = {suffix: (cls.reference / ('sample' + suffix)).read_bytes()
                        for suffix in ('.agp', '.fa', '.Chrom.sizes')}

    @classmethod
    def run_pipeline(cls, name, source, export=False, cpus=1):
        directory = cls.root / name
        directory.mkdir()
        command = [sys.executable, str(ROOT / 'main.py'), '-c', '2', '-m', str(source),
                   '-f', str(cls.fasta), '-p', 'sample', '-s', '100', '-i', '3', '-n', str(cpus)]
        if export:
            # Verify the real export boundary and pair-count preservation, not .hic encoding.
            stub = cls.root / 'juicer_tools_stub'
            stub.write_text('#!' + sys.executable + '\n'
                            'import pathlib, sys\n'
                            'assert sys.argv[1] == "pre"\n'
                            f'assert len(pathlib.Path(sys.argv[2]).read_text().splitlines()) == {len(cls.contacts)}\n'
                            'pathlib.Path(sys.argv[3]).write_text("stub export boundary only")\n')
            stub.chmod(0o755)
            command.extend(['-j', str(stub)])
        else:
            command.append('--skip-hic')
        result = subprocess.run(command, cwd=directory, capture_output=True, text=True, timeout=120,
                                env={**os.environ, 'MPLCONFIGDIR': str(cls.root / 'mpl'),
                                     'XDG_CACHE_HOME': str(cls.root / 'cache'), 'PYTHONHASHSEED': '0'})
        if result.returncode:
            raise AssertionError(result.stdout + '\n' + result.stderr)
        if f'Loaded {len(cls.contacts)} contact pairs.' not in result.stdout:
            raise AssertionError(result.stdout)
        return directory

    def assert_equivalent(self, directory):
        for suffix, expected in self.expected.items():
            self.assertEqual((directory / ('sample' + suffix)).read_bytes(), expected, suffix)
        self.assertEqual(len((directory / 'merged_nodups_short_format.txt').read_text().splitlines()),
                         len(self.contacts))

    def test_juicer_still_assembles_and_skip_hic_omits_export(self):
        self.assertEqual(len(self.expected['.Chrom.sizes'].splitlines()), 2)
        self.assertIn(b'\tU\t100\t', self.expected['.agp'])
        self.assertFalse((self.reference / 'sample.hic').exists())
        self.assertFalse((self.reference / 'tmp' / 'convert.h5').exists())

    def test_gzipped_pairs_equivalent_to_juicer(self):
        self.assert_equivalent(self.run_pipeline('pairs', self.pairs))

    def test_hicpro_equivalent_to_juicer(self):
        self.assert_equivalent(self.run_pipeline('validpairs', self.validpairs))

    @unittest.skipUnless(importlib.util.find_spec('pysam'), 'optional pysam is required')
    def test_bam_and_sam_equivalent_to_juicer(self):
        import pysam
        header = {'HD': {'VN': '1.6', 'SO': 'queryname'},
                  'SQ': [{'SN': f'ctg{i}', 'LN': 10000} for i in range(1, 5)]}
        for extension, mode in [('bam', 'wb'), ('sam', 'w')]:
            source = self.root / ('reads.' + extension)
            with pysam.AlignmentFile(str(source), mode, header=header) as output:
                for i, (a, x, b, y) in enumerate(self.contacts):
                    for side, (chrom, position) in enumerate([(a, x), (b, y)]):
                        record = pysam.AlignedSegment(output.header)
                        record.query_name = f'r{i:06d}'
                        record.flag = 65 if side == 0 else 145
                        record.reference_id = int(chrom[3:]) - 1
                        record.reference_start = position - 1
                        record.query_sequence = 'A'
                        record.cigarstring = '1M'
                        record.mapping_quality = 60
                        record.next_reference_id = int((b if side == 0 else a)[3:]) - 1
                        record.next_reference_start = (y if side == 0 else x) - 1
                        output.write(record)
            with self.subTest(extension=extension):
                self.assert_equivalent(self.run_pipeline(extension, source, cpus=2 if os.cpu_count() > 2 else 1))

    def test_optional_hic_export_receives_all_pairs(self):
        directory = self.run_pipeline('export', self.pairs, export=True)
        self.assert_equivalent(directory)
        self.assertTrue((directory / 'sample.hic').exists())

    def test_export_converter_preserves_first_headerless_contact(self):
        from utils.convert_data import convert_data
        source = self.root / 'headerless.txt'
        source.write_text('0\tctg2\t20\t0\t16\tctg1\t10\t1\n'
                          '0\tctg1\t30\t0\t16\tctg2\t40\t1\n')
        output = self.root / 'headerless.re'
        convert_data({'ctg1': 100, 'ctg2': 100}, source, output)
        self.assertEqual(output.read_text(), '16\tctg1\t10\t1\t0\tctg2\t20\t0\n'
                                             '0\tctg1\t30\t0\t16\tctg2\t40\t1\n')

    def test_cli_requires_explicit_export_choice(self):
        result = subprocess.run([sys.executable, str(ROOT / 'main.py'), '-c', '2', '-m', str(self.pairs),
                                 '-f', str(self.fasta)], capture_output=True, text=True,
                                env={**os.environ, 'MPLCONFIGDIR': str(self.root / 'mpl')})
        self.assertEqual(result.returncode, 2)
        self.assertIn('--skip-hic', result.stderr)


if __name__ == '__main__':
    unittest.main()
