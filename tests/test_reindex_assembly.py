from unittest import TestCase

import os
from opengenomebrowser_tools.reindex_assembly import reindex_assembly

ROOT = os.path.dirname(os.path.dirname(__file__))
INFILE = f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.ffn'  # the assembly ASM2732v1.annotation.nucleotide.1.fasta has only one contig.
TMPFILE = '/tmp/reindexed_assembly.fasta'


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_reindex_assembly(self):
        reindex_assembly(file=INFILE, out=TMPFILE, prefix='TEST')
        with open(TMPFILE) as f:
            firstline = f.readline()
            self.assertEqual(firstline, '>TEST_1\n')
            for line in f:
                if line.startswith('>'):
                    self.assertEqual(line, '>TEST_2\n')
                    break

    def test_reindex_assembly_leading_zeroes(self):
        reindex_assembly(file=INFILE, out=TMPFILE, prefix='TEST', leading_zeroes=5)
        with open(TMPFILE) as f:
            firstline = f.readline()
            self.assertEqual(firstline, '>TEST_00001\n')
            for line in f:
                if line.startswith('>'):
                    self.assertEqual(line, '>TEST_00002\n')
                    break

    def setUp(self) -> None:
        cleanup()

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
