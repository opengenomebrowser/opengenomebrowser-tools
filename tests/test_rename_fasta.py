from unittest import TestCase

import os
from opengenomebrowser_tools.rename_fasta import *

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_fasta.fasta'

fastas = [
    f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.',
    f'{ROOT}/test-data/prokka-good/PROKKA_08112021.',
    f'{ROOT}/test-data/pgap-bad/annot.',
    f'{ROOT}/test-data/pgap-good/annot.'
]
fastas = [f + suffix for suffix in ['faa', 'ffn'] for f in fastas][:-2]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_parse_fasta_header(self):
        for fasta in fastas:
            with open(fasta) as f:
                locus_tag_prefix, gene_id = FastaFile.parse_fasta_header(f.readline())
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])
            self.assertIn(member=gene_id, container=['00001', '000289'])

    def test_detect_locus_tag_prefix(self):
        for fasta in fastas:
            locus_tag_prefix = FastaFile(fasta).detect_locus_tag_prefix()
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])

    def test_rename(self):
        for fasta in fastas:
            cleanup()
            FastaFile(fasta).rename(new_locus_tag_prefix='YOLO_', out=TMPFILE, validate=True)
            with open(fasta) as f_old, open(TMPFILE) as f_new:
                content_old = f_old.read()
                content_new = f_new.read()
            self.assertNotIn(member='tmp', container=content_new)
            self.assertNotIn(member='STRAIN.1', container=content_new)
            self.assertEqual(
                first=content_new.count('YOLO_'),
                second=max(content_old.count('tmp_'), content_old.count('STRAIN.1_'))
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
