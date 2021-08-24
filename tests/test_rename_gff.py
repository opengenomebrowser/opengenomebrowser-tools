from unittest import TestCase

import os
from opengenomebrowser_tools.rename_gff import GffFile

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_gff.gff'

gffs = [
    f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gff',
    f'{ROOT}/test-data/prokka-good/PROKKA_08112021.gff',
    f'{ROOT}/test-data/pgap-bad/annot.gff',
    f'{ROOT}/test-data/pgap-good/annot.gff'
]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_detect_locus_tag_prefix(self):
        for gff in gffs:
            locus_tag_prefix = GffFile(gff).detect_locus_tag_prefix()
            self.assertIn(member=locus_tag_prefix, container=['tmp', 'STRAIN.1'])

    def test_rename(self):
        for gff in gffs:
            cleanup()
            GffFile(gff).rename(new_locus_tag_prefix='YOLO', out=TMPFILE, validate=True)
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('YOLO')
            self.assertNotIn(member='tmp', container=content, msg=f'Found "tmp" in renamed {gff=}!')
            self.assertNotIn(member='STRAIN.1', container=content, msg=f'Found "STRAIN.1" in renamed {gff=}!')
            self.assertGreater(a=count, b=1000)

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
