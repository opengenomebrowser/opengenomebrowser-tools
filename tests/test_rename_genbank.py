from unittest import TestCase

import os
from opengenomebrowser_tools.rename_genbank import *

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_gbk.gbk'

gbks = [
    f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gbk',
    f'{ROOT}/test-data/prokka-good/PROKKA_08112021.gbk',
    f'{ROOT}/test-data/pgap-bad/annot.gbk',
    f'{ROOT}/test-data/pgap-good/annot.gbk'
]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_detect_locus_tag_prefix(self):
        for gbk in gbks:
            strain, locus_tag_prefix = GenBankFile(gbk).detect_strain_locus_tag_prefix()
            self.assertIn(member=strain, container=['replaceme', 'STRAIN'])
            self.assertIn(member=locus_tag_prefix, container=['tmp_', 'STRAIN.1_'])

    def test_get_taxid(self):
        for gbk in gbks:
            self.assertEqual(GenBankFile(gbk).taxid(), 2097)

    def test_rename(self):
        for gbk in gbks:
            cleanup()
            GenBankFile(gbk).rename(new_locus_tag_prefix='YOLO_', out=TMPFILE, validate=True)
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('YOLO_')
            self.assertNotIn(member='tmp_', container=content)
            self.assertNotIn(member='STRAIN.1_', container=content)
            self.assertGreater(a=count, b=1000)

    def test_get_metadata(self):
        for gbk in gbks:
            organism_data, genome_data = GenBankFile(gbk).metadata()
            self.assertIn(genome_data['cds_tool'], container=['PGAP', 'prokka'])
            self.assertIn(genome_data['cds_tool_date'], container=['2021-08-10', '2021-08-11'])
            self.assertIn(genome_data['cds_tool_version'], container=['2021-07-01.build5508', '1.14.5'])

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
