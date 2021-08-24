from unittest import TestCase

import os
from opengenomebrowser_tools.rename_eggnog import *

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_eggnog.eggnog'

eggnogs = [
    f'{ROOT}/test-data/prokka-bad/out.emapper.annotations',
    f'{ROOT}/test-data/pgap-bad/out.emapper.annotations',
]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):
    def test_detect_locus_tag_prefix(self):
        for eggnog in eggnogs:
            locus_tag_prefix = EggnogFile(file=eggnog).detect_locus_tag_prefix()
            self.assertEqual(locus_tag_prefix, 'tmp')

    def test_parse_eggnog_date(self):
        for eggnog in eggnogs:
            dt = EggnogFile(file=eggnog).date_str()
            self.assertIs(type(dt), str)

    def test_guess_eggnog_type(self):
        for eggnog in eggnogs:
            eggnog_type = EggnogFile(file=eggnog).custom_annotation_type
            self.assertTrue(eggnog_type.startswith('eggnog'))

    def test_validate_locus_tags(self):
        for eggnog in eggnogs:
            ef = EggnogFile(file=eggnog)
            ef.validate_locus_tags(locus_tag_prefix=None)
            ef.validate_locus_tags(locus_tag_prefix='tmp')
            with self.assertRaises(AssertionError):
                ef.validate_locus_tags(locus_tag_prefix='xxx')

    def test_rename(self):
        for eggnog in eggnogs:
            cleanup()
            EggnogFile(file=eggnog).rename(new_locus_tag_prefix='YOLO', out=TMPFILE, validate=True)
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('YOLO')
            self.assertNotIn(member='tmp', container=content)
            self.assertGreater(a=count, b=400)

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
