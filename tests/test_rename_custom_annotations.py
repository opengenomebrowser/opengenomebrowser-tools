from unittest import TestCase

import os
from opengenomebrowser_tools.rename_custom_annotations import *

ROOT = os.path.dirname(os.path.dirname(__file__))
TMPFILE = '/tmp/renamed_custom_annotations.KG'

custom_files = [
    f'{ROOT}/test-data/prokka-bad/custom_annotations.KG',
]


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class Test(TestCase):

    def test_detect_locus_tag_prefix(self):
        for custom_file in custom_files:
            locus_tag_prefix = CustomAnnotationFile(custom_file).detect_locus_tag_prefix()
            self.assertEqual(locus_tag_prefix, 'tmp')

    def test_validate_locus_tags(self):
        for custom_file in custom_files:
            CustomAnnotationFile(custom_file).validate_locus_tags(locus_tag_prefix=None)
            CustomAnnotationFile(custom_file).validate_locus_tags(locus_tag_prefix='tmp')
            with self.assertRaises(AssertionError):
                CustomAnnotationFile(custom_file).validate_locus_tags(locus_tag_prefix='xxx')

    def test_rename_custom_annotations(self):
        for custom_file in custom_files:
            cleanup()
            CustomAnnotationFile(custom_file).rename(new_locus_tag_prefix='YOLO', out=TMPFILE, validate=True)
            with open(TMPFILE) as f:
                content = f.read()
            count = content.count('YOLO')
            self.assertNotIn(member='tmp', container=content)
            self.assertEqual(count, 3)

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
