from unittest import TestCase
import os
import shutil
import logging
from opengenomebrowser_tools.init_orthofinder import init_orthofinder

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
GENOMIC_DATABASE = f'{ROOT}/database'
ORTHOFINDER_DIR = f'{GENOMIC_DATABASE}/OrthoFinder'

assert os.path.isdir(GENOMIC_DATABASE)


def cleanup():
    if os.path.isdir(ORTHOFINDER_DIR):
        shutil.rmtree(ORTHOFINDER_DIR)


class Test(TestCase):
    def test_init_orthofinder(self):
        init_orthofinder(database_dir=GENOMIC_DATABASE)

    def setUp(self) -> None:
        cleanup()

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
