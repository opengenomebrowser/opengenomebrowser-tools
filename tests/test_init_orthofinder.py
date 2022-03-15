from unittest import TestCase
import os
import shutil
import logging
from opengenomebrowser_tools.init_orthofinder import init_orthofinder

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
FOLDER_STRUCTURE = f'{ROOT}/folder_structure'
ORTHOFINDER_DIR = f'{FOLDER_STRUCTURE}/OrthoFinder'

assert os.path.isdir(FOLDER_STRUCTURE)


def cleanup():
    if os.path.isdir(ORTHOFINDER_DIR):
        shutil.rmtree(ORTHOFINDER_DIR)


class Test(TestCase):
    def test_init_orthofinder(self):
        init_orthofinder(folder_structure_dir=FOLDER_STRUCTURE)

    def setUp(self) -> None:
        cleanup()

    @classmethod
    def tearDownClass(cls) -> None:
        cleanup()
