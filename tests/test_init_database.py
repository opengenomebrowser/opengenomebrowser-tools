from unittest import TestCase
import os
import shutil
import logging
from opengenomebrowser_tools.init_folder_structure import init_folder_structure

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
FOLDER_STRUCTURE = f'{ROOT}/folder_structure'


def cleanup():
    if os.path.isdir(FOLDER_STRUCTURE):
        shutil.rmtree(FOLDER_STRUCTURE)


class Test(TestCase):
    def test_init_folder_structure(self):
        init_folder_structure(folder_structure_dir=FOLDER_STRUCTURE)

        for dir in [
            FOLDER_STRUCTURE,
            f'{FOLDER_STRUCTURE}/organisms',
            f'{FOLDER_STRUCTURE}/orthologs',
            f'{FOLDER_STRUCTURE}/annotation-descriptions',
            f'{FOLDER_STRUCTURE}/pathway-maps',
            f'{FOLDER_STRUCTURE}/pathway-maps/svg',
        ]:
            self.assertTrue(os.path.isdir(dir))

        for file in [
            f'{FOLDER_STRUCTURE}/annotations.json',
            f'{FOLDER_STRUCTURE}/annotation-descriptions/SL.tsv',
            f'{FOLDER_STRUCTURE}/annotation-descriptions/GO.tsv',
            f'{FOLDER_STRUCTURE}/annotation-descriptions/EC.tsv',
            f'{FOLDER_STRUCTURE}/annotation-descriptions/KG.tsv',
            f'{FOLDER_STRUCTURE}/annotation-descriptions/KR.tsv',
            f'{FOLDER_STRUCTURE}/pathway-maps/type_dictionary.json',
        ]:
            self.assertTrue(os.path.isfile(file), msg=f'File does not exist: {file=}')

    def setUp(self) -> None:
        cleanup()

    # @classmethod
    # def tearDownClass(cls) -> None:
    #     cleanup()
