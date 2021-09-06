from unittest import TestCase
import os
import shutil
import logging
from opengenomebrowser_tools.init_database import init_database

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
GENOMIC_DATABASE = f'{ROOT}/database'


def cleanup():
    if os.path.isdir(GENOMIC_DATABASE):
        shutil.rmtree(GENOMIC_DATABASE)


class Test(TestCase):
    def test_init_orthofinder(self):
        init_database(database_dir=GENOMIC_DATABASE)

        for dir in [
            GENOMIC_DATABASE,
            f'{GENOMIC_DATABASE}/organisms',
            f'{GENOMIC_DATABASE}/orthologs',
            f'{GENOMIC_DATABASE}/annotation-descriptions',
            f'{GENOMIC_DATABASE}/pathway-maps',
            f'{GENOMIC_DATABASE}/pathway-maps/svg',
        ]:
            self.assertTrue(os.path.isdir(dir))

        for file in [
            f'{GENOMIC_DATABASE}/annotations.json',
            f'{GENOMIC_DATABASE}/annotation-descriptions/SL.tsv',
            f'{GENOMIC_DATABASE}/annotation-descriptions/GO.tsv',
            f'{GENOMIC_DATABASE}/annotation-descriptions/EC.tsv',
            f'{GENOMIC_DATABASE}/annotation-descriptions/KG.tsv',
            f'{GENOMIC_DATABASE}/annotation-descriptions/KR.tsv',
            f'{GENOMIC_DATABASE}/pathway-maps/type_dictionary.json',
        ]:
            self.assertTrue(os.path.isfile(file), msg=f'File does not exist: {file=}')

    def setUp(self) -> None:
        cleanup()

    # @classmethod
    # def tearDownClass(cls) -> None:
    #     cleanup()
