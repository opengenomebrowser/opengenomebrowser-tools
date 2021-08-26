import shutil
from unittest import TestCase

import os
import logging
from opengenomebrowser_tools.import_genome import runner, ImportSettings

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
GENOMIC_DATABASE = f'{ROOT}/database'
ORAGNISMS_DIR = f'{GENOMIC_DATABASE}/organisms'
TO_DELETE = [
    f'{ROOT}/test-data/pgap-bad/annot.ffn',
    f'{ROOT}/test-data/pgap-good/annot.ffn',
]

assert os.path.isdir(GENOMIC_DATABASE)


def clean_up():
    if os.path.isdir(ORAGNISMS_DIR):
        shutil.rmtree(ORAGNISMS_DIR)
    os.makedirs(ORAGNISMS_DIR)
    for f in TO_DELETE:
        if os.path.isfile(f):
            os.remove(f)


class Test(TestCase):
    def test_import_settings(self):
        import_settings = ImportSettings()
        res = import_settings.get_path(original_path='bla.faa', genome='XXX', organism='YYY')
        self.assertEqual(res, 'XXX.faa')

        res = import_settings.get_path(original_path='bla.yolo', genome='XXX', organism='YYY')
        self.assertEqual(res, 'rest/bla.yolo')

    def test_import_settings_custom(self):
        import_settings = ImportSettings(dict(
            organism_template=None,
            genome_template=None,
            path_transformer={
                '^.*\.abc$': '{organism}.{genome}.{suffix}',
                '^.*\.ignore$': None,
            }
        ))

        res = import_settings.get_path(original_path='bla.abc', genome='XXX', organism='YYY')
        self.assertEqual(res, 'YYY.XXX.abc')

        res = import_settings.get_path(original_path='bla.ignore', genome='XXX', organism='YYY')
        self.assertEqual(res, None)

        with self.assertRaises(AssertionError):
            res = import_settings.get_path(original_path='bla.yolo', genome='XXX', organism='YYY')

    def test_runner_pgap_good(self):
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/pgap-good')

    def test_runner_pgap_bad(self):
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN', genome='STRAIN.1', rename=True)

    def test_runner_prokka_good(self):
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/prokka-good')

    def test_runner_prokka_bad(self):
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/prokka-bad', organism='STRAIN', genome='STRAIN.1', rename=True)

    def test_runner_conflict(self):
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN', genome='STRAIN.1', rename=True)
        with self.assertRaises(AssertionError):
            runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN', genome='STRAIN.1', rename=True)

    def test_runner_conflict_solved(self):
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN', genome='STRAIN.1', rename=True)
        runner(database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN', genome='STRAIN.2', rename=True)

    def test_runner_advanced_config(self):
        runner(
            database_dir=GENOMIC_DATABASE, import_dir=f'{ROOT}/test-data/prokka-good',
            import_settings=f'{ROOT}/test-data/alternative-config.json',
            organism='STRAIN', genome='STRAIN.1'
        )

    def setUp(self) -> None:
        clean_up()

    @classmethod
    def tearDownClass(cls) -> None:
        clean_up()  # break here to inspect result
