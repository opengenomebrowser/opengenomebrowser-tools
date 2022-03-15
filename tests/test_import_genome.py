import shutil
from unittest import TestCase

import os
import logging
from opengenomebrowser_tools.import_genome import import_genome, ImportSettings

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
FOLDER_STRUCTURE = f'{ROOT}/folder_structure'
ORAGNISMS_DIR = f'{FOLDER_STRUCTURE}/organisms'
TO_DELETE = [
    f'{ROOT}/test-data/pgap-bad/annot.ffn',
    f'{ROOT}/test-data/pgap-good/annot.ffn',
]

assert os.path.isdir(FOLDER_STRUCTURE)


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

    def test_import_pgap_good(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/pgap-good')

    def test_import_pgap_bad(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN',
                      genome='STRAIN.1', rename=True)

    def test_import_prokka_good(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/prokka-good')

    def test_import_prokka_bad(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/prokka-bad',
                      organism='STRAIN', genome='STRAIN.1', rename=True)

    def test_import_no_underline(self):
        import_genome(
            import_dir=f'{ROOT}/test-data/no-underline',
            folder_structure_dir=FOLDER_STRUCTURE,
            organism='BLU',
            genome='BLU.1',
            rename=True,
            check_files=True
        )

    def test_import_conflict(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN',
                      genome='STRAIN.1', rename=True)
        with self.assertRaises(AssertionError):
            import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/pgap-bad',
                          organism='STRAIN', genome='STRAIN.1', rename=True)

    def test_import_conflict_solved(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN',
                      genome='STRAIN.1', rename=True)
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/pgap-bad', organism='STRAIN',
                      genome='STRAIN.2', rename=True)

    def test_import_advanced_config(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/prokka-good',
            import_settings=f'{ROOT}/test-data/alternative-config.json',
            organism='STRAIN', genome='STRAIN.1'
        )

    def test_import_ncbi(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/ncbi-convert',
                      organism='FAM3257', genome='FAM3257-NCBI.1', rename=True)

    def setUp(self) -> None:
        clean_up()

    @classmethod
    def tearDownClass(cls) -> None:
        clean_up()  # break here to inspect result
