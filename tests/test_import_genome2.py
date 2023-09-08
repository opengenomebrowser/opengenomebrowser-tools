import shutil
from unittest import TestCase

import os
import logging
from opengenomebrowser_tools.import_genome2 import import_genome2 as import_genome, ImportSettings2

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
    os.makedirs(ORAGNISMS_DIR, exist_ok=True)
    for entry in os.listdir(ORAGNISMS_DIR):
        shutil.rmtree(os.path.join(ORAGNISMS_DIR, entry))
    for f in TO_DELETE:
        if os.path.isfile(f):
            os.remove(f)


class Test(TestCase):
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
            import_settings=f'{ROOT}/test-data/alternative-config2.json',
            organism='STRAIN', genome='STRAIN.1'
        )

    def test_import_simone_config(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/simone',
            import_settings=f'{ROOT}/test-data/import-config-agroscope2.json',
            organism='FAM24234', genome='FAM24234-i1-2.1'
        )

    def test_import_hatice_config(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/Dog002_mouth-p1-1.1',
            import_settings=f'{ROOT}/test-data/import-config-hatice-staph.json',
            organism='Dog002_mouth', genome='Dog002_mouth-p1-1.1'
        )

    def test_import_hatice_ncbi(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/117',
            import_settings=f'{ROOT}/test-data/import-config-hatice-staph-ncbi.json',
            organism='117', genome='117'
        )

    def test_import_subdir_config(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/subdir',
            import_settings=f'{ROOT}/test-data/import-config-subdir2.json',
            organism='FAM24234', genome='FAM24234-i1-2.1'
        )

    def test_DSM22211(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/DSM22211-i1-1.1',
            import_settings=f'{ROOT}/test-data/import-config-DSM22211-2.json',
            organism='DSM22211', genome='DSM22211-i1-1.1'
        )

    def test_bakta(self):
        import_genome(
            folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/bakta',
            import_settings=f'{ROOT}/test-data/import-settings-bakta.json',
            organism='thomas', genome='thomas-1.1'
        )

    def test_import_ncbi(self):
        import_genome(folder_structure_dir=FOLDER_STRUCTURE, import_dir=f'{ROOT}/test-data/ncbi-convert',
                      organism='FAM3257', genome='FAM3257-NCBI.1', rename=True)

    def setUp(self) -> None:
        clean_up()

    @classmethod
    def tearDownClass(cls) -> None:
        clean_up()  # break here to inspect result
