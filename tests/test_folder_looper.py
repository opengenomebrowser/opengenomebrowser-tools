import os
from unittest import TestCase
from opengenomebrowser_tools.folder_looper import loop, FolderLooper

ROOT = os.path.dirname(os.path.dirname(__file__))
GENOMIC_DATABASE = f'{ROOT}/database'


class Test(TestCase):
    def test_folder_looper_genomes(self):
        for genome in FolderLooper(GENOMIC_DATABASE).genomes():
            print(genome)

    def test_folder_looper_organisms(self):
        for genome in FolderLooper(GENOMIC_DATABASE).organisms():
            print(genome)

    def test_folder_looper_representatives(self):
        for genome in FolderLooper(GENOMIC_DATABASE).genomes(representatives_only=True):
            print(genome)
