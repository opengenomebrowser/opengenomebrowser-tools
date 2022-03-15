import os
from unittest import TestCase
from opengenomebrowser_tools.folder_looper import loop, FolderLooper

ROOT = os.path.dirname(os.path.dirname(__file__))
FOLDER_STRUCTURE = f'{ROOT}/folder_structure'


class Test(TestCase):
    def test_folder_looper_genomes(self):
        for genome in FolderLooper(FOLDER_STRUCTURE).genomes():
            print(genome)

    def test_folder_looper_organisms(self):
        for genome in FolderLooper(FOLDER_STRUCTURE).organisms():
            print(genome)

    def test_folder_looper_representatives(self):
        for genome in FolderLooper(FOLDER_STRUCTURE).genomes(representatives_only=True):
            print(genome)
