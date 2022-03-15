import os
from unittest import TestCase
from opengenomebrowser_tools.import_orthofinder import import_orthofinder

ROOT = os.path.dirname(os.path.dirname(__file__))
FOLDER_STRUCTURE = f'{ROOT}/folder_structure'


class Test(TestCase):
    def test_import_orthofinder_hog(self):
        import_orthofinder(
            folder_structure_dir=FOLDER_STRUCTURE,
            fasta_dir='/home/thomas/PycharmProjects/opengenomebrowser/folder_structure/OrthoFinder/fastas',
            out_annotations=f'{ROOT}/folder_structure/orthologs/orthologs.tsv',
            out_descriptions=f'{ROOT}/folder_structure/annotation-descriptions/OL.tsv',
            which='hog'
        )

    def test_import_orthofinder_og(self):
        import_orthofinder(
            folder_structure_dir=FOLDER_STRUCTURE,
            fasta_dir='/home/thomas/PycharmProjects/opengenomebrowser/folder_structure/OrthoFinder/fastas',
            out_annotations=f'{ROOT}/folder_structure/orthologs/orthologs.tsv',
            out_descriptions=f'{ROOT}/folder_structure/annotation-descriptions/OL.tsv',
            which='og'
        )
