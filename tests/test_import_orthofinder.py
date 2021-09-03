import os
from unittest import TestCase
from opengenomebrowser_tools.import_orthofinder import import_orthofinder

ROOT = os.path.dirname(os.path.dirname(__file__))
GENOMIC_DATABASE = f'{ROOT}/database'


class Test(TestCase):
    def test_import_orthofinder_hog(self):
        import_orthofinder(
            database_dir=GENOMIC_DATABASE,
            fasta_dir='/home/thomas/PycharmProjects/opengenomebrowser/database/OrthoFinder/fastas',
            out_annotations=f'{ROOT}/database/orthologs/orthologs.tsv',
            out_descriptions=f'{ROOT}/database/annotation-descriptions/OL.tsv',
            which='hog'
        )

    def test_import_orthofinder_og(self):
        import_orthofinder(
            database_dir=GENOMIC_DATABASE,
            fasta_dir='/home/thomas/PycharmProjects/opengenomebrowser/database/OrthoFinder/fastas',
            out_annotations=f'{ROOT}/database/orthologs/orthologs.tsv',
            out_descriptions=f'{ROOT}/database/annotation-descriptions/OL.tsv',
            which='og'
        )
