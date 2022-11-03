from unittest import TestCase

import os
import logging
from opengenomebrowser_tools.download_ncbi_genome import get_record_id, download_ncbi_fna_gbk_gff, rename_ncbi_files, \
    download_ncbi_genome

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))
FOLDER_STRUCTURE = f'{ROOT}/folder_structure'
ORAGNISMS_DIR = f'{FOLDER_STRUCTURE}/organisms'

ASSEMBLIES = [('FAM3257', 'GCF_005864195.1', 'FEZ40_RS')]
DOWNLOAD_DIR = f'{ROOT}/test-data/ncbi-download'
CONVERT_DIR = f'{ROOT}/test-data/ncbi-convert'

assert os.path.isdir(FOLDER_STRUCTURE)


def remove_if_exists(file: str):
    if os.path.isfile(file):
        os.remove(file)


def cleanup(files: [str]):
    [remove_if_exists(f) for f in files]


def get_paths(dir: str, prefix: str) -> (str, str, str, str, str, str):
    fna = f'{dir}/{prefix}.fna'
    gbk = f'{dir}/{prefix}.gbk'
    gff = f'{dir}/{prefix}.gff'
    ffn = f'{dir}/{prefix}.ffn'
    faa = f'{dir}/{prefix}.faa'
    return fna, gbk, gff, ffn, faa


class Test(TestCase):
    def test_get_record_id(self):
        self.assertEqual(get_record_id('GCF_005864195.1'), '3205601')
        self.assertEqual(get_record_id('SAMN11653942'), '3205601')
        with self.assertRaises(AssertionError):
            get_record_id('this-is-fake-input')

    def test_download_FAM3257(self):
        cleanup(get_paths(dir=DOWNLOAD_DIR, prefix='GCF_005864195.1'))
        fna, gbk, gff = download_ncbi_fna_gbk_gff(assembly_name='GCF_005864195.1', out_dir=DOWNLOAD_DIR)
        for file in (fna, gbk, gff):
            self.assertTrue(os.path.isfile(file), f'Expected outfile does not exist: {file=}')

    def test_rename_FAM3257(self):
        # run only after previous test!
        raw_fna, raw_gbk, raw_gff, _, _ = get_paths(dir=DOWNLOAD_DIR, prefix='GCF_005864195.1')
        out_fna, out_gbk, out_gff, out_ffn, out_faa = get_paths(dir=CONVERT_DIR, prefix='FAM3257')
        cleanup([out_fna, out_gbk, out_gff, out_ffn, out_faa])
        rename_ncbi_files(
            raw_fna, raw_gbk, raw_gff,
            out_fna, out_gbk, out_gff, out_ffn, out_faa,
            new_locus_tag_prefix='FAM3257_',
            old_locus_tag_prefix='FEZ40_RS',
            scaffold_prefix='FAM3257_scf',
            validate=True
        )

    def test_download_and_rename_FAM3257(self):
        cleanup(get_paths(dir=DOWNLOAD_DIR, prefix='FAM3257'))
        download_ncbi_genome(assembly_name='GCF_005864195.1', new_locus_tag_prefix='FAM3257_', out_dir=DOWNLOAD_DIR)
        for file in get_paths(dir=DOWNLOAD_DIR, prefix='FAM3257'):
            self.assertTrue(os.path.isfile(file), f'Expected outfile does not exist: {file=}')

    def test_download_FAM8105(self):
        cleanup(get_paths(dir=DOWNLOAD_DIR, prefix='CP015496.1'))
        fna, gbk, gff = download_ncbi_fna_gbk_gff(assembly_name='CP015496.1', out_dir=DOWNLOAD_DIR)
        for file in (fna, gbk, gff):
            self.assertTrue(os.path.isfile(file), f'Expected outfile does not exist: {file=}')

    def test_and_rename_biosample(self):
        biosample = 'SAMN11653942'  # FAM3257
        cleanup(get_paths(dir=DOWNLOAD_DIR, prefix=biosample))
        download_ncbi_genome(assembly_name=biosample, new_locus_tag_prefix=f'{biosample}_', out_dir=DOWNLOAD_DIR)
        for file in get_paths(dir=DOWNLOAD_DIR, prefix=biosample):
            self.assertTrue(os.path.isfile(file), f'Expected outfile does not exist: {file=}')

    def test_rename_FAM8105(self):
        raw_fna, raw_gbk, raw_gff, _, _ = get_paths(dir=DOWNLOAD_DIR, prefix='CP015496.1')
        out_fna, out_gbk, out_gff, out_ffn, out_faa = get_paths(dir=CONVERT_DIR, prefix='FAM8105')
        cleanup([out_fna, out_gbk, out_gff, out_ffn, out_faa])
        rename_ncbi_files(
            raw_fna, raw_gbk, raw_gff,
            out_fna, out_gbk, out_gff, out_ffn, out_faa,
            new_locus_tag_prefix='FAM8105_',
            old_locus_tag_prefix='Lh8105_RS',
            scaffold_prefix='FAM8105_scf',
            validate=True
        )

    def test_download_and_rename_AF036485(self):
        # has no record IDs
        with self.assertRaises(AssertionError):
            download_ncbi_genome(assembly_name='AF036485.1', new_locus_tag_prefix='AF036485.2_', out_dir=DOWNLOAD_DIR)

    def test_download_O_1(self):
        cleanup(get_paths(dir=DOWNLOAD_DIR, prefix='O_1_'))
        fna, gbk, gff = download_ncbi_fna_gbk_gff(assembly_name='GCA_900183865.1', out_dir=DOWNLOAD_DIR)
        for file in (fna, gbk, gff):
            self.assertTrue(os.path.isfile(file), f'Expected outfile does not exist: {file=}')

    def test_rename_O_1(self):
        raw_fna, raw_gbk, raw_gff, _, _ = get_paths(dir=DOWNLOAD_DIR, prefix='GCA_900183865.1')
        out_fna, out_gbk, out_gff, out_ffn, out_faa = get_paths(dir=CONVERT_DIR, prefix='GCA_900183865.1')
        cleanup([out_fna, out_gbk, out_gff, out_ffn, out_faa])
        rename_ncbi_files(
            raw_fna, raw_gbk, raw_gff,
            out_fna, out_gbk, out_gff, out_ffn, out_faa,
            new_locus_tag_prefix='O_1_',
            old_locus_tag_prefix='CFS20_RS',
            scaffold_prefix='O_1_scf',
            validate=True
        )

    def test_download_1708(self):
        cleanup(get_paths(dir=DOWNLOAD_DIR, prefix='1708_'))
        fna, gbk, gff = download_ncbi_fna_gbk_gff(assembly_name='GCA_020743785.1', out_dir=DOWNLOAD_DIR)
        for file in (fna, gbk, gff):
            self.assertTrue(os.path.isfile(file), f'Expected outfile does not exist: {file=}')

    def test_download_and_rename_GCA_002994165(self):
        # has no record IDs
        download_ncbi_genome(assembly_name='GCA_002994165.1', new_locus_tag_prefix='33_', out_dir=DOWNLOAD_DIR)
