from unittest import TestCase

import os
import logging
from opengenomebrowser_tools.genbank_to_fasta import GenBankToFasta

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))

GENBANK_FILES = [
    (f'{ROOT}/test-data/pgap-bad/annot.gbk', 'tmp_'),
    (f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gbk', 'tmp_'),
    (f'{ROOT}/test-data/ncbi-download/GCF_005864195.1.gbk', 'FEZ40_RS')
]

FFN_CHARS = set('NATGC')
FAA_CHARS = set('YPMLWDKHEARGTQSNVFCI')

TMPFILE = '/tmp/genbank_to_fasta.fasta'


def cleanup():
    if os.path.isfile(TMPFILE):
        os.remove(TMPFILE)


class TestGenBankToFasta(TestCase):
    def setUp(self) -> None:
        for gbk, locus_tag_prefix in GENBANK_FILES:
            assert os.path.isfile(gbk), f'File does not exist: {gbk}'

    def convert_tester(self, gbk, format):
        cleanup()
        GenBankToFasta.convert(gbk=gbk, out=TMPFILE, format=format, strict=True)
        cleanup()

    def test_convert(self):
        for gbk, locus_tag_prefix in GENBANK_FILES:
            for format in ['faa', 'ffn']:
                self.convert_tester(gbk=gbk, format=format)

    def long_fasta_tester(self, gbk, format, locus_tag_prefix):
        ALLOWED_CHARS = FAA_CHARS if format == 'faa' else FFN_CHARS
        for entry in GenBankToFasta._long_fasta_generator(gbk=gbk, format=format, strict=True):
            self.assertEqual(entry.count('\n'), 1, f'Expect 1 newlines! {entry=}')
            header, sequence = entry.split('\n', 1)
            self.assertTrue(header.startswith(f'>{locus_tag_prefix}'), f'Header does not start with >{locus_tag_prefix}: {entry=}')
            self.assertTrue(set(sequence).issubset(ALLOWED_CHARS),
                            msg=f'FASTA sequences may only contain {ALLOWED_CHARS}. {set(sequence)=} {sequence=}')

    def test_long_fasta_generator(self):
        for gbk, locus_tag_prefix in GENBANK_FILES:
            for format in ['faa', 'ffn']:
                logging.info(f'testing {format=} {gbk=}')
                self.long_fasta_tester(gbk=gbk, format=format, locus_tag_prefix=locus_tag_prefix)

    def short_fasta_tester(self, gbk, format, locus_tag_prefix):
        ALLOWED_CHARS = FAA_CHARS if format == 'faa' else FFN_CHARS
        short_fasta_generator = GenBankToFasta._short_fasta_generator(gbk=gbk, format=format, strict=True)

        for entry in short_fasta_generator:
            self.assertGreater(entry.count('\n'), 0, f'FASTA entry must at least contain 1 newline! {entry=}')
            header, *sequence = entry.split('\n')

            self.assertTrue(header.startswith(f'>{locus_tag_prefix}'), f'Header does not start with >{locus_tag_prefix}: {entry=}')

            for line in sequence[:-3]:
                self.assertEqual(len(line), 60, msg=f'Short FASTA lines must be 60 chars long. {line=}')
                self.assertTrue(set(line).issubset(ALLOWED_CHARS), msg=f'FASTA sequences may only contain {ALLOWED_CHARS}. {set(line)=} {entry=}')

            last_line = sequence[-2]
            self.assertTrue(set(last_line).issubset(ALLOWED_CHARS),
                            msg=f'FASTA sequences may only contain {ALLOWED_CHARS}. {set(last_line)=} {entry=}')
            self.assertTrue(0 < len(last_line) <= 60, msg=f'last FASTA line must be between 1 and 60 chars long. {last_line=}')

            self.assertTrue(sequence[-1] == '', f'FASTA entries must end with a newline. {entry=}')

    def test_short_fasta_generator(self):
        for gbk, locus_tag_prefix in GENBANK_FILES:
            for format in ['faa', 'ffn']:
                logging.info(f'testing {format=} {gbk=}')
                self.short_fasta_tester(gbk=gbk, format=format, locus_tag_prefix=locus_tag_prefix)
