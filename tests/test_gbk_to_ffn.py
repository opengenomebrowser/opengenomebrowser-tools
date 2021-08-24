from unittest import TestCase

import os
from tempfile import NamedTemporaryFile
from opengenomebrowser_tools.gbk_to_ffn import GbkToFfn

ROOT = os.path.dirname(os.path.dirname(__file__))

ORIG_ASSEMBLY = f'{ROOT}/test-data/assembly/ASM2732v1.annotation.nucleotide.1.fasta'
GBK_PGAP_BAD = f'{ROOT}/test-data/pgap-bad/annot.gbk'
GBK_PROKKA_BAD = f'{ROOT}/test-data/prokka-bad/PROKKA_08112021.gbk'

FASTA_CHARS = set('ATGC')


class TestGbkToFfn(TestCase):
    def setUp(self) -> None:
        for file in [GBK_PGAP_BAD, GBK_PROKKA_BAD]:
            assert os.path.isfile(file), f'File does not exist: {file}'

    def test_convert(self):
        with NamedTemporaryFile() as f:
            GbkToFfn.convert(gbk=GBK_PGAP_BAD, ffn=f.name, overwrite=True)
            GbkToFfn.convert(gbk=GBK_PROKKA_BAD, ffn=f.name, overwrite=True)

    def test_long_fasta_generator(self):
        def test(gbk: str):
            identifier_expected = True
            for line in GbkToFfn.long_fasta_generator(gbk=gbk):
                assert line.startswith('>') == identifier_expected, 'lines must alter between identifier and sequence'
                if not identifier_expected:
                    self.assertTrue(len(set(line).difference(FASTA_CHARS)) == 0,
                                    msg=f'fasta sequences may only contain {FASTA_CHARS=}. {set(line)=} {line=}')
                identifier_expected = not identifier_expected
            assert identifier_expected, 'file cannot end with identifier line'

        for gbk in [GBK_PGAP_BAD, GBK_PROKKA_BAD]:
            test(gbk)

    def test_short_fasta_generator(self):
        short_fasta_generator = GbkToFfn.short_fasta_generator(gbk=GBK_PGAP_BAD)
        short_fasta_entry = next(short_fasta_generator).strip()

        self.assertTrue(short_fasta_entry.startswith('>'),
                        msg='fasta must begin with header.')

        trimmed_seq = short_fasta_entry.split('\n')[1:-1]

        for line in trimmed_seq:
            self.assertEquals(len(line), 60,
                              msg=f'short fasta lines must be 60 chars long. {line=}')
            self.assertTrue(len(set(line).difference(FASTA_CHARS)) == 0,
                            msg=f'fasta sequences may only contain {FASTA_CHARS=}. {set(line)=} {line=}')

        last_line = trimmed_seq[-1]
        self.assertTrue(0 < len(last_line) <= 60,
                        msg=f'last fasta line must be between 1 and 60 chars long. {last_line=}')
