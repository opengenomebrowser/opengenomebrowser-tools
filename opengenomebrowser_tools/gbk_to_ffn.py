import io
import os
from Bio import SeqIO


class GbkToFfn:
    @staticmethod
    def convert(gbk, ffn=None, overwrite=False):
        """
        Convert GenBank (gbk) file into nucleotide FASTA (ffn) file.

        :param gbk: path to input GenBank file
        :param ffn: path to output file; if not specified, FASTA will be printed to stdout.
        :param overwrite: True: automatically overwrite output files; False: always ask first (default)
        """
        assert os.path.isfile(gbk)

        if not ffn:
            # Print file to stdout
            for protein in GbkToFfn.short_fasta_generator(gbk):
                print(protein, end='')
            return

        if os.path.isfile(ffn):
            if overwrite or GbkToFfn._query_yes_no(f'Overwrite {ffn}?'):
                os.remove(ffn)
            else:
                print('Aborting.')
                return

        assert not os.path.isfile(ffn)
        with open(ffn, 'w') as results:
            fasta = list(GbkToFfn.long_fasta_generator(gbk))
            fasta = '\n'.join(fasta)
            fasta = io.StringIO(fasta)

            for seq in SeqIO.parse(fasta, "fasta"):
                results.write(seq.format("fasta"))

    @staticmethod
    def short_fasta_generator(gbk):
        """
        Generator that turns gbk-file into standard nucleotide FASTA-format.

        The nucleotide sequence is interrupted by newlines every 50 bp.
        """

        fasta = list(GbkToFfn.long_fasta_generator(gbk))
        fasta = '\n'.join(fasta)
        fasta = io.StringIO(fasta)
        for seq in SeqIO.parse(fasta, "fasta"):
            yield seq.format("fasta")

    @staticmethod
    def long_fasta_generator(gbk):
        """
        Generator that turns gbk-file into long-form nucleotide FASTA-format.

        The nucleotide sequence is not interrupted by newlines.
        """

        with open(gbk) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    if 'locus_tag' in feature.qualifiers and 'product' in feature.qualifiers:
                        gene_identifier = feature.qualifiers['locus_tag'][0]
                        gene_product = feature.qualifiers['product'][0]
                        gene_nucleotides = feature.location.extract(rec).seq
                        yield F">{gene_identifier} {gene_product}"
                        yield F"{gene_nucleotides}"

    @staticmethod
    def _query_yes_no(question: str) -> bool:
        while True:
            reply = input(F'{question} (y/yes or n/no)')
            if reply.lower() in ['y', 'yes']:
                return True
            if reply.lower() in ['n', 'no']:
                return False


def main():
    import fire

    fire.Fire(GbkToFfn.convert)


if __name__ == '__main__':
    main()
