import io
import os
import logging
from typing import Optional

from Bio import SeqIO, SeqRecord, SeqFeature


class GenBankToFasta:
    @classmethod
    def convert(cls, gbk, out: str, format: str, strict: bool = True):
        """
        Convert GenBank (gbk) file into protein FASTA (faa) or nucleotide FASTA (ffn)

        :param gbk: path to input GenBank file
        :param out: path to output file
        :param format: faa for protein FASTA, ffn for nucleotide FASTA
        """
        assert not os.path.isfile(out), f'Output file already exists! {out=}'
        assert format in ('faa', 'ffn'), f'Format must be either faa or ffn! {format=}'

        with open(out, 'w') as f:
            for line in cls._short_fasta_generator(gbk=gbk, format=format, strict=strict):
                f.write(line)

    @classmethod
    def _short_fasta_generator(cls, gbk: str, format: str, strict: bool):
        """
        Generator that turns gbk-file into standard nucleotide FASTA-format.

        The nucleotide sequence is interrupted by newlines every 50 bp.
        """

        fasta = list(cls._long_fasta_generator(gbk=gbk, format=format, strict=strict))
        fasta = '\n'.join(fasta)
        fasta = io.StringIO(fasta)
        for seq in SeqIO.parse(fasta, "fasta"):
            yield seq.format("fasta")

    @classmethod
    def _long_fasta_generator(cls, gbk: str, format: str, strict: bool):
        """
        Generator that turns gbk-file into long-form nucleotide FASTA-format.

        The nucleotide sequence is not interrupted by newlines.
        """

        if format == 'faa':
            expected_entries = cls._get_total_proteins(gbk=gbk)
            expected_type = 'proteins'

            def parse_feature(feature: SeqFeature):
                if 'locus_tag' in feature.qualifiers and 'product' in feature.qualifiers and 'translation' in feature.qualifiers:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    gene_product = feature.qualifiers['product'][0]
                    sequence = feature.qualifiers['translation'][0]
                    return locus_tag, gene_product, sequence
                else:
                    return None, None, None

        elif format == 'ffn':
            expected_entries = cls._get_total_genes(gbk=gbk)
            expected_type = 'genes'

            def parse_feature(feature: SeqFeature):
                if 'locus_tag' in feature.qualifiers and 'product' in feature.qualifiers:
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    gene_product = feature.qualifiers['product'][0]
                    sequence = str(feature.location.extract(rec).seq)
                    return locus_tag, gene_product, sequence
                else:
                    return None, None, None

        else:
            raise AssertionError(f'Format must be either faa or ffn! {format=}')

        locus_tags = set()
        with open(gbk) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    locus_tag, gene_product, sequence = parse_feature(feature=feature)
                    if locus_tag is None:
                        continue

                    # check if locus_tag is unique
                    if locus_tag in locus_tags:
                        logging.warning(f'GenBank is strange: {locus_tag} occurs multiple times!')
                    locus_tags.add(locus_tag)

                    yield f'>{locus_tag} {gene_product}\n{sequence}'

        if expected_entries is not None and expected_entries != len(locus_tags):
            msg = f'GenBank is strange: Claims to have {expected_entries} {expected_type}, but only has {len(locus_tags)} unique locus_tags!'
            if strict:
                raise AssertionError(msg)
            else:
                logging.warning(msg)

    @classmethod
    def _get_total_genes(cls, gbk: str) -> Optional[int]:
        rec, feature = cls._get_first_gbk_rec_feature(gbk=gbk)
        try:
            return int(rec.annotations['structured_comment']['Genome-Annotation-Data']['Genes (total)'])
        except Exception:
            return None

    @classmethod
    def _get_total_proteins(cls, gbk: str) -> Optional[int]:
        rec, feature = cls._get_first_gbk_rec_feature(gbk=gbk)
        try:
            return int(rec.annotations['structured_comment']['Genome-Annotation-Data']['CDSs (with protein)'])
        except Exception:
            return None

    @classmethod
    def _get_first_gbk_rec_feature(cls, gbk: str) -> (SeqRecord, SeqFeature):
        with open(gbk) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    return rec, feature
        raise AssertionError(f'Failed to get rec and feature from {gbk=}')


def main():
    import fire

    fire.Fire(GenBankToFasta.convert)


if __name__ == '__main__':
    main()
