import os
import logging

from Bio import SeqIO, SeqRecord, SeqFeature
from .utils import GenomeFile, query_int, entrez_organism_to_taxid, date_to_string, datetime, create_replace_function, \
    split_locus_tag
from .genbank_to_fasta import GenBankToFasta


class GenBankFile(GenomeFile):
    def rename(
            self,
            out: str,
            new_locus_tag_prefix: str,
            old_locus_tag_prefix: str = None,
            validate: bool = False,
            scf_prefix: str = None,
            scf_leading_zeroes: int = None,
            update_path: bool = True
    ) -> None:
        old_locus_tag_prefix = self._pre_rename_check(out, new_locus_tag_prefix, old_locus_tag_prefix)

        with open(self.path) as in_f:
            content = in_f.readlines()

        if scf_prefix:
            # change scf prefixes
            if type(scf_leading_zeroes) is int and scf_leading_zeroes > 1:
                format = lambda c: f'VERSION     {scf_prefix}{str(c).zfill(scf_leading_zeroes)}\n'
            else:
                format = lambda c: f'VERSION     {scf_prefix}{c}\n'

            counter = 0
            for i, line in enumerate(content):
                if line.startswith('VERSION'):
                    counter += 1
                    content[i] = format(counter)

        content = ''.join(content)
        # change locus tags
        replace_fn = create_replace_function({
            string.format(prefix=old_locus_tag_prefix): string.format(prefix=new_locus_tag_prefix)
            for string in ['/locus_tag="{prefix}', '/protein_id="extdb:{prefix}', ':{prefix}']
        })

        old_hash = hash(content)
        content = replace_fn(content)
        assert old_hash != hash(content), f'The content of {self.path=} has not changed!'

        assert new_locus_tag_prefix in content, f'Something went wrong: did not replace anything!'

        with open(out, 'w') as out_f:
            out_f.write(content)

        if update_path:
            self.path = out

        if validate:
            self.validate_locus_tags(locus_tag_prefix=new_locus_tag_prefix)

    def create_ffn(self, ffn: str):
        GenBankToFasta.convert(gbk=self.path, out=ffn, format='ffn')

    def create_faa(self, faa: str):
        GenBankToFasta.convert(gbk=self.path, out=faa, format='faa')

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        if locus_tag_prefix is None:
            locus_tag_prefix = self.detect_locus_tag_prefix()

        with open(self.path) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    locus_tag = feature.qualifiers.get('locus_tag')
                    if locus_tag is not None:
                        locus_tag = locus_tag[0]
                        real_locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
                        assert real_locus_tag_prefix == locus_tag_prefix, \
                            f'locus_tag_prefix in {self.path=} does not match. expected: {locus_tag_prefix} reality: {real_locus_tag_prefix}'
                        assert gene_id.isdigit(), f'locus_tag in {self.path=} is malformed. expected: {locus_tag_prefix}_[0-9]+ reality: {locus_tag}'

    def metadata(self) -> (dict, dict):
        organism_data, genome_data = {}, {}

        try:
            organism_data['taxid'] = self.taxid(raise_error=True)
        except Exception:
            pass

        rec, feature = self._get_first_gbk_rec_feature(gbk=self.path)
        assert 'comment' in rec.annotations, f'Could not determine annotation information from {self.path=}'
        comment = rec.annotations['comment']
        if 'prokka' in comment:
            genome_data.update(
                cds_tool='prokka',
                cds_tool_date=date_to_string(datetime.strptime(rec.annotations['date'], '%d-%b-%Y')),
                cds_tool_version=comment.split('prokka ')[-1].split(' ')[0]
            )
        elif 'PGAP' in comment:
            pgap_comment = rec.annotations['structured_comment']['Genome-Annotation-Data']
            genome_data.update(
                cds_tool='PGAP',
                cds_tool_date=date_to_string(datetime.strptime(pgap_comment['Annotation Date'], '%m/%d/%Y %H:%M:%S')),
                cds_tool_version=pgap_comment['Annotation Software revision']
            )
        elif 'Bakta' in comment:
            bakta_comment = rec.annotations['comment']
            bakta_version = bakta_comment.split('Software: ',1)[1].split('\n',1)[0]
            bakta_db_verson = bakta_comment.split('Database: ',1)[1].split('\n',1)[0].replace(', ', '')
            genome_data.update(
                cds_tool='Bakta',
                cds_tool_date=date_to_string(datetime.strptime(rec.annotations['date'], '%d-%b-%Y')),  #08-SEP-2023
                cds_tool_version=f'{bakta_version}:{bakta_db_verson}'
            )
        else:
            raise AssertionError(f'Failed to discover annotation information from {self.path=}')

        # get  BioProject / BioSample metadata
        if hasattr(rec, 'dbxrefs'):
            for dbxref in rec.dbxrefs:
                if dbxref.startswith('BioProject:'):
                    genome_data['bioproject_accession'] = dbxref.removeprefix('BioProject:')
                if dbxref.startswith('BioSample:'):
                    genome_data['biosample_accession'] = dbxref.removeprefix('BioSample:')

        return organism_data, genome_data

    def taxid(self, raise_error=True, sample_name: str = None) -> int:
        rec, feature = self._get_first_gbk_rec_feature(gbk=self.path)

        # PGAP files contain taxid in first feature (db_xref qualifier)
        db_xref = feature.qualifiers.get('db_xref')
        if type(db_xref) is list and len(db_xref) == 1:
            taxid = db_xref[0].removeprefix('taxon:')
            assert taxid.isdigit(), f'Failed to extract taxid from db_xref in {self.path=}'
            return int(taxid)

        # prokka files contain organism name, which can be turned into taxid using Entrez
        organism = feature.qualifiers.get('organism')
        if type(organism) is list and len(organism) == 1:
            return entrez_organism_to_taxid(organism[0])

        # ask for input or raise error
        if not raise_error:
            taxid = query_int(
                question=f'What is the taxid of {sample_name if sample_name else self.path}?',
                error_msg=f'The response must be an integer. Hint: taxid=1427524 is "mixed sample".'
            )
            assert taxid != 0, f'0 is not a valid taxid.'
            return taxid
        else:
            raise AssertionError(f'Failed to extract taxid from {self.path=}')

    def detect_locus_tag_prefix(self) -> str:
        strain, locus_tag_prefix = self.detect_strain_locus_tag_prefix()
        return locus_tag_prefix

    def detect_strain_locus_tag_prefix(self) -> (str, str):
        strain, locus_tag = None, None
        with open(self.path) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    if strain is None:
                        strain = feature.qualifiers.get('strain')
                    if locus_tag is None:
                        locus_tag = feature.qualifiers.get('locus_tag')
                    if strain is not None and locus_tag is not None:
                        break

        assert type(locus_tag) is list, f'Could not read genome from .gbk file! {locus_tag=}'

        if type(strain) is list and len(strain) == 1 and type(strain[0]) is str:
            strain = strain[0]
        else:
            strain = os.environ.get('STRAIN', None)
            if strain is None:
                logging.warning(f'Could not read organism from .gbk file! {strain=}')
                strain = input(f'Could not read organism from .gbk file! Please enter it manually and press enter:')
                logging.warning(f'This organism name was manually chosen: {strain}')

        locus_tag_prefix, gene_id = split_locus_tag(locus_tag[0])
        assert type(locus_tag_prefix) is str and type(strain) is str
        return strain, locus_tag_prefix

    @staticmethod
    def _get_first_gbk_rec_feature(gbk: str) -> (SeqRecord, SeqFeature):
        with open(gbk) as f:
            for rec in SeqIO.parse(f, "genbank"):
                for feature in rec.features:
                    return rec, feature
        raise AssertionError(f'Failed to get rec and feature from {gbk=}')


def rename_genbank(
        file: str, out: str,
        new_locus_tag_prefix: str, old_locus_tag_prefix: str = None, validate: bool = False,
        scf_prefix: str = None, scf_leading_zeroes: int = None,
):
    """
    Change the locus tags in a GenBank file

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check on locus tags
    :param scf_prefix: desired scaffold prefix (optional)
    :param scf_leading_zeroes: format scaffold counter with leading zeroes. e.g.: 5 -> >PREFIX_00001 (optional)
    """

    GenBankFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate,
        scf_prefix=scf_prefix,
        scf_leading_zeroes=scf_leading_zeroes
    )


def main():
    import fire

    fire.Fire(rename_genbank)


if __name__ == '__main__':
    main()
