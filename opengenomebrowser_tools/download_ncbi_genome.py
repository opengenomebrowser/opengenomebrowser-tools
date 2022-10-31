import os
import logging
from urllib import request
from tempfile import TemporaryDirectory
from .utils import Entrez, decompress_gz
from .reindex_assembly import reindex_assembly
from .rename_genbank import GenBankFile
from .rename_gff import GffFile
from .genbank_to_fasta import GenBankToFasta


def get_record_id(assembly_name: str) -> str:
    # fetch assembly report
    record = Entrez.read(Entrez.esearch(db='assembly', term=assembly_name, report='full'))
    assert record['Count'] == '1', f'Expected 1 record. Found {record["Count"]}. ' \
                                   f'Ensure a search on https://www.ncbi.nlm.nih.gov/assembly/?term={assembly_name} yields exactlyone result.'

    # extract record id
    ids = record['IdList']
    assert len(ids) == 1, f'Failed to get unique record_id. {ids=}, {assembly_name=}'
    record_id = ids[0]

    logging.info(f'Extracted {record_id=} from {assembly_name=}')
    return record_id


def download_ncbi_file(record_id: str, out: str, file: str):
    assert not os.path.isfile(out), f'Output file already exists! {out=}'

    # fetch summary
    summary = Entrez.read(Entrez.esummary(db='assembly', id=record_id, report='full'))
    document_summary = summary['DocumentSummarySet']['DocumentSummary'][0]

    # extract ftp link
    refseq_url = document_summary['FtpPath_RefSeq']
    genbank_url = document_summary['FtpPath_GenBank']

    assert not (refseq_url == genbank_url == ''), f'Failed to extract FTP path to {file=}; {record_id}'

    # choose refseq if possible
    url = refseq_url if refseq_url != '' else genbank_url

    assert url.startswith('ftp://ftp.ncbi.nlm.nih.gov/genomes'), f'Url does not match: {url=}'

    assembly_label = os.path.basename(url)
    ftp_link = os.path.join(url, assembly_label + file)
    logging.info(f'Extracted {ftp_link=} from {record_id=}')

    with TemporaryDirectory() as tempdir:
        # download
        gz_file = f'{tempdir}/{assembly_label}.{file}.gz'
        request.urlretrieve(ftp_link, gz_file)
        logging.info(f'Downloaded {ftp_link=}')

        # decompress
        decompress_gz(gz=gz_file, out=out)
        logging.info(f'Decompressed {out=}')


def download_ncbi_fna_gbk_gff(assembly_name: str, out_dir: str) -> (str, str, str):
    record_id = get_record_id(assembly_name)

    fna = os.path.join(out_dir, f'{assembly_name}.fna')
    gbk = os.path.join(out_dir, f'{assembly_name}.gbk')
    gff = os.path.join(out_dir, f'{assembly_name}.gff')

    download_ncbi_file(record_id=record_id, out=fna, file=f'_genomic.fna.gz')
    download_ncbi_file(record_id=record_id, out=gbk, file=f'_genomic.gbff.gz')
    download_ncbi_file(record_id=record_id, out=gff, file=f'_genomic.gff.gz')

    return fna, gbk, gff


def rename_ncbi_files(
        raw_fna: str, raw_gbk: str, raw_gff: str,
        out_fna: str, out_gbk: str, out_gff: str, out_ffn: str, out_faa: str,
        new_locus_tag_prefix: str,
        scaffold_prefix: str,
        old_locus_tag_prefix: str = None,
        leading_zeroes: int = None,
        validate: bool = True
):
    if old_locus_tag_prefix is None:
        old_locus_tag_prefix = GenBankFile(raw_gbk).detect_locus_tag_prefix()

    if not new_locus_tag_prefix.endswith('_'):
        logging.warning(f'new_locus_tag_prefix does not begin with _ (underline): {new_locus_tag_prefix}')

    # reindex assembly
    reindex_assembly(file=raw_fna, out=out_fna, prefix=scaffold_prefix, leading_zeroes=leading_zeroes)

    # rename gff
    gff = GffFile(raw_gff)
    gff.rename(new_locus_tag_prefix=new_locus_tag_prefix, old_locus_tag_prefix=old_locus_tag_prefix, out=out_gff,
               validate=validate)

    # rename gbk
    gbk = GenBankFile(raw_gbk)
    gbk.rename(new_locus_tag_prefix=new_locus_tag_prefix, old_locus_tag_prefix=old_locus_tag_prefix, out=out_gbk,
               validate=validate, scf_prefix=scaffold_prefix, scf_leading_zeroes=leading_zeroes)

    # produce ffn
    GenBankToFasta.convert(gbk=out_gbk, out=out_ffn, format='ffn', strict=validate)

    # produce faa
    GenBankToFasta.convert(gbk=out_gbk, out=out_faa, format='faa', strict=validate)


def download_ncbi_genome(
        assembly_name: str,
        out_dir: str,
        new_locus_tag_prefix: str,
        scaffold_prefix: str = None,
        scaffold_leading_zeroes: int = 1,
        validate: bool = True
):
    """
    Download assembly (fna, gbk, gff), rename the locus tags, and generate CDS FASTAs (faa, ffn)

    :param assembly_name: Name of the assembly on NCBI (e.g. 'GCF_005864195.1')
    :param out_dir: Output directory
    :param new_locus_tag_prefix: The desired locus tag prefix
    :param scaffold_prefix: The desired scaffold prefix of the assembly file (fna)
    :param scaffold_leading_zeroes: format scaffold counter with leading zeroes, e.g.: 5 -> >PREFIX_00001
    :param validate: If True, validate the newly created files
    """

    if not new_locus_tag_prefix.endswith('_'):
        logging.warning(f'new_locus_tag_prefix does not begin with _ (underline): {new_locus_tag_prefix}')

    filename_prefix = new_locus_tag_prefix.rstrip('_')
    out_fna = os.path.join(out_dir, f'{filename_prefix}.fna')
    out_gbk = os.path.join(out_dir, f'{filename_prefix}.gbk')
    out_gff = os.path.join(out_dir, f'{filename_prefix}.gff')
    out_ffn = os.path.join(out_dir, f'{filename_prefix}.ffn')
    out_faa = os.path.join(out_dir, f'{filename_prefix}.faa')

    if scaffold_prefix is None:
        scaffold_prefix = f'{filename_prefix}_scf'

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    for file in (out_fna, out_gbk, out_gff, out_ffn, out_faa):
        assert not os.path.isfile(file), f'Output file already exists! {file=}'

    with TemporaryDirectory() as tempdir:
        raw_fna, raw_gbk, raw_gff = download_ncbi_fna_gbk_gff(assembly_name=assembly_name, out_dir=tempdir)

        rename_ncbi_files(
            raw_fna, raw_gbk, raw_gff,
            out_fna, out_gbk, out_gff, out_ffn, out_faa,
            new_locus_tag_prefix=new_locus_tag_prefix,
            scaffold_prefix=scaffold_prefix, leading_zeroes=scaffold_leading_zeroes,
            validate=validate
        )


def main():
    import fire

    fire.Fire(download_ncbi_genome)


if __name__ == '__main__':
    main()
