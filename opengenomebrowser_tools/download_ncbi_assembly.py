import os
import shutil
import logging
from urllib import request
from tempfile import TemporaryDirectory
from .utils import Entrez, decompress_gz
from .rename_genbank import GenBankFile
from .rename_gff import GffFile
from .genbank_to_fasta import GenBankToFasta


def get_record_id(assembly_name: str) -> str:
    # fetch assembly report
    record = Entrez.read(Entrez.esearch(db='assembly', term=assembly_name, report='full'))

    # extract record id
    ids = record['IdList']
    assert len(ids) == 1, f'Failed to get unique record_id. {ids=}, {assembly_name=}'
    record_id = ids[0]

    logging.info(f'Extracted {record_id=} from {assembly_name=}')
    return record_id


def download_ncbi_assembly_file(record_id: str, out: str, file: str):
    assert not os.path.isfile(out), f'Output file already exists! {out=}'

    # fetch summary
    summary = Entrez.read(Entrez.esummary(db='assembly', id=record_id, report='full'))

    # extract ftp link
    url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
    assert url.startswith('ftp://ftp.ncbi.nlm.nih.gov/genomes')
    assembly_label = os.path.basename(url)
    ftp_link = os.path.join(url, assembly_label + file)
    logging.info(f'Extracted {ftp_link=} from {record_id=}')

    with TemporaryDirectory() as tempdir:
        # download
        gz_file = f'{tempdir}/{assembly_label}.{format}.gz'
        request.urlretrieve(ftp_link, gz_file)
        logging.info(f'Downloaded {ftp_link=}')

        # decompress
        decompress_gz(gz=gz_file, out=out)
        logging.info(f'Decompressed {out=}')


def download_ncbi_assembly(assembly_name: str, out_dir: str) -> (str, str, str):
    record_id = get_record_id(assembly_name)

    fna = os.path.join(out_dir, f'{assembly_name}.fna')
    gbk = os.path.join(out_dir, f'{assembly_name}.gbk')
    gff = os.path.join(out_dir, f'{assembly_name}.gff')

    download_ncbi_assembly_file(record_id=record_id, out=fna, file=f'_genomic.fna.gz')
    download_ncbi_assembly_file(record_id=record_id, out=gbk, file=f'_genomic.gbff.gz')
    download_ncbi_assembly_file(record_id=record_id, out=gff, file=f'_genomic.gff.gz')

    return fna, gbk, gff


def convert_ncbi_assembly(
        raw_fna: str, raw_gbk: str, raw_gff: str,
        out_fna: str, out_gbk: str, out_gff: str, out_ffn: str, out_faa: str,
        new_locus_tag_prefix: str,
        old_locus_tag_prefix: str = None,
        validate: bool = True
):
    # copy assembly
    shutil.copy(src=raw_fna, dst=out_fna)

    if old_locus_tag_prefix is None:
        old_locus_tag_prefix = GenBankFile(raw_gbk).detect_locus_tag_prefix()

    # rename gbk
    gbk = GenBankFile(raw_gbk)
    gbk.rename(new_locus_tag_prefix=new_locus_tag_prefix, old_locus_tag_prefix=old_locus_tag_prefix, out=out_gbk, validate=validate)

    # rename gff
    gff = GffFile(raw_gff)
    gff.rename(new_locus_tag_prefix=new_locus_tag_prefix, old_locus_tag_prefix=old_locus_tag_prefix, out=out_gff, validate=validate)

    # produce ffn
    GenBankToFasta.convert(out_gbk, out=out_ffn, format='ffn', strict=validate)

    # produce faa
    GenBankToFasta.convert(out_gbk, out=out_faa, format='ffn', strict=validate)


def download_and_convert_ncbi_assembly(assembly_name: str, out_dir: str, new_locus_tag_prefix: str, validate: bool = True):
    """
    todo: document this, also in readme

    :param assembly_name:
    :param out_dir:
    :param new_locus_tag_prefix:
    :param validate:
    :return:
    """
    filename_prefix = new_locus_tag_prefix.rstrip('_')
    out_fna = os.path.join(out_dir, f'{filename_prefix}.fna')
    out_gbk = os.path.join(out_dir, f'{filename_prefix}.gbk')
    out_gff = os.path.join(out_dir, f'{filename_prefix}.gff')
    out_ffn = os.path.join(out_dir, f'{filename_prefix}.ffn')
    out_faa = os.path.join(out_dir, f'{filename_prefix}.faa')

    with TemporaryDirectory() as tempdir:
        raw_fna, raw_gbk, raw_gff = download_ncbi_assembly(assembly_name=assembly_name, out_dir=tempdir)

        convert_ncbi_assembly(
            raw_fna, raw_gbk, raw_gff,
            out_fna, out_gbk, out_gff, out_ffn, out_faa,
            new_locus_tag_prefix=new_locus_tag_prefix, validate=validate
        )


def main():
    import fire

    fire.Fire(download_and_convert_ncbi_assembly)


if __name__ == '__main__':
    main()
