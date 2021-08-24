import logging
import os
import re
import json
import yaml
import shutil
import tempfile
from textwrap import shorten
from typing import Union, Optional
from schema import SchemaError

from .utils import entrez_organism_to_taxid, GenomeFile
from .rename_genbank import GenBankFile
from .rename_gff import GffFile
from .rename_fasta import FastaFile
from .rename_eggnog import EggnogFile
from .parse_busco import parse_busco
from .rename_custom_annotations import CustomAnnotationFile
from .metadata_schemas import \
    organism_json_schema, genome_json_dummy, genome_json_schema, organism_json_dummy


class ImportSettings:
    settings: dict[re.Pattern:str]
    default_settings = {
        'organism_template': {},
        'genome_template': {},
        'path_transformer': {
            '.*\.fna': '{genome}.{suffix}',
            '.*\.faa': '{genome}.{suffix}',
            '.*\.gbk': '{genome}.{suffix}',
            '.*\.gff': '{genome}.{suffix}',
            '.*\.sqn': '{genome}.{suffix}',
            '.*\.ffn': '{genome}.{suffix}',
            '.*\.emapper.annotations': '{genome}.eggnog',
            '.*\.[A-Z]{2}': '{genome}.{suffix}',
            'genome.md': 'genome.md',
            'organism.md': '../../organism.md',
            'genome.json': None,  # do not copy this file
            'organism.json': None,  # do not copy this file
            '.*': 'rest/{original_path}',
        }
    }

    def __init__(self, settings: Union[dict, str] = None):
        if settings is None:
            settings = os.environ.get('OGB_IMPORT_SETTINGS', None)

        if type(settings) is dict:
            pass
        elif type(settings) is str:
            with open(settings) as f:
                settings = json.load(f)
        else:
            settings = {}

        self.settings = self.default_settings | settings  # PEP-584: union dicts

        assert set(self.settings.keys()) == {'organism_template', 'genome_template', 'path_transformer'}, \
            f'OGB_IMPORT_SETTINGS must contain these JSON keys: organism_template, genome_template and path_transformer! {self.settings.keys()=}'

        self.organism_template = self.settings['organism_template']
        self.genome_template = self.settings['genome_template']
        self.path_transformer: {re.Pattern: str} = {re.compile(pattern=p): n for p, n in self.settings['path_transformer'].items()}

    def get_path(self, original_path: str, genome: str, organism: str) -> Optional[str]:
        pattern: re.Pattern
        for pattern, new_path in self.path_transformer.items():
            if pattern.fullmatch(string=original_path) is not None:
                if new_path is None:
                    return None  # do not copy this file
                return new_path.format(
                    original_path=original_path,
                    suffix=original_path.rsplit('.', 1)[-1],
                    genome=genome,
                    organism=organism,
                )
        raise AssertionError(f'File {original_path} does not match any regex specified in import_settings: {self.settings}')


class OgbImporter:
    genome_json: dict = None
    organism_json: dict = None

    def __init__(self, database_dir: str, import_dir: str, organism: str = None, genome: str = None, import_settings: str = None):
        """
        Easily import files into the OpenGenomeBrowser folder structure.

        :param import_dir: Folder with files to import. Required: [.fna, .faa, .gbk, .gff] Optional: [.ffn, .sqn, custom-annotation-files]
        :param database_dir: Path to the root of the OpenGenomeBrowser folder structure. (Must contain 'organisms' folder.)
        :param organism: Name of the organism.
        :param genome: Identifier of the genome. Must start with organism. May be identical to organism.
        :param rename: Locus tag prefixes must match the genome identifier. If this is not the case, this script can automatically rename relevant files.
        :param import_settings: Path to import settings file. Alternatively, set the environment variable OGB_IMPORT_SETTINGS.
        """
        assert database_dir is not None and os.path.isdir(database_dir), f'Cannot import PGAP files: {database_dir=} does not exist.'
        assert os.path.isdir(import_dir), f'Cannot import PGAP files: {import_dir=} does not exist.'

        self.import_settings = ImportSettings(import_settings)

        self.tempdir = tempfile.TemporaryDirectory()

        self.import_dir = import_dir
        self.database_dir = database_dir
        self.organisms_dir = f'{database_dir}/organisms'
        assert os.path.isdir(self.organisms_dir), f'{database_dir=} does not point to a directory that contains an organisms-folder!'

        files = os.listdir(import_dir)
        self.fna = self.find_file(files, '.fna', FastaFile)  # assembly
        self.faa = self.find_file(files, '.faa', FastaFile)  # protein
        self.gbk = self.find_file(files, '.gbk', GenBankFile)  # genbank
        self.gff = self.find_file(files, '.gff', GffFile)  # general feature format
        self.sqn = self.find_file(files, '.sqn', GenomeFile, raise_error=False)  # GenBank submission file
        self.ffn = self.find_or_create_ffn(files=files)  # nucleic acid sequences
        self.custom_annotations = self.find_custom_annotations(files)  # custom annotation files / eggnog files

        self.genome_md = self.get_genome_file(files, 'genome.md', GenomeFile, raise_error=False)
        self.organism_md = self.get_genome_file(files, 'organism.md', GenomeFile, raise_error=False)

        self.rest_files = files

        auto_organism, auto_genome = self.detect_organism_genome()
        self.organism = organism if organism else auto_organism
        self.genome = genome if genome else auto_genome
        assert self.genome.startswith(self.organism), f'The genome identifier must start with the organism name! {self.genome=} {self.organism=}'

        self.target_dir = f'{self.organisms_dir}/{self.organism}/genomes/{self.genome}'

        assert not os.path.isdir(self.target_dir), f'Could not import {self.organism}:{self.genome}: {self.target_dir=} already exists!'

    def __repr__(self):
        return f'{self.organism}::{self.genome}'

    def _append_new_path_to_file(self, file: GenomeFile) -> str:
        new_path = self.get_new_path(original_path=file.original_path)
        file.new_path = new_path
        file.target_path = os.path.join(self.target_dir, file.new_path)
        return new_path

    def get_new_path(self, original_path: str) -> str:
        return self.import_settings.get_path(original_path=original_path, genome=self.genome, organism=self.organism)

    def add_files_to_json(self, genome_json: dict) -> dict:
        genome_json['cds_tool_faa_file'] = self._append_new_path_to_file(self.faa)
        genome_json['cds_tool_ffn_file'] = self._append_new_path_to_file(self.ffn)
        genome_json['cds_tool_gbk_file'] = self._append_new_path_to_file(self.gbk)
        genome_json['cds_tool_gff_file'] = self._append_new_path_to_file(self.gff)
        genome_json['cds_tool_sqn_file'] = self._append_new_path_to_file(self.sqn)
        genome_json['assembly_fasta_file'] = self._append_new_path_to_file(self.fna)
        genome_json['custom_annotations'] = [
            {'date': ca.date_str(), 'file': self._append_new_path_to_file(ca), 'type': ca.custom_annotation_type}
            for ca in self.custom_annotations
        ]
        return genome_json

    def load_yaml_metadata(self) -> (dict, dict):
        organism_yaml, genome_yaml = {}, {}
        try:
            submol_yaml = self.find_file(files=self.rest_files, suffix='.yaml', file_class=GenomeFile, remove=False)
        except FileNotFoundError:
            return organism_yaml, genome_yaml
        with open(submol_yaml.path) as f:
            submol_yaml = yaml.load(f)

        if 'organism' in submol_yaml and 'genus_species' in submol_yaml['organism']:
            organism_yaml['taxid'] = entrez_organism_to_taxid(submol_yaml['organism']['genus_species'])

        if 'biosample' in submol_yaml:
            genome_yaml['biosample_accession'] = submol_yaml['biosample']
        if 'bioproject' in submol_yaml:
            genome_yaml['bioproject_accession'] = submol_yaml['bioproject']

        if 'publications' in submol_yaml and len(submol_yaml['publications']):
            genome_yaml['literature_references'] = [{
                'url': f"https://pubmed.ncbi.nlm.nih.gov/{p['publication']['pmid']}/",
                'name': shorten(p['publication']['title'], width=30)}
                for p in submol_yaml['publications']
            ]

        return organism_yaml, genome_yaml

    def load_busco_metadata(self) -> dict:
        try:
            busco_file = self.find_file(files=self.rest_files, suffix='_busco.txt', file_class=GenomeFile, remove=False)
        except FileNotFoundError:
            return {}

        return dict(BUSCO=parse_busco(busco_file.path))

    def gather_metadata(self):
        '''
        Load metadata from:
          - pgap_submol.yaml
          - *.gbk
          - *_busco.txt
          - organism.json and genome.json
        :return:
        '''

        # start with dummy jsons
        organism_json = organism_json_dummy.copy()
        genome_json = genome_json_dummy.copy()

        # add import_settings
        organism_json.update(self.import_settings.organism_template)
        genome_json.update(self.import_settings.genome_template)

        # add pgap_submol.yaml
        organism_yaml, genome_yaml = self.load_yaml_metadata()
        organism_json.update(organism_yaml)
        genome_json.update(genome_yaml)

        # add *.gbk
        organism_gbk, genome_gbk = self.gbk.metadata()
        organism_json.update(organism_gbk)
        genome_json.update(genome_gbk)

        # add _busco.txt
        genome_json.update(self.load_busco_metadata())

        # add genome.json / organism.json
        try:
            with open(self.get_file(files=self.rest_files, file='genome.json')) as f:
                genome_json.update(json.load(f))
        except FileNotFoundError:
            pass

        try:
            with open(self.get_file(files=self.rest_files, file='organism.json')) as f:
                organism_json.update(json.load(f))
        except FileNotFoundError:
            pass

        # add elementary identifiers
        organism_json['name'] = self.organism
        organism_json['representative'] = self.genome
        genome_json['identifier'] = self.genome

        # add files
        genome_json = self.add_files_to_json(genome_json)

        try:
            organism_json_schema.validate(organism_json)
        except SchemaError as e:
            logging.warning(f'{self}: FAILED TO CREATE A VALID organism.json!')
            raise e

        try:
            genome_json_schema.validate(genome_json)
        except SchemaError as e:
            logging.warning(f'{self}: FAILED TO CREATE A VALID genome.json!')
            raise e

        self.organism_json = organism_json
        self.genome_json = genome_json

    def detect_organism_genome(self) -> (str, str):
        strain, locus_tag_prefix = self.gbk.detect_strain_locus_tag_prefix()
        organism, genome = strain, locus_tag_prefix
        logging.info(f'autodetected from gbk: {organism=} {genome=}')
        return organism, genome

    def find_or_create_ffn(self, files: [str]) -> GenomeFile:
        try:
            return self.find_file(files, '.ffn', FastaFile)
        except FileNotFoundError:
            self.ffn = self.gbk.path[:-3] + 'ffn'
            logging.info('Creating .ffn based on .gbk...')
            self.gbk.create_ffn(ffn=self.ffn)
            return self.find_file(os.listdir(self.import_dir), '.ffn', FastaFile)

    def check_files(self) -> None:
        self.gbk.validate_locus_tags(locus_tag_prefix=self.genome)
        self.gff.validate_locus_tags(locus_tag_prefix=self.genome)
        self.faa.validate_locus_tags(locus_tag_prefix=self.genome)
        self.ffn.validate_locus_tags(locus_tag_prefix=self.genome)
        for ca in self.custom_annotations:
            ca.validate_locus_tags(locus_tag_prefix=self.genome)

    def _get_temp(self, file: str) -> str:
        return os.path.join(self.tempdir.name, os.path.basename(file))

    def rename_all(self, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None):
        if not old_locus_tag_prefix:
            old_locus_tag_prefix = self.gbk.detect_locus_tag_prefix()
        assert new_locus_tag_prefix != old_locus_tag_prefix, \
            f'old and new locus_tag_prefix are the same! {old_locus_tag_prefix}, {new_locus_tag_prefix=}'

        kwargs = dict(new_locus_tag_prefix=new_locus_tag_prefix, old_locus_tag_prefix=old_locus_tag_prefix)

        self.gbk.rename(out=self._get_temp(self.gbk.path), **kwargs)
        self.gff.rename(out=self._get_temp(self.gff.path), **kwargs)
        self.faa.rename(out=self._get_temp(self.faa.path), **kwargs)
        self.ffn.rename(out=self._get_temp(self.ffn.path), **kwargs)
        for custom_annotation in self.custom_annotations:
            custom_annotation.rename(out=self._get_temp(custom_annotation.path), **kwargs)

    def perform_import(self):
        assert self.genome_json is not None and self.organism_json is not None, \
            f'Cannot perform import yet. Metadata jsons are missing. ' \
            f'(Run gather_metadata or set OgbImporter.genome_json and OgbImporter.organism_json manually)'

        def copy(src: str, dst: str):
            copy_fn = shutil.copy2 if os.path.isfile(src) else shutil.copytree
            os.makedirs(os.path.dirname(dst), exist_ok=True)  # create parent dir if nonexistent
            copy_fn(src=src, dst=dst)

        os.makedirs(self.target_dir)

        for file in [self.fna, self.faa, self.gbk, self.gff, self.sqn, self.ffn, *self.custom_annotations]:
            if file:
                logging.info(f'{file.original_path} >>K>> {file.new_path}')
                copy(src=file.path, dst=file.target_path)

        for original_path in self.rest_files:
            path = os.path.join(self.import_dir, original_path)
            new_path = self.get_new_path(original_path=original_path)
            if new_path is None:
                logging.info(f'not copying: {original_path}')
                continue
            target_path = os.path.join(self.target_dir, new_path)
            logging.info(f'{original_path} >>K>> {new_path}')

            copy(src=path, dst=target_path)

        with open(os.path.join(self.target_dir, 'genome.json'), 'w') as f:
            json.dump(self.genome_json, f, indent=4)

        with open(os.path.join(self.target_dir, '../../organism.json'), 'w') as f:
            json.dump(self.organism_json, f, indent=4)

    def get_file(self, files: [str], file: str, raise_error: bool = True, alternative=None, rel_path=False) -> str:
        if file in files:
            return file if rel_path else os.path.join(self.import_dir, file)
        elif raise_error:
            raise FileNotFoundError(f'Could not find {file}! {files=}')
        else:
            return alternative

    def get_genome_file(self, files: [str], file: str, file_class: type, raise_error: bool = True, alternative=None) -> GenomeFile:
        try:
            file = self.get_file(files=files, file=file, raise_error=True, rel_path=True)
        except FileNotFoundError as e:
            if raise_error:
                raise e
            else:
                return alternative
        genome_file = file_class(os.path.join(self.import_dir, file))
        genome_file.original_path = file
        return genome_file

    def find_file(self, files: [str], suffix: str, file_class: type,
                  raise_error: bool = True, alternative=None, remove: bool = True) -> Union[
        GenomeFile, FastaFile, GffFile, GenBankFile, EggnogFile, CustomAnnotationFile]:
        matches = [f for f in files if f.endswith(suffix)]
        if len(matches) == 0:
            if raise_error:
                raise FileNotFoundError(f'Could not find a {suffix} file! {files=}')
            else:
                return alternative
        if len(matches) > 1:
            raise FileExistsError(f'Found multiple {suffix} files! {matches=}')
        file = matches[0]
        if remove:
            files.remove(file)

        genome_file = file_class(os.path.join(self.import_dir, file))
        genome_file.original_path = file
        return genome_file

    def find_custom_annotations(self, files: [str]) -> [Union[GenomeFile, FastaFile, GffFile, GenBankFile, EggnogFile, CustomAnnotationFile]]:
        custom_annotations = [self.get_genome_file(files, f, CustomAnnotationFile)
                              for f in files
                              if re.fullmatch(pattern=r'.*\.[A-Z]{2}', string=f) is not None]
        eggnog = self.find_file(files, '.emapper.annotations', EggnogFile, raise_error=False)  # GenBank submission file
        if eggnog:
            eggnog: EggnogFile
            logging.info(f'Found eggnog file: {eggnog}')
            custom_annotations.append(eggnog)
        return custom_annotations

    def __del__(self):
        if hasattr(self, 'tempdir'):
            logging.info(f'Deleting {self.tempdir.name}')
            try:
                shutil.rmtree(self.tempdir.name)
            except Exception as e:
                logging.info(f'Could not delete {self.tempdir.name} {str(e)=}')


def runner(
        import_dir: str,
        database_dir: str = None,
        organism: str = None,
        genome: str = None,
        rename: bool = False,
        check_files: bool = True,
        import_settings: str = None
):
    """
    Easily import files into OpenGenomeBrowser folder structure.

    :param import_dir: Folder with files to import. Required: [.fna, .faa, .gbk, .gff] Optional: [.ffn, .sqn, custom-annotation-files]
    :param database_dir: Path to the root of the OpenGenomeBrowser folder structure. (Must contain 'organisms' folder.)
    :param organism: Name of the organism.
    :param genome: Identifier of the genome. Must start with organism. May be identical to organism.
    :param rename: Locus tag prefixes must match the genome identifier. If this is not the case, this script can automatically rename relevant files.
    :param check_files: If true, check if locus tag prefixes match genome identifier.
    :param import_settings: Path to import settings file. Alternatively, set the environment variable OGB_IMPORT_SETTINGS.
    """
    if database_dir is None:
        database_dir = os.environ.get('GENOMIC_DATABASE')

    ogb_importer = OgbImporter(database_dir=database_dir, import_dir=import_dir, organism=organism, genome=genome, import_settings=import_settings)

    if rename:
        ogb_importer.rename_all(new_locus_tag_prefix=ogb_importer.genome)

    ogb_importer.gather_metadata()

    if check_files:
        ogb_importer.check_files()

    ogb_importer.perform_import()


def main():
    import fire

    fire.Fire(runner)


if __name__ == '__main__':
    main()
