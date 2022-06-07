import logging
import os
import json
import yaml
import shutil
import tempfile
from glob import glob
from textwrap import shorten
from typing import Union
from schema import SchemaError

from . import __folder_structure_version__
from .utils import entrez_organism_to_taxid, GenomeFile, merge_json, get_folder_structure_version, WorkingDirectory
from .rename_genbank import GenBankFile
from .rename_gff import GffFile
from .rename_fasta import FastaFile
from .rename_eggnog import EggnogFile
from .parse_busco import parse_busco
from .rename_custom_annotations import CustomAnnotationFile
from .metadata_schemas import \
    organism_json_schema, genome_json_dummy, genome_json_schema, organism_json_dummy


class ImportException(Exception):
    pass


class ImportSettings2:
    settings: dict[str:str]
    default_settings = {
        'organism_template': {},
        'genome_template': {},
        'import_actions': [
            {'type': 'copy', 'from': '*', 'to': '{original_path}', 'expected': True},
        ],
        'file_finder': {
            'fna': {'glob': '*.fna', 'expected': 1},
            'gbk': {'glob': '*.gbk', 'expected': 1},
            'gff': {'glob': '*.gff', 'expected': 1},
            'faa': {'glob': '*.faa', 'expected': False},
            'sqn': {'glob': '*.sqn', 'expected': False},
            'ffn': {'glob': '*.ffn', 'expected': False},
            'eggnog': {'glob': '*.emapper.annotations', 'expected': False},
            'yaml': {'glob': '*.yaml', 'expected': False},
            'busco': {'glob': '*_busco.txt', 'expected': False},
            'custom_annotations': [
                {'glob': f'*.{anno_type}', 'anno_type': anno_type, 'expected': False}
                for anno_type in ('GC', 'GP', 'EP', 'ED', 'EO', 'EC', 'KG', 'KR', 'GO', 'SL', 'OL')]
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

        self.settings = self.default_settings | settings  # PEP-584: dict union: overwrite defaults with new settings
        self.settings['file_finder'] = self.default_settings['file_finder'] | settings.get('file_finder', {})

        assert set(self.settings.keys()) == set(self.default_settings.keys()), \
            f'OGB_IMPORT_SETTINGS must contain these JSON keys: {set(self.default_settings.keys())}! ' \
            f'reality: {self.settings.keys()}'

    @staticmethod
    def _copy(src: str, dst: str):
        if os.path.exists(dst):
            logging.warning(f'Overwriting: {src} -> {dst}')
        copy_fn = shutil.copy2 if os.path.isfile(src) else shutil.copytree
        os.makedirs(os.path.dirname(dst), exist_ok=True)  # create parent dir if nonexistent
        copy_fn(src=src, dst=dst)

    @classmethod
    def copy(cls, source_dir: str, target_dir: str, genome: str, organism: str, action: dict):
        from_ = action['from']
        to = action['to']
        expected = action.get('expected', True)

        with WorkingDirectory(source_dir):
            files = glob(from_)
            cls.check_expected(files, expected, from_)
            for src in files:
                basename = os.path.basename(src)
                rel_dst = to.format(
                    original_path=src,
                    suffix=basename.rsplit('.', 1)[-1] if '.' in basename else '',
                    genome=genome,
                    organism=organism,
                    assembly=genome.rsplit('.', 1)[0]
                )
                dst = os.path.join(target_dir, rel_dst)
                if os.path.isdir(dst):
                    logging.warning(f'Overwriting directory: {src} >>{action}>> {rel_dst}')
                    shutil.rmtree(dst)
                elif os.path.isfile(dst):
                    logging.warning(f'Overwriting file: {src} >>{action}>> {rel_dst}')
                    os.remove(dst)
                else:
                    logging.info(f'{src} >>{action}>> {rel_dst}')
                cls._copy(src=src, dst=dst)

    def execute_actions(self, source_dir: str, target_dir: str, genome: str, organism: str) -> None:
        for action in self.settings['import_actions']:
            action_type = action['type']
            if action_type == 'copy':
                self.copy(source_dir, target_dir, genome, organism, action)
            else:
                raise AssertionError(f'Could not execute action: type must be "copy". {action=}')

    @staticmethod
    def check_expected(files: [str], expected: Union[None, bool, int], glob_pattern:str):
        if type(expected) is int:
            if len(files) != expected:
                raise ImportException(f'Error: {files=} glob={glob_pattern}\n'
                                      f'Found {len(files)}, not {expected} files!')
        elif expected is True:
            if len(files) == 0:
                raise ImportException(f'Error: Found no files using glob={glob_pattern}!')
        else:
            if expected is not False:
                raise ImportException(f'Error: config is bad. type must be integer or boolean. '
                                      f'{expected=} {type(expected)=}')

    def find_files(self, type_: str, root_dir: str) -> [str]:
        settings = self.settings['file_finder'][type_]
        glob_pattern = settings['glob']
        expected = settings.get('expected', False)

        with WorkingDirectory(root_dir):
            files = glob(glob_pattern)

        logging.info(f'Found {len(files)} files of type={type_} using glob={glob_pattern}')
        self.check_expected(files, expected, glob_pattern)
        return files

    def find_file(self, type_: str, root_dir: str, as_class=None, expected: bool = True) \
            -> Union[str, GenomeFile, None]:
        files = self.find_files(type_, root_dir)

        if len(files) == 1:
            if as_class is None:
                return files[0]
            else:
                with WorkingDirectory(root_dir):
                    return as_class(files[0])
        else:
            if expected:
                raise AssertionError(f'Error: found {len(files)} files of {type_=}: {files=}')
            else:
                f'Found no {type_} files.'
                return None

    def find_custom_annotations(self, root_dir: str):
        annotations = []

        eggnog_file = self.find_file(type_='eggnog', root_dir=root_dir, as_class=EggnogFile, expected=False)
        if eggnog_file:
            annotations.append(eggnog_file)

        with WorkingDirectory(root_dir):
            for custom_annotation in self.settings['file_finder']['custom_annotations']:
                glob_pattern = custom_annotation['glob']
                files = glob(glob_pattern)
                expected = custom_annotation.get('expected', False)
                assert len(files) < 2, f'Found multiple {custom_annotation["anno_type"]}: {files=}'

                if files:
                    annotations.append(CustomAnnotationFile(
                        file=files[0],
                        custom_annotation_type=custom_annotation['anno_type'])
                    )
                if not files and expected:
                    raise ImportException(f'Error: Found no custom-file using glob={glob_pattern}!')
        return annotations


def autodetect_organism_genome(root_dir: str) -> (str, str):
    with WorkingDirectory(root_dir):
        gbks = glob('*.gbk')
        for gbk in gbks:
            try:
                strain, locus_tag_prefix = GenBankFile(file=gbk).detect_strain_locus_tag_prefix()
                organism, genome = strain, locus_tag_prefix.rstrip('_')
                logging.info(f'autodetected from gbk: {organism=} {genome=}')
                return organism, genome
            except Exception:
                pass
    raise AssertionError(f'Failed to automatically detect organism and genome name in {gbks=}. '
                         f'Please specify them manually.')


def rename_all(root_dir: str, gbk: GenBankFile, files: [GenomeFile], new_prefix: str, old_prefix: str = None):
    if not old_prefix:
        old_prefix = gbk.detect_locus_tag_prefix()

    assert new_prefix != old_prefix, \
        f'old and new locus_tag_prefix are the same! {old_prefix=} {new_prefix=}'

    kwargs = dict(new_locus_tag_prefix=new_prefix, old_locus_tag_prefix=old_prefix, update_path=False)

    with tempfile.TemporaryDirectory() as rename_tempdir, WorkingDirectory(root_dir):
        for file in files:
            temp_file = os.path.join(rename_tempdir, 'tempfile')
            file.rename(out=temp_file, **kwargs)
            os.remove(file.path)
            os.rename(src=temp_file, dst=file.path)


def load_yaml_metadata(submol_yaml: str) -> (dict, dict):
    organism_yaml, genome_yaml = {}, {}
    with open(submol_yaml) as f:
        submol_yaml = yaml.safe_load(f)

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


def load_cog_metadata(custom_annotations: [GenomeFile]) -> dict:
    for file in custom_annotations:
        if type(file) is EggnogFile:
            try:
                cog = file.cog_categories()
                return {'COG': cog}
            except AssertionError as e:
                logging.info(f'Failed to extract COG information from {file.path}. {str(e)}')
                pass
    return {}  # not eggnog file


def add_files_to_json(genome_json: dict, files: dict, custom_annotations) -> dict:
    def get(key):
        file = files[key]
        return None if file is None else file.path

    genome_json['cds_tool_faa_file'] = get('faa')
    genome_json['cds_tool_ffn_file'] = get('ffn')
    genome_json['cds_tool_gbk_file'] = get('gbk')
    genome_json['cds_tool_gff_file'] = get('gff')
    genome_json['cds_tool_sqn_file'] = get('sqn')
    genome_json['assembly_fasta_file'] = get('fna')
    genome_json['custom_annotations'] = [
        {'date': ca.date_str(), 'file': ca.path, 'type': ca.custom_annotation_type}
        for ca in custom_annotations
    ]
    return genome_json


def gather_metadata(import_settings: ImportSettings2, root_dir: str, files: [GenomeFile],
                    custom_annotations: [GenomeFile], organism_dir: str, import_dir: str,
                    organism: str, genome: str):
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
    organism_json.update(import_settings.settings['organism_template'])
    genome_json.update(import_settings.settings['genome_template'])

    # add pgap_submol.yaml
    try:
        organism_yaml, genome_yaml = load_yaml_metadata(import_settings.find_file(type_='yaml', root_dir=root_dir))
        organism_json.update(organism_yaml)
        genome_json.update(genome_yaml)
    except AssertionError as e:
        logging.info(f'Failed to load metadata from yaml: {e}')

    # add *.gbk
    organism_gbk, genome_gbk = files['gbk'].metadata()
    organism_json.update(organism_gbk)
    genome_json.update(genome_gbk)

    # add _busco.txt
    try:
        busco_file = import_settings.find_file(type_='busco', root_dir=root_dir)
        genome_json['BUSCO'] = parse_busco(busco_file)
    except AssertionError:
        pass

    # add COG from eggnog
    genome_json.update(load_cog_metadata(custom_annotations))

    # add organism.json from folder structure
    organism_json = merge_json(organism_json, os.path.join(organism_dir, 'organism.json'))

    # add organism.json / genome.json from import_dir
    organism_json = merge_json(organism_json, os.path.join(import_dir, 'organism.json'))
    genome_json = merge_json(genome_json, os.path.join(import_dir, 'genome.json'))

    # add elementary identifiers
    organism_json['name'] = organism
    organism_json['representative'] = genome
    genome_json['identifier'] = genome

    # add files
    genome_json = add_files_to_json(genome_json, files, custom_annotations)

    # validate metadata files
    try:
        organism_json_schema.validate(organism_json)
    except SchemaError as e:
        logging.warning(f'FAILED TO CREATE A VALID organism.json! {str(e)}')
        raise e

    try:
        genome_json_schema.validate(genome_json)
    except SchemaError as e:
        logging.warning(f'FAILED TO CREATE A VALID genome.json! {str(e)}')
        raise e

    return organism_json, genome_json


def check_files_(locus_tag_prefix, files: dict, custom_annotations: [GenomeFile]) -> None:
    files['gbk'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    files['gff'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    files['faa'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    files['ffn'].validate_locus_tags(locus_tag_prefix=locus_tag_prefix)
    for ca in custom_annotations:
        ca.validate_locus_tags(locus_tag_prefix=locus_tag_prefix)


def import_genome2(
        import_dir: str,
        folder_structure_dir: str = None,
        organism: str = None,
        genome: str = None,
        rename: bool = False,
        check_files: bool = True,
        import_settings: str = None
):
    """
    Easily import files into OpenGenomeBrowser folder structure.

    :param import_dir: Folder with files to import. Required: [.fna, .faa, .gbk, .gff] Optional: [.ffn, .sqn, custom-annotation-files]
    :param folder_structure_dir: Path to the root of the OpenGenomeBrowser folder structure. (Must contain 'organisms' folder.)
    :param organism: Name of the organism.
    :param genome: Identifier of the genome. Must start with organism. May be identical to organism.
    :param rename: Locus tag prefixes must match the genome identifier. If this is not the case, this script can automatically rename relevant files.
    :param check_files: If true, check if locus tag prefixes match genome identifier.
    :param import_settings: Path to import settings file. Alternatively, set the environment variable OGB_IMPORT_SETTINGS.
    """
    import_dir = os.path.abspath(import_dir)

    if folder_structure_dir is None:
        assert 'FOLDER_STRUCTURE' in os.environ, \
            f'Cannot find the folder_structure. ' \
            f'Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE'
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']

    folder_structure_dir = os.path.abspath(folder_structure_dir)

    organisms_dir = f'{folder_structure_dir}/organisms'
    assert os.path.isdir(organisms_dir), f'Cannot import files: {organisms_dir=} does not exist.'

    current_folder_structure_version = get_folder_structure_version(folder_structure_dir)
    assert current_folder_structure_version == __folder_structure_version__, \
        f'Before importing any genomes, the folder structure needs to be updated to match OpenGenomeBrowser Tools.\n' \
        f'Current version: {current_folder_structure_version}, expected: {__folder_structure_version__}\n' \
        f'Use the script update_folder_structure perform the upgrade!'

    assert os.path.isdir(import_dir), f'Cannot import files: {import_dir=} does not exist.'

    import_settings = ImportSettings2(import_settings)

    if organism is None or genome is None:
        _organism, _genome = autodetect_organism_genome(import_dir)
        if organism is None:
            organism = _organism
        if genome is None:
            genome = _genome

    organism_dir = os.path.join(organisms_dir, organism)
    genome_dir = os.path.join(organism_dir, 'genomes', genome)
    assert not os.path.exists(genome_dir), f'Could not import {organism}:{genome}: {genome_dir=} already exists!'

    work_dir = tempfile.TemporaryDirectory()
    _work_dir = os.getcwd()
    os.chdir(work_dir.name)

    import_settings.execute_actions(import_dir, work_dir.name, genome, organism)

    fna: FastaFile = import_settings.find_file('fna', root_dir=work_dir.name, as_class=FastaFile)  # assembly
    gbk: GenBankFile = import_settings.find_file('gbk', root_dir=work_dir.name, as_class=GenBankFile)  # genbank

    ffn = import_settings.find_file('ffn', root_dir=work_dir.name, as_class=FastaFile,
                                    expected=False)  # nucleic acid sequences
    if ffn is None:
        logging.info(f'Failed to auto-detect ffn.')
        ffn = gbk.path[:-4] + '.ffn'

        gbk.create_ffn(ffn=f'{work_dir.name}/{ffn}')
        ffn = FastaFile(ffn)  # nucleic acid sequences

    faa = import_settings.find_file('faa', root_dir=work_dir.name, as_class=FastaFile, expected=False)  # protein
    if faa is None:
        logging.info(f'Failed to auto-detect faa.')
        faa = gbk.path[:-4] + '.faa'
        gbk.create_faa(faa=f'{work_dir.name}/{faa}')
        faa = FastaFile(faa)  # protein

    gff: GffFile = import_settings.find_file('gff', root_dir=work_dir.name, as_class=GffFile)  # general feature format
    sqn: GenomeFile = import_settings.find_file('sqn', root_dir=work_dir.name,
                                                as_class=GenomeFile, expected=False)  # general feature format

    files = dict(fna=fna, gbk=gbk, ffn=ffn, faa=faa, gff=gff, sqn=sqn)

    custom_annotations = import_settings.find_custom_annotations(
        work_dir.name)  # custom annotation files / eggnog files

    if rename:
        rename_all(
            root_dir=work_dir.name, gbk=gbk,
            files=[gbk, gff, faa, ffn, *custom_annotations],
            new_prefix=f'{genome}_'
        )

    organism_json, genome_json = gather_metadata(import_settings, root_dir=work_dir.name, files=files,
                                                 custom_annotations=custom_annotations,
                                                 organism_dir=organism_dir, import_dir=import_dir, organism=organism,
                                                 genome=genome)

    if check_files:
        check_files_(locus_tag_prefix=f'{genome}_', files=files, custom_annotations=custom_annotations)

    # final movement
    os.makedirs(os.path.dirname(genome_dir), exist_ok=True)
    shutil.copytree(src=work_dir.name, dst=genome_dir)

    os.chdir(_work_dir)
    work_dir.cleanup()

    with open(os.path.join(organism_dir, 'organism.json'), 'w') as f:
        json.dump(organism_json, f, indent=4)
    with open(os.path.join(genome_dir, 'genome.json'), 'w') as f:
        json.dump(genome_json, f, indent=4)


def main():
    import fire

    fire.Fire(import_genome2)


if __name__ == '__main__':
    main()
