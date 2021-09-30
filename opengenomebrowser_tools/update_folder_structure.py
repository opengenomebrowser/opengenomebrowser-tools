import os
import json
import logging

from .folder_looper import FolderLooper, FolderGenome
from .rename_eggnog import EggnogFile
from .utils import query_yes_no, get_folder_structure_version


def _get_database_dir(database_dir: str = None) -> str:
    if database_dir is None:
        assert 'GENOMIC_DATABASE' in os.environ, f'Cannot find the database. Please set --database_dir or environment variable GENOMIC_DATABASE'
        database_dir = os.environ['GENOMIC_DATABASE']
    assert os.path.isdir(database_dir), f'Could not find the database. Folder does not exist: {database_dir}'
    return database_dir


def set_folder_structure_version(new_version: int, database_dir: str) -> None:
    assert type(database_dir) is str
    version_file = f'{database_dir}/version.json'

    with open(version_file) as f:
        version_dict = json.load(f)

    version_dict['folder_structure_version'] = new_version

    with open(version_file, 'w') as f:
        json.dump(version_dict, f, indent=4)

    print()
    print(f'Successfully updated to folder structure version {new_version}!')


def ask(v_from: int, v_to: int, actions: [str], database_dir: str):
    assert type(database_dir) is str

    current_version = get_folder_structure_version(database_dir)
    assert current_version == v_from, \
        f'Cannot proceed: Folder structure version mismatch.\n' \
        f'This script expects version {v_from}, but database/version.json says version {current_version}.'

    question = f'Upgrade folder structure from version {v_from} to {v_to}:'
    for action in actions:
        question += f'\n - {action}'
    question += '\n\nProceed?'
    if not query_yes_no(question=question, default='yes'):
        exit(1)


def loop_genomes(database_dir: str, skip_ignored=False, sanity_check=False, representatives_only=False) -> [FolderGenome]:
    for genome in FolderLooper(database_dir=database_dir).genomes(
            skip_ignored=skip_ignored,
            sanity_check=sanity_check,
            representatives_only=representatives_only
    ):
        if genome.has_json:
            yield genome


def from_1_to_2(database_dir: str = None, skip_ignored=False, sanity_check=False, representatives_only=False):
    """ Upgrade OpenGenomeBrowser folder structure. """
    database_dir = _get_database_dir(database_dir)
    v_from = 1
    v_to = 2

    ask(v_from=v_from, v_to=v_to, actions=['add COG to genome.json'], database_dir=database_dir)

    for genome in loop_genomes(database_dir=database_dir, skip_ignored=skip_ignored, sanity_check=sanity_check,
                               representatives_only=representatives_only):
        genome_json = genome.json
        if 'COG' in genome_json:
            print(f'{genome.identifier}: already has COG in genome.json')
            continue

        COG = {}  # default

        eggnog_files = [f for f in genome_json['custom_annotations'] if f['type'].startswith('eggnog')]
        for file in eggnog_files:
            path = os.path.join(genome.path, file['file'])
            try:
                COG = EggnogFile(file=path).cog_categories()
            except AssertionError as e:
                logging.info(msg=str(e))
                pass

        print(f'{genome.identifier}: adding COG={COG}')
        genome_json['COG'] = COG
        genome.replace_json(genome_json)

    set_folder_structure_version(new_version=v_to, database_dir=database_dir)


def main():
    from fire import Fire

    Fire({
        'get_current_version': get_folder_structure_version,
        '1_to_2': from_1_to_2,
    })


if __name__ == '__main__':
    main()
