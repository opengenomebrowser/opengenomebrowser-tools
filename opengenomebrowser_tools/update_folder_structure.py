import os
import json
import logging

from .folder_looper import FolderLooper, FolderGenome
from .rename_eggnog import EggnogFile
from .utils import query_yes_no, get_folder_structure_version


def _get_folder_structure_dir(folder_structure_dir: str = None) -> str:
    if folder_structure_dir is None:
        assert 'FOLDER_STRUCTURE' in os.environ, f'Cannot find the folder_structure. Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE'
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']
    assert os.path.isdir(folder_structure_dir), f'Could not find the folder_structure. Folder does not exist: {folder_structure_dir}'
    return folder_structure_dir


def set_folder_structure_version(new_version: int, folder_structure_dir: str) -> None:
    assert type(folder_structure_dir) is str
    version_file = f'{folder_structure_dir}/version.json'

    with open(version_file) as f:
        version_dict = json.load(f)

    version_dict['folder_structure_version'] = new_version

    with open(version_file, 'w') as f:
        json.dump(version_dict, f, indent=4)

    print()
    print(f'Successfully updated to folder structure version {new_version}!')


def ask(v_from: int, v_to: int, actions: [str], folder_structure_dir: str):
    assert type(folder_structure_dir) is str

    current_version = get_folder_structure_version(folder_structure_dir)
    assert current_version == v_from, \
        f'Cannot proceed: Folder structure version mismatch.\n' \
        f'This script expects version {v_from}, but folder_structure/version.json says version {current_version}.'

    question = f'Upgrade folder structure from version {v_from} to {v_to}:'
    for action in actions:
        question += f'\n - {action}'
    question += '\n\nProceed?'
    if not query_yes_no(question=question, default='yes'):
        exit(1)


def loop_genomes(folder_structure_dir: str, skip_ignored=False, sanity_check=False, representatives_only=False) -> [FolderGenome]:
    for genome in FolderLooper(folder_structure_dir=folder_structure_dir).genomes(
            skip_ignored=skip_ignored,
            sanity_check=sanity_check,
            representatives_only=representatives_only
    ):
        if genome.has_json:
            yield genome


def from_1_to_2(folder_structure_dir: str = None, skip_ignored=False, sanity_check=False, representatives_only=False):
    """ Upgrade OpenGenomeBrowser folder structure. """
    folder_structure_dir = _get_folder_structure_dir(folder_structure_dir)
    v_from = 1
    v_to = 2

    ask(v_from=v_from, v_to=v_to, actions=['add COG to genome.json'], folder_structure_dir=folder_structure_dir)

    for genome in loop_genomes(folder_structure_dir=folder_structure_dir, skip_ignored=skip_ignored, sanity_check=sanity_check,
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

    set_folder_structure_version(new_version=v_to, folder_structure_dir=folder_structure_dir)


def main():
    from fire import Fire

    Fire({
        'get_current_version': get_folder_structure_version,
        '1_to_2': from_1_to_2,
    })


if __name__ == '__main__':
    main()
