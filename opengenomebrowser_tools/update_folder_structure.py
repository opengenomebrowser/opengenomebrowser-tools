import logging
import os
from .folder_looper import FolderLooper
from .rename_eggnog import EggnogFile


def add_cog(database_dir: str = None, skip_ignored=False, sanity_check=False, representatives_only=False):
    """
    In Sept 2021, a entry was added to the genome.json files: COG.

    This script adds these entries, either by adding a blank, default value ('CDS': {}) or loading the value from eggnog files.
    """
    if database_dir is None:
        assert 'GENOMIC_DATABASE' in os.environ, f'Cannot find the database. Please set --database_dir or environment variable GENOMIC_DATABASE'
        database_dir = os.environ['GENOMIC_DATABASE']

    for genome in FolderLooper(database_dir=database_dir).genomes(
            skip_ignored=skip_ignored,
            sanity_check=sanity_check,
            representatives_only=representatives_only
    ):
        if not genome.has_json:
            continue

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

    print()
    print('Success!')


def main():
    from fire import Fire

    Fire({
        'add_cog': add_cog,
    })


if __name__ == '__main__':
    main()
