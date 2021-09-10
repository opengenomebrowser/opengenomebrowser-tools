import os
import shutil
from urllib import request
from .utils import PACKAGE_ROOT


def download_go_data(out: str) -> None:
    source_url = 'http://purl.obolibrary.org/obo/go.obo'

    print(f'Converting {source_url} -> {out}')

    def go_generator(io) -> [str]:
        go_entry = []

        line = io.readline()
        while line:
            if line == b'[Term]\n':
                yield go_entry
                go_entry.clear()

            go_entry.append(line.decode('utf-8'))
            line = io.readline()

        yield go_entry

    def get_name(entry: list) -> str:
        for line in entry:
            if line.startswith('name: '):
                return line.rstrip()[6:]
        raise TypeError(F'The go.obo file seems to have a wrong format! broken entry: {entry}')

    def get_go(entry: list) -> str:
        entry = entry[1]
        assert entry.startswith('id: GO:') and len(entry) == 15, f'Bad entry in go.obo: {entry}, len={len(entry)}'
        assert entry[7:14].isnumeric()
        return entry[4:14]

    with request.urlopen(source_url) as source_handle, open(out, 'w') as target_handle:
        gos = go_generator(io=source_handle)

        # skip first entry
        file_head = next(gos)
        assert not file_head[0].startswith('[Term]'), F'The go.obo file seems to have a wrong format! file_head looks wrong: {file_head}'

        # save regular entries to file
        for entry in gos:
            target_handle.write(F'{get_go(entry)}\t{get_name(entry)}\n')


def download_kegg_data(src: str, out: str, remove_prefix: str = '', add_prefix: str = '') -> None:
    source_url = f'http://rest.kegg.jp/list/{src}'

    print(f'Converting {source_url} -> {out}')

    with request.urlopen(source_url) as source_handle, open(out, 'w') as target_handle:
        for line in source_handle:
            target_handle.write(f'{add_prefix}{line.decode("utf-8").removeprefix(remove_prefix)}')


def download_sl_data(out: str) -> None:
    source_url = 'https://www.uniprot.org/locations/?query=*&format=tab&force=true&columns=id'

    print(f'Converting {source_url} -> {out}')

    error_msg = 'UniProt must have changed its format. Please contact the developer. error={error}'

    with request.urlopen(source_url) as source_handle, open(out, 'w') as target_handle:
        first_line = source_handle.readline().decode('utf-8')
        assert first_line == 'Subcellular location ID\tDescription\tCategory\tAlias\n', error_msg.format(error=first_line)

        for line in source_handle:
            line = line.decode('utf-8').strip().split('\t')
            assert len(line) == 4, error_msg.format(error=f'{len(line)=}; {line=}')
            sl, description, type, location = line

            target_handle.write(f'{sl}\t{location} ({description})\n')


def init_database(database_dir: str = None) -> None:
    """
    Creates a basic OpenGenomeBrowser folders structure.

    Result:
        database
        ├── organisms
        ├── annotations.json
        ├── annotation-descriptions
        │   ├── SL.tsv
        │   ├── KO.tsv
        │   ├── KR.tsv
        │   ├── EC.tsv
        │   └── GO.tsv
        ├── orthologs
        └── pathway-maps
            ├── type_dictionary.json
            └── svg


    :param database_dir: Path to the root of the OpenGenomeBrowser folder structure. (Will contain 'organisms' folder.)
    """
    if database_dir is None:
        database_dir = os.environ.get('GENOMIC_DATABASE')

    assert os.path.isdir(os.path.dirname(database_dir)), f'Parent dir of {database_dir=} does not exist!'
    assert not os.path.exists(database_dir), f'Error: {database_dir=} already exist!'

    # make main dir
    os.makedirs(database_dir)

    # make organisms dir (empty)
    os.makedirs(f'{database_dir}/organisms')

    # make orthologs dir (empty)
    os.makedirs(f'{database_dir}/orthologs')

    # make pathway maps dir and content
    os.makedirs(f'{database_dir}/pathway-maps')
    os.makedirs(f'{database_dir}/pathway-maps/svg')
    with open(f'{database_dir}/pathway-maps/type_dictionary.json', 'w') as f:
        f.write('{}')

    # Create annotations.json
    shutil.copy(src=f'{PACKAGE_ROOT}/data/annotations.json', dst=f'{database_dir}/annotations.json')

    # download annotation descriptions
    annotation_descriptions_dir = f'{database_dir}/annotation-descriptions'
    os.makedirs(annotation_descriptions_dir)
    download_sl_data(out=f'{annotation_descriptions_dir}/SL.tsv')
    download_kegg_data(src='rn', out=f'{annotation_descriptions_dir}/KR.tsv', remove_prefix='rn:')
    download_kegg_data(src='ko', out=f'{annotation_descriptions_dir}/KG.tsv', remove_prefix='ko:')
    download_kegg_data(src='enzyme', out=f'{annotation_descriptions_dir}/EC.tsv', remove_prefix='ec:', add_prefix='EC:')
    download_go_data(out=f'{annotation_descriptions_dir}/GO.tsv')


def main():
    import fire

    fire.Fire(init_database)


if __name__ == '__main__':
    main()