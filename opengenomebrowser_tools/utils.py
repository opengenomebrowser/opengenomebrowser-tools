import gzip
import json
import logging
import os
import re
from datetime import datetime
from string import digits
from typing import Union, Callable

from Bio import Entrez
from termcolor import colored

DATE_FORMAT = '%Y-%m-%d'
TMPDIR = os.environ.get('TMPDIR', '/tmp')
Entrez.email = os.environ.get('ENTREZ_EMAIL', 'opengenomebrowser@bioinformatics.unibe.ch')

PACKAGE_ROOT = os.path.dirname(__file__)

ANNOTATIONS_JSON = f'{PACKAGE_ROOT}/data/annotations.json'
COG_CATEGORIES_JSON = f'{PACKAGE_ROOT}/data/COG_categories.json'
for f in [ANNOTATIONS_JSON, COG_CATEGORIES_JSON]:
    assert os.path.isfile(f), f'Package is poorly configured: file is missing: {f}'


class WorkingDirectory:
    """
    Context manager to temporarily change the working directory.

    Example:

    with WorkingDirectory('some/path'):
        some_function()
    """
    new: str
    old: str

    def __init__(self, new_working_directory: str):
        assert os.path.isdir(new_working_directory), f'Failed to enter {new_working_directory=}: does not exist!'
        self.new = new_working_directory

    def __enter__(self):
        self.old = os.getcwd()
        os.chdir(self.new)

    def __exit__(self, *args, **kwargs):
        os.chdir(self.old)


class GenomeFile:
    original_path: str
    target_path: str

    def __init__(self, file: str, original_path: str = None):
        self.path = file
        self.original_path = original_path
        assert os.path.isfile(file), f'{self}: {file=} does not exist!'

    def __str__(self):
        return f'{type(self).__name__}: {os.path.basename(self.path)}'

    def metadata(self) -> (dict, dict):
        return {}, {}

    def detect_locus_tag_prefix(self) -> str:
        raise NotImplementedError('This function must be overwritten.')

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        raise NotImplementedError('This function must be overwritten.')

    def rename(self, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None,
               validate: bool = True, update_path: bool = True) -> None:
        raise NotImplementedError('This function must be overwritten.')

    def _pre_rename_check(self, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None) -> str:
        assert not os.path.isfile(out), f'Output file already exists! {out=}'

        if old_locus_tag_prefix is None:
            old_locus_tag_prefix = self.detect_locus_tag_prefix()
            logging.info(f'detected locus tag: {old_locus_tag_prefix}')
        else:
            detected_locus_tag_prefix = self.detect_locus_tag_prefix()
            logging.info(f'detected locus tag: {detected_locus_tag_prefix}')
            assert old_locus_tag_prefix == detected_locus_tag_prefix, f'Locus tag prefix does not match! {detected_locus_tag_prefix=} != {old_locus_tag_prefix=}'

        assert new_locus_tag_prefix != old_locus_tag_prefix, f'Aborting: {new_locus_tag_prefix=} is the same as {old_locus_tag_prefix=}'

        return old_locus_tag_prefix

    def date(self) -> datetime:
        return get_ctime(file=self.path)

    def date_str(self) -> str:
        """
        :returns date in this format: "%Y-%m-%d"
        """
        return date_to_string(self.date())


def query_yes_no(question: str, default: str = None, color='blue') -> bool:
    """
    Ask a yes/no question via raw_input() and return their answer

    :param question: string that is presented to the user
    :param default: the presumed answer if the user just hits <Enter>
    :param color: color in which the question is formatted
    :return: response (bool)
    """
    valid = {'yes': True, 'y': True, 'ye': True, 'no': False, 'n': False}
    if default is None:
        prompt = '[y/n] '
    elif default == 'yes':
        prompt = '[Y/n] '
    elif default == 'no':
        prompt = '[y/N] '
    else:
        raise ValueError(F"invalid default answer: '{default}'")

    while True:
        print(colored(F'{question} {prompt}', color), end='')
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            print(colored("'Please respond with 'y' or 'n'.'", color))


def query_int(question: str, error_msg='Please respond with an integer.', color='blue') -> int:
    while True:
        print(colored(question, color), end='')
        choice = input()
        if choice.isdigit():
            return int(choice)
        else:
            print(colored(error_msg, color))


def date_to_string(dt: datetime) -> str:
    return dt.strftime('%Y-%m-%d')


def get_ctime(file: str) -> datetime:
    ctime = os.path.getctime(file)
    dt = datetime.fromtimestamp(ctime)
    return dt


def is_valid_date(date: str) -> bool:
    try:
        datetime.strptime(date, DATE_FORMAT)
        return True
    except ValueError:
        return False


def get_cog_categories(reload: bool = False) -> dict:
    file = f'{PACKAGE_ROOT}/data/COG_categories.json'

    if not reload and os.path.isfile(file):
        with open(file) as f:
            cog_categories = json.load(f)
        return cog_categories

    else:
        from urllib import request
        with request.urlopen('https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab') as f:
            content = f.read().decode('utf-8').strip()

        def extract_cat(line: str) -> (str, dict):
            cat, color, description = line.split('\t')
            return cat, {'description': description, 'color': f'#{color.lower()}'}

        cog_categories = dict(extract_cat(cat) for cat in content.split('\n'))

        with open(file, 'w') as f:
            json.dump(cog_categories, f, indent=4)

        return cog_categories


def _get_cache_json(cache_file: str) -> dict:
    try:
        with open(cache_file) as f:
            return json.load(f)
    except Exception:
        return {}


def _set_cache_json(cache_file: str, content: Union[dict, list]) -> None:
    try:
        with open(cache_file, 'w') as f:
            json.dump(content, f)
    except Exception:
        logging.warning(f'Failed to write {content=} to {cache_file=}')


def entrez_organism_to_taxid(organism: str) -> int:
    cache_file = os.path.join(TMPDIR, 'organism_to_taxid.json')

    taxid_cache = _get_cache_json(cache_file)

    if organism in taxid_cache:
        logging.info(f'loaded taxid of {organism} from cache')
        return taxid_cache[organism]
    else:
        logging.info(f'loading taxid of {organism} from Entrez')
        handle = Entrez.esearch(db='Taxonomy', term=organism)
        record = Entrez.read(handle)
        assert 'ErrorList' not in record, f'Failed to extract taxid of {organism=} using Entrez!'
        taxids = record["IdList"]
        assert len(taxids) == 1 and taxids[0].isdigit(), f'Failed to extract taxid of {organism=} using Entrez!'
        taxid = int(taxids[0])

        taxid_cache[organism] = taxid
        _set_cache_json(cache_file, taxid_cache)
        return taxid


def clean_locus_tag(locus_tag: str) -> (str):
    return locus_tag.rsplit('|', 1)[-1]


def split_locus_tag(locus_tag: str) -> (str, str):
    locus_tag = clean_locus_tag(locus_tag)
    prefix = locus_tag.rstrip(digits)
    assert len(prefix) < len(locus_tag), f'Failed to detect {prefix=} from {locus_tag=}. Locus tags must end in digits'
    return prefix, locus_tag[len(prefix):]


def create_replace_function(replace_map: {str: str}) -> Callable:
    '''
    Returns a function that will replace all replace_map.keys with their corresponding replace_map.values

    :param replace_map: dictionary that maps strings to be replaced to their desired replacement
    :rtype: Callable
    :return: function that takes str and returns str
    '''
    replace_map = {re.escape(k): v for k, v in replace_map.items()}  # escape key
    pattern = re.compile("|".join(replace_map.keys()))

    def replace_fn(text: str) -> str:
        return pattern.sub(lambda m: replace_map[re.escape(m.group(0))], text)

    return replace_fn


def decompress_gz(gz: str, out: str):
    with gzip.open(gz, 'rb') as f_in, open(out, 'w') as f_out:
        for line in f_in:
            f_out.write(line.decode('utf-8'))


def merge_json(dict_: dict, new: str = None) -> dict:
    if new is not None and os.path.isfile(new):
        with open(new) as f:
            new_dict = json.load(f)

            for key, value in new_dict.items():
                if value is None:
                    # do not overwrite with None
                    continue
                if type(value) in [list, dict] and len(value)==0:
                    # do not overwrite with empty
                    continue

                dict_[key]=value

    return dict_


def get_folder_structure_version(folder_structure_dir: str) -> int:
    """
    Determine current folder structure version. (Read folder_structure/version.json)

    :param folder_structure_dir: Path to the root of the OpenGenomeBrowser folder structure. (Must contain 'organisms' folder.)
    :return: version (integer)
    """
    assert type(folder_structure_dir) is str
    version_file = f'{folder_structure_dir}/version.json'

    if not os.path.isfile(version_file):
        with open(version_file, 'w') as f:
            json.dump({'folder_structure_version': 1}, f, indent=4)

    try:
        with open(version_file) as f:
            version_dict = json.load(f)
    except Exception as e:
        raise AssertionError(f'Failed to {version_file} as json: {str(e)}')

    assert 'folder_structure_version' in version_dict, f"Key 'folder_structure_version' missing in {version_file}"
    version = version_dict['folder_structure_version']
    assert type(version) is int, f"Key 'folder_structure_version' in {version_file} must be an integer!"

    return version
