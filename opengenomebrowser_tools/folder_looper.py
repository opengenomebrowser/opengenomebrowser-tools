import json
from functools import cached_property
from json.decoder import JSONDecodeError
import os
import shutil
from datetime import datetime
from .metadata_schemas import organism_json_schema, genome_json_schema


def set_to_list(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError(f'Could not serialize {obj}')


class FolderEntity:
    path: str
    json_path: str

    def __init__(self, path: str):
        assert os.path.isdir(path), path
        self.path = path

    @property
    def is_ignored(self):
        return os.path.isfile(f'{self.path}/ignore')

    @property
    def has_json(self):
        return self.json is not None

    @cached_property
    def json(self):
        try:
            with open(self.json_path) as f:
                return json.load(f)
        except (FileNotFoundError, JSONDecodeError):
            return None

    def get_json_attr(self, attr: str):
        assert self.has_json, f'Could not extract {attr=} from {self.json_path}.'
        return self.json[attr]

    def sanity_check(self):
        if self.is_ignored:
            return True
        else:
            assert self.has_json, f'{self} :: does not have a (valid) json'

    def replace_json(self, data):
        assert type(data) is dict
        # ensure data is serializable
        try:
            json.dumps(data, default=set_to_list)
        except json.JSONDecodeError as e:
            raise AssertionError(f'Could not save dictionary as json: {e}')

        date = datetime.now().strftime("%Y_%b_%d_%H_%M_%S")

        # create backup
        bkp_dir = F'{self.path}/.bkp'
        os.makedirs(bkp_dir, exist_ok=True)
        bkp_file = F'{bkp_dir}/{date}_folderlooper_organism.json'
        shutil.move(src=self.json_path, dst=bkp_file)

        # write new file
        assert not os.path.isfile(self.json_path)
        with open(self.json_path, 'w') as f:
            json.dump(data, f, sort_keys=True, indent=4, default=set_to_list)


class FolderOrganism(FolderEntity):
    def __init__(self, path: str):
        super().__init__(path)
        self.name = os.path.basename(self.path)
        self.json_path = self.path + "/organism.json"

    def __str__(self):
        return f'<Organism {self.name}>'

    def sanity_check(self):
        super().sanity_check()

        assert self.name == self.get_json_attr('name'), \
            f"{self} :: 'name' in organism.json doesn't match folder name: {self.path}"

        organism_json_schema.validate(self.json)

        representative = self.representative()

        assert representative.identifier.startswith(self.name), \
            f'{self} :: representative identifer must start with organism name'
        assert os.path.isdir(representative.path), \
            F"{self} :: Representative doesn't exist! Representative: {representative}"
        assert not representative.is_ignored, \
            F"{self} :: Representatives may not be ignored! Representative: {representative}, {representative.path}, {representative.is_ignored}"

    @property
    def genomes_path(self):
        return f'{self.path}/genomes'

    @property
    def representative_path(self):
        return f'{self.path}/genomes/{self.get_json_attr("representative")}'

    def representative(self, sanity_check=True):
        """returns: Genome"""
        representative = FolderGenome(path=self.representative_path, organism=self)
        if sanity_check:
            representative.sanity_check()
        return representative

    def genomes(self, skip_ignored: bool, sanity_check=True) -> []:
        """generator, yields [Genome]"""
        genomes_folders = os.scandir(self.genomes_path)
        for genome_folder in genomes_folders:
            genome = FolderGenome(path=genome_folder.path, organism=self)
            if skip_ignored and genome.is_ignored:
                continue
            if sanity_check:
                genome.sanity_check()
            yield genome


class FolderGenome(FolderEntity):
    def __init__(self, path: str, organism: FolderOrganism):
        super().__init__(path)
        self.identifier = os.path.basename(self.path)
        self.organism = organism
        self.json_path = self.path + "/genome.json"

    def __str__(self):
        return f'<Genome {self.identifier}>'

    def sanity_check(self):
        super().sanity_check()

        assert self.identifier == self.get_json_attr('identifier'), \
            F"{self} :: 'identifier' in genome.json doesn't match folder name: {self.path}"

        assert self.identifier.startswith(self.organism.name), f'{self} :: identifer must start with organism name'

        genome_json_schema.validate(self.json)


class FolderLooper:
    def __init__(self, folder_structure_dir: str):
        self.folder_structure_dir = folder_structure_dir
        self.organism_path = os.path.join(self.folder_structure_dir, 'organisms')

        for dir in (self.folder_structure_dir, self.organism_path):
            assert os.path.isdir(dir), f'Folder does not exist: {dir=}'

    def organisms(self, skip_ignored: bool = True, sanity_check: bool = True) -> [FolderOrganism]:
        """generator"""
        for organism_folder in os.scandir(self.organism_path):
            organism = FolderOrganism(path=organism_folder.path)
            if skip_ignored and organism.is_ignored:
                continue
            if sanity_check:
                organism.sanity_check()
            yield organism

    def genomes(self, skip_ignored: bool = True, sanity_check: bool = True, representatives_only: bool = False) -> [FolderGenome]:
        """generator"""
        for organism in self.organisms(skip_ignored=skip_ignored, sanity_check=sanity_check):
            if representatives_only:
                yield organism.representative(sanity_check=sanity_check)
            else:
                for genome in organism.genomes(skip_ignored=skip_ignored, sanity_check=sanity_check):
                    yield genome


def loop(folder_structure_dir: str = None, what: str = 'genomes', skip_ignored: bool = True, sanity_check: bool = True, representatives_only: bool = False):
    if folder_structure_dir is None:
        folder_structure_dir = os.environ.get('FOLDER_STRUCTURE')

    folder_looper = FolderLooper(folder_structure_dir)
    if what == 'genomes':
        for organism in folder_looper.organisms(skip_ignored=skip_ignored, sanity_check=sanity_check):
            print(organism.path)
    elif what == 'organisms':
        for genome in folder_looper.genomes(skip_ignored=skip_ignored, sanity_check=sanity_check, representatives_only=representatives_only):
            print(genome.path)
    else:
        raise AssertionError(f"{what=} must be either 'genomes' or 'organisms'.")


def main():
    import fire

    fire.Fire(loop)


if __name__ == '__main__':
    main()
