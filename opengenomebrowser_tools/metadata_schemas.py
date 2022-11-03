from schema import Schema, And, Optional, Or
from .utils import is_valid_date, get_cog_categories


class Dummy:
    def __repr__(self):
        return '<Dummy>'


dummy = Dummy()

organism_json_schema = Schema({
    'name': str,
    'alternative_name': Or(str, None),
    'taxid': int,
    'restricted': bool,
    'tags': [str],
    'representative': str
})

organism_json_dummy = {
    'name': dummy,
    'alternative_name': None,
    'taxid': dummy,
    'restricted': False,
    'tags': [],
    'representative': dummy
}

genome_json_schema = Schema({
    'identifier': And(str, lambda x: all(char not in x for char in ' :/\\')),
    'contaminated': bool,
    'old_identifier': Or(str, None),
    'isolation_date': Or(is_valid_date, None),
    'env_broad_scale': [str],
    'env_local_scale': [str],
    'env_medium': [str],
    'growth_condition': Or(str, None),
    'geographical_coordinates': Or(str, None),
    'geographical_name': Or(str, None),
    'library_preparation': Or(str, None),
    'sequencing_tech': Or(str, None),
    'sequencing_tech_version': Or(str, None),
    'sequencing_date': Or(is_valid_date, None),
    'sequencing_coverage': Or(str, None),
    'read_length': Or(str, None),
    'assembly_tool': Or(str, None),
    'assembly_version': Or(str, None),
    'assembly_date': Or(is_valid_date, None),
    'nr_replicons': Or(int, None),
    'cds_tool': Or(str, None),
    'cds_tool_date': Or(is_valid_date, None),
    'cds_tool_version': Or(str, None),
    'cds_tool_faa_file': str,
    'cds_tool_ffn_file': str,
    'cds_tool_gbk_file': str,
    'cds_tool_gff_file': str,
    'cds_tool_sqn_file': Or(str, None),
    'assembly_fasta_file': str,
    'custom_annotations': [{
        "date": is_valid_date,
        "file": str,
        "type": str
    }],
    'BUSCO': Or({}, {
        'C': int,
        'D': int,
        'F': int,
        'M': int,
        'S': int,
        'T': int,
        Optional('dataset'): str,
        Optional('dataset_creation_date'): is_valid_date,
    }),
    'COG': {
        Optional(cog): Or(int, float)
        for cog in get_cog_categories().keys()
    },
    'bioproject_accession': Or(str, None),
    'biosample_accession': Or(str, None),
    'genome_accession': Or(str, None),
    'literature_references': [{
        "url": str,
        "name": str
    }],
    'custom_tables': dict,
    'tags': [str]
})

genome_json_dummy = {
    "identifier": dummy,
    "contaminated": False,
    "old_identifier": None,
    "isolation_date": None,
    "env_broad_scale": [],
    "env_local_scale": [],
    "env_medium": [],
    "growth_condition": None,
    "geographical_coordinates": None,
    "geographical_name": None,
    "library_preparation": None,
    "sequencing_tech": None,
    "sequencing_tech_version": None,
    "sequencing_date": None,
    "sequencing_coverage": None,
    "read_length": None,
    "assembly_tool": None,
    "assembly_version": None,
    "assembly_date": None,
    "nr_replicons": None,
    "cds_tool": None,
    "cds_tool_date": None,
    "cds_tool_version": None,
    "cds_tool_faa_file": dummy,
    "cds_tool_ffn_file": dummy,
    "cds_tool_gbk_file": dummy,
    "cds_tool_gff_file": dummy,
    "cds_tool_sqn_file": None,
    "assembly_fasta_file": dummy,
    "custom_annotations": [],
    "BUSCO": {},
    "COG": {},
    "bioproject_accession": None,
    "biosample_accession": None,
    "genome_accession": None,
    "literature_references": [],
    "custom_tables": {},
    "tags": []
}
