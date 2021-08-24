# OpenGenomeBrowser Tools

A set of scripts that helps to import genome data into the OpenGenomeBrowser folder structure.

## Installation

This package requires at least `Python 3.9`.

```bash
pip install git+https://github.com/opengenomebrowser/opengenomebrowser-tools.git
```

## Usage

```bash
import_genome --help
```

### Import genome files that require no renaming

If the annotation is performed using the proper organism name, genome identifier and taxonomic information (recommended), the import is
straightforward because no files need to be renamed.

Suppose the desired organism name is `STRAIN`, the genome identifier is `STRAIN.1`, this is how to run prokka:

```bash
prokka \
  --strain STRAIN \ 
  --locustag STRAIN.1 \
  --genus Mycoplasma --species genitalium \  # Optional. If set, this script can automatically detect the taxid.
  --out /prokka/out/dir \
  assembly.fasta
```

These are the lines in PGAPs `submol.yaml` that are relevant to this script:

```yaml
organism:
  genus_species: 'Mycoplasma genitalium'  # Optional. If set, this script can automatically detect the taxid.
  strain: 'STRAIN'
locus_tag_prefix: 'STRAIN.1'
bioproject: 'PRJNA9999999'  # Optional. If set, this script can automatically add it to bioproject_accession in genome.json.
biosample: 'SAMN99999999'  # Optional. If set, this script can automatically add it to biosample_accession in genome.json.
publications: # Optional. If set, this script can automatically add it to the literature_references in genome.json.
  - publication:
      pmid: 16397293
```

This is how to import the resulting files:

```bash
export GENOMIC_DATABASE=/path/to/database  # this directory contains the 'organisms' folder
import_genome --import_dir=/prokka/out/dir  # optional: add "--organism STRAIN --genome STRAIN.1" as sanity check
```

Required files in `import_dir`:

- `.fna`: assembly (FASTA)
- `.faa`: protein sequences (FASTA)
- `.gbk`: GenBank file
- `.gff`: General feature format file

Optional files:

- `.ffn`: nucleotides file (FASTA). If non-existent, it will automatically be generated from the `.gbk` file
- `.sqn`: required for submission to GenBank, not really used by OpenGenomeBrowser
- `.emapper.annotations`: Eggnog annotation file
- `.XX`: custom annotation file (e.g. `EC`, `.GO`, etc.; any files with a suffix of two upper case letters are detected as custom annotations)
- `_busco.txt`: BUSCO output file, content will be added to `genome.json`
- `genome.json`: content will be added to final `genome.json`, may be as simple as `{"restricted": true}`
- `organism.json`: content will be added to final `organism.json`, may be as simple as `{"assembly_tool": "PGAP"}`

This will result in the following result:

```text
#### folder structure ####
database
└── organisms
    └── STRAIN
        ├── genomes
        │    └── STRAIN.1
        │	     ├── genome.json
        │	     ├── rest
        │	     │	 ├── PROKKA_08112021.err
        │	     │	 ├── PROKKA_08112021.fsa
        │
        	     │	 ├── PROKKA_08112021.log
        │	     │	 ├── PROKKA_08112021.tbl
        │	     │	 ├── PROKKA_08112021.tsv
        │	     │	 ├── PROKKA_08112021.txt
        │	     │	 └── short_summary.specific.lactobacillales_odb10.FAM3228-i1-1_busco.txt
        │	     ├── STRAIN.1.faa
        │	     ├── STRAIN.1.ffn
        │	     ├── STRAIN.1.fna
        │	     ├── STRAIN.1.gbk
        │	     ├── STRAIN.1.gff
        │	     └── STRAIN.1.sqn
        └── organism.json
```

### Import genome files that require renaming

```bash
export GENOMIC_DATABASE=/path/to/database  # this directory contains the 'organisms' folder
import_genome --import_dir=/prokka/out/dir --organism STRAIN --genome STRAIN.1 --rename
```

The renaming is provided as-is, and was only tested on files produced by certain versions of prokka and PGAP. If there is an error, you must rename
the files manually (with or without the help of my [renaming scripts](#rename-single-files)) and then import them as described in the previous
section.

### Add more metadata to genome.json or organism.json

There are two ways to achieve this.

### Modify how files are imported

It is possible to change where files end up in the folder structure. The behaviour is determined by a config file in json format that can be specified
with the --import_settings parameter or the OGB_IMPORT_SETTINGS environment variable.

```bash
export OGB_IMPORT_SETTINGS=/path/to/import_config.json
```

This is the default setting:

```text
{
    "organism_template": {},                           # use this to add metadata to all imported organism.json files, e.g. {"restricted": true}
    "genome_template": {},                             # use this to add metadata to all imported genome.json files, e.g. {"assembly_tool": "PGAP"}
    "path_transformer": {
        ".*\\.fna": "{genome}.{suffix}",               # all files that match the regex will end up in organisms/STRAIN/genomes/STRAIN.1/STRAIN.1.fna
        ".*\\.faa": "{genome}.{suffix}",
        ".*\\.gbk": "{genome}.{suffix}",
        ".*\\.gff": "{genome}.{suffix}",
        ".*\\.sqn": "{genome}.{suffix}",
        ".*\\.ffn": "{genome}.{suffix}",
        ".*\\.emapper.annotations": "{genome}.eggnog",
        ".*\\.[A-Z]{2}": "{genome}.{suffix}",
        "genome.md": "genome.md", 
        "organism.md": "../../organism.md",            # this file will end up in /organisms/STRAIN/organism.md
        "genome.json": null,                           # this file will not be copied
        "organism.json": null,                         # this file will not be copied
        ".*": "rest/{original_path}"                   # this regex matches all files, thus all files that did not match any previous regex will
                                                       #   will end up in .../STRAIN.1/rest/
    }
}

```

## Additional scripts

### Rename single files

This package contains a set of scripts to change the locus tags of genome-associated files:

| Script                         | Purpose                                             |
|--------------------------------|-----------------------------------------------------|
| `rename_faa_fnn`            | FASTA files with locus tags (i.e. `STRAIN.1_00001`) |
| `rename_genbank`            | GenBank files (tested with prokka and PGAP files)   |
| `rename_gff`                | General feature format files                        |
| `rename_eggnog`             | Eggnog `.emapper.annotations`-file                  |
| `rename_custom_annotations` | OpenGenomeBrowser custom annotations file           |

The syntax is always the same:

```bash
rename_??? \
  --file /path/to/input.file \
  --out /path/to/output.file \
  --new_locus_tag_prefix STRAIN.2 \
  --old_locus_tag_prefix STRAIN.1  # optional, good as sanity check
```

There is another script, `reindex_assembly`, which changes the header of FASTA files. Usage:

```bash
reindex_assembly \
  --file /path/to/input.file \
  --out /path/to/output.file \
  --prefix STRAIN_scf \
  --leading_zeroes 5  # optional
```

This would transform a FASTA header like this `>anything here` into `>STRAIN_scf_00001`.
