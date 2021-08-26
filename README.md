# OpenGenomeBrowser Tools

A set of scripts that helps to import genome data into the OpenGenomeBrowser folder structure.

## Installation

This package requires at least `Python 3.9`.

```bash
pip install git+https://github.com/opengenomebrowser/opengenomebrowser-tools.git
```

## Scripts

This package contains the following scripts

| Script                      | Purpose                                                                |
|-----------------------------|------------------------------------------------------------------------|
| `import_genome`             | Import genome-associated files into OpenGenomeBrowser folder structure, automatically generate metadata files |
| `rename_fasta`              | Change locus tags of FASTA files                                       |
| `rename_genbank`            | Change locus tags of GenBank files (tested with prokka and PGAP files) |
| `rename_gff`                | Change locus tags of gff (general feature format) files                |
| `rename_eggnog`             | Change locus tags of Eggnog files (`.emapper.annotations`)             |
| `rename_custom_annotations` | Change locus tags of custom annotations files                          |
| `reindex_assembly`          | Change FASTA headers of assembliy files                                |
| `gbk_to_ffn`                | Convert GenBank (`.gbk`) to nucleotide FASTA (`.ffn`)                  |

All of these scripts have help functions, for example:

```bash
import_genome --help
```

## import_genome

If the annotation was performed using the proper organism name, genome identifier and taxonomic information (recommended), the import is
straightforward because no files need to be renamed.

```bash
export GENOMIC_DATABASE=/path/to/database   # this directory contains the 'organisms' folder
import_genome --import_dir=/prokka/out/dir  # optional: add "--organism STRAIN --genome STRAIN.1" as sanity check
```

### `import_genome`: Renaming files

```bash
export GENOMIC_DATABASE=/path/to/database   # this directory contains the 'organisms' folder
import_genome --import_dir=/prokka/out/dir --organism STRAIN --genome STRAIN.1 --rename
```

The renaming is provided as-is, and was only tested on files produced by certain versions of prokka and PGAP. If there is an error, you must rename
the files manually (with or without the help of my [renaming scripts](#rename_-rename-locus-tags-in-genome-associated-files)) and then import them as
described in the previous section.

### `import_genome`: Required files

These files need to be in `import_dir`:

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

<details>
  <summary>This will result in the following result:</summary>

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

</details>

### `import_genome`: How to run annotation software to get correct locus tags

- **Prokka**

<details>
  <summary>Suppose the desired organism name is `STRAIN`, the genome identifier is `STRAIN.1`, this is how to run prokka:</summary>

```bash
prokka \
  --strain STRAIN \ 
  --locustag STRAIN.1 \
  --genus Mycoplasma --species genitalium \  # Optional. If set, this script can automatically detect the taxid.
  --out /prokka/out/dir \
  assembly.fasta
```

</details>

- **PGAP**

<details>
  <summary>These are the lines in PGAPs `submol.yaml` that are relevant to this script:</summary>

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

</details>

### `import_genome`: Modify where files are moved to

It is possible to change where files end up in the folder structure. The behaviour is determined by a config file in json format that can be specified
with the --import_settings parameter or the OGB_IMPORT_SETTINGS environment variable.

```bash
export OGB_IMPORT_SETTINGS=/path/to/import_config.json
```

<details>
  <summary>These are the default settings:</summary>

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

</details>

<details>
  <summary>This is an example of an alternative configuration:</summary>

```text
{
    "organism_template": {},
    "genome_template": {},
    "path_transformer": {
        
        # raw reads
        ".*fastqc?\\..*": "0_raw_reads/{original_path}",
        
        # assembly
        ".*\\.fna": "1_assembly/{genome}.{suffix}",
        
        # coding sequence (CDS) calling
        ".*\\.faa": "2_cds/{genome}.{suffix}",
        ".*\\.gbk": "2_cds/{genome}.{suffix}",
        ".*\\.gff": "2_cds/{genome}.{suffix}",
        ".*\\.ffn": "2_cds/{genome}.{suffix}",
        ".*\\.sqn": "2_cds/{genome}.{suffix}",
        "PROKKA_.*": "2_cds/{original_path}",
        
        # functional annotations
        ".*\\.emapper.annotations": "3_annotation/{genome}.eggnog",
        ".*\\.[A-Z]{2}": "3_annotation/{genome}.{suffix}",
        ".*_busco\\.txt": "3_annotation/{original_path}",
        
        # special files
        "genome.md": "genome.md",
        "organism.md": "../../organism.md",
        "genome.json": null,
        "organism.json": null,
        
        # rest
        ".*": "rest/{original_path}"
    }
}
```

This is what the result looks like:

```text
#### folder structure ####
database
└── organisms
    └── STRAIN
       ├── genomes
       │     └── STRAIN.1
       │         ├── 1_assembly
       │         │     └── STRAIN.1.fna
       │         ├── 2_cds
       │         │     ├── PROKKA_08112021.err
       │         │     ├── PROKKA_08112021.fsa
       │         │     ├── PROKKA_08112021.log
       │         │     ├── PROKKA_08112021.tbl
       │         │     ├── PROKKA_08112021.tsv
       │         │     ├── PROKKA_08112021.txt
       │         │     ├── STRAIN.1.faa
       │         │     ├── STRAIN.1.ffn
       │         │     ├── STRAIN.1.gbk
       │         │     ├── STRAIN.1.gff
       │         │     └── STRAIN.1.sqn
       │         ├── 3_annotation
       │         │     └── short_summary_busco.txt
       │         └── genome.json
       └── organism.json
```

</details>

### `import_genome`: Add custom metadata

There are two ways to achieve this:

1) Add a `organism.json` and/or `genome.json` file into `import_dir` (see [import_genome: Required files](#import_genome-required-files))
2) Set a global `organism.json` and/or `genome.json` file that is used as a basis for all future imports (
   see [import_genome: Modify where files are moved to](#import_genome-modify-where-files-are-moved-to))

## `rename_*`: Rename locus tags in genome-associated files

All rename-scripts (`rename_fasta`, `rename_genbank`, `rename_gff`, `rename_eggnog`, `rename_custom_annotations`) have the same syntax:

```bash
rename_fasta \
  --file /path/to/input.file \
  --out /path/to/output.file \
  --new_locus_tag_prefix STRAIN.2 \
  --old_locus_tag_prefix STRAIN.1  # optional, good as sanity check
```

## `reindex_assembly`: Change the header line of FASTA files

This script changes the header of FASTA files.

```bash
reindex_assembly \
  --file /path/to/input.file \
  --out /path/to/output.file \
  --prefix STRAIN_scf \
  --leading_zeroes 5  # optional
```

This would transform a FASTA header like this `>anything here` into `>STRAIN_scf_00001`.

## `gbk_to_ffn`: Convert GenBank to nucleotide FASTA

Usage:

```bash
gbk_to_ffn \
  --gbk /path/to/input.gbk \
  --ffn /path/to/output.ffn
```
