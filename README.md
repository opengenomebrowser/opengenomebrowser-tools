# OpenGenomeBrowser Tools

A set of scripts that helps to import genome data into the OpenGenomeBrowser folder structure.

## Installation

This package requires at least `Python 3.9`.

```bash
pip install opengenomebrowser-tools
```

## Help function

All scripts have a help function, for example:

```bash
import_genome --help
```

## `init_folder_structure`

Creates a basic OpenGenomeBrowser folders structure.

<details>
  <summary>More details:</summary>

Once the folder structure has been initiated...

- use [`import_genome`](#import_genome) to add genomes to the folder structure
- use [`download_ncbi_genome`](#download_ncbi_genome) and [`import_genome`](#import_genome) to download and add genomes from NCBI
- when all genomes have been added, use [`init_orthofinder`](#init_orthofinder) and [`import_orthofinder`](#import_orthofinder) to calculate
  orthologs (optional)

Usage:

```shell
export FOLDER_STRUCTURE=/path/to/folder_structure
init_folder_structure  # or --folder_structure_dir=/path/to/folder_structure
```

<details>
  <summary>Result:</summary>

```
  folder_structure
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
```

</details>

</details>

## `import_genome`

Import genome-associated files into OpenGenomeBrowser folder structure, automatically generate metadata files.

<details>
  <summary>More details:</summary>

If the annotation was performed using the proper organism name, genome identifier and taxonomic information (recommended), the import is
straightforward because no files need to be renamed.

```shell
export FOLDER_STRUCTURE=/path/to/folder_structure   # this directory contains the 'organisms' folder
import_genome --import_dir=/prokka/out/dir  # optional: add "--organism STRAIN --genome STRAIN.1" as sanity check
```

<details>
  <summary>How to run prokka to get correct locus tags</summary>

Suppose the desired organism name is `STRAIN`, the genome identifier is `STRAIN.1`, this is how to run prokka:

```shell
prokka \
  --strain STRAIN \ 
  --locustag STRAIN.1 \
  --prefix STRAIN.1 \
  --genus Mycoplasma --species genitalium \  # Optional. If set, this script can automatically detect the taxid.
  --out /prokka/out/dir \
  assembly.fasta
```

</details>

<details>
  <summary>How to run PGAP to get correct locus tags</summary>

Suppose the desired organism name is `STRAIN`, the genome identifier is `STRAIN.1`, these are the lines in PGAPs `submol.yaml` that are relevant to
this script:

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

### Rename files during import

Should the locus tags not start with the genome identifier, the files need to be renamed accordingly. The `import_genome` command can do this
automatically sung the `--rename` flag.

```shell
export FOLDER_STRUCTURE=/path/to/folder_structure   # this directory contains the 'organisms' folder
import_genome --import_dir=/prokka/out/dir --organism STRAIN --genome STRAIN.1 --rename
```

The renaming is provided as-is, and was only tested on files produced by certain versions of prokka and PGAP. If there is an error, you must rename
the files manually (with or without the help of my [renaming scripts](#rename_-rename-locus-tags-in-genome-associated-files)) and then import them as
described in the previous section.

### Required files

These files need to be in `import_dir`:

- `.fna`: assembly (FASTA)
- `.gbk`: GenBank file
- `.gff`: General feature format file

Optional files:

- `.faa`: protein sequences (FASTA). If non-existent, it will automatically be generated from the `.gbk` file
- `.ffn`: nucleotides file (FASTA). If non-existent, it will automatically be generated from the `.gbk` file
- `.sqn`: required for submission to GenBank, not really used by OpenGenomeBrowser
- `.emapper.annotations`: Eggnog annotation file
- `.XX`: custom annotation file (e.g. `EC`, `.GO`, etc.; any files with a suffix of two upper case letters are detected as custom annotations)
- `_busco.txt`: BUSCO output file, content will be added to `genome.json`
- `genome.json`: content will be added to final `genome.json`, may be as simple as `{"restricted": true}`
- `organism.json`: content will be added to final `organism.json`, may be as simple as `{"assembly_tool": "PGAP"}`

<details>
  <summary>Example result:</summary>

```text
#### folder structure ####
folder_structure
└── organisms
    └── STRAIN
        ├── organism.json
        └── genomes
             └── STRAIN.1
         	     ├── genome.json
         	     ├── STRAIN.1.faa
         	     ├── STRAIN.1.ffn
         	     ├── STRAIN.1.fna
         	     ├── STRAIN.1.gbk
         	     ├── STRAIN.1.gff
         	     ├── STRAIN.1.sqn
         	     └── rest
         	      	 ├── PROKKA_08112021.err
         	      	 ├── PROKKA_08112021.fsa
        	      	 ├── PROKKA_08112021.log
         	      	 ├── PROKKA_08112021.tbl
         	      	 ├── PROKKA_08112021.tsv
         	      	 ├── PROKKA_08112021.txt
         	      	 └── short_summary.specific.lactobacillales_odb10.FAM3228-i1-1_busco.txt
```

</details>

### Modify where files are moved to

It is possible to change where files end up in the folder structure. The behaviour is determined by a config file in json format that can be specified
with the --import_settings parameter or the `OGB_IMPORT_SETTINGS` environment variable.

```shell
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

Result:

```text
#### folder structure ####
folder_structure
└── organisms
    └── STRAIN
       ├── organism.json
       └── genomes
             └── STRAIN.1
                 ├── genome.json
                 ├── 1_assembly
                 │     └── STRAIN.1.fna
                 ├── 2_cds
                 │     ├── PROKKA_08112021.err
                 │     ├── PROKKA_08112021.fsa
                 │     ├── PROKKA_08112021.log
                 │     ├── PROKKA_08112021.tbl
                 │     ├── PROKKA_08112021.tsv
                 │     ├── PROKKA_08112021.txt
                 │     ├── STRAIN.1.faa
                 │     ├── STRAIN.1.ffn
                 │     ├── STRAIN.1.gbk
                 │     ├── STRAIN.1.gff
                 │     └── STRAIN.1.sqn
                 └── 3_annotation
                       └── short_summary_busco.txt
```

</details>

### Add custom metadata

There are two ways to achieve this:

1) Add a `organism.json` and/or `genome.json` file into `import_dir` (see [import_genome: Required files](#required-files))
2) Set a global `organism.json` and/or `genome.json` file that is used as a basis for all future imports (
   see [import_genome: Modify where files are moved to](#modify-where-files-are-moved-to))

</details>

## `rename_*`

The following scripts change the locus tags in the respective file formats.

| Script                      | Purpose                                                                  |
|-----------------------------|--------------------------------------------------------------------------|
| `rename_fasta`              | Change locus tags of protein or nucleotide FASTA files                   |
| `rename_genbank`            | Change locus tags of GenBank files (tested with prokka and PGAP files)   |
| `rename_gff`                | Change locus tags of gff (general feature format) files                  |
| `rename_eggnog`             | Change locus tags of Eggnog files (`.emapper.annotations`)               |
| `rename_custom_annotations` | Change locus tags of custom annotations files                            |

<details>
  <summary>More details:</summary>

The syntax is always the same.

```shell
rename_fasta \
  --file /path/to/input.file \
  --out /path/to/output.file \
  --new_locus_tag_prefix STRAIN.2 \
  --old_locus_tag_prefix STRAIN.1  # optional, good as sanity check
```

</details>

## `reindex_assembly`

This script changes the header of assembly FASTA (`.fna`) files.

<details>
  <summary>More details:</summary>

```shell
reindex_assembly \
  --file /path/to/input.file \
  --out /path/to/output.file \
  --prefix STRAIN_scf \
  --leading_zeroes 5  # optional
```

This would transform a FASTA header like this `>anything here` into `>STRAIN_scf_00001`.

</details>

## `genbank_to_fasta`

Convert GenBank to nucleotide (`.ffn`) or protein FASTA (`.faa`).

<details>
  <summary>More details:</summary>

Usage:

```shell
genbank_to_fasta \
  --gbk /path/to/input.gbk \
  --out /path/to/output.fasta \
  --format faa  # or ffn
```

</details>

## `download_ncbi_genome`

Download genome-associated files (`.fna`, `.gbk`, `.gff`) from NCBI, rename the locus tags, and generate `.ffn` and `faa` files.

<details>
  <summary>More details:</summary>

Usage:

```shell
download_ncbi_genome \
  --assembly_name GCF_005864195.1 \
  --out_dir /path/to/outdir \
  --new_locus_tag_prefix FAM3257_ 
```

Result:

```text
outdir
├── FAM3257.faa
├── FAM3257.ffn
├── FAM3257.fna
├── FAM3257.gbk
└── FAM3257.gff
```

The next step might be to import these genomes into the OpenGenomeBrowser folder structure like this:

```shell
import_genome --import_dir=/path/to/outdir --organism FAM3257 --genome FAM3257
```

</details>

## `init_orthofinder`

This script collects the protein FASTAs in `folder_structure/OrthoFinder/fastas` and prints the command to run OrthoFinder.

<details>
  <summary>More details:</summary>

Usage:

```shell
export FOLDER_STRUCTURE=/path/to/folder_structure
init_orthofinder --representatives_only
```

Result:

```
  folder_structure
  ├── ...
  └── OrthoFinder
      └── fastas
          ├── GENOME1.faa
          ├── GENOME2.faa
          └── ...
```

</details>

## `import_orthofinder`

The output of OrthoFinder needs to be processed for OpenGenomeBrowser. This script creates two files:

- `annotation-descriptions/OL.tsv`: maps orthologs to the most common gene name, i.e. `OG0000005` -> `MFS transporter`
- `orthologs/orthologs.tsv`: maps orthologs to genes, i.e. `OG0000005` -> `STRAIN1_000069, STRAIN2_000128, STRAIN2_000137`

<details>
  <summary>More details:</summary>

Usage:

```shell
export FOLDER_STRUCTURE=/path/to/folder_structure
import_orthofinder --which hog  # 'hog' for hierarchical orthogroups and 'og' for regular orthogroups
```

Once these files exist, run the following command from within the OpenGenomeBrowser docker container:

```shell
python db_setup/manage_ogb.py import-orthologs
```

</details>

## `update_folder_structure`

From time to time, changes are made to the OpenGenomeBrowser folder structure. The current version of your folder structure is denoted
in `version.json`. Use this script to upgrade to a new version.

<details>
  <summary>More details:</summary>

- `1_to_2`: add `COG` to genome.json

</details>

