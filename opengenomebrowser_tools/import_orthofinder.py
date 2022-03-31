import os
import pandas as pd
from Bio import SeqIO
from collections import Counter
from .utils import clean_locus_tag


class OrthogroupToGeneName:
    def __init__(self, fasta_dir: str, file_endings='fasta'):
        assert os.path.isdir(fasta_dir), F'fasta_dir does not exist: "{fasta_dir}"'
        self.fasta_dir = fasta_dir
        self.file_endings = file_endings

    def majority_dict(self):
        assert hasattr(self, 'majority_df'), F'Load Orthogroups.tsv or N0.tsv first!'
        return {
            orthogroup: majority_name
            for orthogroup, majority_name
            in self.majority_df['Best Gene Name'].iteritems()
        }

    def save_majority_df(self, out_file: str):
        """
        Writes the following file:

        HOG Best Gene Name  Gene Name Occurrences
        N0.HOG0000000   amino acid ABC transporter Counter({'amino acid ABC transporter': 43})
        ...
        """
        assert hasattr(self, 'majority_df'), F'Load Orthogroups.tsv or N0.tsv first!'
        self.majority_df.to_csv(path_or_buf=out_file, sep='\t')

    def save_orthogroup_to_gene_ids(self, out_file: str):
        """
        Writes the following file (no header):

        N0.HOG0000000   gene_1, gene_2
        N0.HOG0000001   gene_3, gene_4, gene_5
        ...
        """
        with open(out_file, 'w') as f:
            for orthogroup, row in self.gene_ids_df.iterrows():
                gene_ids = [gene_id for gene_ids in row for gene_id in gene_ids]
                gene_ids = ', '.join(gene_ids)
                f.write(F'{orthogroup}\t{gene_ids}\n')

    def save_orthogroup_to_best_name(self, out_file: str):
        """
        Writes the following file (no header):

        N0.HOG0000000	amino acid ABC transporter ATP-binding protein
        N0.HOG0000001	ATP-binding cassette domain-containing protein
        ...
        """

        self.majority_df.drop(
            columns=['Gene Name Occurrences'],
            inplace=False
        ).to_csv(
            path_or_buf=out_file,
            sep='\t',
            header=False
        )

    def load_og(self, og_tsv: str):
        assert os.path.isfile(og_tsv), F'og_tsv does not exist: "{og_tsv}"'

        # read "Orthogroups.tsv"
        gene_ids_df = pd.read_csv(og_tsv, sep='\t', dtype=str)
        gene_ids_df.set_index('Orthogroup', inplace=True)
        gene_ids_df = gene_ids_df.applymap(
            lambda x: [] if pd.isnull(x) else [clean_locus_tag(locus_tag) for locus_tag in x.split(', ')])

        self.gene_ids_df = gene_ids_df
        self.__load_gene_names(self.gene_ids_df.__deepcopy__())

    def load_hog(self, n0_tsv: str):
        assert os.path.isfile(n0_tsv), F'n0_tsv does not exist: "{n0_tsv}"'

        # read "N0.tsv" and drop extra columns
        gene_ids_df = pd.read_csv(n0_tsv, sep='\t', dtype=str)
        gene_ids_df.set_index('HOG', inplace=True)
        gene_ids_df.drop(columns=['OG', 'Gene Tree Parent Clade'], inplace=True)
        gene_ids_df = gene_ids_df.applymap(
            lambda x: [] if pd.isnull(x) else [clean_locus_tag(locus_tag) for locus_tag in x.split(', ')])

        self.gene_ids_df = gene_ids_df
        self.__load_gene_names(self.gene_ids_df.__deepcopy__())

    def __load_gene_names(self, gene_names_df: pd.DataFrame):
        self.strains = gene_names_df.columns

        for strain in self.strains:
            gene_id_to_name = self.__get_gene_id_to_name_dict(strain)
            gene_names_df[strain] = gene_names_df[strain].apply(lambda ids: [gene_id_to_name[id] for id in ids])

        self.majority_df = pd.DataFrame(index=gene_names_df.index)
        self.majority_df[['Best Gene Name', 'Gene Name Occurrences']] = gene_names_df.apply(
            lambda row: pd.Series(self.__majority_vote(row)), axis=1
        )

    def __majority_vote(self, row):
        """ Get gene names set per HOG and count them. Disregards names with eAED in description,
        typically useless maker annotations"""

        all_names = [name for cell in row for name in cell if "eAED" not in name]
        if len(all_names) == 0:
            all_names = ["NO DESCRIPTION"]  # if all names contain "eAED", all_names is an empty list

        names_to_count = Counter(all_names)

        return names_to_count.most_common(1)[0][0], str(
            names_to_count)  # return best gene name and gene name occurrences as string

    def __get_gene_id_to_name_dict(self, strain):
        fasta_file_path = os.path.join(self.fasta_dir + f'/{strain}.{self.file_endings}')
        assert os.path.isfile(fasta_file_path), F'{self.file_endings} file "{fasta_file_path}" is missing!'
        genes = SeqIO.parse(fasta_file_path, "fasta")

        def extract_description(gene):
            description = gene.description.split(' ', maxsplit=1)
            assert len(description) == 2, \
                F'Failed to extract description for strain={strain}:\n' \
                F'gene.id={gene.id}\n' \
                F'gene.description={gene.description}\n' \
                F'description={description}'
            description = description[1]
            if description.endswith(']'):
                # remove species description
                return description.rsplit(' [', maxsplit=1)[0]
            else:
                return description

        return {clean_locus_tag(gene.id): extract_description(gene) for gene in genes}


def import_orthofinder(
        folder_structure_dir: str = None,
        fasta_dir: str = None,
        out_annotations: str = None,
        out_descriptions: str = None,
        which: str = 'hog'
):
    """
    Creates two files:
      - folder_structure/orthologs/orthologs.tsv
      - folder_structure/annotation-descriptions/OL.tsv

    :param folder_structure_dir: Path to the root of the OpenGenomeBrowser folder structure. (Must contain 'organisms' folder.)
    :param fasta_dir: Path to the OrthoFinder FASTA dir
    :param out_annotations: Path to output annotation file (maps OG -> genes)
    :param out_descriptions: Path to output descriptions file (maps OG -> description)
    :param which: either 'hog' for N0.tsv or 'og' for 'Orthogroups.tsv
    """
    if folder_structure_dir is None:
        assert 'FOLDER_STRUCTURE' in os.environ, f'Cannot find the folder_structure. Please set --folder_structure_dir or environment variable FOLDER_STRUCTURE'
        folder_structure_dir = os.environ['FOLDER_STRUCTURE']

    if fasta_dir is None:
        fasta_dir = f'{folder_structure_dir}/OrthoFinder/fastas'

    if out_descriptions is None:
        out_descriptions = f'{folder_structure_dir}/annotation-descriptions/OL.tsv'  # is AnnotationDescriptionFile, maps OG -> description

    if out_annotations is None:
        out_annotations = f'{folder_structure_dir}/orthologs/orthologs.tsv'  # is CustomAnnotationsFile, maps HOG -> [gene]

    assert os.path.isdir(fasta_dir), f'Error: directory does not exist: {fasta_dir=}'
    assert os.path.isdir(
        f'{fasta_dir}/OrthoFinder'), f'Error: {fasta_dir=} contains no OrthoFinder directory. You must run OrthoFinder first!'
    for file in (out_descriptions, out_annotations):
        assert not os.path.isfile(file), f'Output file already exists! {file=}'

    results_dir = os.listdir(f'{fasta_dir}/OrthoFinder')
    assert len(
        results_dir) == 1, f'Error: {fasta_dir}/OrthoFinder must contain exactly one directory! dirs={results_dir}'
    results_dir = results_dir[0]

    og_tsv = f'{fasta_dir}/OrthoFinder/{results_dir}/Orthogroups/Orthogroups.tsv'
    hog_tsv = f'{fasta_dir}/OrthoFinder/{results_dir}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv'

    otg = OrthogroupToGeneName(fasta_dir, file_endings='faa')

    if which == 'hog':
        assert os.path.isfile(hog_tsv), f'Could not find {hog_tsv=}'
        otg.load_hog(n0_tsv=hog_tsv)
    elif which == 'og':
        assert os.path.isfile(og_tsv), f'Could not find {og_tsv=}'
        otg.load_og(og_tsv=og_tsv)
    else:
        raise AssertionError(f"which must be 'og' or 'hog'")

    otg.save_orthogroup_to_gene_ids(out_annotations)
    otg.save_orthogroup_to_best_name(out_descriptions)

    print(f'Sucessfully generated {out_annotations} and {out_descriptions}')


def main():
    import fire

    fire.Fire(import_orthofinder)


if __name__ == '__main__':
    main()
