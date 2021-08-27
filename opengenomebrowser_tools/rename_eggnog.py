from datetime import datetime
from functools import cached_property

from .utils import GenomeFile, split_locus_tag

EGGNOG_VERSIONS = {
    'eggnog-2.1.2':
        '#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs\n',
    'eggnog':
        '#query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	best_tax_level	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction\n',
}


class EggnogFile(GenomeFile):
    def rename(self, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None, validate: bool = False) -> None:
        old_locus_tag_prefix = self._pre_rename_check(out, new_locus_tag_prefix, old_locus_tag_prefix)

        with open(self.path) as in_f:
            content = in_f.readlines()

        def rename_line(line: str):
            if line.startswith('#'):
                return line

            locus_tag, rest = line.split('\t', 1)
            locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
            assert locus_tag_prefix == old_locus_tag_prefix, f'Eggnog line does not contain old_locus_tag_prefix!' \
                                                             f'{old_locus_tag_prefix=}, {line=}, {self.path=}'
            return locus_tag.replace(old_locus_tag_prefix, new_locus_tag_prefix, 1) + '\t' + rest

        content = [rename_line(line) for line in content]

        with open(out, 'w') as out_f:
            out_f.writelines(content)

        self.path = out

        if validate:
            self.validate_locus_tags(locus_tag_prefix=new_locus_tag_prefix)

    def detect_locus_tag_prefix(self) -> str:
        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                locus_tag = line.split('\t', 1)[0]
                locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
                return locus_tag_prefix

        raise KeyError(f'Could not extract locus_tag from {self.path=}, it does not appear to contain annotations!')

    def date(self) -> datetime:
        with open(self.path) as f:
            head = [next(f) for x in range(4)]

        if head[0].startswith('##'):
            date_line = head[0].strip()
            dt = datetime.strptime(date_line[3:], '%c')
        elif head[0].startswith('#'):
            date_line = head[2].strip()
            dt = datetime.strptime(date_line[8:], '%c')
        else:
            dt = super().date()  # use creation date of file

        return dt

    @cached_property
    def custom_annotation_type(self) -> str:
        with open(self.path) as f:
            head = '\n'.join(next(f) for x in range(5))  # read 5 lines

        for type, columns_header in EGGNOG_VERSIONS.items():
            if columns_header in head:
                return type

        raise KeyError(f'Could not discover eggnog type! {self.path=}')

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        if locus_tag_prefix is None:
            locus_tag_prefix = self.detect_locus_tag_prefix()

        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                locus_tag = line.split('\t', 1)[0]
                real_locus_tag_prefix, gene_id = split_locus_tag(locus_tag)

                assert real_locus_tag_prefix == locus_tag_prefix, \
                    f'locus_tag_prefix in {self.path=} does not match. expected: {locus_tag_prefix} reality: {real_locus_tag_prefix}'
                assert gene_id.isdigit(), f'locus_tag in {self.path=} is malformed. expected: {locus_tag_prefix}_[0-9]+ reality: {locus_tag}'


def rename_eggnog(file: str, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None, validate: bool = False):
    """
    Change the locus tags in a eggnog output file (.emapper.annotations)

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check
    """

    EggnogFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate
    )


def main():
    import fire

    fire.Fire(rename_eggnog)


if __name__ == '__main__':
    main()
