from re import compile
from datetime import datetime
from functools import cached_property

from .utils import GenomeFile, split_locus_tag, get_cog_categories

EGGNOG_VERSIONS = {
    'eggnog-2.1.2':
        '#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\tPFAMs\n',
    'eggnog':
        '#query_name\tseed_eggNOG_ortholog\tseed_ortholog_evalue\tseed_ortholog_score\tbest_tax_level\tPreferred_name\tGOs\tEC\tKEGG_ko\tKEGG_Pathway\tKEGG_Module\tKEGG_Reaction\tKEGG_rclass\tBRITE\tKEGG_TC\tCAZy\tBiGG_Reaction\n',
}


class EggnogFile(GenomeFile):
    def rename(
            self,
            out: str,
            new_locus_tag_prefix: str,
            old_locus_tag_prefix: str = None,
            validate: bool = False,
               update_path: bool = True
    ) -> None:
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

        if update_path:
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

    def cog_categories(self) -> dict:
        cog_categories = get_cog_categories()

        re_categories = compile(pattern=f"-|[{''.join(cog_categories.keys())}]+")

        cog_to_count = {cat: 0 for cat in cog_categories.keys()}
        cog_to_count['-'] = 0  # eggnog assigns '-' sometimes; surely it means the same as 'S': Function unknown
        n_genes = 0

        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                n_genes += 1

                cog_cats = line.split('\t')[6]

                assert re_categories.fullmatch(cog_cats) is not None, f'Failed to interpret "{cog_cats}" as COG categories. {self.path=} {line=}'

                weight = 1 / len(cog_cats)
                for cog_cat in cog_cats:
                    cog_to_count[cog_cat] += weight

        # merge '-' and 'S'
        cog_to_count['S'] += cog_to_count['-']
        del cog_to_count['-']

        # remove empty categories
        cog_to_count = {cog_cat: count for cog_cat, count in cog_to_count.items() if count != 0}

        assert n_genes == round(sum(cog_to_count.values()))

        return cog_to_count


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
