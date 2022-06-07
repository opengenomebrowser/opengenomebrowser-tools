from .utils import GenomeFile, create_replace_function, split_locus_tag


class NoLocusTagInGffLine(KeyError):
    pass


class GffFile(GenomeFile):
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
            content = in_f.read()

        replace_fn = create_replace_function({
            string.format(prefix=old_locus_tag_prefix): string.format(prefix=new_locus_tag_prefix)
            for string in ['-{prefix}', '={prefix}', ':{prefix}']
        })

        old_hash = hash(content)
        content = replace_fn(content)
        assert old_hash != hash(content), f'The content of {self.path=} has not changed!'

        assert new_locus_tag_prefix in content, f'Something went wrong: did not replace anything!'

        with open(out, 'w') as out_f:
            out_f.write(content)

        if update_path:
            self.path = out

        if validate:
            self.validate_locus_tags(locus_tag_prefix=new_locus_tag_prefix)

    def detect_locus_tag_prefix(self) -> str:
        with open(self.path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                try:
                    locus_tag_prefix, gene_id = self._extract_gff_locus_tag(line)
                except NoLocusTagInGffLine:
                    continue
                return locus_tag_prefix

        raise KeyError(f'Could not extract locus_tag from {self.path=}')

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        with open(self.path) as f:
            for line in f:
                if line == '##FASTA\n':
                    break  # prokka
                if line.startswith('#'):
                    continue
                try:
                    real_locus_tag_prefix, gene_id = self._extract_gff_locus_tag(line)
                except NoLocusTagInGffLine:
                    continue  # in PGAP gffs, some lines contain no locus_tag
                assert real_locus_tag_prefix == locus_tag_prefix, \
                    f'locus_tag_prefix in {self.path=} does not match. expected: {locus_tag_prefix} reality: {real_locus_tag_prefix}'

    @staticmethod
    def _extract_gff_data(line: str) -> dict:
        line = line.rstrip('\n').split('\t')
        assert len(line) == 9, f'gff line is malformed! {len(line)=} {line=}'
        return dict(info.split('=', 1) for info in line[8].split(';') if '=' in info)

    @classmethod
    def _extract_gff_locus_tag(cls, line: str) -> (str, str):
        data = cls._extract_gff_data(line)
        if 'locus_tag' not in data:
            raise NoLocusTagInGffLine(f'gff data contains no locus_tag! {line}')
        locus_tag = data['locus_tag']
        locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
        assert ' ' not in locus_tag_prefix, f'The locus_tag may not contain blanks! {locus_tag=}'
        return locus_tag_prefix, gene_id


def rename_gff(
        file: str, out: str,
        new_locus_tag_prefix: str,
        old_locus_tag_prefix: str = None,
        validate: bool = False
):
    """
    Change the locus tags in a general feature format (gff) file

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check
    """
    GffFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate
    )


def main():
    import fire

    fire.Fire(rename_gff)


if __name__ == '__main__':
    main()
