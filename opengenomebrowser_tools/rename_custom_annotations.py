from .utils import GenomeFile, split_locus_tag


class CustomAnnotationFile(GenomeFile):
    def __init__(self, file: str, original_path: str = None, custom_annotation_type: str = None):
        if custom_annotation_type:
            self.custom_annotation_type = custom_annotation_type
        else:
            self.custom_annotation_type = file.rsplit('.', 1)[-1]
        super().__init__(file=file, original_path=original_path)

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
            assert line.startswith(
                old_locus_tag_prefix), f'custom_annotations_file line does not contain old_locus_tag_prefix!' \
                                       f'{old_locus_tag_prefix=}, {line=}, {self.path=}'
            return line.replace(old_locus_tag_prefix, new_locus_tag_prefix, 1)

        content = [rename_line(line) for line in content]

        with open(out, 'w') as out_f:
            out_f.writelines(content)

        if update_path:
            self.path = out

        if validate:
            self.validate_locus_tags(locus_tag_prefix=new_locus_tag_prefix)

    def detect_locus_tag_prefix(self) -> str:
        with open(self.path) as f:
            line = f.readline()
            locus_tag = line.split('\t', 1)[0]
            locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
            assert gene_id.isdigit(), f'locus_tag in {self.path=} is malformed. expected: {locus_tag_prefix}_[0-9]+ reality: {locus_tag}'
            return locus_tag_prefix

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        if locus_tag_prefix is None:
            locus_tag_prefix = self.detect_locus_tag_prefix()

        with open(self.path) as f:
            for line in f:
                locus_tag = line.split('\t')[0]
                real_locus_tag_prefix, gene_id = split_locus_tag(locus_tag)

                assert real_locus_tag_prefix == locus_tag_prefix, \
                    f'locus_tag_prefix in {self.path=} does not match. expected: {locus_tag_prefix} reality: {real_locus_tag_prefix}'
                assert gene_id.isdigit(), f'locus_tag in {self.path=} is malformed. expected: {locus_tag_prefix}_[0-9]+ reality: {locus_tag}'


def rename_custom_annotations(file: str, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None,
                              validate: bool = False):
    """
    Change the locus tags in a custom annotations file

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check
    """

    CustomAnnotationFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate
    )


def main():
    import fire

    fire.Fire(rename_custom_annotations)


if __name__ == '__main__':
    main()
