import logging

from .utils import GenomeFile, split_locus_tag


class FastaFile(GenomeFile):
    def rename(self, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None,
               validate: bool = False, update_path: bool = True) -> None:
        old_locus_tag_prefix = self._pre_rename_check(out, new_locus_tag_prefix, old_locus_tag_prefix)

        with open(self.path) as in_f:
            content = in_f.readlines()

        def rename_line(line: str):
            if line.startswith('>'):
                assert old_locus_tag_prefix in line, \
                    f'Fasta header does not contain old_locus_tag_prefix! {old_locus_tag_prefix=}, {line=}, fasta={self.path}'
                return line \
                    .replace(old_locus_tag_prefix, new_locus_tag_prefix, 1) \
                    .replace(f'hypothetical protein {old_locus_tag_prefix}',
                             f'hypothetical protein {new_locus_tag_prefix}')
            else:
                return line

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
                if not line.startswith('>'):
                    assert line.strip() == '', f'Could not extract locus_tag from {self.path=}, it does not start with a header line!'
                    continue
                locus_tag_prefix, gene_id = self.parse_fasta_header(line)
                return locus_tag_prefix

        raise KeyError(
            f'Could not extract locus_tag from {self.path=}, it does not appear to contain a header line (>)!')

    def validate_locus_tags(self, locus_tag_prefix: str = None):
        with open(self.path) as f:
            for line in f:
                if line.startswith('>'):
                    real_locus_tag_prefix, gene_id = self.parse_fasta_header(header=line)
                assert real_locus_tag_prefix == locus_tag_prefix, \
                    f'locus_tag_prefix in {self.path=} does not match. expected: {locus_tag_prefix} reality: {real_locus_tag_prefix}'
                assert gene_id.isdigit(), f'locus_tag in {self.path=} is malformed. gene_id is expected to be: [0-9]+ reality: {gene_id}'

    @staticmethod
    def parse_fasta_header(header: str) -> (str, str):
        header = header.rstrip()
        error_message = f'This fasta file does not start with gene identifiers (>gene-identifier_00001)! {header=}'
        assert header.startswith('>'), error_message
        if not '_' in header:
            logging.warning(error_message)
        locus_tag = header[1:].split(' ', 1)[0]
        locus_tag_prefix, gene_id = split_locus_tag(locus_tag)
        assert ' ' not in locus_tag_prefix and len(locus_tag_prefix) > 0, error_message
        assert gene_id.isdigit(), error_message
        return locus_tag_prefix, gene_id


def rename_fasta(file: str, out: str, new_locus_tag_prefix: str, old_locus_tag_prefix: str = None,
                 validate: bool = False):
    """
    Change the locus tags in a protein/nucleotide FASTA file

    :param file: input file
    :param out: output file
    :param new_locus_tag_prefix: desired locus tag
    :param old_locus_tag_prefix: locus tag to replace
    :param validate: if true, perform sanity check
    """
    FastaFile(
        file=file
    ).rename(
        out=out,
        new_locus_tag_prefix=new_locus_tag_prefix,
        old_locus_tag_prefix=old_locus_tag_prefix,
        validate=validate
    )


def main():
    import fire

    fire.Fire(rename_fasta)


if __name__ == '__main__':
    main()
