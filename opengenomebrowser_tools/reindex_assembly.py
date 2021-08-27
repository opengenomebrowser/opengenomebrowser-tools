import os


def reindex_assembly(file: str, out: str, prefix: str, leading_zeroes: int = None):
    """
    Change the header line of FASTA to: f'>{prefix}_{counter}'.

    :param file: input file
    :param out: output file
    :param prefix: desired prefix
    :param leading_zeroes: format counter with leading zeroes (optional). e.g.: 5 -> >PREFIX_00001
    """
    assert not os.path.isfile(out), f'Output file already exists! {out=}'

    if type(leading_zeroes) is int and leading_zeroes > 1:
        format = lambda c: f'>{prefix}{str(c).zfill(leading_zeroes)}\n'
    else:
        format = lambda c: f'>{prefix}{c}\n'

    counter = 0
    with open(file) as in_f, open(out, 'w') as out_f:
        for line in in_f:
            if line.startswith('>'):
                counter += 1
                out_f.write(format(counter))
            else:
                out_f.write(line)


def main():
    import fire

    fire.Fire(reindex_assembly)


if __name__ == '__main__':
    main()
