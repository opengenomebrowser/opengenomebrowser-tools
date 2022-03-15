from re import compile
from datetime import datetime

extr_digits = compile('[0-9]+')


def parse_busco(file: str) -> dict:
    entries = [
        ('Complete BUSCOs (C)', 'C'),
        ('Complete and single-copy BUSCOs (S)', 'S'),
        ('Complete and duplicated BUSCOs (D)', 'D'),
        ('Fragmented BUSCOs (F)', 'F'),
        ('Missing BUSCOs (M)', 'M'),
        ('Total BUSCO groups searched', 'T')
    ]

    busco_dict = {"C": None, "D": None, "F": None, "M": None, "S": None, "T": None}

    error_msg = F'Error while parsing BUSCO: {file} - {busco_dict}'

    with open(file) as f:
        busco_file = f.readlines()

    for line in busco_file:
        for descr, abbr in entries:
            if descr in line:
                x = extr_digits.findall(line)
                assert len(x) == 1, error_msg
                busco_dict[abbr] = int(x[0])

    assert len(busco_dict.keys()) == 6 and None not in busco_dict.values(), error_msg
    assert busco_dict['C'] + busco_dict['M'] + busco_dict['F'] == busco_dict['T'], error_msg
    assert busco_dict['S'] + busco_dict['D'] == busco_dict['C'], error_msg

    dataset_line = [l for l in busco_file if l.startswith('# The lineage dataset is:')][0]
    dataset, rest = dataset_line[26:].split(' (Creation date: ', maxsplit=1)
    dataset_date = rest[:10]
    try:
        datetime.strptime(dataset_date, '%Y-%m-%d')
    except Exception as e:
        raise AssertionError(F'Bad busco creation date ! {file}n{str(e)}')

    busco_dict['dataset'] = dataset
    busco_dict['dataset_creation_date'] = dataset_date

    return busco_dict


def main():
    import fire
    import json

    def runner(file: str):
        print(json.dumps(parse_busco(file), indent=4))

    fire.Fire(runner)


if __name__ == '__main__':
    main()
