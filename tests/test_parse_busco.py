from unittest import TestCase
import os
import logging
from opengenomebrowser_tools.parse_busco import parse_busco

logging.basicConfig(level=logging.INFO)

ROOT = os.path.dirname(os.path.dirname(__file__))

buscos = [
    f'{ROOT}/test-data/prokka-good/short_summary.specific.lactobacillales_odb10.FAM3228-i1-1_busco.txt',
]


class Test(TestCase):
    def test_init_orthofinder(self):
        for busco in buscos:
            busco_dict = parse_busco(busco)
            self.assertEqual(busco_dict['C'] + busco_dict['M'] + busco_dict['F'], busco_dict['T'])
            self.assertEqual(busco_dict['S'] + busco_dict['D'], busco_dict['C'])
