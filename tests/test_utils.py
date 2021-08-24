from unittest import TestCase

from opengenomebrowser_tools.utils import *

logging.basicConfig(level=logging.INFO)


class Test(TestCase):
    def test_entrez_organism_to_taxid(self):
        self.assertEqual(entrez_organism_to_taxid(organism='root'), 1)
        self.assertEqual(entrez_organism_to_taxid(organism='Mycoplasma genitalium'), 2097)
        with self.assertRaises(AssertionError):
            entrez_organism_to_taxid(organism='Organism fake-ensis')

    def test_create_replace_function(self):
        original = 'One two three four five six _six_ eight nine ten.'
        expected = 'Zero one two three four five six seven eight nine ten eleven.'

        replace_function = create_replace_function(replace_map={
            'One': 'Zero one',
            '_six_': 'seven',
            'ten': 'ten eleven',
        })

        result = replace_function(original)

        self.assertEqual(result, expected)
