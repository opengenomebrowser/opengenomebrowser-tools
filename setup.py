# based on https://realpython.com/pypi-publish-python-package
import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / 'README.md').read_text()

# This call to setup() does all the work
setup(
    name='opengenomebrowser-tools',
    version='0.0.1',
    description=' Set of scripts to aid OpenGenomeBrowser administrators import data',
    long_description=README,
    long_description_content_type='text/markdown',
    url='https://github.com/opengenomebrowser/opengenomebrowser-tools',
    author='Thomas Roder',
    author_email='roder.thomas@gmail.com',
    license='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ],
    packages=['opengenomebrowser_tools'],
    include_package_data=True,
    install_requires=['schema', 'biopython', 'termcolor', 'fire', 'pyyaml'],
    entry_points={
        'console_scripts': [
            'import_genome=opengenomebrowser_tools.import_genome:main',
            'gbk_to_ffn=opengenomebrowser_tools.gbk_to_ffn:main',
            'reindex_assembly=opengenomebrowser_tools.reindex_assembly:main',
            'rename_custom_annotations=opengenomebrowser_tools.rename_custom_annotations:main',
            'rename_eggnog=opengenomebrowser_tools.rename_eggnog:main',
            'rename_fasta=opengenomebrowser_tools.rename_fasta:main',
            'rename_genbank=opengenomebrowser_tools.rename_genbank:main',
            'rename_gff=opengenomebrowser_tools.rename_gff:main',
        ]
    },
)
