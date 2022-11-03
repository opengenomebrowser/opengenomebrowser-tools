# based on https://realpython.com/pypi-publish-python-package
# How to upload:
#  - change package version in `setup.py` and `__init__.py`
#  - `python setup.py sdist`
#  - `twine upload dist/opengenomebrowser-tools-?.tar.gz`
import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / 'README.md').read_text()

# This call to setup() does all the work
setup(
    name='opengenomebrowser-tools',
    version='0.0.8',
    description='Set of scripts to aid OpenGenomeBrowser administrators import data',
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
    include_package_data=True,  # see MANIFEST.in
    install_requires=['schema', 'biopython', 'termcolor', 'fire', 'pyyaml'],
    entry_points={
        'console_scripts': [
            'init_folder_structure=opengenomebrowser_tools.init_folder_structure:main',
            'import_genome=opengenomebrowser_tools.import_genome:main',
            'import_genome2=opengenomebrowser_tools.import_genome2:main',
            'download_ncbi_genome=opengenomebrowser_tools.download_ncbi_genome:main',
            'genbank_to_fasta=opengenomebrowser_tools.genbank_to_fasta:main',
            'reindex_assembly=opengenomebrowser_tools.reindex_assembly:main',
            'rename_custom_annotations=opengenomebrowser_tools.rename_custom_annotations:main',
            'rename_eggnog=opengenomebrowser_tools.rename_eggnog:main',
            'rename_fasta=opengenomebrowser_tools.rename_fasta:main',
            'rename_genbank=opengenomebrowser_tools.rename_genbank:main',
            'rename_gff=opengenomebrowser_tools.rename_gff:main',
            'init_orthofinder=opengenomebrowser_tools.init_orthofinder:main',
            'import_orthofinder=opengenomebrowser_tools.import_orthofinder:main',
            'folder_looper=opengenomebrowser_tools.folder_looper:main',
            'update_folder_structure=opengenomebrowser_tools.update_folder_structure:main',
        ]
    },
)
