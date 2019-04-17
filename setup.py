from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='scanpy-scripts',
    version='0.2.0',
    author='nh3',
    author_email='nh3@users.noreply.github.com',
    description='Scripts for using scanpy from the command line',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ebi-gene-expression-group/scanpy-scripts',
    packages=find_packages(),
    scripts=[
        'scanpy-scripts-tests.sh',
        'scanpy-scripts-tests.bats',
    ],
    entry_points=dict(
        console_scripts='scanpy-cli=scanpy_scripts.cli:cli',
    ),
    install_requires=[
        'matplotlib',
        'pandas',
        'scanpy>=1.4.0',
        'louvain',
        'leidenalg',
        'MulticoreTSNE',
        'Click'
    ],
)
