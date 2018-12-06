from setuptools import setup, find_namespace_packages
from pathlib import Path

HERE = Path(__file__).parent

script_names = [
    path.with_suffix("").name
    for path in (HERE / 'scanpy_scripts').glob('scanpy_*.py')
]
scripts = [
    f'{name.replace("_", "-")}=scanpy_scripts.{name}:main'
    for name in script_names
]

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='scanpy-scripts',
    version='0.1.0.2',
    author='nh3',
    author_email='',
    description='Scripts for using scanpy from the command line',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ebi-gene-expression-group/scanpy-scripts',
    packages=find_namespace_packages(),
    scripts=[
        'scanpy-scripts-tests.sh',
        'scanpy-scripts-tests.bats',
    ],
    entry_points=dict(
        console_scripts=scripts,
    ),
    install_requires=[
        'matplotlib',
        'pandas',
        'scanpy==1.3.2',
        'louvain'
    ],
)
