from setuptools import setup
from pathlib import Path

HERE = Path(__file__).parent

script_files = (HERE / 'scanpy_cli').glob('scanpy_*.py')
scripts = [
    f'{file.with_suffix("").name.replace("_", "-")}=scanpy_cli.{file.name}:main'
    for file in script_files
]

setup(
    name='scanpy-cli',
    version='1.0',
    description='Scripts for using scanpy from the command line',
    author='nh3',
    #author_email='',
    url='https://github.com/ebi-gene-expression-group/scanpy-scripts',
    packages=['scanpy_cli'],
    scripts=[
        'scanpy-cli-tests.sh',
        'scanpy-cli-tests.bats',
    ],
    entry_points=dict( 
        console_scripts=scripts,
    ),
    install_requires=[
        'matplotlib',
        'pandas',
        'scanpy',
    ],
)
