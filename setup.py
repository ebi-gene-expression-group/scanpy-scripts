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
        'scanpy-scripts-tests.bats',
    ],
    entry_points=dict(
        console_scripts=[
            'scanpy-cli=scanpy_scripts.cli:cli',
            'scanpy-read-10x=scanpy_scripts.cmds:READ_CMD',
            'scanpy-filter-cells=scanpy_scripts.cmds:FILTER_CMD',
            'scanpy-filter-genes=scanpy_scripts.cmds:FILTER_CMD',
            'scanpy-normalise-data=scanpy_scripts.cmds:NORM_CMD',
            'scanpy-find-variable-genes=scanpy_scripts.cmds:HVG_CMD',
            'scanpy-scale-data=scanpy_scripts.cmds:SCALE_CMD',
            'scanpy-regress=scanpy_scripts.cmds:REGRESS_CMD',
            'scanpy-run-pca=scanpy_scripts.cmds:PCA_CMD',
            'scanpy-neighbors=scanpy_scripts.cmds:NEIGHBOR_CMD',
            'scanpy-run-tsne=scanpy_scripts.cmds:TSNE_CMD',
            'scanpy-run-umap=scanpy_scripts.cmds:UMAP_CMD',
            'scanpy-find-cluster=scanpy_scripts.cli:cluster',
            'scanpy-find-markers=scanpy_scripts.cmds:DIFFEXP_CMD',
        ]
    ),
    install_requires=[
        'scipy>=1.2.0,<1.3.0',
        'matplotlib',
        'pandas',
        'scanpy>=1.4.2',
        'louvain',
        'leidenalg',
        'loompy',
        'MulticoreTSNE',
        'Click'
    ],
)
