"""
Provides version, author and exports
"""
import pkg_resources

__version__ = pkg_resources.get_distribution('scanpy-scripts').version

__author__ = ', '.join([
    'Ni Huang',
    'Pablo Moreno',
    'Jonathan Manning',
    'Philipp Angerer',
])

from . import lib
