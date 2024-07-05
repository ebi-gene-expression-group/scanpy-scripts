"""
Provides version, author and exports
"""
import importlib.metadata

__version__ = importlib.metadata.version("scanpy-scripts")

__author__ = ", ".join(
    [
        "Ni Huang",
        "Pablo Moreno",
        "Jonathan Manning",
        "Philipp Angerer",
    ]
)

from . import lib
