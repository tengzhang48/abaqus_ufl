"""Make the in-tree package importable under the ``pytest`` console script.

``python -m pytest`` inserts the current directory into ``sys.path``;
the bare ``pytest`` entry point does not, so an uninstalled checkout
would fail to import ``abaqus_ufl``. Tests must exercise THIS checkout
even when another copy of the package is installed.
"""
import sys
from pathlib import Path

ROOT = str(Path(__file__).resolve().parent)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
