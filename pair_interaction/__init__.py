"""Top‑level package initialisation."""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

from .pair_interaction import PairInteraction

__all__ = ["PairInteraction"]

try:
    __version__ = _pkg_version(__name__)
except PackageNotFoundError:
    from ._version import __version__


# File: pae_ee_interaction/_version.py
"""Fallback static version (used in editable installs)."""

__version__ = "0.1.0"