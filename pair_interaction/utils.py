"""Utility geometry routines for linker‑mediated interactions."""

from __future__ import annotations
import numpy as np


def shell_volume(radius: float, rc: float) -> float:
    """Volume (nm³) of a spherical shell of thickness *rc* around a sphere.

    Parameters
    ----------
    radius : float
        Core radius of the particle (nm).
    rc : float
        Contour length of a grafted linker; effectively the shell thickness (nm).

    Returns
    -------
    float
        Volume of the shell region in nm³.
    """
    return 4.0 / 3.0 * np.pi * ((radius + rc) ** 3 - radius ** 3)


def overlap_sphere(r1: float, r2: float, d: float) -> float:
    """Intersection volume of two spheres.

    Parameters
    ----------
    r1, r2 : float
        Radii of the two spheres (nm).
    d : float
        Centre-to-centre separation between the spheres (nm).

    Returns
    -------
    float
        Overlap volume in nm³. Zero when the spheres are disjoint.
    """
    if d >= r1 + r2:
        return 0.0

    return (
        np.pi
        * (r1 + r2 - d) ** 2
        * (d ** 2 + 2.0 * d * (r1 + r2) - 3.0 * (r1 - r2) ** 2)
        / (12.0 * d)
    )
