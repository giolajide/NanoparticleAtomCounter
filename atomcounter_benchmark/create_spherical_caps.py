"""
Create spheres with their specified radii,
then cut them to yield spherical caps
of specified contact angles
"""

from ase.io import write
from os import system
import math
import numpy as np
from ase import Atoms
from ase.build import bulk
import warnings
from ase.data import chemical_symbols
from typing import Union, Optional, List, Tuple
from atomcounter_benchmark.atomistic_utils import (
    scaler,
    create_unit_support,
)

# from functools import lru_cache


# @lru_cache(maxsize=100)
def _calculate_bulk_density(
    element: str,
) -> float:
    """
    Calculate bulk density in atoms/A^3,
    given the element symbol
    """
    atoms = bulk(element)
    volume = atoms.get_volume()
    return len(atoms) / volume


def _atoms_to_radius(
    n_atoms: int,
    element: str,
) -> float:
    """
    Convert number of atoms to radius (of curvature)
    using the formula for a sphere
    returns radius in Ang
    """
    bulk_density = _calculate_bulk_density(
        element=element,
    )
    volume = n_atoms / bulk_density  # A^3
    return ((3 * volume) / (np.pi * 4)) ** (1 / 3)  # A


def _radius_to_atoms(
    radius: float,
    element: str,
) -> int:
    """
    Reverse of atoms_to_radius()
    """
    volume = (4 * np.pi * (radius**3)) / 3  # A^3
    bulk_density = _calculate_bulk_density(
        element=element,
    )  # atoms/A^3
    return volume * bulk_density  # n_atoms


def create_sphere(
    *,  # must use this function by specifying the keywords too, not just values
    element: str,
    n_atoms: Optional[int] = None,
    radius: Optional[float] = None,
) -> Atoms:
    """
    Create a spherical nanoparticle by finding the center of mass of a bulk structure
    then keeping atoms within a specified distance from that center of mass

    Requires:
        n_atoms (int):  required number of atoms in the cluster.
                        will try to fulfil this but may not
        radius (float): required radius of cluster, in Ang.
        element (str):  element type. Mandatory

    Returns:
        atoms_to_keep (Atoms): required spherical nanoparticle/cluster
    """
    element = element.capitalize()
    if radius <= 0.0:
        raise ValueError(f"Invalid radius {radius} supplied")
    if element not in chemical_symbols:
        raise ValueError(f"Invalid element {element} specified")
    if n_atoms is None and radius is None:
        raise ValueError("Provide at least one of radius and n_atoms")
    if n_atoms and radius:
        warnings.warn(
            f"Both radius and n_atoms given; will use only radius", category=UserWarning
        )
        n_atoms = None

    if n_atoms is not None:
        radius = _atoms_to_radius(
            n_atoms=n_atoms,
            element=element,
        )
    diameter = radius * 2
    ##build a bulk cell of safely larger dimensions than the diameter
    base_atoms = bulk(element)
    required_length = 2.5 * diameter
    current_length = min(base_atoms.cell.lengths())
    multiplying_factor = math.ceil(required_length / current_length)
    base_atoms *= (multiplying_factor, multiplying_factor, multiplying_factor)
    ##get the center of mass and draw a circle of the specified diameter around it
    com = base_atoms.get_center_of_mass()
    atoms_to_keep = []
    for index, atom in enumerate(base_atoms):
        position = atom.position
        distance = np.linalg.norm(position - com, ord=2)
        if distance <= radius:
            atoms_to_keep.append(atom)
    atoms_to_keep = Atoms(atoms_to_keep)

    return atoms_to_keep


def cut_particle(atoms: Atoms, angle: float) -> Atoms:
    """
    Cut symetric particle into an asymetric one of a specified contact angle
    """
    atoms = atoms.copy()
    positions = atoms.get_positions()
    z_positions = positions[:, 2]

    ##calculate the min and max position along the z-axis
    z_min = np.min(z_positions)
    z_max = np.max(z_positions)
    z_center = (z_min + z_max) / 2

    angle_radians = np.radians(angle)
    cos_angle = np.cos(
        angle_radians
    )  # calculate the normal vector of the cutting plane

    cut_plane_z = z_center + (z_max - z_center) * cos_angle  # define the cutting plane

    ##select atoms below the cutting plane
    mask = z_positions >= cut_plane_z
    cut_atoms = atoms[mask]

    return cut_atoms
