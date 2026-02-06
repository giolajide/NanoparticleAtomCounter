"""
Supports the discrimination script
"""

from ase.io import read, write, iread
import numpy as np
from tqdm import tqdm
from ase import Atoms, Atom
import math
from typing import Tuple, Literal, List, Union
from ase.constraints import FixAtoms
from ase.neighborlist import NeighborList, natural_cutoffs
import numpy as np
from ase.geometry.analysis import Analysis

LATERAL_SPACING = (
    14.5  # spacing between between each nanoparticles repeated image. Note, though that
)
# this is a little lower than the spacing, due to the ceil function being used
Z_SPACING = 13.5  # spacing in z direction
ADSORPTION_HEIGHT = 2.2  # adsorption height
INTERFACE_SCALING = 1.3  # scaling factor for interfacial NP atoms' covalent radii
NANOPARTICLE_ELEMENT = "Ag"  # what kind of atom is the NP made of?
SCALING_FACTOR = 1.0
SKIN = 0.3  # ASE's default value=0.3A


def create_unit_support(support: Literal["graphene", "au", "mgo"]) -> Atoms:
    """Creates a unit cell of the support: graphene, mgo, or gold"""
    support = support.lower()

    graphene = Atoms(
        "C4",
        pbc=np.array([True, True, True]),
        cell=np.array(
            [
                [4.2687822100000004, 0.0000000000000000, 0.0000000000000000],
                [0.0000000000000000, 2.4645825600000002, 0.0000000000000000],
                [0.0000000000000000, 0.0000000000000000, 9.5824915799999992],
            ]
        ),
        positions=np.array(
            [
                [1.4229274000000001, 0.0000000000000000, 0.0000000000000000],
                [2.8458548100000001, 0.0000000000000000, 0.0000000000000000],
                [0.7114637000000000, 1.2322912800000001, 0.0000000000000000],
                [3.5573185100000000, 1.2322912800000001, 0.0000000000000000],
            ]
        ),
    )

    gold = Atoms(
        "Au4",
        pbc=np.array([True, True, True]),
        cell=np.array(
            [
                [2.8849956672411143, 0.0000000000000000, 0.0000000000000000],
                [0.0000000000000000, 4.9969590752776831, 0.0000000000000000],
                [0.0000000000000000, 0.0000000000000000, 6.3555890982936738],
            ]
        ),
        positions=np.array(
            [
                [0.0000000000000000, 1.6656530250925610, 2.0000000000000000],
                [1.4424978336205572, 4.1641325627314023, 2.0000000000000000],
                [0.0000000000000000, 0.0000000000000000, 4.3555890982936738],
                [1.4424978336205572, 2.4984795376388416, 4.3555890982936738],
            ]
        ),
    )

    mg = Atoms(
        "Mg4",
        pbc=graphene.pbc,
        cell=np.array(
            [
                [4.2602238819630642, 0.0000000000000000, 0.0000000000000000],
                [0.0000000000000003, 4.2602238819630642, 0.0000000000000000],
                [0.0000000000000000, 0.0000000000000000, 24.9107835868707248],
            ]
        ),
        positions=np.array(
            [
                [0.0000000000000021, 0.0000000000000009, 17.7806716458891927],
                [2.1301119409815343, 2.1301119409815330, 17.7806716458891927],
                [0.0000000000000025, 2.1301119409815330, 19.9107835868707248],
                [2.1301119409815348, 0.0000000000000011, 19.9107835868707248],
            ]
        ),
    )

    o = Atoms(
        "O4",
        pbc=graphene.pbc,
        cell=mg.cell,
        positions=np.array(
            [
                [0.0000000000000023, 2.1301119409815330, 17.7806716458891927],
                [2.1301119409815343, 0.0000000000000009, 17.7806716458891927],
                [0.0000000000000025, 0.0000000000000011, 19.9107835868707248],
                [2.1301119409815348, 2.1301119409815334, 19.9107835868707248],
            ]
        ),
    )

    mg.extend(o)

    requested = {"graphene": graphene, "au": gold, "mgo": mg}
    unit_support = requested[support]

    return unit_support


def setup_neighborlist(
    atoms: Atoms,
    scaling_factor: float = SCALING_FACTOR,
    radial_cutoff: float = None,
    skin: float = SKIN,
):
    """
    Setup neighborlist object,
    given the Atoms object and the scaling factor for covalent radii

    if radial_cutoff is given, then we will consider only neighbors within
    the lower of that cutoff and the bond length
    """
    if not radial_cutoff:
        cutoffs = natural_cutoffs(atoms, mult=scaling_factor)
    else:
        cutoffs = np.minimum(radial_cutoff / 2, np.array(natural_cutoffs(atoms)))

    neighbors = NeighborList(
        cutoffs=cutoffs, self_interaction=False, bothways=True, skin=skin
    )
    neighbors.update(atoms)
    return neighbors


def setup_analyzer(atoms: Atoms, neighborlist: NeighborList = None) -> Analysis:
    """
    Setup the analyzer object, given an optional NeighborList object
    If not given a NeighborList object, use defaults
    """
    analyzer = Analysis(atoms, nl=neighborlist)
    return analyzer


def divider(
    atoms: Atoms, element: str = NANOPARTICLE_ELEMENT
) -> Union[None, Tuple[Atoms]]:
    """
    Divide the structure into NP and surface
    Preserve the constraints of the support; assume the NP has none

    Inputs:
                Supported_NP
                element         (what element is the NP?)
                                default = 'Ag'
    Returns:    NP
                Bare_Surface
    """
    if not isinstance(atoms, Atoms):
        return None

    element = element.capitalize()

    constrained_indices = (
        atoms.constraints[0].get_indices() if atoms.constraints else list()
    )
    silver_indices = [
        index for index, atom in enumerate(atoms) if atom.symbol == element
    ]
    support_indices = [
        index for index, atom in enumerate(atoms) if atom.symbol != element
    ]

    silvers, support = atoms[silver_indices], atoms[support_indices]
    silvers.set_cell(atoms.get_cell())
    support.set_cell(atoms.get_cell())
    silvers.pbc, support.pbc = True, True

    old_to_new_support_indices = dict()
    for new_index, old_index in enumerate(support_indices):
        old_to_new_support_indices[old_index] = new_index

    new_constrained_support_indices = [
        old_to_new_support_indices[old_i]
        for old_i in constrained_indices
        if old_i in old_to_new_support_indices
    ]

    if new_constrained_support_indices:
        support.set_constraint(FixAtoms(indices=new_constrained_support_indices))

    return silvers, support


def centralize(atoms: Atoms, element: str = NANOPARTICLE_ELEMENT) -> Atoms:
    """
    Centralize the NP to make viewing easier
    Inputs:     Supported_NP
                element the NP is made of
    Returns:    Supported_NP (with NP centralized)
    """
    silvers, support = divider(atoms, element=element)

    center_x = support.cell.lengths()[0] / 2
    center_y = support.cell.lengths()[1] / 2

    silvers_com_x = silvers.get_center_of_mass()[0]
    silvers_com_y = silvers.get_center_of_mass()[1]

    x_move = center_x - silvers_com_x
    y_move = center_y - silvers_com_y
    silvers.translate(displacement=(x_move, y_move, 0))

    support.extend(silvers)
    support.pbc = atoms.pbc
    support.set_constraint(atoms.constraints)

    support.pbc = True

    return support


def scaler(
    image: Atoms,
    adsorption_height: float = ADSORPTION_HEIGHT,
    z_spacing: float = Z_SPACING,
    lateral_spacing: float = LATERAL_SPACING,
    unit_support: Union[Atoms, None] = None,
) -> Atoms:
    """
    Return an NP supported upon a support of a good size
    Inputs:
        Image:              NP or Supported NP
        element:            What element the NP is made of. default = 'Ag'
        adsorption_height:  Desired adsorption height (Optional)
        z_spacing:          Desired Z-spacing of periodic images (Optional)
        lateral_spacing:    Desired X- and Y-spacing of periodic NP images (Optional)
        unit_support:       Atoms object of the unit support.
    Returns:
                            Supported_NP
    """
    if not unit_support:
        unit_support = create_unit_support(support="mgo")
    unit_cell = unit_support.cell
    unit_cell_x = unit_cell[0, 0]
    unit_cell_y = unit_cell[1, 1]
    unit_cell_z = unit_cell[2, 2]
    unit_cell_max_z = max(unit_support.positions[:, 2])

    silvers, support = divider(image, element=image.get_chemical_symbols()[0])
    silvers = Atoms(silvers)
    silvers.center(vacuum=10)

    min_x = min(silvers.positions[:, 0])
    min_y = min(silvers.positions[:, 1])
    max_x = max(silvers.positions[:, 0])
    max_y = max(silvers.positions[:, 1])
    min_z = min(silvers.positions[:, 2])

    x_diameter = max_x - min_x
    y_diameter = max_y - min_y

    required_x = x_diameter + lateral_spacing
    required_y = y_diameter + lateral_spacing

    ratio_x = math.ceil(required_x / unit_cell_x)
    ratio_y = math.ceil(required_y / unit_cell_y)

    adsorption_height = adsorption_height + unit_cell_max_z
    new_support = unit_support * (ratio_x, ratio_y, 1)
    silvers_displacement = adsorption_height - min_z
    silvers.translate(displacement=(0, 0, silvers_displacement))
    # I think the logic here is tedious, but never mind
    new_support.cell[2, 2] = (
        10  # deliberately set too low so that the logic within the if statement (below) will play out fine
    )
    new_support.extend(silvers)

    cell_bottom = 0
    min_z = min(new_support.positions[:, 2])
    new_cell_displacement = cell_bottom - min_z
    new_support.translate(displacement=(0, 0, new_cell_displacement))

    max_height = max(new_support.positions[:, 2])
    cell_top = new_support.cell[2, 2]
    distance = cell_top - max_height

    if distance < z_spacing:
        new_support.cell[2, 2] += lateral_spacing - cell_top + max_height

    new_support = centralize(new_support)
    new_support.info.update(image.info)

    return new_support
