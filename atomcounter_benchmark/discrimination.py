"""
In memory of Korosensei

This is the atomistic method

Use:                    Discrimination of kinds of atoms in supported NPs
                        The interface absolutely MUST be flat or nearly so

Functions:              get_interface
                        get_np_surface_by_CN
                        get_perimeter
                        discriminate


To remove complications regarding alphashape and its very specific set of dependencies
(and on how slow optimizing alpha is for large nanoparticles)
I have eliminated it, commenting it out entirely

Just to note it here, in case a future version has me putting back alphashape, here's the process:
    pip uninstall alphashape
    pip install git+https://github.com/bellockk/alphashape.git
    pip install trimesh==4.2.0 numpy==1.26.4

"""

import pandas as pd
from argparse import ArgumentParser
from copy import deepcopy
from ase.build import bulk
from tqdm import tqdm
import statistics
from ase.io import read, write
import numpy as np
from ase import Atoms, Atom
import math
from sys import argv, exit

# import alphashape
# from shapely.geometry import Polygon, MultiPolygon
from ase.geometry.analysis import Analysis
from typing import Tuple, Literal, List, Union
from ase.neighborlist import natural_cutoffs, NeighborList
import warnings
from os import system, environ
from atomcounter_benchmark.atomistic_utils import (
    divider,
    setup_neighborlist,
    setup_analyzer,
    NANOPARTICLE_ELEMENT,  # SKIN
)
from joblib import Parallel, delayed
from ascii_colors import ASCIIColors


INTERFACE_SCALING = 1.3  # 1.3 scaling factor for interfacial NP atoms' covalent radii
CN_SCALING = 1.05  # 1.2 scaling factor for discriminating bulk from NP surface atoms
FCC_AND_HCP_COORD_CUTOFF = 12  # any atom with coord < COORD_CUTOFF is a surface atom
BCC_COORD_CUTOFF = 8  # any atom with coord < COORD_CUTOFF is a surface atom
INTERFACIAL_FACET = (1, 1, 1)  # facet of NP that faces the support;
# also that which is on the outer surface of the NP
# which means this is a misnomer
DO_ALPHA_SHAPE = False  # don't use alpha shape to get perimeter atoms
PROCESSES = -1  # How many processes to run in parallel.
NL_SKIN = 0.05  # Ang; nl skin for distinguishing surface NP atoms;
# set low because we normally use this on perfect crystals,
# not structures from MD; this value seems to work well


def get_interface(
    atoms: Atoms,
    nl: NeighborList = None,
    analyzer: Analysis = None,
    support_elements: Union[str, List[str]] = None,
    np_element: str = NANOPARTICLE_ELEMENT,
    scaling_factor: float = INTERFACE_SCALING,
) -> Union[List[None], Tuple[Atoms, Atoms, List[int], List[int]]]:
    """
    Get indices of interfacial atoms between NP and Support

    Requires:
        atoms (Atoms):                      atoms object of NP + Support
        nl (NeighborList):                  neighborlist object of whole system (i.e. NP + Support).
                                            Optional but strongly recommended to provide analyzer or nl
        analyzer (Analysis):                Analysis object  of whole system (i.e. NP + Support).
                                            Optional but strongly recommended to provide analyzer or nl
        support_elements (Union[str, List[str]]):
                                            what element(s) the support is made of. If not supplied, we
                                            say it is whatever np_element is not
        np_element (str):                   what element the NP is made of. default = 'Ag'
        scaling_factor (float):             scaling factor for covalent radii. Optional
                                            Advised to be > 1.0 to fully capture the interface
    Returns:
        np_interfacial_indices (List[int]):         indices of interfacial (+ perimeter) NP atoms
        support_interfacial_indices (List[int]):    indices of interfacial support atoms

    If no NP-Support bonds found,
        Returns:    None, None
    """
    elements = np.unique(atoms.get_chemical_symbols())
    np_element = np_element.capitalize()
    if not support_elements:
        support_elements = [i for i in elements if i != np_element]
    if isinstance(support_elements, str):
        support_elements = [support_elements]

    if not analyzer:
        if not nl:  # reuse neighborlist if given
            natural_cutoff = (
                np.max([natural_cutoffs(Atoms(element)) for element in elements]) * 2
            )  # no need to search beyond this distance
            nl = setup_neighborlist(atoms, scaling_factor, radial_cutoff=natural_cutoff)
            nl.update(atoms)

        analyzer = setup_analyzer(
            atoms, neighborlist=nl
        )  # reuse Analysis object if given

    NP_Support_bonds = {
        e.capitalize(): analyzer.get_bonds(e.capitalize(), np_element, unique=True)[0]
        for e in support_elements
    }
    support_interfacial_indices, np_interfacial_indices = [], []
    for support_e, bonds in NP_Support_bonds.items():
        if bonds:
            np_interfacial_indice, support_interface_indice = zip(
                *[(j, i) for (i, j) in bonds]
            )
            np_interfacial_indice = list(np_interfacial_indice)
            support_interface_indice = list(support_interface_indice)
        else:
            warnings.warn(
                f"""There are no {support_e}-{np_element} bonds in your system!
            This may be a very small cluster""",
                category=RuntimeWarning,
            )
            np_interfacial_indice = list()
            support_interface_indice = list()

        support_interfacial_indices.extend(support_interface_indice)
        np_interfacial_indices.extend(np_interfacial_indice)

    if not (support_interfacial_indices and np_interfacial_indices):
        print(f"No NP-Support bonds found! Failed!")
        return None, None

    support_interfacial_indices = list(np.unique(support_interfacial_indices))
    np_interfacial_indices = list(np.unique(np_interfacial_indices))

    # CAUTION!! np_interfacial_indices includes perimeter atoms.
    # To separate out the perimeter, use get_perimeter on results from this function
    return np_interfacial_indices, support_interfacial_indices


def get_perimeter(
    atoms: Atoms,
    interfacial_results: [List[int], List[int]],
    do_alphashape: bool = DO_ALPHA_SHAPE,
) -> Union[List[None], Tuple[Atoms, Atoms, List[int], List[int], float]]:
    f"""
    Get indices of perimeter atoms of the NP
    By:
        1. Getting the interface
        2. Convex/concave hull of interfacial atoms to trace the perimeter (misses concave regions);
        from experience, the bigger the NP, the more it'll undercount the perimeter, especially the convex
        hull.
        3. To correct (2), also gets CNmax (the most-coordinated atom amongst those flagged by (2)
        as perimeter) minus 1, and defines that any atom of coordination less than the max of CNmax and
        CNmedian (median coordination of all atoms flagged by (2) as perimeter).
        i.e.:
            threshold_CN = max(CNmedian, CNmax - 1)

        Note: simply setting threshold_CN = CNmax would be good but a little risky in rare cases

    Requires:
        atoms (Atoms):                      atoms object of NP + Support
        interfacial_results: Tuple[List[int], List[int]]. Required; results of the get_interface function,
        do_alphashape (bool):               whether or not to get perimeter atoms by alphashape rather than only by CN
                                            default = False
                                            It's very expensive to optimize alpha for large NPs

    Returns:
        perimeter_indices (List[int])      indices of perimeter atoms

    If get_interface function failed,
        Returns:    None

    Warning:
        This threshold_CN method will only work when the atomic density
        at the interface is homogeneous
    """
    np_interfacial_indices, support_interfacial_indices = interfacial_results
    np_interface = atoms[np_interfacial_indices]  # NPinterface (i.e. no support)

    # get perimeter atoms
    perimeter_indices = list()

    # criterion 1: concave hull vertices
    # August 1st, 2025: hardcoded to never run
    if False:  # do_alphashape:
        positions_2D = np_interface.positions[:, :2]
        optimum_alpha = alphashape.optimizealpha(positions_2D)
        print(f"Optimal alpha = {optimum_alpha}")
        alpha_shape = alphashape.alphashape(
            positions_2D, optimum_alpha * 0.88
        )  # scale down optimum_alpha because it'll tend to overcount
        perim_set = {tuple(pt) for pt in alpha_shape.exterior.coords[:-1]}
        perimeter_indices = [
            i for (i, p) in enumerate(positions_2D) if tuple(p) in perim_set
        ]

    # criterion 2: coordination threshold
    ##we have to create a new nl because we are now using a subset of the provided Atoms object
    neighbors = NeighborList(
        cutoffs=natural_cutoffs(np_interface, mult=INTERFACE_SCALING),
        self_interaction=False,
        bothways=True,
    )
    neighbors.update(np_interface)
    CNs = [
        (atom_index, len(neighbors.get_neighbors(atom_index)[0]))
        for (atom_index, _) in enumerate(np_interface)
    ]  # get number of neighbors for each interfacial NP atom
    # remember: only interface-interface bonds are considered
    neighbor_counts = [num_neighbors for (_, num_neighbors) in CNs]
    medianCN, maxCN = np.median(neighbor_counts), max(neighbor_counts)
    thresholdCN = max(medianCN, maxCN)  # - 1)
    missed_perimeter_indices = [
        atom_index
        for (atom_index, num_neighbors) in CNs
        if (num_neighbors < thresholdCN and atom_index not in perimeter_indices)
    ]
    perimeter_indices.extend(missed_perimeter_indices)

    # convert the perimeter indices to the indices in the Atoms object
    perimeter_indices = [np_interfacial_indices[i] for i in perimeter_indices]

    return perimeter_indices


def get_np_surface_by_CN(
    atoms: Atoms,
    nl: NeighborList = None,
    support_elements: Union[str, List[str]] = None,
    np_element: str = NANOPARTICLE_ELEMENT,
    scaling_factor: float = CN_SCALING,
    coord_cutoff: int = None,
) -> Tuple[List[int], List[int]]:
    """
    Separate NP surface's from bulk based on CN (coordination numbers)

    Requires:
        atoms (Atoms):                      atoms object of NP + Support
        nl (NeighborList):                  neighborlist object of whole system (i.e. NP + Support).
                                            Optional but strongly recommended to provide analyzer or nl
        support_elements (Union[str, List[str]]):
                                            what element(s) the support is made of.
        np_element (str):                   what element the NP is made of. Default = 'Ag'
        scaling_factor (float):             scaling factor for covalent radii.
                                            optional. default = 1.05
        coord_cutoff (int):                 coordination cutoff for discriminating surface from bulk atoms
                                            optional. default = 12 for FCC and 8 for BCC
    Returns
        surface_plus_interface_indices (List[int]):
                                            indices of nanoparticle surface atoms.
                                            Note: This INCLUDES THE INTERFACE + PERIMETER ATOMS
    """
    if not support_elements:
        support_elements = [i for i in elements if i != np_element]
    if isinstance(support_elements, str):
        support_elements = [support_elements]
    if not nl:
        # create neighborlist object
        cutoff_kwargs = {element.capitalize(): 0 for element in support_elements}
        scaling_factor = 1.0
        cutoffs = natural_cutoffs(
            atoms, mult=scaling_factor, **cutoff_kwargs
        )  # set non NP atoms to zero cutoffs
        nl = NeighborList(
            cutoffs=cutoffs, self_interaction=False, bothways=True, skin=NL_SKIN
        )
        nl.update(atoms)

    # apply criterion
    if not coord_cutoff:
        lattice = bulk(np_element).cell.get_bravais_lattice().__class__.__name__
        if lattice in ["FCC", "HCP", "HEX"]:  # should confirm later, but it seems ASE
            # has only hex not hcp (they mean the same, I suppose)
            coord_cutoff = FCC_AND_HCP_COORD_CUTOFF
        elif lattice == "BCC":
            coord_cutoff = BCC_COORD_CUTOFF
        else:
            warnings.warn(
                f"""NP element's crystal system is neither fcc, hcp, nor bcc!
            Setting coordination cutoff to that of fcc/hcp ({FCC_AND_HCP_COORD_CUTOFF})!
            """,
                category=UserWarning,
            )

    CN = [
        (index, len(nl.get_neighbors(index)[0]))
        for index, i in enumerate(atoms)
        if i.symbol == np_element
    ]
    bulk_indices = [i for i, j in CN if j == coord_cutoff]
    surface_plus_interface_indices = [i for i, j in CN if i not in bulk_indices]

    return surface_plus_interface_indices


def discriminate(
    atoms: Atoms,
    nl: NeighborList,
    analyzer: Analysis,
    support_elements: Union[str, List[str]],
    np_element: str = NANOPARTICLE_ELEMENT,
    interface_scaling: float = INTERFACE_SCALING,
    surface_scaling: float = CN_SCALING,
    coord_cutoff: int = None,
    do_alphashape: bool = DO_ALPHA_SHAPE,
) -> Tuple[
    List[int],
    List[int],
    List[int],
    List[int],
    List[int],
]:
    """Main function for separating a supported NP into NP surface, NP Bulk, NP interface (non-perimeter),
    NP perimeter, and Support

    Requires:
        atoms (Atoms):                      atoms object of NP + Support
        nl (NeighborList):                  neighborlist object of whole system (i.e. NP + Support).
                                            Optional but strongly recommended to provide analyzer or nl
        analyzer (Analysis):                Analysis object  of whole system (i.e. NP + Support).
                                            Optional but strongly recommended to provide analyzer or nl
        support_elements (Union[str, List[str]]):
                                            what element(s) the support is made of.
        np_element (str):                   what element the NP is made of. default = 'Ag'
        interface_scaling (float):          scaling factor for covalent radii to discriminate interface from rest of system.
                                            default = 1.3
        surface_scaling (float):             scaling factor for covalent radii. Default = 1.05
        coord_cutoff (int):                 coordination cutoff for discriminating surface from bulk atoms
                                            optional. default = 12 for FCC and 8 for BCC
        do_alphashape (bool):               whether or not to get perimeter atoms by alphashape rather than only by CN
                                            default = False
                                            !!!It's very expensive to optimize alpha for large NPs!!!
                                            So, I hope you know what you are doing!

    Returns:
        np_surface (List[int]):             indices of NP surface atoms (excluding interface)
        np_bulk (List[int]):                indices of bulk NP atoms
        np_interface (List[int]):           indices of NP interfacial atoms (excluding perimeter)
        np_perimeter (List[int]):           indices of NP perimeter atoms
        substrate_interface (List[int]):     indices of support interfacial atoms
        substrate (List[int]):              indices of support (including interfacial support) atoms

    Returns [0,0,0,0,0,0] if get_interface didnt get any NP-support bonds
    """
    do_alphashape = False  # hardcoded to False (August 2025)

    elements = np.unique(atoms.get_chemical_symbols())
    np_element = np_element.capitalize()
    if not support_elements:
        support_elements = [i for i in elements if i != np_element]
    if isinstance(support_elements, str):
        support_elements = [support_elements.capitalize()]
    else:
        support_elements = [i.capitalize() for i in support_elements]
    np_and_support_elements = support_elements + [np_element]
    check_elements = [i in elements for i in np_and_support_elements]
    if not np.all(check_elements):
        raise ValueError(
            f"Some of the given elements ({np_and_support_elements}) not in Atoms object!"
        )

    if not analyzer:
        if not nl:  # reuse neighborlist if given
            natural_cutoff = (
                np.max([natural_cutoffs(Atoms(element)) for element in elements]) * 2
            )  # no need to search beyond this distance
            nl = setup_neighborlist(
                atoms, scaling_factor=interface_scaling, radial_cutoff=natural_cutoff
            )  # , skin = NL_SKIN) adding this ruins things!
            nl.update(atoms)

        analyzer = setup_analyzer(
            atoms, neighborlist=nl
        )  # reuse Analysis object if given

    nps = [index for index, _ in enumerate(atoms) if atoms[index].symbol == np_element]
    substrate = [
        index
        for index, _ in enumerate(atoms)
        if atoms[index].symbol in support_elements
    ]

    try:
        np_interface_plus_perimeter, substrate_interface = get_interface(
            atoms,
            nl=nl,
            analyzer=analyzer,
            np_element=np_element,
            support_elements=support_elements,
            scaling_factor=interface_scaling,
        )

        np_perimeter = get_perimeter(
            atoms,
            interfacial_results=(np_interface_plus_perimeter, substrate_interface),
            do_alphashape=do_alphashape,
        )

        np_interface = [i for i in np_interface_plus_perimeter if i not in np_perimeter]
        nl = None  # unfortunately we have to reset in order to set skin to zero for the surface discriminaion
        np_surface_plus_interface = get_np_surface_by_CN(
            atoms,
            nl=nl,
            np_element=np_element,
            support_elements=support_elements,
            scaling_factor=surface_scaling,
            coord_cutoff=coord_cutoff,
        )

        np_surface = [
            i for i in np_surface_plus_interface if i not in np_interface_plus_perimeter
        ]
        np_bulk = [i for i in nps if i not in np_surface_plus_interface]

    except TypeError:  # get_interface returned None, None
        ASCIIColors.red("get_interface probably failed")
        (
            np_surface,
            np_bulk,
            np_interface,
            np_perimeter,
            substrate_interface,
            substrate,
        ) = [0] * 6
    except Exception as e:
        ASCIIColors.red(f"Unexpected error occurred:\t{e}")
        exit(1)

    return (
        np_surface,
        np_bulk,
        np_interface,
        np_perimeter,
        substrate_interface,
        substrate,
    )


if __name__ == "__main__":
    parser = ArgumentParser(
        description="""Splits up a supported NP into surface, substrate, bulk, perimeter, and interface
    Also calculates interfacial radius by any/both of two methods if requested"""
    )

    parser.add_argument(
        "--traj",
        "-t",
        type=str,
        required=True,
        help="Traj file on which to run the script. Required",
    )
    parser.add_argument(
        "--traj_output",
        "-to",
        type=str,
        required=True,
        help="Traj file to create in which atoms have been discriminated. Required",
    )
    parser.add_argument(
        "--processes",
        "-p",
        type=int,
        default=PROCESSES,
        help=f"How many processes to run in parallel. default = {PROCESSES}",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        default="output_atomistic.csv",
        help="Output filename.",
    )
    parser.add_argument(
        "--np_element",
        "-ne",
        type=str,
        default=None,  # NANOPARTICLE_ELEMENT,
        help=f"Element of which the NP is composed. by default we'll read atoms.info['np_element']",
    )

    parser.add_argument(
        "--support_elements",
        "-se",
        type=str,
        nargs="+",
        default=None,
        help="""List of elements of which the substrate is composed.
    In the format: a b c d. Defaults to elements absent from the NP""",
    )

    parser.add_argument(
        "--interface_scaling",
        "-is",
        type=float,
        default=INTERFACE_SCALING,
        help=f"""scaling factor for interfacial atoms' covalent radii.
            Default = {INTERFACE_SCALING}""",
    )
    parser.add_argument(
        "--surface_scaling",
        "-ss",
        type=float,
        default=CN_SCALING,
        help=f"""scaling factor for discriminating NP surface from bulk.
            Default = {CN_SCALING}""",
    )

    do_alphashape = False

    args = parser.parse_args()
    output_file = args.output
    output_traj = args.traj_output
    atoms = read(args.traj, ":")
    atoms = atoms if isinstance(atoms, list) else [atoms]
    processes = environ.get("SLURM_NTASKS_PER_NODE", args.processes)
    processes = int(processes)
    np_element_list = [
        image.info["np_element"] if not args.np_element else args.np_element
        for image in atoms
    ]

    with Parallel(n_jobs=processes) as parallel:
        results = parallel(
            delayed(discriminate)(
                image,
                np_element=np_element_list[index],
                support_elements=args.support_elements,
                interface_scaling=args.interface_scaling,
                surface_scaling=args.surface_scaling,
                do_alphashape=do_alphashape,  # args.do_alphashape,
                nl=None,
                analyzer=None,
            )
            for index, image in enumerate(
                tqdm(atoms, desc="Discriminating..", total=len(atoms))
            )
        )

    np_surface, np_bulk, np_interface, np_perimeter, substrate_interface, substrate = (
        zip(*results)
    )
    np_total = [
        (j + np_bulk[i] + np_perimeter[i] + np_interface[i])
        for i, j in enumerate(np_surface)
    ]

    np_interface_plus_peri = [j + np_perimeter[i] for i, j in enumerate(np_interface)]
    interfacial_results_list = [
        (j, substrate_interface[i]) for i, j in enumerate(np_interface_plus_peri)
    ]

    # return Traj but with atoms' identities switched
    for index, image in enumerate(
        tqdm(atoms, total=len(atoms), desc="Changing identities..")
    ):
        np_surf, np_int, np_peri = (
            np_surface[index],
            np_interface[index],
            np_perimeter[index],
        )

        # hardcoded; should include in args sometime but this is sufficient for tests
        for i in np_int:
            image[i].symbol = "Cu"
        for i in np_peri:
            image[i].symbol = "Mo"
        for i in np_surf:
            image[i].symbol = "Pd"

    print("Saving results")
    write(output_traj, atoms)

    ##time to write results
    # write the output for the atomistic_discriminator
    np_perimeter_len = [len(i) for i in np_perimeter]
    np_interface_len = [len(i) for i in np_interface]
    np_surface_len = [len(i) for i in np_surface]
    np_total_len = [len(i) for i in np_total]

    df = pd.DataFrame(
        {
            "Perimeter": np_perimeter_len,
            "Interface": np_interface_len,
            "Surface": np_surface_len,
            "Total": np_total_len,
        }
    )
    df.to_csv(output_file, index=False)
