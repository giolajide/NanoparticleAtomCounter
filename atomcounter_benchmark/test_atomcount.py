"""
Here, I'll compare the results given by the npatomcounter to those
from making an approximate atomistic model for nanoparticles of an FCC
lattice and flat interface
"""

import subprocess
from ase.io import write
import shlex
import sys
from os import system, environ, makedirs
import numpy as np
import pandas as pd
from itertools import product
from ase import Atoms
from tqdm import tqdm
from typing import List, Union, Tuple, Literal
from ase.visualize import view
import warnings
from argparse import ArgumentParser
from pathlib import Path
from ascii_colors import ASCIIColors
from time import perf_counter
from atomcounter_benchmark.create_spherical_caps import (
    create_sphere,
    cut_particle,
)
from atomcounter_benchmark.atomistic_utils import (
    scaler,
    create_unit_support,
)
import pytest

MIN_ANGLE = 30
MAX_ANGLE = 160
MIN_RADIUS = 10  # Ang
MAX_RADIUS = 30  # Ang
N_ANGLES = 5
N_RADII = 5
PROCESSES = -1
OUTPUT_TRAJECTORY = "atoms.traj"
NP_ELEMENTS = ["Ag"]
SUPPORT_ELEMENTS = ["graphene"]  # , "au"]

script_dir = Path(__file__).resolve().parent


def create_outputdir() -> str:
    home_dir = environ["HOME"]
    output_dir = f"{home_dir}/nanoparticle_atom_counter_benchmark/"
    makedirs(output_dir, exist_ok=True)

    return output_dir


def create_trajectory(
    min_angle: float,
    max_angle: float,
    n_angles: int,
    min_radius: float,
    max_radius: float,
    n_radii: int,
    output_trajectory: str,
    np_elements: Union[str, List[str]],
    support_element: Union[str, List[str]],
) -> Tuple[List[float], List[float], List[str], List[str], int]:
    """
    Creates the atomistic model

    Requires:

        min_angle:          lowest contact angle desired
        max_angle:          max contact angle desired
        n_angles:           how many contact angles desired
        min_radius:         lowest radius OF CURVATURE (NOT footprint radius) desired
        max_radius:         figure it out
        n_radii:            figure it out
        output_trajectory:  name of file in which to write the models
        np_elements:        the type of atom of which the nanoparticle is composed,
                            as a list, one for each nanoparticle/support combination
        support_element:   the type of atom of which the support is composed,
                            as a list, one for each nanoparticle/support combination.
                            Should be 'graphene', 'mgo', or 'au'

    Returns:

        contact_angles:     list of contact angles
        radii_angstrom:     list of curvature radii
        nanoparticles:      list of nanoparticle elements
        supports:           list of support elements
        n_calculations:     number of calculations to run
    """

    if max_radius > 50:
        warnings.warn(
            """This will take some time.
                I hope you can wait 5-10 minutes!""",
            category=UserWarning,
        )

    nanoparticles = [np_elements] if isinstance(np_elements, str) else np_elements
    supports = (
        [support_element] if isinstance(support_element, str) else support_element
    )

    n_calculations = int(n_radii * n_angles * len(nanoparticles) * len(supports))
    requested_system = (
        "Requested system: \n"
        f"All possible combinations of {np_elements} nanoparticles on {support_element}\n"
        f"Nanoparticle curvature radii are from {min_radius} to {max_radius} A,\n"
        f"with contact angles from {min_angle} to {max_angle}\n\n"
        f"Total number of calculations: {n_calculations}"
    )
    print(requested_system)

    radii_angstrom = np.linspace(min_radius, max_radius, n_radii)
    contact_angles = np.linspace(min_angle, max_angle, n_angles)

    atoms_list: List[Atoms] = []

    for r in tqdm(
        radii_angstrom,
        total=len(radii_angstrom),
        desc="creating atomistic models of the requested systems",
    ):
        for theta in contact_angles:
            for nanoparticle in nanoparticles:
                for support in supports:
                    atoms = create_sphere(element=nanoparticle, radius=r)
                    atoms = cut_particle(atoms, theta)
                    unit_support = create_unit_support(support)
                    atoms = scaler(image=atoms, unit_support=unit_support)
                    atoms.info["np_element"] = nanoparticle
                    atoms_list.append(atoms)

    write(output_trajectory, atoms_list)
    print("created and written atoms objects")

    return contact_angles, radii_angstrom, nanoparticles, supports, n_calculations


def run_atomistic(
    processes: int,
    trajectory_file: str,
    radii_angstrom: List[float],
    contact_angles: List[float],
    nanoparticles: List[str],
    supports: List[str],
    input_to_atomcounter: str,
    atomistic_output: str,
    new_atoms_output: str,
    output_dir: str,
) -> None:
    """
    Use the atomistic model to get the number of each kind of atom in each nanoparticle

    Requires:

        processes:              number of parallel processes to use.
        trajectory_file:        the file containing the atomistic models
        radii_angstrom:         list of curvature radii
        contact_angles:         list of contact angles
        nanoparticles:          list of nanoparticle elements
        supports:               list of support elements
        input_to_atomcounter:   name of input file this will generate;
                                will contain the curvature radius, element, and contact angle,
                                which the atomcounter will need
        atomistic_output:       name of output file to be generated, containing the number
                                of each kind of atom in each nanoparticle
        new_atoms_output:       name of trajectory file to be created, in which the perimeter,
                                interface, and surface atoms are colored differently than the bulk
        output_dir:             name of output directory in which all files will be created

    """
    print("discriminating...")
    command = shlex.split(
        f"python {script_dir}/discrimination.py"
        f" -t {trajectory_file} -p {processes} -o {atomistic_output}"
        f" -to {new_atoms_output}"
    )

    discrimination_out = output_dir + "discrimination.out"
    discrimination_err = output_dir + "discrimination.err"

    with open(discrimination_out, "w") as out_f, open(discrimination_err, "w") as err_f:
        subprocess.run(command, stdout=out_f, stderr=err_f, check=True)

    rows = []
    for R, theta, element, _ in product(
        radii_angstrom, contact_angles, nanoparticles, supports
    ):
        rows.append(
            {
                "r (A)": "",
                "R (A)": R,
                "Theta": theta,
                "Element": element.capitalize(),
                "Interface Facet": "(1,0,0)",
                "Surface Facet": "",  # "(1,1,1)",
            }
        )

    df = pd.DataFrame(rows)
    df.to_csv(input_to_atomcounter, index=False)
    print("Finished the atomistic modelling")


def run_atomcounter(input_file: str, output_file: str, output_dir: str) -> float:
    """
    Runs the atomcounter to get the number of each kind of atom in each nanoparticle

    Requires:

        input_file:     name of input file, having the element, curvature radius, etc
        output_file:    name of output file to be generated
        output_dir:     name of output directory in which all files will be created

    Returns:

        timing:          timing in seconds
    """

    print("running NanoparticleAtomCounter ...")

    atomcounter_out = output_dir + "atomcounter.out"
    atomcounter_err = output_dir + "atomcounter.err"

    command = shlex.split(f"nanoparticle-atom-count -i {input_file} -o {output_file}")
    with open(atomcounter_out, "w") as out_f, open(atomcounter_err, "w") as err_f:
        subprocess.run(command, stdout=out_f, stderr=err_f, check=True)

    timing = None
    with open(atomcounter_out, "r") as t:
        lines = t.readlines()
        last = lines[-1].strip()
        if last.startswith("Calculation took"):
            parts = last.split()
            timing = float(parts[2])

    print("Finished running the nanoparticle atom counter")

    return timing  # seconds


def plot_parities(
    atomistic_output: str, atomcounter_output: str, input_file: str, output_dir: str
):
    """
    Creates parity plots of the results, comparing the atomistic to the atomcounter

    Requires:
        atomistic_output:       name of output file from the atomistic
        atomcounter_output:     name of output file from the atomcounter
        input_file:             name of input file to both
    """
    print("comparing . . .")
    command = shlex.split(
        f"python {script_dir}/plot-parity.py {atomistic_output} "
        f"{atomcounter_output} {input_file} --output_dir {output_dir}"
    )

    parity_out = output_dir + "parity.out"
    parity_err = output_dir + "parity.err"

    with open(parity_out, "w") as out_f, open(parity_err, "w") as err_f:
        subprocess.run(command, stdout=out_f, stderr=err_f, check=True)


##ikimashou
def main() -> None:
    """
    if any function fails, we exist with a non-zero exit code
    details on the failure will be written in the .err file (by the above functions)
    """
    output_dir = create_outputdir()
    print(f"\n\nWriting all results to {output_dir}\n\n")

    traj_file = output_dir + OUTPUT_TRAJECTORY
    input_to_atomcounter = output_dir + "input.csv"
    atomcounter_output = output_dir + "counter.csv"
    atomistic_output = output_dir + "atomistic.csv"
    new_atoms_output = output_dir + "identified.traj"

    theory = (
        "- Calculating the total number of atoms by assuming a spherical cap\n"
        "- Calculating perimeter atoms by assuming the interface is an annular ring\n"
        "  (this might introduce some errors)\n"
        "- Calculating surface atoms by assuming the nanoparticle surface is an annulus\n"
        "  (this might also introduce some errors)\n"
    )

    ASCIIColors.print(
        theory,
        color=ASCIIColors.color_yellow,
        style=ASCIIColors.style_bold,
        background=ASCIIColors.color_black,
        end="\n\n",
        flush=True,
        file=sys.stdout,
    )

    exit_code = 0
    try:
        contact_angles, radii_angstrom, nanoparticles, supports, n_calculations = (
            create_trajectory(
                min_angle=MIN_ANGLE,
                max_angle=MAX_ANGLE,
                n_angles=N_ANGLES,
                min_radius=MIN_RADIUS,
                max_radius=MAX_RADIUS,
                n_radii=N_RADII,
                output_trajectory=traj_file,
                np_elements=NP_ELEMENTS,
                support_element=[SUPPORT_ELEMENTS[0]],
            )
        )

        run_atomistic(
            processes=PROCESSES,
            trajectory_file=traj_file,
            contact_angles=contact_angles,
            radii_angstrom=radii_angstrom,
            nanoparticles=nanoparticles,
            supports=supports,
            input_to_atomcounter=input_to_atomcounter,
            atomistic_output=atomistic_output,
            new_atoms_output=new_atoms_output,
            output_dir=output_dir,
        )

        timing = run_atomcounter(
            input_file=input_to_atomcounter,
            output_file=atomcounter_output,
            output_dir=output_dir,
        )

        plot_parities(
            atomistic_output=atomistic_output,
            atomcounter_output=atomcounter_output,
            input_file=input_to_atomcounter,
            output_dir=output_dir,
        )

        print(f"\n\nAll output written to {output_dir}")

        explanation = """\
        - input.csv: input file for the NanoparticleAtomCounter
        - counter.csv: atom counts according to the NanoparticleAtomCounter
        - atomistic.csv: atom counts by using an atomistic model
        - atoms.traj: atomistic model
        - identified.traj: atomistic model, with different types of atoms distinguished
        - parity*.png: parity plots comparing the NanoparticleAtomCounter to the atomistic model
        - heatmap*.png: heatmaps comparing the NanoparticleAtomCounter to the atomistic model
        """

        speed_benchmark = f"""\
        nanoparticleatomcounter took {timing * 1000} milliseconds to run {n_calculations} calculations
        """

        readme = output_dir + "README.md"
        speed = output_dir + "timing.log"

        with open(readme, "w") as explain_file, open(speed, "w") as speed_file:
            explain_file.write(explanation)
            speed_file.write(speed_benchmark)
    except Exception as e:
        print(
            f"benchmarking failed with error ({e}). see .err files in {output_dir} for details"
        )
        exit_code = 1

    return exit_code
