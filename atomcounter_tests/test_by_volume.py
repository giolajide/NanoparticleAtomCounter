"""
Tests all functions in the by_volume script
"""

import pytest
import numpy as np
import warnings
from NanoparticleAtomCounter.by_volume import (
    calculate_by_volume,
    calculate_volumes,
)


def test_calculate_volumes():
    ##test angles
    for angle in [0, 180, -20, 200]:
        with pytest.raises(ValueError):
            calculate_volumes("Pt", 3.2, angle)
    ##test invalid elements
    with pytest.raises(KeyError):
        calculate_volumes("Fol", 3.2, 33)
    ##test negative, low, and zero radii
    with pytest.raises(ValueError):
        calculate_volumes("Mn", -3.2, 33, (1, 0, 0))
    with pytest.raises(ValueError):
        calculate_volumes("V", 0, 33, (1, 0, 0))
    with pytest.warns(UserWarning):
        calculate_volumes("Ti", 1.2, 134, (1, 0, 0))
    ##test valid values
    interface_V, peri_V, total_V = calculate_volumes("Fe", 66, 90, (1, 1, 1))
    expected = [10442, 888, 602130]
    for index, region_volume in enumerate([interface_V, peri_V, total_V]):
        assert region_volume == pytest.approx(expected[index], abs=2)


def test_calculate_by_volume():
    """
    No need, here, to test theta since a downstream function does so;
    Also no need to test element for the same reason
    """
    ##test valid values, using a hemispherical cap,
    # whose equations are simpler to calculate by hand
    peri_atoms, interface_atoms, surf_atoms, total_atoms = calculate_by_volume(
        "Ni", 150, 90, None, None
    )
    expected = [432, 12740, 25994, 648283]
    for index, region_atoms in enumerate(
        [peri_atoms, interface_atoms, surf_atoms, total_atoms]
    ):
        assert region_atoms == pytest.approx(expected[index], abs=2)
