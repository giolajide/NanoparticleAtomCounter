"""
Tests all functions in the by_area script
"""

import pytest
import numpy as np
import warnings
from NanoparticleAtomCounter.by_area import (
    calculate_by_area,
    calculate_areas,
)


def test_calculate_areas():
    """
    No need, here, to test theta since it'll
    be tested when calculate_areas calls calculate_surface_area
    """
    ##test invalid elements
    with pytest.raises(KeyError):
        calculate_areas("Jc", 6.0, 99, None, None)
    ##test negative, low, and zero radii
    with pytest.raises(ValueError):
        calculate_areas("Co", -3.2, 33, None, None)
    with pytest.raises(ValueError):
        calculate_areas("Co", 0.0, 155, (1, 0, 0))
    with pytest.warns(UserWarning):
        calculate_areas("Cu", 1.2, 67, (1, 0, 0))
    ##test valid values
    interface_A, peri_A, surf_A = calculate_areas("Cr", 55, 90, None, (1, 1, 1))
    expected = [8567, 936, 18328]
    for index, region_area in enumerate([interface_A, peri_A, surf_A]):
        assert region_area == pytest.approx(expected[index], abs=2)


def test_calculate_by_area():
    """
    No need, here, to test theta since a downstream function does so;
    Also no need to test element for the same reason
    """
    ##test valid values, using a hemispherical cap,
    # whose equations are simpler to calculate by hand
    peri_atoms, interface_atoms, surf_atoms, total_atoms = calculate_by_area(
        "Ag", 50, 90, (1, 0, 0), (1, 1, 1)
    )
    expected = [106, 833, 2083, 15306]
    for index, region_atoms in enumerate(
        [peri_atoms, interface_atoms, surf_atoms, total_atoms]
    ):
        assert region_atoms == pytest.approx(expected[index], abs=2)
