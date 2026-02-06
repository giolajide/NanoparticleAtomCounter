"""
Tests all functions in the count_utilities script
"""

import pytest
import numpy as np
import warnings
from NanoparticleAtomCounter.count_utilities import (
    volume_to_atoms,
    beta,
    alpha,
    area_to_atoms,
    calculate_atomic_density,
    calculate_total_volume,
    calculate_surface_area,
    calculate_constants,
    DEFAULT_FACETS,
)

TOLERANCE = 1e-4


def test_default_facets():
    assert DEFAULT_FACETS["hcp"] == "(0, 0, 1)"


def test_calculate_constants():
    pd_molarvol, pd100_planarspacing, pd_dia = calculate_constants("Pd", (1, 0, 0))
    assert pd_molarvol == pytest.approx(8.86216261990501e24, abs=TOLERANCE)
    assert pd100_planarspacing == pytest.approx(1.945, abs=TOLERANCE)
    assert pd_dia == pytest.approx(2.78, abs=TOLERANCE)


def test_alpha():
    assert alpha(180) == np.inf
    assert alpha(90) == pytest.approx(1.0, abs=TOLERANCE)


def test_beta():
    for angle in [0, 180]:
        with pytest.raises(ValueError):
            beta(angle)
    assert beta(90) == pytest.approx(2.0, abs=TOLERANCE)


def test_calculate_surface_area():
    ##test angles
    for angle in [0, 180, -20, 200]:
        with pytest.raises(ValueError):
            calculate_surface_area("Pt", 3.2, angle, (1, 0, 0))
    ##test invalid elements
    with pytest.raises(KeyError):
        calculate_surface_area("Px", 3.2, 33, (1, 0, 0))
    ##test negative, low, and zero radii
    with pytest.raises(ValueError):
        calculate_surface_area("Pd", -3.2, 33, (1, 0, 0))
    with pytest.raises(ValueError):
        calculate_surface_area("Pd", 0, 33, (1, 0, 0))
    with pytest.warns(UserWarning):
        calculate_surface_area("Pd", 1.2, 134, (1, 0, 0))
    ##test valid values
    area = calculate_surface_area("Pd", 70, 122, None)
    assert area > 0


def test_calculate_total_volume():
    vol = calculate_total_volume(9, 99)
    assert vol > 0


def test_calculate_atomic_density():
    assert calculate_atomic_density("Ag", (1, 1, 1)) == pytest.approx(
        0.1380551931635095, abs=TOLERANCE
    )


def test_area_to_atoms():
    assert area_to_atoms(200, "Au", (1, 1, 1)) == pytest.approx(27, abs=1)


def test_volume_to_atoms():
    assert volume_to_atoms(1000, "Cr", None) == pytest.approx(83, abs=1)
