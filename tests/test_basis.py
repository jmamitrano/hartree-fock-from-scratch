"""Tests for src/basis.py: function count and normalization."""

import numpy as np
import pytest
from src.basis import build_basis, BasisFunction
from tests.basis_reference import contracted_norm_reference

Z = np.array([8,1,1]) #Z values for H2O
COORDS = np.array([
    [ 0.000000, 0.0, 0.000000],  # O (Bohr)
    [ 1.432337, 0.0, 1.108069],  # H (Bohr)
    [-1.432337, 0.0, 1.108069]   # H (Bohr)
], dtype=float)
EXPECTED_N_FUNCTIONS = {"sto-3g": 7, "6-31g": 13}

@pytest.mark.parametrize('basis_name,basis_value', EXPECTED_N_FUNCTIONS.items())
def test_number_of_functions(basis_name, basis_value):
    """H2O has the expected number of basis functions in each basis set."""    
    assert len(build_basis(Z, COORDS, basis_name)) == basis_value

@pytest.mark.parametrize('basis_name', EXPECTED_N_FUNCTIONS.keys())
def test_normalization(basis_name):
    """Every basis function has unit norm according to the quad reference.

    rtol=1e-11 sits above the quad accuracy (~1e-13) and below the smallest
    deviation (~2e-11) that appears if the contraction normalization is skipped.
    """    
    expected = 1.0
    for basis_functions in build_basis(Z, COORDS, basis_name):
        L = sum(basis_functions.lmn)
        np.testing.assert_allclose(contracted_norm_reference(L, basis_functions.exponents, basis_functions.coefficients),
                                    expected, rtol=1e-11, err_msg = f'{basis_functions} failed')