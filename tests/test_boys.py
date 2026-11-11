"""Tests for src/boys.py against two independent references."""

import numpy as np
import pytest
from tests.boys_reference import boys_reference, boys_hyp1f1
from src.boys import boys_array, SMALL_X, LARGE_X

X_VALUES = [0.0, 1e-12, SMALL_X-0.1, SMALL_X, SMALL_X+0.1, 1.0, 5.0, 10.0, LARGE_X/2, LARGE_X-0.1, LARGE_X, LARGE_X+0.1, LARGE_X*10, LARGE_X*100, LARGE_X*1000]
N_MAX_VALUES = [0, 1, 2, 3, 4]

@pytest.mark.parametrize("n_max", N_MAX_VALUES)
def test_value_at_zero(n_max):
    """F_n(0) = 1 / (2n + 1), exact analytic value at x = 0."""
    expected = [1 / (2 * n + 1) for n in range(n_max + 1)]
    np.testing.assert_allclose(boys_array(n_max, 0.0), expected, rtol=1e-12)

@pytest.mark.parametrize("n_max", N_MAX_VALUES)
@pytest.mark.parametrize("x", X_VALUES)
def test_oracle_agree(n_max, x):
    """The quad and hyp1f1 references agree with each other.

    hyp1f1 alone cannot validate the Taylor branch, since with these
    parameters its series is the same Taylor series. Agreement with quad
    settles this independently.
    """
    expected = [boys_hyp1f1(n,x) for n in range(n_max + 1)]
    np.testing.assert_allclose([boys_reference(n, x) for n in range(n_max + 1)], expected, rtol=1e-12)

@pytest.mark.parametrize("n_max", N_MAX_VALUES)
@pytest.mark.parametrize("x", X_VALUES)
def test_against_oracle(n_max, x):
    """boys_array matches the hyp1f1 reference for every branch."""
    expected = [boys_hyp1f1(n,x) for n in range(n_max + 1)]
    np.testing.assert_allclose(boys_array(n_max, x), expected, rtol=1e-12)

@pytest.mark.parametrize("n_max", N_MAX_VALUES)
@pytest.mark.parametrize("x", X_VALUES)
def test_properties(n_max, x):
    """Length, positivity and strict decrease in n of the F_n array."""
    result = boys_array(n_max, x)
    assert len(result) == n_max + 1 #the F_n array must have the length of n_max + 1
    assert np.all(result > 0) #all elements in the F_n array should be positive
    assert np.all(np.diff(result) < 0) #the F_n array should decrease when n increase
