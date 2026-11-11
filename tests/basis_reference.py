"""Slow reference implementations for testing src/basis.py.

The functions integrate numerically with scipy.integrate.quad 
and share no formulas with src/basis.py,
so they can serve as an independent oracle.
"""

import numpy as np
from scipy.integrate import quad

def contracted_norm_reference(L, exponents, coefficients):

    """Norm of a contracted Cartesian Gaussian by numerical integration.

    Computes the integral of phi**2 over all space as

        A * integral_0^inf (sum_k coefficients[k] * exp(-exponents[k] * r**2))**2
                           * r**(2L + 2) dr

    where A is the angular factor: 4*pi for s (L=0) and 4*pi/3 for p (L=1).
    The result does not depend on the center of the function.

    Parameters
    ----------
    L : int
        Angular momentum, l + m + n. Only 0 and 1 are supported.
    exponents : np.ndarray, shape (K,)
        Primitive exponents alpha_k, in bohr**-2.
    coefficients : np.ndarray, shape (K,)
        Final coefficients multiplying each primitive, including all
        normalization factors.

    Returns
    -------
    float
        The integral of phi**2. Equals 1 for a normalized function.
    """

    exponents = np.asarray(exponents, dtype=float)
    coefficients = np.asarray(coefficients, dtype=float)

    if L == 0:
        A = 4*np.pi
    elif L == 1:
        A = 4/3 * np.pi 
    else:
        raise ValueError("L must be either 0 or 1")
    
    f = lambda r: A*(np.sum(coefficients*np.exp(-exponents * r**2)))**2 * r**(2*L + 2)

    return quad(f, 0, np.inf, epsabs=0, epsrel=1.5e-13)[0]