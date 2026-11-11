"""Contracted Cartesian Gaussian basis functions.

Builds the basis set for a molecule from Basis Set Exchange data. Only s and p
shells (L <= 1) are supported. All lengths are in bohr.
"""

from dataclasses import dataclass
import basis_set_exchange as bse
import numpy as np

@dataclass(frozen=True)
class BasisFunction:
    center: np.ndarray
    lmn: tuple[int, int, int]
    exponents: np.ndarray
    coefficients: np.ndarray

    def __post_init__(self):
        self.center.setflags(write=False)
        self.exponents.setflags(write=False)
        self.coefficients.setflags(write=False)

def _primitive_norm(exponents,L):
    """N_k = (2 alpha_k / pi)**(3/4) * (4 alpha_k)**(L/2). Valid only for L <= 1."""    
    nk = (2*np.sqrt(exponents))**L * (2*exponents/np.pi)**(3/4)
    return nk

def _prim_overlap(exponents,L):
    """Overlap matrix S_kk' between normalized primitives with the same center and (l, m, n)."""
    skk = (4 * exponents[:, None] * exponents / (exponents[:, None] + exponents)**2) ** ((3+2*L)/4)
    return skk

def _contraction_norm(dk, skk):
    """c = (d^T S d)**(-1/2), with dk given for normalized primitives (without N_k)."""
    c = (dk @ skk @ dk)**(-1/2)
    return c

def build_basis(z_list, coords_list, basis_name):
    """One contracted Cartesian Gaussian basis function.

    phi(r) = sum_k coefficients[k] * x_A**l * y_A**m * z_A**n
             * exp(-exponents[k] * |r - A|**2)

    with x_A = x - center[0], and likewise for y and z.

    Attributes
    ----------
    center : np.ndarray, shape (3,)
        Position A of the atom the function sits on, in bohr.
    lmn : tuple[int, int, int]
        Cartesian powers (l, m, n); L = l + m + n.
    exponents : np.ndarray, shape (K,)
        Primitive exponents alpha_k, in bohr**-2.
    coefficients : np.ndarray, shape (K,)
        Final coefficients c * d_k * N_k. They already include the primitive
        normalization N_k and the contraction normalization c, so phi is
        normalized as written above.
    """
    elements = list(set(z_list))
    basis = bse.get_basis(basis_name, elements=elements)
    basis_functions = []
    lmnp = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1)
    ]
    for z,coords in zip(z_list,coords_list):
        center = np.array(coords, dtype=float)
        element = basis['elements'][str(z)]
        shells = element['electron_shells']
        for shell in shells:
            angular_momentum = shell['angular_momentum']
            exponents = np.array(shell['exponents'], dtype=float)
            if len(angular_momentum) != len(shell['coefficients']):
                raise ValueError( f"General contractions are not supported: basis '{basis_name}', Z={z} "
    f"has a shell with {len(angular_momentum)} angular momenta and "
    f"{len(shell['coefficients'])} coefficient sets. "
    "Only segmented and SP shells are handled.")
            for i,L in enumerate(angular_momentum):
                dk = np.array(shell['coefficients'][i], dtype=float)
                nk = _primitive_norm(exponents,L)
                skk = _prim_overlap(exponents,L)
                c = _contraction_norm(dk, skk)
                coefficients = nk * c * dk
                if L == 1:
                    for lmn in lmnp:
                        basis_functions.append(BasisFunction(center, lmn, exponents, coefficients))
                elif L == 0:
                    lmn = (0,0,0)
                    basis_functions.append(BasisFunction(center, lmn, exponents, coefficients))
                else:
                    raise ValueError(f"Unsupported angular momentum L={L} in basis '{basis_name}', Z={z}. "
    "Only s (L=0) and p (L=1) shells are supported.")
    return basis_functions

Z = np.array([8,1,1]) #Z values for H2O
COORDS = np.array([
    [ 0.000000, 0.0, 0.000000],  # O (Bohr)
    [ 1.432337, 0.0, 1.108069],  # H (Bohr)
    [-1.432337, 0.0, 1.108069]   # H (Bohr)
], dtype=float)
