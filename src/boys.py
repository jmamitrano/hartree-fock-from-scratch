"""Boys function F_n(x), the auxiliary function of Gaussian integrals.

F_n(x) = integral_0^1 t**(2n) * exp(-x * t**2) dt

It appears in the nuclear attraction and electron repulsion integrals.
"""

import numpy as np
from scipy.special import erf
from scipy.special import factorial2

# Below SMALL_X the Taylor series is used; above LARGE_X, the asymptotic formula.
# In between, F_0 from erf plus upward recurrence.
SMALL_X = 2.0
LARGE_X = 55

eps = np.finfo(float).eps

def boys_array(n_max, x, small_x=SMALL_X, large_x=LARGE_X):
    """Boys function F_n(x) for all orders n = 0, ..., n_max.

    Three branches, chosen by the value of x:

    - x < small_x: Taylor series around x = 0,
      F_n(x) = sum_k (-x)**k / (k! * (2n + 2k + 1)),
      truncated when a term falls below machine epsilon relative to the sum
      (at most 40 terms).
    - small_x <= x <= large_x: F_0(x) = 0.5 * sqrt(pi / x) * erf(sqrt(x)),
      then upward recurrence
      F_{n+1}(x) = ((2n + 1) * F_n(x) - exp(-x)) / (2x).
    - x > large_x: asymptotic formula
      F_n(x) = (2n - 1)!! / 2**(n + 1) * sqrt(pi / x**(2n + 1)).

    Parameters
    ----------
    n_max : int
        Highest order to compute, n_max >= 0.
    x : float
        Argument of the Boys function, x >= 0.
    small_x : float, optional
        Threshold between the Taylor series and the recurrence branch.
    large_x : float, optional
        Threshold between the recurrence branch and the asymptotic formula.

    Returns
    -------
    np.ndarray, shape (n_max + 1,)
        The values F_0(x), ..., F_{n_max}(x).
    """
    f_array = np.array([])
    if x < small_x:
        n = 0
        while n < n_max + 1:
            tk = 1
            summ = 0
            for k in range(40):
                nk = 1/(2*n + 2*k +1)
                old_summ = tk * nk
                summ += old_summ
                tk = -x/(k+1) * tk

                if abs(old_summ) < eps * abs(summ):
                    break
            n += 1
            f_array = np.append(f_array, summ)
                
        return f_array

    if x >= small_x and x <= large_x:
        f_array = np.array([0.5 * np.sqrt(np.pi/x) * erf(np.sqrt(x))])
        f_sum = f_array[0]
        n = 0
        while n < n_max:
            f_sum = ((2*n + 1)* f_sum - np.exp(-x))/(2*x) 
            n += 1
            f_array = np.append(f_array, f_sum)

        return f_array
    
    if x > large_x:
        f_array = np.array([0.5 * np.sqrt(np.pi/x)])
        f_sum = f_array[0]
        n = 0
        while n < n_max:
            n += 1
            f_sum = factorial2(2*n - 1) / 2**(n + 1) * np.sqrt(np.pi/x**(2*n+1))
            f_array = np.append(f_array, f_sum)

        return f_array