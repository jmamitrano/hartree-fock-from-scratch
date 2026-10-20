import numpy as np
from scipy.special import erf
from scipy.special import factorial2

SMALL_X = 2.0
LARGE_X = 55

eps = np.finfo(float).eps

def boys_array(n_max, x, small_x=SMALL_X, large_x=LARGE_X):
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