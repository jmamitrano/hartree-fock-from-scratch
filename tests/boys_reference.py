"""Slow, independent reference implementations used only by the tests"""
import numpy as np
from scipy.integrate import quad
from scipy.special import hyp1f1

def boys_reference(n, x):
    """Boys function from its definition, integrated numerically"""
    f = lambda t: t ** (2 * n) * np.exp(-x * t * t)
    return quad(f, 0, 1, epsabs=0, epsrel=1.5e-13)[0]

def boys_hyp1f1(n,x):
    """Boys function via the confluent hypergeometric function"""
    return hyp1f1(n + 0.5, n + 1.5, -x)/ (2 * n + 1)
