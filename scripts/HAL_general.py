#!/usr/bin/env python3 
import numpy as np
import scipy, warnings

def upper_gamma(a, x):
    """ Incomplete gamma function, defined as 
    \int_x^{\infty} t^{s-1} e^{-t} dt 
    https://en.wikipedia.org/wiki/Incomplete_gamma_function
    
    Args:
        a (float): positive parameter
        x (float): nonnegative argument

    Returns:
        float
    """
    return scipy.special.gammaincc(a, x) * scipy.special.gamma(a)

def sel_coeff(d_w, K, stirling_approx = True):
    """ Selection coefficient from Eq. 11 of Gavrilets 1999.

    Args:
        d_w (float >= 0): within population mean distance
        K (float > 0): threshold for interfertility  
        stirling_approx (bool, optional): use Stirling approximation for the calculation
            Default True.

    Returns:
        float: selection coefficient of the purifying selection.
    """
    d_w, K = np.float64(d_w), np.float64(K)
    if stirling_approx:
        return np.exp(K - d_w) * np.power(d_w / K, K, dtype = np.longdouble) / (np.sqrt(2*np.pi*K) * scipy.special.gammaincc(K+1, d_w))
    elif K * np.log10(d_w) >= 150.0:
        warnings.warn("Trying to call the function with large values of dw and K could lead to numerical errors: consider using the Stirling approximation (stirling_approx = True)")
    return np.exp(-d_w) * np.power(d_w, K, dtype = np.longdouble) / upper_gamma(K+1, d_w)


def R_neutral(s, N, taylor_approx = True):
    """ Calculates the probability of fixation of a mutation relative to the case 
    without outbreeding depression, in the neutral case. This term is called 
    'rate of divergence relative to the neutral case' in Gavrilets 1999 paper, 
    equation 13b.

    Args:
        s (float): coefficient of selection induced by outbreeding depression, >0
        N (float): population size, > 0
        taylor_approx (bool, optional): use a Taylor approximation for the 
            calculation. Defaults to True.

    Returns:
        float
    """
    S = N*s / 2
    if taylor_approx and S < 1e-5: # approximate by Taylor expansion at order 1 of erf
        return np.exp(-S)
    elif S < 1e-5:
        warnings.warn("Trying to calculate R with low value for selection coefficient: consider using Taylor approximation by setting taylor_approx = True.")
    rS = np.sqrt(S, dtype = np.float64)
    return 2 * np.exp(-S) * rS / (np.sqrt(np.pi) * scipy.special.erf(rS))

def mean_within_fitness(dw, K):
    """ Calculates the probability of compatibility between two individuals 
    at random within one population. 

    Args:
        dw (float): within population mean distance
        K (float): threshold for interfertility

    Returns:
        float
    """
    return scipy.special.gammaincc(K+1, dw)

def mean_between_fitness(db, k, K):
    """ Calculates the probability of compatibility between two individuals 
    at random between two populations. 

    Args:
        d_b (float): between population mean distance
        k (float) : nb of loci close to fixation in 1 subpop (but not in the other subpop)
        K (float): threshold for interfertility
        
    Returns:
        float
    """
    if k > K:
        return 0.0 
    return scipy.special.gammaincc(K+1-k, db - k)


def _dichotomy_decreasing(f, t0, t1, dt):
    """ Finds the root of a decreasing function given two bounds. 
    Internal function, used to find the speciation time
    """
    # 
    while t1 - t0 > dt:
        t = (t1+t0)/2
        if f(t) < 0:
            t1 = t
        else:
            t0 = t
    return (t1+t0)/2