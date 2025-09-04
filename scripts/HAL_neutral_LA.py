#!/usr/bin/env python3 
# Scripts for the resolution of the HAL model without migration.
# Author: Pierre Veron
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import scipy.special
from scripts.HAL_general import mean_between_fitness, mean_within_fitness, sel_coeff, R_neutral

def Rlimit_LA(N, sLA):
    """ Calculates the R coefficient in the case with local adaptation when the 
    coefficient of selection is small. 
    
    Args:
        N (float): population size
        sLA (float): coefficient of local adaptation

    Returns:
        float
    """
    # limit of R with local adaptation when s-->0
    return 2*N*sLA / (1-np.exp(-2*N*sLA))

def R_LA(s, sLA, N, taylor_approx = True):
    """ Calculates the probability of fixation of a mutation relative to the case 
    without outbreeding depression, in the case with local adatpation. 
    See equation 19, Gavrilets 1999.

    Args:
        s (float): coefficient of selection induced by outbreeding depression, >0
        sLA (float): coefficient of selection for local adaptation
        N (float): population size, > 0
        taylor_approx (bool, optional): use a Taylor approximation for the 
            calculation. Defaults to True.

    Returns:
        float
    """
    S = N * s / 2
    if sLA > 0: # with local adaptation
        if s == 0.0:
            return Rlimit_LA(N, sLA)
        alpha = sLA / s
        rS = np.sqrt(S)
        erf1 = scipy.special.erf(rS * (1 + alpha))
        erf2 = scipy.special.erf(rS * (1 - alpha))
        if erf1 == 1.0 or erf2 == -1.0:
            return Rlimit_LA(N, sLA)  
        sumerf = erf1 + erf2
        return 4 * np.exp(-S*(1-alpha)**2) * rS / (np.sqrt(np.pi) * sumerf)

    # neutral case (without local adaptation)
    return R_neutral(s, N, taylor_approx)

def R(d, N, K, sLA = 0.0, stirling_approx = True, taylor_approx = True):
    """ Rate of divergence relative to the neutral case. Between 0 (strong selection) and 
    1 (neutral). Defined in Eq. 13b of Gavrilets 1999 and 19 for the case with 
    local adaptation. 

    Args:
        d (float): within population mean distance
        N (float): population size
        K (float): threshold for interfertility
        sLA (float, optional): average strenght of selection per locus for local adaptation.
            Default 0.0.
        stirling_approx (bool, optional): use Stirling approximation for
            the calculation of sel_coeff. Default True.  
        taylor_approx (bool, optional): use Taylor approximation of order 1 for 
            erf for low selection coefficients. 

    Returns:
        float
    """
    s = sel_coeff(d, K, stirling_approx=stirling_approx)
    return R_LA(s, sLA, N, taylor_approx=taylor_approx)


def odes_LA(t, X, nu, N, K, t_g, sLA, n, stirling_approx = True):
    """ Returns 3-dimensional vector of derivative
    
    Args:
        t (float): time (for scipy ivp)
        X (3 dimensional vector) : (d_w, d_b, k) : the initial conditions of the odes
        nu, N, K, t_g, n, sLA : cf. the functions in which these parameters are used
        
    Returns:
        3 dimensional vector : (dD_w_pred, dD_b_pred, dk_pred) : 
            the predictions of dDw, dDb and dk over time
    """
    d_w, d_b, k = X
    s = sel_coeff(d_w, K, stirling_approx=stirling_approx)
    if n == 1.0:
        dd_b, dk = 0.0, 0.0
    else:
        R_ = R_LA(s, sLA, N)
        dd_b = 2*nu*R_ / t_g
        dk = dd_b

    

    dd_w = ((-s - 1/N) * d_w + 2*nu)/t_g
    return np.array([dd_w, dd_b, dk], dtype = np.float64)

def solve_D_w_star(nu, N, K, sLA, t_g = 1.0, stirling_approx = True, **kwargs):
    """ Find the mutation-selection-drift equilibrium solution D_w* of the within-population
    mean distance D_w.

    Args:
        nu (float): per generation per individual probability of mutation
        N (float): population size
        K (float): threshold for interfertility
        sLA (float): average strenght of selection per locus for local adaptation.
        t_g (float, optional): Generation time, useless here. Defaults to 1.0.
        stirling_approx (bool, optional): use Stirling approximation for
            the calculation of sel_coeff. Default True.  

    Returns:
        float: 
    """
    # For high values of 2 N nu, we can not calculate dD_w so choose another initial value
    if K * (np.log(2 * nu * N) - np.log(K)) > 91:
        x0 = K
    else:
        x0 = 2 * nu * N
    dD_w = lambda d: odes_LA(0, [d[0],0,0], nu, N, K, t_g, sLA, 1, stirling_approx)[0]
    return scipy.optimize.root(dD_w, x0 = x0).x[0]


def solve_ODE_HAL_LA(t, nu, N, K, t_g, burnin, sLA = 0.0,
                       stirling_approx = True, 
                       solver_kwargs = dict(), wb_tol = 0.0, 
                       error_burnin_convergence = True, dw_burnin_tol = 0.005,
                       extend_burnin = 1):
    """ Solves the whole system of ODE for the Holey Adaptive Landscape with 
    migration, including a burnin phase with one population. 

    Args:
        t (float): duration of the prediction (not including burnin)
        nu (float): mutation rate
        N (float): population size (for one patch). Initial population is 2*N.
        K (float): threshold for interfertility
        t_g (float): generation time
        burnin (float): duration of the burnin phase (= splitting time)
        sLA (float, optional): average strength of selection per locus for local adaptation.
            Default 0.0.
        stirling_approx (bool, optional): use Stirling approximation for
            the calculation of sel_coeff. Defaults to True.
        solver_kwargs (dict, optional): kwargs passed to the solver. See arguments
            accepted by scipy.integrate.solve_ivp. Defaults to dict().
        wb_tol (float or str, optional): tolerance for speciation in the value of
            w_b. Must be a number >= 0 and < 1, or "1/N". In the latter case, the
            value is set to 1/N. If 0.0, speciation occurs when the probability of 
            a viable offspring between the two populations is exactly 0. Otherwise,
            speciation occurs when this probability is <= wb_tol. Defaults to 0.0.
        error_burnin_convergence (bool, optional): if True, an error is raised
            when the convergence is not reached after the burnin phase. 
            Defaults to True.
        dw_burnin_tol (float > 0, optional): tolerance used for the convergence
            test on D_w after the burnin phase. The test is:
            d (D_w) / dt < dw_burnin_tol * 2 * nu
            Defaults to 0.005.
        extend_burnin (int > 0, optional): maximum number of times to repeat the
            burnin phase if the convergence is not reached. If 1, the burnin is 
            not repeated.
            Default to 1. 

    Returns:
        dict: solutions to the burnin phase and the split phase. Objects:
            T_burnin, Dw_burnin,T, Dw, Db, k (arrays): solutions of the burnin 
                and the split phase.
            m_e, s, wb, ww (arrays): additional quantities calculated on the 
                split phase (effective migration rate, selection coefficient, 
                inter-pop mean fitness, intra-pop mean fitness).
            burnin_converge (float): result on the convergence test on Dw_burnin.
            speciation (float): True if speciation is reached
            t_spec (float): duration of speciation of np.inf if not reached. 
    """
    if wb_tol == "1/N":
        wb_tol = 1.0/N    
    elif not(isinstance(wb_tol, float) or isinstance(wb_tol, int)) or wb_tol < 0 or wb_tol >= 1:
        raise(ValueError("wb_tol should be '1/N' or a number between 0 (included) and 1 (excluded)."))
    
    # Burnin 
    step, start = 0, 0.0 # number of times of burnin extension until convergence 
    Dw_convergence = False
    T_burnin, Dw_burnin = [0.0], [0.0]
    while not(Dw_convergence) and step < extend_burnin:
        sol_burnin = solve_ivp(odes_LA, t_span = (start, start + burnin), 
                               y0 = [Dw_burnin[-1],0.0,0.0],
                               args = (nu, 2*N, K, t_g, 0.0, 1, stirling_approx),
                               t_eval = np.linspace(start, start + burnin, 10000),
                               **solver_kwargs)
        T_burnin += list(sol_burnin['t'])
        Dw_burnin += list(sol_burnin['y'][0,:])
        diff_Dw_burnin = (Dw_burnin[-1] - Dw_burnin[-2]) / (T_burnin[-1] - T_burnin[-2])
        Dw_convergence = (diff_Dw_burnin /(2*nu) < dw_burnin_tol)
        step += 1
        start += burnin 
    if error_burnin_convergence and not(Dw_convergence):
        raise Exception("Burnin phase did not satisfy convergence criterion.")
    
    # Split
    sol_split = solve_ivp(odes_LA, t_span = (start, start + t), 
                          y0 = [Dw_burnin[-1], Dw_burnin[-1], 0.0],
                          args = (nu, N, K, t_g, sLA, 2, stirling_approx), 
                          t_eval = np.linspace(start, start + t, 10000),
                          **solver_kwargs)  
    T = sol_split['t']
    Dw = sol_split['y'][0,:]
    Db = sol_split['y'][1,:]
    k = sol_split['y'][2,:]

    # Calculate additional quantities
    R_ = np.array([R(dw, N, K, sLA, stirling_approx) for dw in Dw])
    s = np.array([sel_coeff(dw, K, stirling_approx) for dw in Dw])
    wb = np.array([mean_between_fitness(Db[i], k[i], K) for i in range(len(T))])
    ww = np.array([mean_within_fitness(dw, K) for dw in Dw])

    # Calculate speciation time
    speciation, t_spec, i = False, np.inf, 0
    while i < len(T) and not(speciation):
        if wb[i] <= wb_tol:
            speciation = True
            t_spec = T[i] - start
        i += 1
    
    return dict(T_burnin = np.array(T_burnin), Dw_burnin = np.array(Dw_burnin), 
                burnin_converge = Dw_convergence,
                T = T, Dw = Dw,  Db = Db, k = k, R_ = R_, s = s,
                wb = wb, ww = ww, speciation = speciation, t_spec = t_spec)