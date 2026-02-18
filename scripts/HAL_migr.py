#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Anaïs Spire, Agathe Chave-Lucas and Pierre Veron
Date: 2026-02-18
Description: implements the deterministic prediction tool for the scenario 
    with migration (parapatric).
License: MIT License
"""

import numpy as np
import scipy
from scipy.optimize import elementwise
from scripts.HAL_general import mean_between_fitness, mean_within_fitness, sel_coeff, R_neutral, _dichotomy_decreasing

def eff_migr_rate(m, d_w, d_b, K, k):
    """ Effective migration rate, ie nb of migrants that actually mate and produce
    viable offspring (takes into account the Wallace effect), from Eq. 15 of 
    Gavrilets 1999. 
    Args:
        m (float): observed proportion of migrants (proportion of individuals of
            a given subpop replaced by migrants from another subpop). Per generation
        d_w (float): within population mean distance
        d_b (float): between population mean distance
        K (float): threshold for interfertility
        k (float) : nb of loci close to fixation in 1 subpop (but not in the other subpop)
    
    Returns:
        float: m_e : effective migration rate, ie proportion of migrants that 
            actually survive, by replacing individuals of the resident pop
    """
    mean_ww = scipy.special.gammaincc(K+1, d_w) 
    mean_wb = scipy.special.gammaincc(K+1-k, d_b-k)
    if k <= K:
        return m*mean_wb/mean_ww
    else:
        return 0.0


def R(d, N, K, stirling_approx = True, taylor_approx = True):
    """ Rate of divergence relative to the neutral case. Between 0 (strong selection) and 
    1 (neutral). Defined in Eq. 13b of Gavrilets 1999.

    Args:
        d (float): within population mean distance
        N (float): population size
        K (float): threshold for interfertility
        stirling_approx (bool, optional): use Stirling approximation for
            the calculation of sel_coeff. Default True.  
        taylor_approx (bool, optional): use Taylor approximation of order 1 for 
            erf for low selection coefficients. 

    Returns:
        float
    """
    if d == 0.0:
        return 1.0
    s = sel_coeff(d, K, stirling_approx=stirling_approx)
    return R_neutral(s, N, taylor_approx)


def odes_migr(t, X, nu1, nu2, N1, N2, K, m12, m21, t_g, contact, stirling_approx = True,
              N1_N2_are_funct = False, migr_rates_are_funct = False, tsplit = None):
    """ Returns 3-dimensional vector of derivative
    
    Args:
        t (float): time (for scipy ivp)
        X (4 dimensional vector) : (d_w1, dw2, d_b, k) : the initial conditions of the odes
        nu, N, K, m, t_g: cf. the functions in which these parameters are used
        contact (bool): if True, the populations are considered as one unique 
            population with size N1 and mutation rate nu1
        stirling_approx: if True, uses the stirling approximation for the calculation
            of the Gamma function.
        
    Returns:
        3 dimensional vector : (dD_w1_migr_pred, dD_w2_migr_pred, dD_b_migr_pred, 
            dk_pred) : the predictions of dDw1, dDw2, dDb and dk over time
    """
    d_w1, d_w2, d_b, k = X
    if N1_N2_are_funct:
        N1_,N2_ = N1(t-tsplit), N2(t-tsplit)
    else:
        N1_,N2_ = N1,N2
    if migr_rates_are_funct:
        m12_, m21_ = m12(t-tsplit), m21(t-tsplit)
    else:
        m12_, m21_ = m12, m21
    s1 = sel_coeff(d_w1, K, stirling_approx=stirling_approx)
    s2 = sel_coeff(d_w2, K, stirling_approx=stirling_approx)
    if contact:
        dd_w2, dd_b, dk, m_e21 = 0.0, 0.0, 0.0, 0.0
    else:
        m_e12 = eff_migr_rate(m12_, d_w2, d_b, K, k)
        m_e21 = eff_migr_rate(m21_, d_w1, d_b, K, k)
        R1 = R(d_w1, N1_, K)
        R2 = R(d_w2, N2_, K)
        dd_w2 = (-s2*d_w2 + 2*nu2 - d_w2/N2_ + 2*m_e12*(d_b - d_w2))/t_g
        dd_b = (nu1 + nu2 + 
                m_e12 * (d_w1-d_b) +
                m_e21 * (d_w2-d_b) -
                s1/2 * (d_b-k) - 
                s2/2 * (d_b-k)) / t_g 
        dk = (nu1 * R1 * 2**(-2*N1_*m_e21) + 
              nu2 * R2 * 2**(-2*N2_*m_e12) - 
              k * (m_e21 * R1 * (np.e/2)**(2*N1_*m_e21) + 
                   m_e12 * R2 * (np.e/2)**(2*N2_*m_e12)))

    
    dd_w1 = (-s1*d_w1 + 2*nu1 - d_w1/N1_ + 2*m_e21*(d_b - d_w1))/t_g

    return np.array([dd_w1, dd_w2, dd_b, dk])


def _bracketable_wb(X, K):
    """Internal function, used to find the speciation time"""
    d_w1, d_w2, d_b, k = X
    wb = mean_between_fitness(d_b, k, K)
    if wb == 0.0:
        return -1
    else:
        return wb 


def find_speciation_time(solution, tmin, tmax, K, dt = 1.0):
    """Finds the speciation time of a solved system with a precision dt.

    Args:
        solution (func): function f(t) -> [d_w1, d_w2, d_b, k]
        tmin (float): time such that w_b(t) > 0
        tmax (float): time such that w_b(t) = 0
        K (float): Threshold for outbreeding depression. 
        dt (float, optional): precision for the speciation time. Defaults to 1.0.

    Returns:
        float: speciation time (not duration)
    """
    f = lambda t: _bracketable_wb(solution(t), K)
    if f(tmin) < 0 or f(tmax) > 0:
        raise Exception("Provide valid bounds for the speciation time. ")
    br = elementwise.bracket_root(f, xl0 = tmin, xmin = tmin, xmax = tmax)
    if br.success:
        tspec =  _dichotomy_decreasing(f, *br.bracket, dt = dt)
        return tspec
    else:
        raise Exception(f"Unable to find speciation time within the given time frame: solver error {br.status}")

def solve_ODE_HAL_migr(t, nu1, nu2, N1, N2, K, t_g, burnin, m12 = 0.0, m21 = 0.0, Na = "sum",
                       stirling_approx = True, 
                       solver_kwargs = dict(),
                       error_burnin_convergence = True, dw_burnin_tol = 0.005,
                       extend_burnin = 1, N1_N2_are_funct = False,
                       migr_rates_are_funct = False, dt = 10):
    """ Solves the whole system of ODE for the Holey Adaptive Landscape with 
    migration, including a burnin phase with one population. 

    Args:
        t (float): duration of the prediction (not including burnin)
        nu1, nu2 (float): mutation rates of population 1 and 2 (initial 
            population has mutation rate nu1)
        N1, N2 (float or functions): population size (for one patch).
        K (float): threshold for interfertility
        t_g (float): generation time
        burnin (float): duration of the burnin phase (= splitting time)
        m12 (float or function, optional): migration rate defined as the proportion of population 2
            that was in population 1 in the previous generation. Delaults to 0.0.
        m21 (float or function, optional): migration rate defined as the proportion of population 1
            that was in population 2 in the previous generation. Delaults to 0.0.
        Na (float or string, optional): ancestral population size. If "sum", is 
            set to N1+N2. Defaults to "sum".
        stirling_approx (bool, optional): use Stirling approximation for
            the calculation of sel_coeff. Defaults to True.
        solver_kwargs (dict, optional): kwargs passed to the solver. See arguments
            accepted by scipy.integrate.solve_ivp. Defaults to dict().
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
        N1_N2_are_funct (bool, optional): specify if N1_N2 are function of the 
            time (time after split). In this case, Na has to be specified. 
        migr_rates_are_funct (bool, optional): specify if the migration rates are 
            function of time (after split). 
        dt (float > 0, optional): resolution for the solver. Default to 10. 

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
    # Burnin 
    step, start = 0, 0.0 # number of times of burnin extension until convergence 
    Dw_convergence = False
    T_burnin, Dw_burnin = [0.0], [0.0]
    if Na == "sum":
        Na = N1+N2
    elif type(Na) == str:
        raise ValueError("Na must be a float or 'sum'.")
    while not(Dw_convergence) and step < extend_burnin:
        sol_burnin = scipy.integrate.solve_ivp(odes_migr, t_span = (start, start + burnin), 
                               y0 = [Dw_burnin[-1],0.0,0.0,0.0],
                               args = (nu1, 0.0, Na, 0.0, K, 0.0, 0.0, 1.0, True,    stirling_approx), 
                                       #nu1,nu2, N1,     N2, K, m12, m21, t_g, contact, stirling_approx 
                               **solver_kwargs)
        T_burnin += list(sol_burnin['t'])
        Dw_burnin += list(sol_burnin['y'][0,:])
        diff_Dw_burnin = (Dw_burnin[-1] - Dw_burnin[-2]) / (T_burnin[-1] - T_burnin[-2])
        Dw_convergence = (diff_Dw_burnin /(2*nu1) < dw_burnin_tol)
        step += 1
        start += burnin 
    if error_burnin_convergence and not(Dw_convergence):
        raise Exception("Burnin phase did not satisfy convergence criterion.")
    # Split
    sol_split = scipy.integrate.solve_ivp(odes_migr, t_span = (start, start + t), 
                         t_eval = np.arange(start, start + t, dt),
                         dense_output = True, 
                         y0 = [Dw_burnin[-1], Dw_burnin[-1], Dw_burnin[-1], 0.0],
                         args = (nu1, nu2, N1, N2, K, m12, m21, t_g, False,    stirling_approx, N1_N2_are_funct, migr_rates_are_funct, start), 
                                 #nu1,nu2, N1, N2, K, m12, m21, t_g, contact, stirling_approx,  N1_N2_are_funct, migr_rates_are_funct, tsplit
                          **solver_kwargs)  
    T = sol_split['t']
    Dw1 = sol_split['y'][0,:]
    Dw2 = sol_split['y'][1,:]
    Db = sol_split['y'][2,:]
    k = sol_split['y'][3,:]

    # Calculate additional quantities
    tsplit = T[0]
    if N1_N2_are_funct:
        N1t = lambda t: N1(t-tsplit)
        N2t = lambda t: N2(t-tsplit)
    else:
        N1t = lambda t: N1
        N2t = lambda t: N2
    if migr_rates_are_funct:
        m12t = lambda t: m12(t-tsplit)
        m21t = lambda t: m21(t-tsplit)
    else:
        m12t = lambda t: m12
        m21t = lambda t: m21
    m_e12 = np.array([eff_migr_rate(m12t(t), Dw2[i], Db[i], K, k[i]) for i,t in enumerate(T)])
    m_e21 = np.array([eff_migr_rate(m21t(t), Dw1[i], Db[i], K, k[i]) for i,t in enumerate(T)])
    R1 = np.array([R(dw, N1t(t), K, stirling_approx) for dw, t in zip(Dw1, T)])
    R2 = np.array([R(dw, N2t(t), K, stirling_approx) for dw, t in zip(Dw2, T)])
    s1 = np.array([sel_coeff(dw, K, stirling_approx) for dw in Dw1])
    s2 = np.array([sel_coeff(dw, K, stirling_approx) for dw in Dw1])
    wb = np.array([mean_between_fitness(Db[i], k[i], K) for i in range(len(T))])
    ww1 = np.array([mean_within_fitness(dw, K) for dw in Dw1])
    ww2 = np.array([mean_within_fitness(dw, K) for dw in Dw2])
    
    # Calculate speciation time
    if wb[-1] > 0:
        speciation = False
        t_spec = np.inf
    else:
        speciation = True
        t_spec = find_speciation_time(sol_split.sol, start, start + t, K) - start
    
    return dict(T_burnin = np.array(T_burnin), Dw_burnin = np.array(Dw_burnin), 
                burnin_converge = Dw_convergence,
                T = T, Dw1 = Dw1, Dw2 = Dw2, Db = Db, k = k, m_e12 = m_e12, 
                m_e21 = m_e21, R1 = R1, R2 = R2, s1 = s1, s2 = s2, 
                wb = wb, ww1 = ww1, ww2 = ww2, speciation = speciation, t_spec = t_spec,
                message_burnin = sol_burnin.message, message_split = sol_split.message)