#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pierre Veron
Date: 2026-02-18
Description: Wrapper script to run the HAL deterministic prediction from command-line.
License: MIT License
"""
import argparse
import os
import json
import matplotlib.pyplot as plt
import numpy as np
from scripts import HAL_neutral_LA
from scripts import HAL_migr
from scripts import HAL_plot
from scripts import utils




def HAL_prediction(time, popsize, Na, K, nu, timeburnin, m, sla, precision, migr_outwards = False, **kwargs):
    """ Computes the resolution of the Holey Adaptive Lanscape model.

    Args:
        time (float): Number of generations after split.
        popsize (float or list of floats with size 1 or 2): Population size(s)
        Na (float or None): Ancestral population size (if None, set to N1+N2)
        K (float): Threshold for outbreeding depression
        nu (float or list[1 or 2] floats): Mutation rate(s)
        timeburnin (float or None): Burnin time, if None set automatically.
        m (float or list[1 or 2] floats or None): Migration rate(s). If None,
            assumed zero. 
        sla (float or None): Coefficient for local adaptation. 
        precision (float): Precision for the solver.
        migr_outwards (bool, optional): Specify the convention for the 
            definition of the migration rates. Defaults to False. 

    Returns:
        dict: Solution of the model.
    """
    if not(hasattr(popsize, "__len__")):
        popsize = [popsize]
    if not(hasattr(nu, "__len__")):
        nu = [nu]
    if not(m is None):
        if not(hasattr(m, "__len__")):
            m = [m, m]
        elif len(m) == 1:
            m = [m[0],m[0]]
    if Na is None:
        Na = "sum"
    
    check_requirements(time, popsize, K, nu, timeburnin, m, sla, precision)
    
    # Prepare the resolution
    solver_kwargs = dict(atol = precision, rtol = precision)
    par = dict(solver_kwargs = solver_kwargs, t = time, K = K, t_g = 1.0)

    if timeburnin is None:
        par["error_burnin_convergence"] = True
        par["burnin"] = 0.1 * par["t"] # increment the burnin time by 10% of the resolution time until convergence 
        par["extend_burnin"] = 250 # until we reach convergence or 25*time
    else:
        par["error_burnin_convergence"] = False
        par["burnin"] = timeburnin
        par["extend_burnin"] = 1



    if not(sla is None):
        par["nu"] = nu[0]
        par["N"] = popsize[0]
        par["sLA"] = sla
        sol = HAL_neutral_LA.solve_ODE_HAL_LA(**par)
    else:
        par["N1"] = popsize[0]
        if len(popsize) == 1:
            par["N2"] = popsize[0]
        else:
            par["N2"] = popsize[1]
        par["Na"] = Na
        par["nu1"] = nu[0]
        if len(nu) == 1:
            par["nu2"] = nu[0]
        else:
            par["nu2"] = nu[1]
        if not(m is None):
            if migr_outwards:
                par["m12"] = m[0] * popsize[0] / popsize[1]
                par["m21"] = m[1] * popsize[1] / popsize[0]
            else:
                par["m12"] = m[0]
                par["m21"] = m[1]
        sol = HAL_migr.solve_ODE_HAL_migr(**par)

    return sol
 
def save_pred(sol, outdir):
    # Save the output
    fp = outdir

    output_files = []
    summary = dict()
    for k in sol.keys():

        if isinstance(sol[k], np.ndarray):
            np.savetxt("{}/{}.txt".format(fp, k), sol[k])
            output_files.append("{}.txt".format(k))
        else:
            summary[k] = sol[k]
    with open("{}/_output_files.txt".format(fp), "w") as f:
        f.writelines(output_files)
    with open("{}/SUMMARY.json".format(fp), "w") as f:
        json.dump(summary, f, cls = utils.NpEncoder)
    
def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog = "HAL_predictions", 
        description = "Numerically resolves the equations of the holey adaptive landscape (HAL) of speciation.",
        epilog = "Once resolved, the results are saved in the output directory and can be plotted.",
        formatter_class=argparse.MetavarTypeHelpFormatter
        )


    parser.add_argument(
        "--time", 
        help = "Number of generations after split, > 0", 
        type = float, required = True)

    parser.add_argument(
        "--popsize", 
        help = "Population size(s), one or two values accepted, > 0. "
        "If one value given, the two populations are considered to have the same size.",
        nargs="+", type=float, required = True)
    
    parser.add_argument(
        "--K", 
        help = "Threshold for outbreeding depression, > 0", 
        type = float, 
        required = True)
    
    parser.add_argument(
        "--nu", 
        help = "Mutation rate(s), > 0", 
        type = float, 
        required = True, 
        nargs="+")

    parser.add_argument(
        "--plot",
        help = "Plot the solution", 
        action = "store_true")
    
    parser.add_argument(
        "--output", 
        help = "Path to store the solution, default to ./output", 
        type = str, 
        default = "./output")
    
    parser.add_argument(
        "--timeburnin", 
        help = "Manually specify the burnin time (before split), > 0. If not "
        "specified, the burnin time is automatically adjusted to reach the "
        "equilibrium (limited to 25*time before an error occurs). If specified, "
        "the equilibrium might not be reached.", 
        type = float)
    
    parser.add_argument(
        "--m", 
        help = "Specify (a) migration rate(s), >= 0. If two migration rates are "
        "specified, they are in order m12, m21. If one is specified, it is "
        "assumed that m12 = m21. If none specified, assumed zero.", 
        default = None, 
        nargs = "+", 
        type = float)
    
    parser.add_argument(
        "--Na", 
        help = "Ancestral population size. If not specified, automatically set to N1+N2.", 
        default = None, 
        type = float)
    
    parser.add_argument(
        "--sla", 
        help = "For the simulation with local adaptation, specify a coefficient of selection for local adaptation, >= 0", 
        default = None, 
        type = float)
    
    parser.add_argument(
        "--precision", 
        help = "Solver precision", 
        type = float, 
        default = 1e-6)
    
    parser.add_argument(
        "--migr_outwards", 
        help = "Specify if migration rates are defined as the probability of an "
        "individual to move to the other population (emigration). If not specified, " 
        "it is assumed that the migration rates are the probability that an "
        "individual was in the other population at the previous "
        "generation (the default).", 
        action = "store_true")
    
    return parser

def check_requirements(time, popsize, K, nu, timeburnin, m, sla, precision):
    if len(popsize) < 1 or len(popsize) > 2:
        raise ValueError("popsize must have one or two values, and {} were provided.".format(len(popsize)))
    elif len(popsize) == 1:
        N1, N2 = popsize[0], popsize[0]
    else:
        N1, N2 = popsize[0], popsize[1]
    if len(nu) < 1 or len(nu) > 2:
        raise ValueError("nu must have one or two values, and {} were provided.".format(len(nu)))
    elif len(nu) == 1:
        nu1, nu2 = nu[0], nu[0]
    else:
        nu1, nu2 = nu[0], nu[1]
    if not(m is None) and ((len(m) > 2) or (m[0] < 0) or (m[1] < 0)):
        raise ValueError("If specified --m must be one ore two positive or zero values")

    if not(sla is None):
        if not(m is None):
            raise NotImplementedError("Model with local adaptation not compatible with the model with migration.") 
        if (N1 != N2) or (nu1 != nu2):
            raise NotImplementedError("Model with local adaptation and different population sizes and/or different mutation rates is not implemented.")
        if sla < 0:
            raise ValueError("For the model with local adaptation, please specify a non-negative coefficient of selection for local adaptation with --sla")
    if not(timeburnin is None) and timeburnin <= 0:
        raise ValueError("If specified, timeburnin must be > 0.")

    if time <= 0 or K <= 0 or nu1 <= 0 or nu2 <= 0 or N1 <= 0 or N2 <= 0:
        raise ValueError("time, K, nu and popsize(s) must be > 0.")
    min_prec = 2.220446049250313e-14
    if precision < min_prec:
        raise ValueError("Precision too small, must be > " + str(min_prec))


def main():
    parser = create_parser()
    args = parser.parse_args()

    
    # Compute solution
    sol = HAL_prediction(**vars(args))

    # Save solution
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    with open(args.output + "/_call.json", "w") as f:
        json.dump(vars(args), f)
    save_pred(sol, args.output)

    # Plot the output
    if args.plot:
        fig, ax = plt.subplots()
        HAL_plot.plot_sol(ax, sol)
        plt.show()

if __name__ == '__main__':
    main()
    