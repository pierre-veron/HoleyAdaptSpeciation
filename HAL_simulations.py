#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pierre Veron
Date: 2026-02-18
Description: Wrapper script to run the HAL stochastic simulations from command-line.
License: MIT License
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from scripts import simulation_migr, simulation_neutral_LA, HAL_plot
import HAL_predictions


def main():
    parser = create_parser()  
    args = parser.parse_args()
    
    # Validations
    if args.time <= 0 or args.timeburnin <= 0 or args.K <= 0 or args.nu <= 0:
        raise ValueError("time, timeburnin, K, nu must be > 0")
    if args.recombrate < 0:
        raise ValueError("recombrate must be >= 0")
    if len(args.popsize) < 1 or len(args.popsize) > 2:
        raise ValueError("popsize must have one or two values, and {} were provided.".format(len(args.popsize)))
    for n in args.popsize:
        if n <= 0:
            raise ValueError("popsize(s) must be > 0 integer(s)")
    
    if not(args.m is None):
        if not(args.sla is None):
            raise NotImplementedError("Models with migration and with local adaptation are incompatible.")
    
    migr = not(args.m is None)
    la = not(args.sla is None)
        
    
    

    if len(args.popsize) == 1:
        popsize = [args.popsize[0],args.popsize[0]]
    else:
        popsize = [args.popsize[0],args.popsize[1]]
    
    # Validations specific to the predictions
    if args.pred:
        if la and popsize[0] != popsize[1]:
            raise NotImplementedError("For the predictions, the model with local adaptation and N1 != N2 is not implemented. It is still possible to run the simulations in this configurations. For this, remove --pred to run only the simulations.")
        min_prec = 2.220446049250313e-14
        if args.precision < min_prec:
            raise ValueError("Precision too small, must be > " + str(min_prec))
        
    par = dict(nu=args.nu, K=args.K, n_gen=args.timeburnin+args.time, 
               n_gen_burnin=args.timeburnin, recomb_rate = args.recombrate, 
               t_g = 1.0, od = args.output)
    
    # Prepare output 
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Migration
    if migr:
        if args.migrmodel is None:
            migr_model = "island"
        else:
            migr_model = migr_model
        par = {**par, **dict(n = 2, N = popsize, m = args.m, 
                             migr_model = migr_model, savedist = False)}
        simulation_migr.simul(**par)
    else: 
        # Neutral or LA
        if la:
            sla = args.sla
        else:
            sla = 0.0
        par = {**par, **dict(N1 = popsize[0], N2 = popsize[1], sLA = sla)}
        simulation_neutral_LA.simul(**par)
        
    if args.pred:
        par_predict = vars(args)
        par_predict["Na"] = np.sum(par_predict["popsize"])
        sol = HAL_predictions.HAL_prediction(**par_predict)
        
        od = args.output + "/predictions"
        if not os.path.exists(od):
            os.makedirs(od)

        HAL_predictions.save_pred(sol, od)
    
    if args.plot:
        fig, ax = plt.subplots()
        
        if migr:
            HAL_plot.plot_sim_migr_from_files(ax, args.output)
        else:
            HAL_plot.plot_sim_neutal_la_from_files(ax, args.output)
        
        if args.pred:
            HAL_plot.plot_sol(ax, sol, plot_compatibility=False)
        ax.legend()
        plt.show()

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog = "HAL_simulations", 
                                     description = "Run simulation of the holey adaptive landscape (HAL) of speciation.",
                                     epilog = "",
                                     formatter_class=argparse.MetavarTypeHelpFormatter)


    
    parser.add_argument("--time", help = "Time of the resolution (after split), int > 0", 
                        type = int, required = True)
    parser.add_argument("--timeburnin", help = "Burnin time (before split), int > 0.",
                        type = int, required=True)
    parser.add_argument("--popsize", help = "Population size(s), one or two values accetped, int > 0. If one value given, the two populations are considered to have the same size.",
                        nargs="+", type=int, required = True)
    parser.add_argument("--K", help = "Threshold for outbreeding depression, > 0", 
                        type = float, required = True)
    parser.add_argument("--nu", help = "Mutation rate, > 0", type = float, required = True)
    parser.add_argument("--recombrate", help = "Recombination rate, >= 0", type = float, required = True)
    parser.add_argument("--plot", help = "Plot the solution", action = "store_true")
    parser.add_argument("--pred", help = "Run also the deterministic prediction. In this case, the results are stored in <output>/predictions/", 
                     action = "store_true")
    parser.add_argument("--output", help = "Path to store the solution, default to ./output", 
                        type = str, default = "./output")
    
    parser.add_argument("--m", help = "In the model with migration, specify a migration rate, >= 0", 
                        default = None, type = float)
    parser.add_argument("--sla", help = "In the model with local adaptation, specify a coefficient of selection for local adaptation, >= 0", 
                        default = None, type = float)
    parser.add_argument("--migrmodel", help = "Migration model: island, directional_stepping_stone or stepping_stone. If unspecified, island is assumed.", 
                     default = None, type = str)
    
    parser.add_argument("--precision", help = "Solver precision (only if --pred is used)", type = float, 
                        default = 1e-6)
    return parser

if __name__ == '__main__':
    main()