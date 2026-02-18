#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Anaïs Spire, Agathe Chave-Lucas and Pierre Veron
Date: 2026-02-18
Description: funnctions to plot the results of simulations or predictions of the 
    HAL model
License: MIT License
"""

import numpy as np
import json

def plot_sol(ax, sol, legend = True, xlabel = True, ylabel = True, 
             ylabelright = True, legendkwargs = dict(), plot_compatibility = True):
    """ Plot the time dynamics of the solution of a HAL model.

    Args:
        ax (matplotlib.axes.Axes): axes on which to plot the solution
        sol (dict): solution of the HAL model, as returned by the functions 
            HAL_neutral_LA.solve_ODE_HAL_LA or HAL_migr.solve_ODE_HAL_migr
        legend (bool, optional): add a legend. Defaults to True.
        xlabel (bool, optional): add a xlabel. Defaults to True.
        ylabel (bool, optional): add a y label. Defaults to True.
        ylabelright (bool, optional): add a y label for the secondary axis. 
            Defaults to True.
        legendkwargs (dict, optional): kwargs passed to ax.legend(). 
            Defaults to dict().
        plot_compatibility (bool, optional): plot the between-populations 
            compatibility. Defaults to True. 

    Returns:
        dict:
            * 'ax' contains the main axis
            * 'ax2' contains the secondary axis.
    """
    if plot_compatibility:
        ax2 = ax.twinx()
    
    
    if "Dw1" in sol.keys():
        ax.plot(sol['T_burnin'], sol['Dw_burnin'], label = "$D_{w,\\text{anc}}$", 
                color = "#98D4E2")
        ax.plot(sol["T"], sol["Dw1"], color = "#98D4E2", ls = ":", label = "$D_{w1}$")
        ax.plot(sol["T"], sol["Dw2"], color = "#98D4E2", ls = "--", label = "$D_{w2}$")
    else:
        ax.plot(sol['T_burnin'], sol['Dw_burnin'], label = "$D_{w}$",color = "#98D4E2")
        ax.plot(sol["T"], sol["Dw"], color = "#98D4E2")
    ax.plot(sol["T"], sol["Db"], label = "$D_b$", color = "#B90845")
    
    if "k" in sol.keys():
        ax.plot(sol["T"], sol["k"], label = "$k$", color = "#8D5E2A")
    
    
    
    if plot_compatibility:
        ax2.plot(sol["T"], sol["wb"], ls = "-.", color = "k")
    
    if sol["speciation"]:
        ax.fill_between(sol["T"], 0, 1, where= sol["wb"]  > 0.0,
                        color='gray', alpha=0.2, transform=ax.get_xaxis_transform(),
                        ec = "none", label = "speciation")
    ax.plot([],[], ls = "-.", color = "k", label = "$w_b$ (right axis)")
    ax.set_zorder(1)
    
    if plot_compatibility:
        ax.set_frame_on(False) # make it transparent
        ax2.set_frame_on(True) 
    if xlabel:
        ax.set_xlabel("Time")
    if ylabel:
        ax.set_ylabel("Genetic distance")
    if ylabelright and plot_compatibility:
        ax2.set_ylabel("Fitness")
    if legend:
        ax.legend(**legendkwargs)
    if plot_compatibility:
        return dict(ax = ax, ax2 = ax2)
    return dict(ax = ax)


def load_from_files(fp):
    """ Load the solution of a HAL resolution from the folder

    Args:
        fp (str): path to the output directory

    Returns:
        dict: solution of the HAL model. 
    """
    sol = dict()
    with open(f"{fp}/_output_files.txt", "r") as f:
        output_files = f.readlines()
    for file in output_files:
        fname = file.rstrip()
        values = np.loadtxt(f"{fp}/{fname}")
        sol[fname.replace(".txt","")] = values
    with open(f"{fp}/SUMMARY.json", "r") as f:
        summary = json.load(f)
    sol.update(summary)
    
    return sol 

def plot_from_files(ax, fp, legend = True, xlabel = True, ylabel = True, 
             ylabelright = True, legendkwargs = dict()):
    """ Plot the time dynamics of the solution of a HAL model from the folder 
    containing the solution.

    Args:
        ax (matplotlib.axes.Axes): axes on which to plot the solution
        fp (str): path of the folder containing the solution.
        legend (bool, optional): add a legend. Defaults to True.
        xlabel (bool, optional): add a xlabel. Defaults to True.
        ylabel (bool, optional): add a y label. Defaults to True.
        ylabelright (bool, optional): add a y label for the secondary axis. 
            Defaults to True.
        legendkwargs (dict, optional): kwargs passed to ax.legend(). 
            Defaults to dict().

    Returns:
        dict:
            * 'ax' contains the main axis
            * 'ax2' contains the secondary axis.
    """
    
    sol = load_from_files(fp)
    return plot_sol(ax, sol, legend, xlabel, ylabel, ylabelright, legendkwargs)


def plot_sim_migr_from_files(ax, dir):
    Db = np.loadtxt("{}/D_b01_mean.txt".format(dir))
    Dw0 = np.loadtxt("{}/D_w0_mean.txt".format(dir))
    Dw1 = np.loadtxt("{}/D_w1_mean.txt".format(dir))
    Dw_burnin = np.loadtxt("{}/D_w_mean.txt".format(dir))
    
    with open("{}/results.json".format(dir), "r") as f:
        results = json.load(f)
    
    time_spec = results["time_spec"]
    
    T_burnin = np.arange(len(Dw_burnin))
    T = len(Dw_burnin) + np.arange(len(Db))
    
    ax.plot(T_burnin, Dw_burnin, label = "sim $D_w$", color = "#7CAEBA")
    ax.plot(T, Dw0, label = "sim $D_{w1}$", color = "#7CAEBA", ls = "-.")
    ax.plot(T, Dw1, label = "sim $D_{w2}$", color = "#7CAEBA", ls = "--")
    ax.plot(T, Db, label = "sim $D_b$", color = "#F8528C")
    
    if np.isfinite(time_spec):
        ax.plot([time_spec],[0], ls = "", marker = 7, label = "sim speciation", color = "k")

def plot_sim_neutal_la_from_files(ax, dir):
    Db = np.loadtxt("{}/D_b_mean.txt".format(dir))
    Dw0 = np.loadtxt("{}/D_w0_mean.txt".format(dir))
    Dw1 = np.loadtxt("{}/D_w1_mean.txt".format(dir))
    Dw_burnin = np.loadtxt("{}/D_w_mean.txt".format(dir))
    
    with open("{}/results.json".format(dir), "r") as f:
        results = json.load(f)
    
    time_spec = results["time_spec"]
    
    T_burnin = np.arange(len(Dw_burnin))
    T = len(Dw_burnin) + np.arange(len(Db))
    
    ax.plot(T_burnin, Dw_burnin, label = "sim $D_w$", color = "#7CAEBA")
    ax.plot(T, Dw0, label = "sim $D_{w1}$", color = "#7CAEBA", ls = "-.")
    ax.plot(T, Dw1, label = "sim $D_{w2}$", color = "#7CAEBA", ls = "--")
    ax.plot(T, Db, label = "sim $D_b$", color = "#F8528C")
    
    if np.isfinite(time_spec):
        ax.plot([time_spec],[0], ls = "", marker = 7, label = "sim speciation", color = "k")