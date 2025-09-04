#!/usr/bin/env python3 
# Scripts for the plotting the solutions of the HAL model.
# Author: Pierre Veron

def plot_sol(ax, sol, legend = True, xlabel = True, ylabel = True, 
             ylabelright = True, legendkwargs = dict()):
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

    Returns:
        dict:
            * 'ax' contains the main axis
            * 'ax2' contains the secondary axis.
    """
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
    
    
    
    
    ax2.plot(sol["T"], sol["wb"], ls = "-.", color = "k")
    
    if sol["speciation"]:
        ax.fill_between(sol["T"], 0, 1, where= sol["wb"]  > 0.0,
                        color='gray', alpha=0.2, transform=ax.get_xaxis_transform(),
                        ec = "none", label = "speciation")
    ax.plot([],[], ls = "-.", color = "k", label = "$w_b$ (right axis)")
    ax.set_zorder(1)
    ax.set_frame_on(False) # make it transparent
    ax2.set_frame_on(True) # make sure there is any background
    if xlabel:
        ax.set_xlabel("Time")
    if ylabel:
        ax.set_ylabel("Genetic distance")
    if ylabelright:
        ax2.set_ylabel("Fitness")
    if legend:
        ax.legend(**legendkwargs)
    return dict(ax = ax, ax2 = ax2)