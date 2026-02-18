Main scripts to run the simulations and the predictions. 
* `Genotype.py`: defines the class `Genotype` that represents the individuals used in the simulation framework.
* `HAL_general.py`: implements functions used in the general HAL framework.
* `HAL_migr.py`: implements the deterministic prediction tool for the scenario with migration (parapatric). See in particular the function `solve_ODE_HAL_migr`. 
* `HAL_neutral_LA.py` implements the deterministic prediction tool for the allopatric scenarios (neutral or with local adaptation). See in particular the function `solve_ODE_HAL_LA`.
* `HAL_plot.py`: implements useful functions to plot the results of the simulations and the predictions.
* `simulation_migr.py`: implements the stochastic simulation tool for the scenario with migration. The simulation is done with the function `simul`.
* `simulation_neutral_LA.py`: implements the stochastic simulation tool for the allopatric scenarios (neutral or with local adaptation). The simulation is done with the function `simul`. 
* `utils.py`: some useful functions that are used multiple times throughout the project. 