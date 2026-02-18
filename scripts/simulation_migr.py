#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Anaïs Spire and Pierre Veron
Date: 2026-02-18
Description: implements the stochastic simulation tool for the scenario with 
    migration
License: MIT License
"""

import argparse
import numpy as np
import json
from tqdm import tqdm
import math
import sys 
sys.path.append("scripts")
import utils
import Genotype


def main(): # Function used for running on the cluster with a set of parameters
    # load global parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indexsim", type = int)
    parser.add_argument("-p", "--paramfile", type = str)
    parser.add_argument("-o", "--outdir", type = str)
    args = parser.parse_args()

    savedist = False # save the distribution of genetic distances step by step (takes time+space)
    # load simulation parameters
    import pandas as pd
    param = pd.read_csv(args.paramfile, sep = ";", encoding = "utf-8")
    param = dict(param.loc[args.indexsim])
    np.random.seed(int(param['seed']))
    od = args.outdir
    print(param)
    n = int(param['n'])
    N = json.loads(param['N'])
    nu = float(param['nu'])
    m = float(param['m'])
    K = int(param['K'])
    n_gen = int(param['n_gen'])
    n_gen_burnin = int(param['n_gen_split'])
    recomb_rate = int(param['recomb_rate'])
    migr_model = param['migr_model']
    t_g = float(param['t_g'])

    simul(n, N, nu, m, K, n_gen, n_gen_burnin, recomb_rate, migr_model, t_g, savedist, od)
        
def simul(n, N, nu, m, K, n_gen, n_gen_burnin, recomb_rate, migr_model, t_g, savedist, od):
    """Simulate the HAL process with migration. 

    Args:
        n (int >0): number of populations
        N (list): population sizes
        nu (float >0): mutation rate
        m (float): migration rate
        K (float): threshold of incompatibility
        n_gen (int >0): number of generations in total (including burnin)
        n_gen_burnin (int >0): number of generations for the burnin phase (before split)
        recomb_rate (float >=0): recombination rate
        migr_model (str): type of migration, between 'island', 'stepping_stone' or
            'directional_stepping_stone'
        t_g (float >0): scaling factor for generations, set to 1 for counting time
            in generations
        savedist (bool): whether to save the distributions of genetic distances 
            through time (takes time and space) 
        od (str): directory to save the output of the simulations
    """
    if savedist:
        fp_dist_w = od + "/D_w_dist.txt" # dist_w before split
        meta_fp_dist_w = [] # dist_w after split
        for p in range(n):
            meta_fp_dist_w.append(od + "/D_w{}_dist.txt".format(p))

        fp_dist_b_global = od + "/D_b_dist_global.txt" # dist_b over all subpops (not very relevant)
        fp_dist_b = [] #dist_b_ij (two subpop by two subpop)
        for p1 in range(n):
            for p2 in range(p1+1, n):
                #index = int((p2 * (p2 - 1)) / 2 + p1)
                fp_dist_b.append(od + "/D_b{}{}".format(p1,p2))

        for fp in [fp_dist_w, *meta_fp_dist_w, fp_dist_b_global, *fp_dist_b]:
            file = open(fp, 'x')
            file.close()

    # initiate population
    gen = Genotype.Genotype()
    population = [gen.inherit() for i in range(sum(N))]

    D_w_mean, D_w_sd = np.zeros(n_gen_burnin), np.zeros(n_gen_burnin)

    # burnin
    for g in tqdm(range(n_gen_burnin), desc = "Burnin"):
        if savedist:
            d_dist = dict()
        # generate next generation
        population = simulate_generation(population, sum(N), nu, K, recomb_rate)


        # calculate within distance
        D = []
        for i in range(sum(N)):
            for j in range(i):
                d = population[i].nbDiff(population[j])
                D.append(d)
                if savedist:
                    utils.increment_dist(d_dist, d)
        D_w_mean[g] = np.mean(D)
        D_w_sd[g] = np.std(D)

        # store the distribution for this step
        if savedist:
            with open(fp_dist_w, 'a') as f:
                f.write(str(d_dist) + "\n")

    print("Burnin done")



    # split populations 
    npairs = [N[p] * (N[p]-1) // 2 for p in range(n)] # nb of potential couples inside a subpop

    metapop = [population[:N[0]]]
    for p in range(n-2):
        metapop.append(population[N[p]:N[p]+N[p+1]])
    metapop.append(population[sum(N)-N[n-2]:])

    meta_D_w_mean, meta_D_w_sd = [[np.zeros(n_gen - n_gen_burnin) for p in range(n)] for q in range(2)]
    D_b_mean, D_b_sd = [np.zeros(((n * (n - 1)) // 2, n_gen-n_gen_burnin)) for q in range(2)] # D_b_ij
    D_b_global_mean, D_b_global_sd = [np.zeros(n_gen - n_gen_burnin) for q in range(2)] # D_b_global

    Fertile_w = [np.zeros(n_gen - n_gen_burnin) for p in range(n)]
    Fertile_b = np.zeros(((n * (n - 1)) // 2, n_gen-n_gen_burnin))
    Fertile_b_global = np.zeros(n_gen - n_gen_burnin)

    # simulate each population
    g = n_gen_burnin
    speciation = False
    count_migrants = 0
    
    for g in tqdm(range(n_gen_burnin, n_gen), desc = "Simulation"):
        # migration step
        migrants = [[] for p in range(n)]

        if migr_model == 'island': # same model as in the predictions
            # loop over all individuals to see who's migrating and deleting them from their ancient pop
            for p in range(n):
                for gen in list(metapop[p]): # a genotype = an individual / list() function so that when we remove a genotype, we don't skip the next one
                    # Island model:  migration to any other subpop with probabilty m
                    p_iter = 0
                    migr = 0  # variable that allows us to know whether an individual has already migrated or not
                    while p_iter < n and migr == 0 : # Bernouilli : n=1
                        if p_iter != p :
                            if np.random.binomial(1, m): # probability m of migrating to any other subpop
                                metapop[p].remove(gen)
                                migrants[p_iter].append(gen)
                                migr = 1
                        p_iter += 1

        elif migr_model == 'stepping_stone':
             # on the borders (p = 0 or p=n-1)
            for gen in list(metapop[0]):
                if np.random.binomial(1, m):
                    metapop[0].remove(gen)
                    migrants[1].append(gen)
            for gen in list(metapop[n-1]):
                if np.random.binomial(1, m):
                    metapop[n-1].remove(gen)
                    migrants[n-2].append(gen)
            # not on the borders
            for p in range(1, n-1):
                for gen in list(metapop[p]): # a genotype = an individual / list() function so that when we remove a genotype, we don't skip the next one
                    # 1D stepping stone model
                    if np.random.binomial(1, m): # probability m of migrating to a neighboring subpop
                        metapop[p].remove(gen)
                        migrants[p-1].append(gen)
                    elif np.random.binomial(1, m): # probability m of migrating to a neighboring subpop
                        metapop[p].remove(gen)
                        migrants[p+1].append(gen)

        elif migr_model == 'directional_stepping_stone':
            migrants = [[] for p in range(n)]

            # on the borders (p = 0 or p=n-1)
            for gen in list(metapop[0]):
                if np.random.binomial(1, m):
                    metapop[0].remove(gen)
                    migrants[1].append(gen)

            for gen in list(metapop[n-1]):
                if np.random.binomial(1, m/2):
                    metapop[n-1].remove(gen)
                    migrants[n-2].append(gen)

            # not on the borders
            for p in range(1, n-1):
                for gen in list(metapop[p]): # a genotype = an individual / list() function so that when we remove a genotype, we don't skip the next one
                    # 1D stepping stone model
                        if np.random.binomial(1, m): # probability m/2 of migrating to a neighboring subpop
                            metapop[p].remove(gen)
                            migrants[p+1].append(gen)
                        elif np.random.binomial(1, m/2): # probability m/2 of migrating to a neighboring subpop
                            metapop[p].remove(gen)
                            migrants[p-1].append(gen)
        else:
            raise ValueError("migr_model must be 'island', 'stepping_stone' or 'directional_stepping_stone'")        
        # append the migrants to their new pops
        for p in range(n): # number of subpops
            metapop[p] += migrants[p]
            count_migrants += len(migrants[p])




        # reproduction step
        for p in range(n):
            if savedist:
                d_dist = dict()
            metapop[p] = simulate_generation(metapop[p], N[p], nu, K, recomb_rate)

            # calculate D_w for this population
            D = []
            count_fertile = 0
            for i in range(N[p]):
                for j in range(i):
                    d = metapop[p][i].nbDiff(metapop[p][j])
                    D.append(d)
                    count_fertile += int(d <= K)
                    if savedist:
                        utils.increment_dist(d_dist, d)
            meta_D_w_mean[p][g - n_gen_burnin] = np.mean(D)
            meta_D_w_sd[p][g - n_gen_burnin] = np.std(D)
            Fertile_w[p][g - n_gen_burnin] = count_fertile / npairs[p] # fraction of total fertile couples

            # store the distribution
            if savedist:
                with open(meta_fp_dist_w[p], 'a') as f:
                    f.write(str(d_dist) + "\n")


        # compute Db, Db_global (n>=2)
        D = []
        if savedist:
            d_dist = dict()
        count_fertile = 0
        weights = []

        for p1 in range(n):
            for i in range(N[p1]):
                x = metapop[p1][i]
                for p2 in range(p1+1, n):
                    for j in range(N[p2]):
                        y = metapop[p2][j]
                        # Count pairwise differences for D_b
                        d = x.nbDiff(y)
                        D.append(d)
                        count_fertile += int(d <= K)
                        if savedist:
                            utils.increment_dist(d_dist, d)
                    # Compute D_b_ij
                    index = int((p2 * (p2 - 1)) / 2 + p1)
                    D_b_mean[index][g - n_gen_burnin] = np.mean(D)
                    D_b_sd[index][g - n_gen_burnin] = np.std(D)
                    Fertile_b[index][g - n_gen_burnin] = count_fertile / (N[p1]*N[p2])
                    # store the distribution
                    if savedist:
                        with open(fp_dist_b[index], 'a') as f:
                            f.write(str(d_dist) + "\n")

                    # Compute D_b global
                    D_b_global_mean[g - n_gen_burnin] += D_b_mean[index][g - n_gen_burnin]*(N[p1]+N[p2])
                    weights.append(N[p1]+N[p2])

        # Compute D_b global
        D_b_global_mean[g - n_gen_burnin] /= sum(N)
        Fertile_b_global[g - n_gen_burnin] = count_fertile / (math.prod(N))

        if not(speciation) and (count_fertile == 0):
            gen_spec = g
            speciation = True 


    print("Simulation done.")
    
    # find speciation time 
    if not(speciation):
        gen_spec = np.inf

    time_spec = t_g * gen_spec
    iend = g - n_gen_burnin 

    # store results of the simulation
    res = dict()
    res['time_spec'] = float(time_spec)

    # conversion from numpy format to standard Python format (for json formatting)
    for key, value in res.items():
        if isinstance(value, np.generic):
            res[key] = value.item()

    # save the results 
    np.savetxt(od + "/D_w_mean.txt", D_w_mean)
    np.savetxt(od + "/D_w_sd.txt", D_w_sd)
    for p in range(n):
        np.savetxt(od + "/D_w{}_mean.txt".format(p), meta_D_w_mean[p][:iend])
        np.savetxt(od + "/D_w{}_sd.txt".format(p), meta_D_w_sd[p][:iend])
        np.savetxt(od + "/Fertile_w{}.txt".format(p), Fertile_w[p][:iend])
    for p1 in range(n):
        for p2 in range(p1+1, n):
            index = int((p2 * (p2 - 1)) / 2 + p1)
            np.savetxt(od + "/D_b{}{}_mean.txt".format(p1,p2), D_b_mean[index][:iend])
            np.savetxt(od + "/D_b{}{}_sd.txt".format(p1,p2), D_b_sd[index][:iend])
            np.savetxt(od + "/Fertile_b{}{}.txt".format(p1,p2), Fertile_b[index][:iend])
    np.savetxt(od + "/D_b_global_mean.txt", D_b_global_mean[:iend])
    np.savetxt(od + "/Fertile_b_global.txt", Fertile_b_global[:iend])
    
    with open(od + "/results.json", "w") as f: 
        json.dump(res, f, indent = 4, cls = utils.NpEncoder)


def simulate_generation(population, N, nu, K, recomb_rate, mutation_track_path = None,
                        generation = None, failed_cross_path = None):
    """ Simulate one generation of haploid reproduction under the holey adaptive landscape.

    Args:
        population (list of Genotype.Genotype objects): list of individuals in 
            the previous generation
        nu (float): per generation per individual probability of mutation
        K (float): threshold for interfertility 
        recomb_rate (float): recombination rate in number of events per generations 
        mutation_track_path (str or None): path for the file to store the track 
            of the mutations or None if the mutations are not tracked (default).
            File must exist. 
        generation (int or None): generation, only used if mutation_track_path is
            not None.
        failed_cross_path (str or None): path to store the tracker of failed crossings. 
            If None (default), the failed crossing are not tracked. 
            File is created. 

    Returns:
        list of Genotype.Genotype: the next generation.
    """
    newPop = []
    if not(failed_cross_path is None):
        with open(failed_cross_path, "w") as f:
            pass
    while len(newPop) < N:
        g1, g2 = np.random.choice(population, replace = False, size = 2)
        if g1.nbDiff(g2) < K:
            ch1, ch2 = g1.inherit(), g2.inherit()
            for ch in [ch1, ch2]:
                if np.random.binomial(1, nu):
                    id_mut = ch.addMutation()
                    if not(mutation_track_path is None):
                        with open(mutation_track_path, 'a') as f:
                            f.write("\n{};{}".format(id_mut, generation))
            if recomb_rate > 0: 
                nb_recomb = np.random.poisson(recomb_rate)
                for i in range(nb_recomb):
                    ch1.recombine(ch2)

            newPop.append(ch1)
            if len(newPop) < len(population): # for the last reproduction avoid having one more child
                newPop.append(ch2)
        elif not(failed_cross_path is None): #failed cross
            with open(failed_cross_path, "a") as f:
                f.write(g1.binaryRepr() + "\n")
                f.write(g2.binaryRepr() + "\n")
    return newPop


if __name__ == '__main__':
    main()