#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pierre Veron
Date: 2026-02-18
Description: implements the stochastic simulation tool for the neutral and 
    adaptive scenarios
License: MIT License
"""
import argparse
import numpy as np
from tqdm import tqdm
import json
import sys
sys.path.append("scripts")
import utils
import Genotype

def main():
    import pandas as pd
    # load global parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indexsim", type = int)
    parser.add_argument("-p", "--paramfile", type = str)
    parser.add_argument("-o", "--outdir", type = str)
    args = parser.parse_args()

    # load simulation parameters
    param = pd.read_csv(args.paramfile, sep = ";", encoding = "utf-8")
    param = dict(param.loc[args.indexsim])
    np.random.seed(int(param['seed']))
    od = args.outdir
    N1 = int(param['N1'])
    N2 = int(param['N2'])
    nu = param['nu']
    K = param['K']
    n_gen = int(param['n_gen'])
    n_gen_burnin = int(param['n_gen_split'])
    recomb_rate = param['recomb_rate']
    t_g = param['t_g']
    sLA = param['sLA']
    
    simul(N1, N2, nu, K, n_gen, n_gen_burnin, recomb_rate, sLA, t_g, od)

def simul(N1, N2, nu, K, n_gen, n_gen_burnin, recomb_rate, sLA, t_g, od):
    """Simulate the HAL process in allopatry with or without local adaptation. 

        Args:
            N, N2 (int): population sizes
            nu (float >0): mutation rate
            K (float): threshold of incompatibility
            n_gen (int >0): number of generations in total (including burnin)
            n_gen_burnin (int >0): number of generations for the burnin phase (before split)
            recomb_rate (float >=0): recombination rate
            sLA (float >=0): coefficient of local adaptation
            t_g (float >0): scaling factor for generations, set to 1 for counting time
                in generations
            od (str): directory to save the output of the simulations
        """
    # Trackers of the distribution of distance 
    fp_dist_w = od + "/D_w_dist.txt"
    meta_fp_dist_w = [od + "/D_w0_dist.txt", od + "/D_w1_dist.txt"]
    fp_dist_b = od + "/D_b_dist.txt"

    for fp in [fp_dist_w, *meta_fp_dist_w, fp_dist_b]:
        with open(fp, "w") as f: # create or delete the file
            pass

    # initiate population
    gen = Genotype.Genotype()
    population = [gen.inherit() for i in range(N1+N2)]

    D_w_mean, D_w_sd = np.zeros(n_gen_burnin), np.zeros(n_gen_burnin)

    # burnin
    fertile_pairs = None
    for k in tqdm(range(n_gen_burnin), desc = "Burnin"):
        # generate next generation
        population = simulate_generation(population, nu, K, recomb_rate, 
                                          mutation_track_path=None,
                                          generation=k, 
                                          fertile_pairs=fertile_pairs,
                                          sLA=sLA)

        parsing = parse_intra_pop_get_dist_fertile(population, K)
        fertile_pairs = parsing["fertile_pairs"]

        D_w_mean[k] = parsing["meanD"]
        D_w_sd[k] = parsing["sdD"]

        # store the distribution for this step
        with open(fp_dist_w, 'a') as f:
            f.write(str(parsing["d_dist"]) + "\n")

    print("Burnin done")

    # split populations 
    N = [N1, N2]
    npairs = [N[p] * (N[p]-1) / 2 for p in range(2)]
    metapop = [population[:N1], population[N1:]]
    meta_D_w_mean, meta_D_w_sd = [[np.zeros(n_gen - n_gen_burnin) for p in range(2)] for q in range(2)]
    D_b_mean, D_b_sd = [np.zeros(n_gen - n_gen_burnin) for p in range(2)]

    Fertile_w = [np.zeros(n_gen - n_gen_burnin) for p in range(2)]
    Fertile_b = np.zeros(n_gen - n_gen_burnin)

    # simulate each populations 
    k = n_gen_burnin
    speciation = False
    fertile_pairs = [None, None]
    pbar = tqdm(total = n_gen - n_gen_burnin, desc = "Simulation")
    while k < n_gen:
        for p in range(2):
            metapop[p] = simulate_generation(metapop[p], nu, K, recomb_rate,
                                             mutation_track_path=None,
                                             generation=k, 
                                             fertile_pairs=fertile_pairs[p], 
                                             sLA=sLA)
            # calculate D_w for this population
            parsing = parse_intra_pop_get_dist_fertile(metapop[p], K)
            fertile_pairs[p] = parsing["fertile_pairs"]

            meta_D_w_mean[p][k - n_gen_burnin] = parsing["meanD"]
            meta_D_w_sd[p][k - n_gen_burnin] = parsing["sdD"]
            Fertile_w[p][k - n_gen_burnin] = parsing["count_fertile"] / npairs[p] # fraction of total fertiles couples

            # store the distribution
            with open(meta_fp_dist_w[p], 'a') as f:
                f.write(str(parsing["d_dist"]) + "\n")

        # pairs between the two populations
        parsing = parse_extra_pop_get_dist_fertile(*metapop,K)
        D_b_mean[k - n_gen_burnin] = parsing["meanD"]
        D_b_sd[k - n_gen_burnin] = parsing["sdD"]
        Fertile_b[k - n_gen_burnin] = parsing["count_fertile"] / (N1 * N2)

        # store the distribution
        with open(fp_dist_b, 'a') as f:
            f.write(str(parsing["d_dist"]) + "\n")

        if not(speciation) and (parsing["count_fertile"] == 0):
            gen_spec = k
            speciation = True

        k+= 1
        pbar.update(1)
    pbar.close()

    print("Simulation done.")

    if not(speciation):
        gen_spec = np.inf

    time_spec = t_g * gen_spec
    iend = k - n_gen_burnin 

    # store results of the simulation
    res = dict()
    res['time_spec'] = time_spec

    # save the results 
    np.savetxt(od + "/D_w_mean.txt", D_w_mean)
    np.savetxt(od + "/D_w_sd.txt", D_w_sd)
    for p in range(2):
        np.savetxt(od + "/D_w{}_mean.txt".format(p), meta_D_w_mean[p][:iend])
        np.savetxt(od + "/D_w{}_sd.txt".format(p), meta_D_w_sd[p][:iend])
        np.savetxt(od + "/Fertile_w{}.txt".format(p), Fertile_w[p][:iend])
    np.savetxt(od + "/D_b_mean.txt", D_b_mean[:iend])
    np.savetxt(od + "/D_b_sd.txt", D_b_sd[:iend])
    np.savetxt(od + "/Fertile_b.txt", Fertile_b[:iend])

    with open(od + "/results.json", "w") as f:
        json.dump(res, f, cls = utils.NpEncoder)


def simulate_generation(population, nu, K, recomb_rate, mutation_track_path = None,
                        generation = None, fertile_pairs = None, sLA = 0.0):
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
        fertile_pairs (list of tuples or None): list of indices of fertile pairs
            that will be used to build the new population. If None it will be 
            built before (but takes more time). 

    Returns:
        list of Genotype.Genotype: the next generation.
    """
    newPop = []
    n = len(population)
    # Build list of interfertile parents 
    if fertile_pairs is None:
        fertile_pairs = parse_intra_pop_get_dist_fertile(population, K)['fertile_pairs']
    p = None # default  chance of picking pairs, without fitness = uniform
    if sLA > 0:
        p = get_fitness_pair(population, fertile_pairs, sLA)
    mating_pairs = np.random.choice(len(fertile_pairs), replace=True, size = (n+1)//2, p = p)
    for i_pair in mating_pairs:
        pair = fertile_pairs[i_pair]
        g1, g2 = population[pair[0]], population[pair[1]]
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
        if len(newPop) < n: # for the last reproduction avoid having one more child
            newPop.append(ch2)

    return newPop

def parse_intra_pop_get_dist_fertile(population, K):
    fertile_pairs = []
    D = []
    d_dist = dict()
    count_fertile = 0
    for i1, x1 in enumerate(population):
        for i2 in range(i1):
            x2 = population[i2]
            d = x1.nbDiff(x2)
            utils.increment_dist(d_dist, d)
            D.append(d)
            if d <= K:
                fertile_pairs.append((i1, i2))
                count_fertile += 1
    
    # compute some statistics
    meanD = np.mean(D)
    sdD = np.std(D)

    return dict(fertile_pairs = fertile_pairs, count_fertile = count_fertile, 
                d_dist = d_dist, meanD = meanD, sdD = sdD)

def get_fitness_pair(population, fertile_pairs, sLA):
    """ Calculate fitness of each fertile pairs in case with local adapatation.

    Args:
        population (list of Genotype.Genotype objects): population 
        fertile_pairs (list of tuples): list of fertile pairs in the population.
            Elements of the list are (i,j) with i and j indices of individuals in
            the population.
        single_mut_log_fitness (float): exp(sLA), log fitness induced by a 
            single adapted locus.

    Returns:
        np.array, same size as fertile_pairs: normalized fitness of each pair. 
    """
    nb_locus_adap = [g.countMutations() for g in population]
    fit_pair = np.zeros(len(fertile_pairs))
    for i_pair, pair in enumerate(fertile_pairs):
        na1, na2 = [nb_locus_adap[pair[i]] for i in [0,1]]
        fit_pair[i_pair] = np.power(1+sLA, na1+na2) #np.exp(sLA * (na1+na2))
    return fit_pair / np.sum(fit_pair)

def parse_extra_pop_get_dist_fertile(pop1, pop2, K):
    D = []
    d_dist = dict()
    count_fertile = 0
    for x1 in pop1:
        for x2 in pop2:
            d = x1.nbDiff(x2)
            utils.increment_dist(d_dist, d)
            D.append(d)
            if d <= K:
                count_fertile += 1
    
    # compute some statistics
    meanD = np.mean(D)
    sdD = np.std(D)

    return dict(count_fertile = count_fertile,  d_dist = d_dist, meanD = meanD, sdD = sdD)

if __name__ == "__main__":
    main()