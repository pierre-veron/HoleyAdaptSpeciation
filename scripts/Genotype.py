#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pierre Veron
Date: 2026-02-18
Description: Implements the class Genotype, used to represent the haploid individuals
    in stochastic simulations
License: MIT License
"""
import numpy as np
import matplotlib.pyplot as plt

class InfiniteIter:
    # Iterable of "infinite" size for int, starting at 0.
    def __iter__(self):
        self.a = 0
        return self

    def __next__(self):
        x = self.a
        self.a += 1
        return x
    
    def size(self):
        return self.a

def readBinary(file):
    """Reads a file of genotypes and stores it to a dictionnary. 
    The file must be structured as such:
           a:010001010
           c:010001110
           f:000100100
           ... 
    The return would be:
            {'a':'010001010', 'c':'010001110', 'f':'000100100'}
    Args:
        file (str): file path. 

    Returns:
        dict: the genotypes
    """
    with open(file, "r", newline="\r\n") as f:
        gtypes = [line.rstrip() for line in f]
    out = {}
    for line in gtypes:
        x = line.split(":")
        out[x[0]] = x[1]
    return out 

alphabet = [chr(value) for value in range(97, 123)]
def alphaName(i):
    """
    Converts an nonnegative integer to a string with the order :
    a, b, c, ..., z, aa, ab, ..., az, ba, ..., zz, aaa...

    Parameters
    ----------
    i : nonnegative integer

    Returns
    -------
    str
        The corresponding string.

    """
    n = len(alphabet)
    if i < n:
        return alphabet[i]
    else: 
        return alphaName(i // n - 1) + alphabet[i % n]

def nbDiffBinary(bin1, bin2):
    """Calculates the number of pairwise differences between strings with same sizes.

    Args:
        bin1, bin2 (str): String with same sizes

    Returns:
        int: The number of pairwise differences between characters of bin1 and bin2.
    """
    n = 0
    for c in range(len(bin1)):
        n+= int(bin1[c] != bin2[c])
    return n

#%% Class Genotype
class Genotype:
    def __init__(self, mutiter = None):
        """An object representing a genotype. The mutiter argument is used to generate new mutations.

        Args:
            mutiter (iterator or None, optional): an iterator popping integers. 
                If None, a new one is created. Defaults to None.
        """
        self.mutations = set()
        self.loci = dict()
        
        if mutiter is None:
            mutiter = iter(InfiniteIter())
        self.mutiter = mutiter
    
    def __repr__(self):
        return "Genotype with {} mutations: {}".format(len(self.mutations), self.mutations)
    
    def addMutation(self):
        """Adds a mutation to the genotype on a locus that has not been mutated already. 
        """
        id_mut = next(self.mutiter)
        self.mutations.add(id_mut)
        self.loci[id_mut] = np.random.uniform()
        return id_mut

    def countMutations(self):
        return len(self.mutations)
    
    def recombine(self, other):
        crossover = np.random.uniform() # position of the crossover
        diff1 = self.mutations.difference(other.mutations)
        diff2 = other.mutations.difference(self.mutations)
        for id_mut in diff1:
            if self.loci[id_mut] > crossover:
                self.mutations.remove(id_mut)
                other.mutations.add(id_mut)
        for id_mut in diff2:
            if self.loci[id_mut] > crossover:
                other.mutations.remove(id_mut)
                self.mutations.add(id_mut)
    
    def inherit(self):
        """Create a new genotype with the same mutations.

        Returns:
            object of type Genotype
        """
        child = Genotype(mutiter = self.mutiter)
        child.mutations = self.mutations.copy()
        child.loci = self.loci
        return child
        
    def nbDiff(self, other):
        """Computes the number of differences among locii of two genotypes

        Args:
            other (object of class Genotype)

        Returns:
            int: the number of mutations that are different to the two genotypes.
        """
        return len(self.mutations.symmetric_difference(other.mutations))
    
    def binaryRepr(self):
        """Binary representation of a genotype, each locus is represented by a 0 if unmuted and 1 if it carries a mutation. 
        Position is random, unrelated to locus.

        Returns:
            str
        """
        out = ['0' for i in range(max(1, self.mutiter.size()))]
        for k in self.mutations:
            out[k] = '1'
        chrom = ""
        for loc in out:
            chrom += loc 
        return chrom

    def plot(self, ax = None, y = 0):
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot([0,1], [y,y], color = "k")
        for id_mut in self.mutations:
            ax.plot([self.loci[id_mut]], [y], ls = "", marker = "|")
            ax.text(x = self.loci[id_mut], y = y + 0.1, s = str(id_mut), ha = "center", size = "small")

    def size(self):
        """Returns the size of a genotype.

        Returns:
            int
        """
        return self.mutiter.size()
    
    def __str__(self):
        out = [alphaName(i) for i in range(max(1, self.mutiter.size()))]
        for k in self.mutations:
            out[k] = out[k].upper()
        chrom = out[0]
        for loc in out[1:]:
            chrom += "-" + loc
        return chrom

