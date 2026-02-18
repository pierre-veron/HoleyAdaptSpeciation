#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pierre Veron
Date: 2026-02-18
Description: useful functions used multiple times in the project
License: MIT License
"""


import numpy as np
import json
    
class NpEncoder(json.JSONEncoder): # used to properly save numpy object in json format
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def increment_dist(dist, x):
    """ Adds an occurence of the value x to an existing distribution. 

    Args:
        dist (dict): a dictionnary representing the distribution
        x (int): the value to increment
    """
    x = int(x)
    try:
        dist[str(x)] += 1
    except:
        dist[str(x)] = 1