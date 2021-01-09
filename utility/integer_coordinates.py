#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 12:30:26 2021

@author: leonard
"""
from Parameters.constants import BITS_FOR_POSITIONS
from Parameters.parameter import NDIM
from numpy import zeros 

def convert_to_int_position(position, Problem):
    """
    converts a coordinate with double accuracy to an integer position
    The left boundary corresponds to 0 and the right boundary to 2^32
    """
    intpos = zeros(NDIM, dtype = int)
    for i in range(NDIM):
        intpos[i] += int(position[i]/Problem.FacIntToCoord[i])
    return intpos

def convert_to_phys_position(intpos, Problem):
    "converts an integer coordinate vector to a physical position vector"
    pos = zeros(NDIM, dtype = float)
    for i in range(NDIM):
        pos[i] += intpos[i] * Problem.FacIntToCoord[i]
    return pos
        
        
def get_distance_vector(x, y, Problem):
    "returns the minimum distance vector taking into account periodic boundaries"
    dx = zeros(NDIM, dtype = int)
    for i in range(NDIM):
        if Problem.Periodic[i]:
            if x[i] <= y[i]:
                if abs(x[i]-y[i]) < (1 << BITS_FOR_POSITIONS) + x[i] - y[i]:
                    dx[i] += x[i] - y[i]
                else:
                    dx[i] += (1 << BITS_FOR_POSITIONS) + x[i] - y[i]
            else:
                if x[i] - y[i] < abs(x[i] - y[i] - (1 << BITS_FOR_POSITIONS)):
                    dx[i] += x[i] - y[i]
                else:
                    dx[i] += x[i] - y[i] - (1 << BITS_FOR_POSITIONS) 
        else:
            dx[i] += x[i] - y[i]
    return dx