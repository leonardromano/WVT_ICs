#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:33:16 2021

@author: leonard
"""
from numpy import sqrt
from time import time

from Parameters.constants import LARGE_NUM
from utility.utility import relative_density_error
from Parameters.parameter import Npart


def compute_l1_error(Particles, Problem, Functions):
    t0 = time()
    
    err_min = LARGE_NUM
    err_max = 0.
    err_mean = 0.
    err_sigma = 0.
    for particle in Particles:
        err        = relative_density_error(particle, Problem, Functions)
        err_min    = min(err, err_min)
        err_max    = max(err, err_max)
        err_mean  += err
        err_sigma += err * err
    err_mean /= Npart
    err_sigma = sqrt(err_sigma/Npart - err_mean * err_mean)
    
    t1 = time()
    Problem.Timer["L1-ERROR"] += t1 - t0
    return err_min, err_max, err_mean, err_sigma