#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:33:16 2021

@author: leonard
"""
from numpy import sqrt, array
from numpy.linalg import norm

from Parameters.constants import LARGE_NUM
from utility.utility import relative_density_error
from Parameters.parameter import Npart


def compute_l1_error(Particles, Problem, Functions):
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
    return err_min, err_max, err_mean, err_sigma

def compute_field_stats(field):
    s_min = LARGE_NUM
    s_max = 0.
    s_mean = 0.
    s_sigma = 0.
    for i in range(Npart):
        v        = norm(field[i])
        s_min    = min(v, s_min)
        s_max    = max(v, s_max)
        s_mean  += v
        s_sigma += v * v
    s_mean /= Npart
    s_sigma = sqrt(s_sigma/Npart - s_mean * s_mean)
    return array([s_min, s_max, s_mean, s_sigma])