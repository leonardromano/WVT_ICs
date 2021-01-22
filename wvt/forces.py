#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:39:15 2021

@author: leonard
"""
#public libraries
from numpy import zeros, log
from time import time

#custom libraries
from Parameters.parameter import BiasCorrection, NDIM
from Parameters.constants import DESNNGB, NORM_COEFF
from tree.forcewalk import find_neighbors
from utility.integer_coordinates import convert_to_phys_position
from utility.utility import relative_density_error, norm

def compute_wvt_weights(Particles, Problem, Functions, avg_boxsize):
    "In this function the weights of the WVT are computed"
    t0 = time()
    
    v_sph_sum = 0
    for particle in Particles:
        rho = Functions.Density_func(particle, Problem, BiasCorrection)
        particle.Rho_Model = rho
        particle.Hwvt  = (DESNNGB * Problem.Mpart / rho / NORM_COEFF)**(1/NDIM)
        v_sph_sum     += particle.Hwvt**NDIM
    norm_hwvt = avg_boxsize * (DESNNGB / v_sph_sum / NORM_COEFF)**(1/NDIM)
    #now normalise
    for particle in Particles:
        particle.Hwvt *= norm_hwvt
        
    t1 = time()
    Problem.Timer["WVT"] += t1 - t0

def compute_wvt_forces(Particles, Problem, Functions, NgbTree, step):
    "In this function the WVT-forces are computed"
    t0 = time()
    
    for particle in Particles:
        p = particle.ID
        particle.delta = zeros(NDIM)
        err = relative_density_error(particle, Problem, Functions)
        delta_fac = err/(1 + err)
        #see if we even need to update the list of neighbors
        if particle.Hsml < particle.Hwvt:
            find_neighbors(particle, Problem, NgbTree)
        #now add contributions for all neighbors
        for [neighbor, intDistance] in particle.neighbors:
            n = neighbor.ID
            if p == n:
                continue
            if neighbor.IsGhost:
                #if we have a ghost need to make sure Hwvt has the right value
                neighbor.Hwvt = Particles[n].Hwvt
            dist = convert_to_phys_position(intDistance, Problem)
            r = norm(dist)
            h    = 0.5 * (particle.Hwvt + neighbor.Hwvt)
            if r > h:
                continue
            if NDIM == 1:
                wk = log(r/h + 1e-3)
            else:
                wk = (r/h + 1e-3)**(-(NDIM-1))
            particle.delta += h * wk * dist / r
        #reduce stepsize for particles with small error
        particle.delta *= step * delta_fac
        
    t1 = time()
    Problem.Timer["WVT"] += t1 - t0