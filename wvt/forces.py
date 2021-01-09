#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 11:39:15 2021

@author: leonard
"""
#public libraries
from numpy import zeros
from numpy.linalg import norm
#custom libraries
from Parameters.parameter import BiasCorrection, NDIM, Npart
from Parameters.constants import DESNNGB, NORM_COEFF
from tree.forcewalk import find_neighbors
from utility.integer_coordinates import convert_to_phys_position
from utility.utility import relative_density_error
from sph.Kernel import kernel

def compute_wvt_weights(Particles, Problem, Functions, avg_boxsize):
    "In this function the weights of the WVT are computed"
    v_sph_sum = 0
    hsml  = zeros(Npart)
    for particle in Particles:
        rho = Functions.Density_func(particle, Problem, BiasCorrection)
        particle.Rho_Model = rho
        hsml[particle.ID]  = (DESNNGB * Problem.Mpart / rho / NORM_COEFF)**(1/NDIM)
        v_sph_sum         += hsml[particle.ID]**NDIM
    norm_hsml = avg_boxsize * (DESNNGB / v_sph_sum / NORM_COEFF)**(1/NDIM)
    hsml *= norm_hsml
    return hsml

def compute_wvt_forces(Particles, Problem, Functions, NgbTree, hsml, step):
    "In this function the WVT-forces are computed"
    delta = zeros((Npart, NDIM))
    for particle in Particles:
        p   = particle.ID
        err = relative_density_error(particle, Problem, Functions)
        delta_fac = err/(1 + err)
        #see if we even need to update the list of neighbors
        if particle.Hsml < hsml[p]:
            find_neighbors(particle, Problem, NgbTree, hsml[p])
        #now add contributions for all neighbors
        for [neighbor, intDistance] in particle.neighbors:
            n = neighbor.ID
            if p == n or neighbor.IsGhost:
                continue
            dist = convert_to_phys_position(intDistance, Problem)
            r = norm(dist)
            h    = 0.5 * (hsml[p] + hsml[n])
            if r > h:
                continue
            wk = kernel(r/h, h) * h**NDIM
            delta[p] += h * wk * dist / r
        #reduce stepsize for particles with small error
        delta[p] *= step * delta_fac
    return delta