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
from Parameters.parameter import NDIM
from tree.treewalk import get_minimum_distance_from_wall, add_ghost
from utility.integer_coordinates import get_distance_vector
from utility.utility import relative_density_error, norm


def compute_wvt_forces(Particles, Problem, NgbTree, step):
    "In this function the WVT-forces are computed"
    t0 = time()
    
    for particle in Particles:
        p = particle.ID
        particle.delta = zeros(NDIM)
        err = relative_density_error(particle)
        delta_fac = err/(1 + err)
        if particle.CloseToWall:
            min_dist_from_wall = get_minimum_distance_from_wall(particle, NgbTree)
        #now add contributions for all neighbors
        for no in particle.neighbors:
            ngb = NgbTree.Tp[no]
            n = ngb.ID
            if p == n:
                continue
            
            dist = get_distance_vector(particle.position, ngb.position, NgbTree)
            if particle.CloseToWall and add_ghost(particle, ngb, dist, \
                                                  min_dist_from_wall, NgbTree):
                #particle and ghost contributions cancel out
                continue
            r = norm(dist)
            h    = 0.5 * (particle.Hsml + ngb.Hsml)
            if NDIM == 1:
                wk = log(r/h + 1e-3)
            else:
                wk = (r/h + 1e-3)**(-(NDIM-1))
            particle.delta += h * wk * dist / r
        #reduce stepsize for particles with small error
        particle.delta *= step * delta_fac
        #reduce memory load
        particle.neighbors = list()
        
    t1 = time()
    Problem.Timer["WVT"] += t1 - t0