#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 12:59:39 2021

@author: leonard
"""
from time import time

from tree.treewalk import density
from Parameters.constants import DESNNGB, BITS_FOR_POSITIONS, NORM_COEFF
from Parameters.parameter import NDIM

def find_sph_quantities(Particles, Problem, NgbTree):
    t0 = time()
    #now walk the tree to determine the weights and sph-quantities
    density(Particles, NgbTree)
    
    t1 = time()
    Problem.Timer["DENSITY"] += t1 - t0

def initial_guess_hsml(Particles, NgbTree):
    "computes an initial guess for the smoothing lengths"
    for particle in Particles:
        no = NgbTree.Father[particle.ID]
        while(10 * DESNNGB * NgbTree.Mpart > NgbTree.get_nodep(no).Mass):
            p = NgbTree.get_nodep(no).father
            if p < 0:
                break
            no = p
        if(NgbTree.get_nodep(no).level > 0):
            length = (1 << (BITS_FOR_POSITIONS - NgbTree.get_nodep(no).level)) * \
                     NgbTree.FacIntToCoord.max()
        else:
            length = NgbTree.Boxsize.max()

        particle.Hsml = (DESNNGB * NgbTree.Mpart / NgbTree.get_nodep(no).Mass \
                         / NORM_COEFF)**(1/NDIM) * length