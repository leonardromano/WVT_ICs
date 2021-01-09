#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 12:59:39 2021

@author: leonard
"""
from tree.treewalk import density
from Parameters.constants import DESNNGB, BITS_FOR_POSITIONS, NORM_COEFF
from Parameters.parameter import NDIM

def find_sph_quantities(Particles, Problem, Functions, NgbTree, niter):
    if niter == 0:
        initial_guess_hsml(Particles, NgbTree, Problem)
    #now walk the tree to determine the weights and sph-quantities
    density(Particles, NgbTree, Problem)

def initial_guess_hsml(Particles, NgbTree, Problem):
    "computes an initial guess for the smoothing lengths"
    for i in range(NgbTree.MaxPart):
        no = NgbTree.Father[i]
        while(10 * DESNNGB * Problem.Mpart > NgbTree.get_nodep(no).Mass):
            p = NgbTree.get_nodep(no).father
            if p < 0:
                break
            no = p
        if(NgbTree.get_nodep(no).level > 0):
            length = (1 << (BITS_FOR_POSITIONS - NgbTree.get_nodep(no).level)) * \
                     Problem.FacIntToCoord.max()
        else:
            length = Problem.Boxsize.max()

        Particles[i].Hsml =  length * \
            (DESNNGB * Problem.Mpart / NgbTree.get_nodep(no).Mass \
             / NORM_COEFF)**(1/NDIM)