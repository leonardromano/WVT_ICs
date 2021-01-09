#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:36:40 2021

@author: leonard

    Settle SPH particle with weighted Voronoi tesselations (Diehl+ 2012).
    Here hsml is not the SPH smoothing length, but is related to a local
    metric defined ultimately by the density model.
    Relaxation is done in units of the boxsize, hence the box volume is 1
    Return code true means that rerun could be usefull
"""
#public libraries
from numpy import array
from time import time
from sys import exit
#custom libraries
from output.diagnostics import write_diagnostics, write_step_file
from Parameters.constants import LARGE_NUM
from Parameters.parameter import Npart, NDIM, Maxiter, MpsFraction, \
    StepReduction, LimitMps, LimitMps10, LimitMps100, LimitMps1000, \
    LastMoveStep, RedistributionFrequency, SAVE_WVT_STEPS
from sph.sph import find_sph_quantities
from tree.tree import ngbtree
from utility.errors import compute_l1_error
from wvt.drift import drift_particles
from wvt.forces import compute_wvt_weights, compute_wvt_forces
from wvt.redistribution import redistribute

def regularise_particles(Particles, Problem, Functions):
    "This is the main loop of the IC making"
    print("Starting iterative SPH regularisation \n" + \
          "   Maxiter=%d, MpsFraction=%g StepReduction=%g "\
              %(Maxiter, MpsFraction, StepReduction) + \
          "LimitMps=(%g,%g,%g,%g)\n\n"\
              %(LimitMps, LimitMps10, LimitMps100, LimitMps1000))
    t0 = time()
    #compute some necessary variables
    avg_boxsize = 0
    for length in Problem.Boxsize:
        avg_boxsize += length
    avg_boxsize /= NDIM
    
    npart_1D = Npart**(1/NDIM)
    step     = 1./(npart_1D * MpsFraction)
    last_cnt = LARGE_NUM
    err_last = LARGE_NUM
    err_diff = LARGE_NUM
    niter    = 0
    
    while(niter <= Maxiter):
        #build the search tree and update SPH quantities
        
        
        
        NgbTree = ngbtree(Particles, Problem)
        find_sph_quantities(Particles, Problem, Functions, NgbTree, niter)
        niter += 1
        if SAVE_WVT_STEPS:
            write_step_file(Particles, Problem, niter)
        
        #redistribute particles
        NgbTree = redistribute(Particles, Problem, Functions, NgbTree, niter)
            
        #next find minimum, maximum and average error and their variance
        err_min, err_max, \
        err_mean, err_sigma = compute_l1_error(Particles, Problem, Functions)
        
        #update err_diff and err_last
        err_diff  = ( err_last - err_mean ) / err_mean
        print("   #%02d: Err min=%3g max=%3g mean=%03g sigma=%03g diff=%03g step=%g\n"\
              %(niter, err_min, err_max, err_mean, err_sigma, err_diff, step))
        err_last  = err_mean
        
        #Now compute the sph volume and hsml weights
        compute_wvt_weights(Particles, Problem, Functions, avg_boxsize)
        
        #now compute the forces
        compute_wvt_forces(Particles, Problem, Functions, \
                                   NgbTree, step)
        #now drift all particles
        cnts = drift_particles(Particles, Problem)
        #now some diagnostics
        move_Mps = cnts * 100/Npart
        print("        Del %g > Dmps; %g > Dmps/10; %g > Dmps/100; %g > Dmps/1000\n"\
              %(move_Mps[0], move_Mps[1], move_Mps[2], move_Mps[3]))
        if niter == 1:
            if move_Mps[0] > 80:
                print("WARNING: A lot of initial movement detected." )
                print("Consider increasing MpsFraction in the parameter file!\n")
                exit()
        write_diagnostics(niter, err_diff, move_Mps, \
                          array([err_min, err_max, err_mean, err_sigma]))
        
        #check whether we're done
        if move_Mps[0] < LimitMps or move_Mps[1] < LimitMps10 \
            or move_Mps[2] < LimitMps100 or move_Mps[3] < LimitMps1000:
                break
        
        #enforce convergence if distribution doesnt tighten
        if cnts[1] > last_cnt and ( niter > LastMoveStep or niter % RedistributionFrequency != 0 ):
            step *= StepReduction
        
        last_cnt = cnts[1]
    
    if niter > Maxiter:
        print("Max iterations reached, result might not be converged properly.")
    t1 = time()
    print("Finished WVT relaxation after %03d iterations. Took %g seconds.\n"\
          %(niter, t1 - t0))
                      
            
        
        
        
    
