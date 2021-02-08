#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 13:29:49 2021

@author: leonard
"""
from numpy import zeros, asarray
from math import ceil
from sys import exit

from Parameters.constants import DESNNGB, NNGBDEV, NORM_COEFF, MAX_INT
from Parameters.parameter import NDIM
from sph.Kernel import kernel, wendland_bias_correction
from utility.integer_coordinates import get_distance_vector
from utility.utility import norm


###############################################################################
#Treewalk related functions

def evaluate_particle_node_opening_criterion(particle, node, NgbTree):
    """
    This function checks whether there is a spatial overlap between the 
    (rectangular) enclosing box of the particles contained in a node, 
    and the search region.
    """
    if node.level <= 0:
        return 1
    
    #compute the offset from the particle position
    part_offset  = zeros(NDIM, dtype = int)
    part_offset += particle.position
    for i in range(NDIM):
        part_offset[i] -= ceil(particle.Hsml/NgbTree.FacIntToCoord[i])
        
    left  = get_distance_vector(node.center_offset_min + node.center, \
                                part_offset, NgbTree)
    right = get_distance_vector(node.center_offset_max + node.center, \
                                part_offset, NgbTree)
    for i in range(NDIM):
        if left[i] > 2 * particle.Hsml and right[i] > left[i]:
            return 0
    return 1

def sph_density_open_node(particle, nop, NgbTree):
    "Continues to walk the tree for the particle by opening a node."
    p = nop.nextnode
    while p != nop.sibling:
        if p < 0:
            print("p=%d < 0  node.sibling=%d node.nextnode=%d" \
                   %(p, nop.sibling, nop.nextnode))
            exit()
        nextp = 0
        typep = ""
        if p < NgbTree.MaxPart:
            nextp = NgbTree.Nextnode[p]
            typep = "particle"
        else:
            node = NgbTree.get_nodep(p)
            nextp = node.sibling
            typep = "node"
        sph_density_interact(particle, p, typep, NgbTree)
        
        p = nextp

def sph_density_interact(particle, no, no_type, NgbTree):
    """
    Take care of SPH density interaction between the particle, and the node
    referenced through no. The node can either be a node or a particle.
    """
    if no_type == "particle":
        #we have a particle check whether it's a neighbor
        ngb = NgbTree.Tp[no]
        r = norm(get_distance_vector(particle.position, ngb.position, NgbTree))
        if r > particle.Hsml:
            return
        particle.neighbors.append(no)
    else:
        node = NgbTree.get_nodep(no)
        if not node.notEmpty:
            return
        if evaluate_particle_node_opening_criterion(particle, node, NgbTree):
            sph_density_open_node(particle, node, NgbTree)

def densities_determine(NgbTree, Workstack):
    """
    for each target walk the tree to determine the neighbors and then 
    compute density and thermodynamic quantities
    """
    for particle in Workstack:
        particle.Rho            = 0
        particle.neighbors      = list()
        sph_density_interact(particle, NgbTree.MaxPart, "node", NgbTree)
        evaluate_kernel(particle, NgbTree)
###############################################################################
#main loop

def density(Particles, NgbTree):
    "For each particle compute density, smoothing length and thermodynamic variables"
    Left = zeros(len(Particles))
    Right = zeros(len(Particles))
    Workstack = copy_list(Particles)
    Donestack = list()
    niter = 0
    while True:
        # now do the primary work with this call
        densities_determine(NgbTree, Workstack)
        # do final operations on results
        nstack = len(Workstack)
        npleft = 0
        for _ in range(nstack):
            particle = Workstack.pop(0)
            #do some postprocessing on density
            finish_density_update(particle, NgbTree)
            NumNgb = NORM_COEFF * particle.Hsml**(NDIM) \
                                * particle.Rho/NgbTree.Mpart
            
            if (abs(NumNgb-DESNNGB) > NNGBDEV):    
                #check whether we're done                
                i = particle.ID
                if Left[i] > 0 and Right[i] > 0 and \
                   Right[i]-Left[i] < 1e-3 * Left[i]:
                    #this one should be ok
                    Donestack.append(particle)
                    continue
                #need to redo this one
                npleft += 1
                Left[i], Right[i] = update_bounds(Left[i], Right[i], \
                                                  NumNgb, particle.Hsml)
                update_smoothing_length(Left[i], Right[i], particle)
                Workstack.append(particle)
            else:
                Donestack.append(particle)
            
        niter += 1
        if npleft <= 0:
            break
    
    print("DENSITY: Finished after %d iterations."%niter)
    Particles = Donestack

###############################################################################
#Bisection algorithm related functions

def update_bounds(lowerBound, upperBound, NumNgb, h):
    "Update the bounds for the smoothing length in the bisection algorithm"
    if NumNgb < DESNNGB - NNGBDEV:
        lowerBound = max(lowerBound, h)
    else:
        if upperBound != 0:
            upperBound = min(h, upperBound)
        else:
            upperBound = h
    return lowerBound, upperBound

def update_smoothing_length(lowerBound, upperBound, particle):
    "Perform the bisection part of the bisection algorithm"
    if lowerBound > 0 and upperBound > 0:
        particle.Hsml = ((lowerBound**3 + upperBound**3)/2)**(1/3)
    else:
        if upperBound == 0 and lowerBound == 0:
            print("Upper and Lower bounds not updated!")
            exit()
            
        if upperBound == 0 and lowerBound > 0:
            particle.Hsml *= 1.26
                    
        if upperBound > 0 and lowerBound  == 0:
            particle.Hsml /= 1.26

###############################################################################
#ghost particle related functions    

def close_to_wall(particle, NgbTree):
    for i in range(NDIM):
        if not NgbTree.Periodic[i]:
            dist_from_wall = get_distance_from_wall(particle, NgbTree, i)
            if particle.Hsml > dist_from_wall:
                particle.CloseToWall = 1
                return
    particle.CloseToWall = 0

def add_ghost(particle, ngb, dist, min_dist_from_wall, NgbTree):
    for i in range(NDIM):
        if not NgbTree.Periodic[i]:
            dist_from_wall = get_distance_from_wall(particle, NgbTree, i)
            if dist_from_wall < get_distance_from_wall(ngb, NgbTree, i) and \
               dist_from_wall < abs(dist[i]) + 0.75 * min_dist_from_wall[i]:
                   return 1
    return 0

def get_distance_from_wall(particle, NgbTree, axis):
    "returns the distance from the nearest non-periodic boundary"
    dx = min(particle.position[axis], MAX_INT - particle.position[axis]) 
    return NgbTree.FacIntToCoord[axis] * dx

def get_minimum_distance_from_wall(particle, NgbTree):
    "Determines the minimum distance from the wall among a list of particles"
    min_dist_from_wall = [*NgbTree.Boxsize]
    for i in range(NDIM):
        if not NgbTree.Periodic[i]:
            for n in particle.neighbors:
                dist = get_distance_from_wall(NgbTree.Tp[n], NgbTree, i)
                if dist < min_dist_from_wall[i]:
                    min_dist_from_wall[i] = dist
    return asarray(min_dist_from_wall)

def copy_list(oldList):
    "creates a hard copy of a list"
    return [item for item in oldList]

###############################################################################
#density calculation

def evaluate_kernel(particle, NgbTree):
    "Perform the neighbor sum to compute density and SPH correction factor"
    close_to_wall(particle, NgbTree)
    if particle.CloseToWall:
        min_dist_from_wall = get_minimum_distance_from_wall(particle, NgbTree)
    for n in particle.neighbors:
        ngb = NgbTree.Tp[n]
        dist = get_distance_vector(particle.position, ngb.position, NgbTree)
        h = particle.Hsml
        r = norm(dist)
        wk = kernel(r/h, h)
        if particle.CloseToWall and add_ghost(particle, ngb, dist, min_dist_from_wall, NgbTree):
            #if we have a ghost, double the neighbors contribution instead
            wk *= 2.0
        particle.Rho += wk
    
def finish_density_update(particle, NgbTree):
    "some final postprocessing steps in density calculation"
    if particle.Rho > 0:
        wendland_bias_correction(particle)
        particle.Rho *= NgbTree.Mpart
