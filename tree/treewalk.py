#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 13:29:49 2021

@author: leonard
"""
from Parameters.constants import DESNNGB, NNGBDEV, BITS_FOR_POSITIONS, NORM_COEFF
from Parameters.parameter import NDIM
from sph.Kernel import kernel, wendland_bias_correction
from data.structures import particle_data
from utility.integer_coordinates import get_distance_vector, \
    convert_to_phys_position
from numpy.linalg import norm
from numpy import zeros
from math import ceil
from sys import exit

###############################################################################
#Treewalk related functions

def evaluate_particle_node_opening_criterion(particle, node, Problem):
    """
    This function checks whether there is a spatial overlap between the 
    (rectangular) enclosing box of the particles contained in a node, 
    and the search region.
    """
    if node.level <= 0:
        return 1
    part_offset  = zeros(NDIM, dtype = int)
    part_offset += particle.position
    for i in range(NDIM):
        part_offset[i] -= ceil(particle.Hsml/Problem.FacIntToCoord[i])
        
    left  = get_distance_vector(node.center_offset_min + node.center, \
                                part_offset, Problem)
    right = get_distance_vector(node.center_offset_max + node.center, \
                                part_offset, Problem)
    for i in range(NDIM):
        if left[i] > 2 * ceil(particle.Hsml/Problem.FacIntToCoord[i]) and \
           right[i] > left[i]:
               return 0
    return 1

def sph_density_open_node(particle, nop, NgbTree, Problem):
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
        sph_density_interact(particle, p, typep, NgbTree, Problem)
        
        p = nextp

def sph_density_interact(particle, no, no_type, NgbTree, Problem):
    """
    Take care of SPH density interaction between the particle, and the node
    referenced through no. The node can either be a node or a particle.
    """
    if no_type == "particle":
        #we have a particle check whether it's a neighbor
        ngb = NgbTree.Tp[no]
        dx = get_distance_vector(particle.position, ngb.position, Problem)
        if norm(convert_to_phys_position(dx, Problem)) > particle.Hsml:
            return
        particle.neighbors.append([ngb, dx])
    else:
        node = NgbTree.get_nodep(no)
        if not node.notEmpty:
            return
        if evaluate_particle_node_opening_criterion(particle, node, Problem):
            sph_density_open_node(particle, node, NgbTree, Problem)

def densities_determine(targetlist, NgbTree, Particles, Problem):
    """
    for each target walk the tree to determine the neighbors and then 
    compute density and thermodynamic quantities
    """
    for i in targetlist:
        particle                = Particles[i]
        particle.Rho            = 0
        particle.VarHsmlFac     = 0
        particle.neighbors      = list()
        sph_density_interact(particle, NgbTree.MaxPart, "node", NgbTree, Problem)
        # If the particle is close to the wall we need to add ghosts
        add_ghosts(particle, Problem)
        evaluate_kernel(particle, Problem)
###############################################################################
#main loop

def density(Particles, NgbTree, Problem):
    "For each particle compute density, smoothing length and thermodynamic variables"
    Left = zeros(len(Particles))
    Right = zeros(len(Particles))
    targetlist = list()
    for particle in Particles:
        targetlist.append(particle.ID)
    while True:
        # now do the primary work with this call
        densities_determine(targetlist, NgbTree, Particles, Problem)
        # do final operations on results
        npleft = 0
        for p in range(len(targetlist)):
            i = targetlist[p]
            particle = Particles[i]
            #do some postprocessing on density
            finish_density_update(particle, Problem)
            numberOfNeighbors = NORM_COEFF * particle.Hsml**(NDIM) \
                                * particle.Rho/Problem.Mpart
            if (abs(numberOfNeighbors-DESNNGB)>NNGBDEV):    
                #check whether we're done
                if Left[i] > 0 and Right[i] > 0 and \
                   Right[i]-Left[i] < 1e-3 * Left[i]:
                    #this one should be ok
                    continue
                #need to redo this one
                targetlist[npleft] = i
                npleft += 1
                Left[i], Right[i] = update_bounds(Left[i], Right[i], \
                                                  numberOfNeighbors, \
                                                  particle.Hsml)
                update_smoothing_length(Left[i], Right[i], \
                                        numberOfNeighbors, particle)
        targetlist = targetlist[:npleft]
        if len(targetlist) <=0:
            break

###############################################################################
#Bisection algorithm related functions

def update_bounds(lowerBound, upperBound, numberOfNeighbors, h):
    "Update the bounds for the smoothing length in the bisection algorithm"
    if numberOfNeighbors < DESNNGB - NNGBDEV:
        lowerBound = max(lowerBound, h)
    else:
        if upperBound != 0:
            if h < upperBound:
                upperBound = h
        else:
            upperBound = h
    return lowerBound, upperBound

def update_smoothing_length(lowerBound, upperBound, numberOfNeighbors, particle):
    "Perform the bisection part of the bisection algorithm"
    if lowerBound > 0 and upperBound > 0:
        particle.Hsml = ((lowerBound**3 + upperBound**3)/2)**(1/3)
    else:
        if upperBound == 0 and lowerBound == 0:
            print("Upper and Lower bounds not updated!")
            exit()
            
        if upperBound == 0 and lowerBound>0:
            if abs(numberOfNeighbors - DESNNGB) < 0.5 * DESNNGB:
                fac = 1-(numberOfNeighbors - DESNNGB)/numberOfNeighbors/NDIM*\
                    particle.VarHsmlFac
                if fac < 1.26:
                    particle.Hsml *= fac
                else:
                    particle.Hsml *= 1.26
            else:
                particle.Hsml *= 1.26
                    
        if upperBound > 0 and lowerBound  == 0:
            if abs(numberOfNeighbors - DESNNGB) < 0.5*DESNNGB:
                fac = 1-(numberOfNeighbors - DESNNGB)/numberOfNeighbors/NDIM*\
                    particle.VarHsmlFac
                if fac > 1/1.26:
                    particle.Hsml *= fac
                else:
                    particle.Hsml /= 1.26
            else:
                particle.Hsml /= 1.26

###############################################################################
#ghost particle related functions           

def add_ghosts(particle, Problem):
    for i in range(NDIM):
        if not Problem.Periodic[i]:
            intHsml = ceil(particle.Hsml/Problem.FacIntToCoord[i])
            dist_from_wall = get_distance_from_wall(particle, i)
            #check if we are touching the wall
            if dist_from_wall < intHsml:
                min_dist_from_wall = \
                    get_minimum_distance_from_wall(particle.neighbors, i)
                listOfNeighbors = copy_list(particle.neighbors)
                for [neighbor, intDistance] in listOfNeighbors:
                    ngb_dist_from_wall = get_distance_from_wall(neighbor, i)
                    #check if it makes sense to mirror neighbor
                    if ngb_dist_from_wall > dist_from_wall and \
                       dist_from_wall < abs(intDistance[i]) + \
                       0.75 * min_dist_from_wall:
                           particle.neighbors.\
                               append(mirror_particle(particle, neighbor, \
                                                      i, Problem))
                                   
def mirror_particle(particle, neighbor, axis, Problem):
    """
    creates an instance of a ghost particle mirroring a neighbor 
    of a particle close to the wall
    """
    mirrorParticle = particle_data(neighbor.ID)
    mirrorParticle.position +=neighbor.position
    mirrorParticle.position += 2*(particle.position[axis] \
                                  - neighbor.position[axis])
    mirrorParticle.Hsml     = neighbor.Hsml
    mirrorParticle.IsGhost  = 1
    dx = get_distance_vector(particle.position, mirrorParticle.position, \
                                                Problem)
    return [mirrorParticle, dx]

def get_distance_from_wall(particle, axis):
    "returns the distance from the nearest non-periodic boundary"
    return int(min(particle.position[axis], \
                   (1 << BITS_FOR_POSITIONS) - particle.position[axis]))

def get_minimum_distance_from_wall(particles, axis):
    "Determines the minimum distance from the wall among a list of particles"
    minimumDistanceFromWall = 1 << BITS_FOR_POSITIONS
    for particle in particles:
        dist = get_distance_from_wall(particle[0], axis)
        if dist < minimumDistanceFromWall:
            minimumDistanceFromWall = dist
    return minimumDistanceFromWall

def copy_list(oldList):
    "creates a hard copy of a list"
    newList = list()
    for item in oldList:
        newList.append(item)
    return newList

###############################################################################
#density calculation

def evaluate_kernel(particle, Problem):
    "Perform the neighbor sum to compute density and SPH correction factor"
    for [neighbor, intDistance] in particle.neighbors:
        dist = norm(convert_to_phys_position(intDistance, Problem))
        particle.Rho +=  kernel(dist/particle.Hsml, particle.Hsml)
        particle.VarHsmlFac -= dist * \
            kernel(dist/particle.Hsml, particle.Hsml, True)
            
def finish_density_update(particle, Problem):
    "some final postprocessing steps in density calculation"
    if particle.Rho > 0:
        wendland_bias_correction(particle)
        #now if we have more than one neighbor update the dhsml factor
        if particle.VarHsmlFac > 0:
            particle.VarHsmlFac = NDIM * particle.Rho/particle.VarHsmlFac
        else:
            particle.VarHsmlFac = 1
        particle.Rho *= Problem.Mpart

