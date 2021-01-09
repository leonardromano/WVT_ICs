#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:00:39 2021

@author: leonard
"""
from Parameters.constants import BITS_FOR_POSITIONS
from Parameters.parameter import NDIM
from data.structures import particle_data
from utility.integer_coordinates import get_distance_vector, \
    convert_to_phys_position
from numpy.linalg import norm
from math import ceil
from sys import exit
    
def sph_density_evaluate_particle_node_opening_criterion(particle, node, \
                                                         Problem, h):
    """
    This function checks whether there is a spatial overlap between the 
    (rectangular) enclosing box of the particles contained in a node, 
    and the search region.
    """
    if node.level <= 0:
        return 1
    part_offset = particle.position
    for i in range(NDIM):
        part_offset[i] -= ceil(h/Problem.FacIntToCoord[i])
    left  = get_distance_vector(node.center_offset_min + node.center, \
                                part_offset, Problem)
    right = get_distance_vector(node.center_offset_max + node.center, \
                                part_offset, Problem)
    for i in range(NDIM):
        if left[i] > 2 * ceil(h/Problem.FacIntToCoord[i]) and \
           right[i] > left[i]:
               return 0
    return 1

def sph_density_open_node(particle, nop, NgbTree, Problem, h):
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
        sph_density_interact(particle, p, typep, NgbTree, Problem, h)
        
        p = nextp

def sph_density_interact(particle, no, no_type, NgbTree, Problem, h):
    """
    Take care of SPH density interaction between the particle, and the node
    referenced through no. The node can either be a node or a particle.
    """
    if no_type == "particle":
        #we have a particle check whether it's a neighbor
        ngb = NgbTree.Tp[no]
        dx = get_distance_vector(particle.position, ngb.position, Problem)
        if norm(convert_to_phys_position(dx, Problem)) > h:
            return
        particle.neighbors.append([ngb, dx])
    else:
        node = NgbTree.get_nodep(no)
        if not node.notEmpty:
            return
        if sph_density_evaluate_particle_node_opening_criterion(particle, node, \
                                                                Problem, h):
            sph_density_open_node(particle, node, NgbTree, Problem, h)
                

def add_ghosts(particle, Problem, h):
    for i in range(NDIM):
        intHsml = ceil(h/Problem.FacIntToCoord[i])
        if not Problem.Periodic[i]:
            dist_from_wall = get_distance_from_wall(particle, i)
            if dist_from_wall < intHsml:
                #this particle is close to the wall
                min_dist_from_wall = \
                    get_minimum_distance_from_wall(particle.neighbors, i)
                listOfNeighbors = copy_list(particle.neighbors)
                for [neighbor, intDistance] in listOfNeighbors:
                    ngb_dist_from_wall = get_distance_from_wall(neighbor, i)
                    if ngb_dist_from_wall > dist_from_wall and \
                       dist_from_wall < abs(intDistance[i]) + \
                       0.75 * min_dist_from_wall:
                           particle.neighbors.\
                               append(mirror_particle(particle, neighbor, \
                                                      i, Problem))

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

def copy_list(oldList):
    "creates a hard copy of a list"
    newList = list()
    for item in oldList:
        newList.append(item)
    return newList

def find_neighbors(particle, Problem, NgbTree, h):
    particle.neighbors = list()
    sph_density_interact(particle, NgbTree.MaxPart, "node", NgbTree, Problem, h)
    # If the particle is close to the wall we need to add ghosts
    add_ghosts(particle, Problem, h)