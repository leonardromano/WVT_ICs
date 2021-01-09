#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:15:01 2021

@author: leonard
"""
from numpy import zeros
from numpy.linalg import norm

from Parameters.parameter import NDIM
from Parameters.constants import DESNNGB, NORM_COEFF, BITS_FOR_POSITIONS, MAX_INT
from utility.integer_coordinates import convert_to_int_position

def drift_particles(Particles, Problem):
    "In this function the particles are drifted according to their WVT forces"
    cnts = zeros(4, dtype = int)
    for particle in Particles:
        d = norm(particle.delta)
        d_mps = (NORM_COEFF/DESNNGB)**(1/NDIM) * particle.Hsml
        for i in range(cnts.shape[0]):
            if d > (0.1)**(i) * d_mps:
                cnts[i] += 1
        particle.position += convert_to_int_position(particle.delta, Problem)
        keep_inside_box(particle, Problem)
    return cnts
        

def keep_inside_box(particle, Problem):
    "Makes sure the particle is within the domain"
    for axis in range(NDIM):
        if Problem.Periodic[axis]:
            while particle.position[axis] < 0:
                particle.position[axis] += (1 << BITS_FOR_POSITIONS)
            while particle.position[axis] > MAX_INT:
                particle.position[axis] -= (1 << BITS_FOR_POSITIONS)
        else:
            while particle.position[axis] < 0 or \
                particle.position[axis] > MAX_INT:
                    if particle.position[axis] < 0:
                        if particle.position[axis] < -(1 << (BITS_FOR_POSITIONS-1)):
                            particle.position += (1 << BITS_FOR_POSITIONS)
                        else:
                            particle.position[axis] *= -1
                    else:
                        if particle.position[axis] > MAX_INT + (1 << (BITS_FOR_POSITIONS-1)):
                            particle.position[axis] -= (1 << BITS_FOR_POSITIONS)
                        else:
                            particle.position[axis] += 2 * (MAX_INT - \
                                                        particle.position[axis])
            
        