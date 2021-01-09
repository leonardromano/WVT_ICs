#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 14:48:49 2021

@author: leonard
"""
from time import time

from output.output import write_output

from init.positions import make_Positions, make_Velocities, make_Entropies
from init.setup import setup
from sph.bias import calculate_Bias
from wvt.relax import regularise_particles

def main():
    "Set up a box with SPH particles, then WVT relax the particles"
    t0 = time()
    Particles, Problem, Functions = setup()
    make_Positions(Particles, Functions)
    regularise_particles(Particles, Problem, Functions)
    make_Velocities(Particles, Problem, Functions)
    make_Entropies(Particles, Problem, Functions)
    #bias as measure of convergence
    calculate_Bias(Particles)
    write_output(Particles, Problem, Functions)
    t1 = time()
    print("Successfully created ICs! Took %g seconds"%(t1 - t0))
    
main()