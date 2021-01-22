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
from wvt.relax import regularise_particles

def main():
    "Set up a box with SPH particles, then WVT relax the particles"
    t0 = time()
    Particles, Problem, Functions = setup()
    make_Positions(Particles, Problem, Functions)
    regularise_particles(Particles, Problem, Functions)
    make_Velocities(Particles, Problem, Functions)
    make_Entropies(Particles, Problem, Functions)
    write_output(Particles, Problem, Functions)
    t1 = time()
    T = t1 - t0
    print("Successfully created ICs! Took %g seconds.\n"%T)
    print("Compuational cost of individual parts:")
    print("INIT: %g (%g)"%(Problem.Timer["INIT"], Problem.Timer["INIT"]/T))
    print("TREE: %g (%g)"%(Problem.Timer["TREE"], Problem.Timer["TREE"]/T))
    print("DENSITY: %g (%g)"%(Problem.Timer["DENSITY"], Problem.Timer["DENSITY"]/T))
    print("L1-ERROR: %g (%g)"%(Problem.Timer["L1-ERROR"], Problem.Timer["L1-ERROR"]/T))
    print("WVT: %g (%g)"%(Problem.Timer["WVT"], Problem.Timer["WVT"]/T))
    print("REDIST: %g (%g)"%(Problem.Timer["REDIST"], Problem.Timer["REDIST"]/T))
    print("OUTPUT: %g (%g)"%(Problem.Timer["OUTPUT"], Problem.Timer["OUTPUT"]/T))
main()