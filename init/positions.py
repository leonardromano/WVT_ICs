#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:16:40 2021

@author: leonard
"""
from time import time

def make_Positions(Particles, Problem, Functions):
    "Sample positions on random uniform grid"
    print("INIT: Sampling positions...")
    t0 = time()
    for particle in Particles:
        particle.position = Functions.Position_func(particle)
    t1 = time()
    Problem.Timer["INIT"] += t1-t0
    
def make_Velocities(Particles, Problem, Functions):
    "Sample velocities according to the analytical model"
    print("INIT: Sampling velocities...")
    t0 = time()
    for particle in Particles:
        particle.velocity = Functions.Velocity_func(particle, Problem)
    t1 = time()
    Problem.Timer["INIT"] += t1-t0
    
def make_Entropies(Particles, Problem, Functions):
    "Sample internal energy according to the analytical model"
    print("INIT: Sampling Entropy...")
    t0 = time()
    for particle in Particles:
        particle.Entropy = Functions.Entropy_func(particle, Problem)
    t1 = time()
    Problem.Timer["INIT"] += t1-t0

        
        
    