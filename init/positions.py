#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:16:40 2021

@author: leonard
"""

def make_Positions(Particles, Functions):
    "Sample positions on random uniform grid"
    print("Sampling positions...")
    for particle in Particles:
        particle.position = Functions.Position_func(particle)
    print("done.")
    
def make_Velocities(Particles, Problem, Functions):
    "Sample velocities according to the analytical model"
    print("Sampling velocities...")
    for particle in Particles:
        particle.velocity = Functions.Velocity_func(particle, Problem)
    print("done.")
    
def make_Entropies(Particles, Problem, Functions):
    "Sample internal energy according to the analytical model"
    print("Sampling Entropy...")
    for particle in Particles:
        particle.Entropy = Functions.Entropy_func(particle, Problem)
    print("done.")

        
        
    