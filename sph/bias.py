#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:41:18 2021

@author: leonard
"""
from Parameters.constants import SMALL_NUM
from Parameters.parameter import Npart

def calculate_Bias(Particles):
    "As suggested by Klaus Dolag"
    #compute average density
    rho_mean = 0
    for particle in Particles:
        rho_mean += particle.Rho
    rho_mean /= Npart
    #compute bias
    bias     = 0
    for particle in Particles:
        bias += (particle.Rho_Model - particle.Rho)\
               / (  particle.Rho_Model - rho_mean \
                  + sign(particle.Rho_Model - rho_mean) * rho_mean * SMALL_NUM)
    bias /= Npart
    print("\nDensity bias: %g\n"%bias)
                

def sign(x):
    return x/abs(x)