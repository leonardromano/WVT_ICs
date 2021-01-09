#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:47:28 2021

@author: leonard
"""
from Parameters.parameter import BiasCorrection

def volume(vector):
    dV = 1
    for dx in vector:
        dV *= dx
    return dV

def factorial(n):
    result = 1
    for i in range(2, n+1):
        result *= i
    return result

def relative_density_error_with_sign(particle, Problem, Functions):
    rhoModel = Functions.Density_func(particle, Problem, BiasCorrection)
    return (particle.Rho - rhoModel)/rhoModel

def relative_density_error(particle, Problem, Functions):
    return abs(relative_density_error_with_sign(particle, Problem, Functions))