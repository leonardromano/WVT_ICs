#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:16:40 2021

@author: leonard
"""
from time import time

def sample(Particles, Problem, funcs):
    "Sample the quantities specified by funcs"
    t0 = time()
    
    for particle in Particles:
        for func in funcs:
            func(particle)
            
    t1 = time()
    Problem.Timer["INIT"] += t1-t0

        
        
    