#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 12:43:29 2021

@author: leonard
"""
import h5py
from numpy import zeros

from Parameters.parameter import output, Npart, NDIM
from utility.errors import compute_field_stats
from utility.integer_coordinates import convert_to_phys_position

def write_diagnostics(niter, delta, err_diff, move_Mps, err_quad):
    "iteration err_min err_max err_mean err_sigma err_diff move_Mps[4] s_delta[4]"
    s_delta = compute_field_stats(delta)
    #write to file
    if niter == 1:
        file = open(output + "diagnostics.log", "w")
    else:
        file = open(output + "diagnostics.log", "a")
    file.write("%03d %g %g %g %g %g %g %g %g %g %g %g %g %g\n"\
               %(niter, *err_quad, err_diff, *move_Mps, *s_delta))
    file.close()

def write_step_file(Particles, Problem, niter):
    "writes all the particle data in a hdf5 file."
    f = h5py.File("%s/stepfiles/step_%03d.hdf5"%(output, niter), "w")
    #first dump all the header info
    header = f.create_group("Header")
    h_att = header.attrs
    h_att.create("NumPart", Npart)
    h_att.create("Mpart", Problem.Mpart)
    h_att.create("Boxsize", Problem.Boxsize)
    h_att.create("NDIM", NDIM)
    #now make the data sets for the particle data
    IDs        = zeros((Npart), dtype = int)
    positions  = zeros((Npart, NDIM), dtype = float)
    velocities = zeros((Npart, NDIM), dtype = float)
    densities  = zeros((Npart), dtype = float)
    hsml       = zeros((Npart), dtype = float)
    for i in range(Npart):
        IDs[i]        += Particles[i].ID
        positions[i]  += convert_to_phys_position(Particles[i].position, Problem)
        velocities[i] += Particles[i].velocity
        densities[i]  += Particles[i].Rho
        hsml[i]       += Particles[i].Hsml
    f.create_dataset("PartData/IDs", data = IDs,  dtype = "u4")
    f.create_dataset("PartData/Coordinates", data = positions, dtype = "f4")
    f.create_dataset("PartData/Velocity", data = velocities, dtype = "f4")
    f.create_dataset("PartData/Density", data = densities, dtype = "f4")
    f.create_dataset("PartData/SmoothingLength", data = hsml, dtype = "f4")
    f.close()