#!/usr/bin/env python
# coding: utf-8

import argparse
import json
import sys
from pathlib import Path
from time import time

from ngsolve import TaskManager, pi, H1,GridFunction, HDiv, grad, curl, Integrate, dx, sqrt, Norm, Conj, InnerProduct
from ngsolve.webgui import Draw
from src import (
    create_spherical_geometry,
    create_ellipsoid_scatterer_geometry,
    create_box_scatterer_geometry,
    create_two_box_scatterer_geometry,
    create_two_ellipsoid_scatterer_geometry,
    create_cylinder_scatterer_geometry,
    create_incident_wave,
    MaxwellProblem,
    solve_gmres,
    solve_cg,
    solve_direct,
    solve_bvp,
)

# PARAMETERS
wavelength = 1.0
prop_dir = [0,0,1]
polarization = [1,0,0]
scatterer_radius = 0.1
outer_radius = 1.0
pml_width = 0.25
mesh_size = 0.25
order = 5
solver = "gmres"  # options: gmres, bvp, cg, direct
num_threads = 4

# Compute wavenumber
k = 2 * pi / wavelength

print(f"Physical Parameters:")
print(f"  Wavelength: {wavelength}")
print(f"  Wavenumber: {k:.4f}")
print(f"  Propagation: {prop_dir}")
print(f"  Polarization: {polarization}")

# GEOMETRY
print(f"\nCreating geometry...", flush=True)
mesh = None

if False:
    mesh = create_spherical_geometry(
        R=outer_radius,
        PMLw=pml_width,
        r=scatterer_radius,
        h_max=mesh_size
    )
if True:
    mesh = create_ellipsoid_scatterer_geometry(
        wavelength=wavelength,
        semi_axis_a=0.125,
        semi_axis_b=0.25,
        semi_axis_c=0.25,
        domain_radius=outer_radius,
        pml_width=pml_width,
        max_mesh_size=mesh_size,
        scatterer_mesh_size=mesh_size,
        curve_order=order
    )
if False:
    mesh = create_box_scatterer_geometry(
        wavelength=wavelength,
        axis_a=0.7,
        axis_b=0.54,
        axis_c=0.65,
        domain_radius=outer_radius,
        pml_width=pml_width,
        max_mesh_size=mesh_size,
        scatterer_mesh_size=mesh_size,
        curve_order=order,
        box_radius=0.1
    )
if False:
    mesh = create_two_box_scatterer_geometry(
        wavelength=wavelength,
        dist=0.5,
        b1_axis_a=0.1,
        b1_axis_b=0.4,
        b1_axis_c=0.4,
        b2_axis_a=0.1,
        b2_axis_b=0.4,
        b2_axis_c=0.4,
        box_radius=0.049,
        domain_radius=outer_radius,
        pml_width=pml_width,
        max_mesh_size=mesh_size,
        curve_order=order
    )
if False:
    mesh = create_two_ellipsoid_scatterer_geometry(
        wavelength=wavelength,
        dist=0.4,
        b1_axis_a=0.1,
        b1_axis_b=0.25,
        b1_axis_c=0.1,
        b2_axis_a=0.1,
        b2_axis_b=0.25,
        b2_axis_c=0.1,
        box_radius=0.04,
        domain_radius=outer_radius,
        pml_width=pml_width,
        max_mesh_size=mesh_size,
        curve_order=order
    )
if False:
    mesh = create_cylinder_scatterer_geometry(
        wavelength=wavelength,
        origin=(0,0,0),
        height=0.5,
        radius=0.2,
        radius_2=0.15,
        box_radius=0.05,
        domain_radius=outer_radius,
        pml_width=pml_width,
        max_mesh_size=mesh_size,
        curve_order=order
    ) 

print(f"Mesh generated:")
print(f"  Elements: {mesh.ne}")
print(f"  Vertices: {mesh.nv}")

clipping = { "function" : False,  "pnt" : (0,0,0.25), "vec" : (0,0,-1) }
Draw(mesh, clipping=clipping, filename="bin/imgs/scatterer_mesh_test.html")


# """
# l = []    # l = list of estimated total error
# space_flux = HDiv(mesh,  complex=True, autoupdate=True)
# gf_flux = GridFunction(space_flux, "flux", autoupdate=True)

# def CalcError(gfu):
#     fes = problem.fes

#     flux = curl(gfu)   # the FEM-flux      
#     gf_flux.Set(flux)        # interpolate into H(div)

#     # compute estimator:
#     err = InnerProduct(flux-gf_flux, flux-gf_flux)
#     eta2 = Integrate(err*dx, mesh, element_wise=True)
#     l.append ((fes.ndof, sqrt(sum(eta2))))
#     print("ndof =", fes.ndof, " toterr =", sqrt(sum(eta2)))

#     # mark for refinement based on error:
#     maxerr = max(eta2.real) 
#     for el in mesh.Elements():
#        if eta2[el.nr].real > 0.25*maxerr:
#            mesh.SetRefinementFlag(el, True)

#     # additional refinement at material interfaces (box edges):
#     for el in mesh.Elements():
#         for f in el.faces:
#             if len(f.elements) == 2:  # internal face
#                 mat1 = f.elements[0].mat
#                 mat2 = f.elements[1].mat
#                 if mat1 != mat2:
#                     mesh.SetRefinementFlag(el, True)
#                     break

#     clipping = { "function" : False,  "pnt" : (0,0,0.25), "vec" : (0,0,-1) }
#     Draw(err, mesh, clipping=clipping);
# """

# INCIDENT WAVE
E_inc = create_incident_wave(
    k=k,
    propagation_dir=tuple(prop_dir),
    polarization=tuple(polarization)
)

Draw(E_inc,mesh, filename="bin/imgs/scatterer_Einc_test.html")

# SETUP PROBLEM
print(f"\nSetting up scattering problem...", flush=True)
problem = MaxwellProblem(mesh, k, E_inc, fes_order=order)

print(f"\nAssembling system...", flush=True)
problem.assemble_system()

# SOLVE
solution = None

with TaskManager():#pajetrace=10**8):
    # 'block_jacobi', 'bddc', 'hcurlamg'
    if solver == "gmres":
        solution = solve_gmres(problem.a, problem.l, problem.fes, preconditioner="block_jacobi", maxsteps=400)#, restart=100)
    elif solver == "bvp":
        solution = solve_bvp(problem.a, problem.l, problem.fes, preconditioner="block_jacobi", maxsteps=1000)
    elif solver == "cg":
        solution = solve_cg(problem.a, problem.l, problem.fes, preconditioner="block_jacobi", maxsteps=1000)
    else:
        solution = solve_direct(problem.a, problem.l, problem.fes)

problem.set_solution(solution)

clipping = { "function" : False,  "pnt" : (0,0,0.26), "vec" : (0,0,-1) }
vectors = {"grid_size" : 50, "offset" : 0.5 }

line_1 = { "type": "lines", "position": [1,1,0, 1+0.5*prop_dir[0], 1+0.5*prop_dir[1], 0.5*prop_dir[2]], "name": "propagation direction", "color": "red",}
line_2 = { "type": "lines", "position": [1,1,0, 1+0.5*polarization[0], 1+0.5*polarization[1], 0.5*polarization[2]], "name": "polarization direction", "color": "blue"}
points = { "type": "points", "position": [1,1,0], "size":20, "color": "black", "name": "origin"}
text_1 = { "type": "text", "name": "info1", "text": f" wavelength = {wavelength}, outer radius = {outer_radius}, PML width = {pml_width}, mesh size = {mesh_size}", "position": [-1,-1.2,0]}
text_2 = { "type": "text", "name": "info2", "text": f" elements = {mesh.ne}, vertices = {mesh.nv}, free DOFs = {sum(problem.fes.FreeDofs())}", "position": [-1,-1.3,0]}

Draw(solution, mesh, objects=[line_1,line_2,points,text_1,text_2], clipping=clipping, filename="bin/imgs/scatterer_solution_test.html");

# CalcError(solution)
# mesh.Refine()
