#!/usr/bin/env python3
"""
Run electromagnetic scattering simulation.

This script solves the Maxwell-Helmholtz scattering problem for a spherical
scatterer illuminated by a plane wave.
"""

import argparse
import json
import sys
from pathlib import Path
from time import time

from ngsolve import TaskManager, pi
from src import (
    create_spherical_geometry,
    create_incident_wave,
    ScattererProblem,
    solve_gmres,
    solve_direct,
)


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Electromagnetic scattering simulation using NGSolve"
    )

    # Physical parameters
    parser.add_argument(
        "--wavelength", type=float, default=0.65,
        help="Wavelength in meters (default: 0.65)"
    )
    parser.add_argument(
        "--prop-dir", type=float, nargs=3, default=[0, 0, 1],
        metavar=("X", "Y", "Z"),
        help="Propagation direction of incident wave (default: 0 0 1)"
    )
    parser.add_argument(
        "--polarization", type=float, nargs=3, default=[1, 0, 0],
        metavar=("X", "Y", "Z"),
        help="Polarization direction of incident wave (default: 1 0 0)"
    )

    # Geometry parameters
    parser.add_argument(
        "--outer-radius", type=float, default=1.0,
        help="Outer computational domain radius (default: 1.0)"
    )
    parser.add_argument(
        "--scatterer-radius", type=float, default=0.1,
        help="Scatterer radius (default: 0.1)"
    )
    parser.add_argument(
        "--pml-width", type=float, default=0.25,
        help="PML layer width (default: 0.25)"
    )
    parser.add_argument(
        "--mesh-size", type=float, default=0.5,
        help="Maximum mesh element size (default: 0.5)"
    )

    # FEM parameters
    parser.add_argument(
        "--order", type=int, default=5,
        help="Finite element order (default: 5)"
    )

    # Solver parameters
    parser.add_argument(
        "--solver", type=str, default="gmres", choices=["gmres", "direct"],
        help="Solver type (default: gmres)"
    )
    parser.add_argument(
        "--num-threads", type=int, default=4,
        help="Number of threads for parallel computation (default: 4)"
    )

    # Output parameters
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output file path for results (JSON format)"
    )
    parser.add_argument(
        "--quiet", action="store_true",
        help="Suppress output messages"
    )

    return parser.parse_args()


def main():
    """Main execution function."""
    args = parse_arguments()

    if not args.quiet:
        print("=" * 60)
        print("Electromagnetic Scattering Simulation")
        print("=" * 60)

    # Start timing
    t_start = time()

    # Compute wavenumber
    k = 2 * pi / args.wavelength

    if not args.quiet:
        print(f"\nPhysical Parameters:")
        print(f"  Wavelength: {args.wavelength}")
        print(f"  Wavenumber: {k:.4f}")
        print(f"  Propagation: {args.prop_dir}")
        print(f"  Polarization: {args.polarization}")

    # Create geometry
    if not args.quiet:
        print(f"\nCreating geometry...")
    t_mesh = time()
    mesh = create_spherical_geometry(
        R=args.outer_radius,
        PMLw=args.pml_width,
        r=args.scatterer_radius,
        h_max=args.mesh_size
    )
    t_mesh = time() - t_mesh

    if not args.quiet:
        print(f"  Mesh generated in {t_mesh:.3f}s")
        print(f"  Elements: {mesh.ne}")
        print(f"  Vertices: {mesh.nv}")

    # Create incident wave
    E_inc = create_incident_wave(
        k=k,
        propagation_dir=tuple(args.prop_dir),
        polarization=tuple(args.polarization)
    )

    # Setup problem
    if not args.quiet:
        print(f"\nSetting up scattering problem...")
    problem = ScattererProblem(mesh, k, E_inc, fes_order=args.order)

    if not args.quiet:
        print(f"\nAssembling system...")
    t_assembly = time()
    problem.assemble_system()
    t_assembly = time() - t_assembly

    if not args.quiet:
        print(f"  Assembly time: {t_assembly:.3f}s")

    # Solve
    if not args.quiet:
        print(f"\nSolving with {args.solver}...")
    t_solve = time()

    with TaskManager(pajetrace=10**8):
        if args.solver == "gmres":
            solution = solve_gmres(problem.a, problem.l, problem.fes)
        else:
            solution = solve_direct(problem.a, problem.l, problem.fes)

    t_solve = time() - t_solve

    problem.set_solution(solution)

    t_total = time() - t_start

    if not args.quiet:
        print(f"  Solve time: {t_solve:.3f}s")
        print(f"\nTotal time: {t_total:.3f}s")
        print("=" * 60)
        print("Simulation completed successfully!")
        print("=" * 60)

    # Save results if output specified
    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        results = {
            "problem_type": "scattering",
            "parameters": {
                "wavelength": args.wavelength,
                "wavenumber": k,
                "propagation_direction": args.prop_dir,
                "polarization": args.polarization,
                "outer_radius": args.outer_radius,
                "scatterer_radius": args.scatterer_radius,
                "pml_width": args.pml_width,
                "mesh_size": args.mesh_size,
                "fem_order": args.order,
                "solver": args.solver,
            },
            "mesh_info": {
                "elements": mesh.ne,
                "vertices": mesh.nv,
            },
            "fem_info": {
                "ndof": problem.fes.ndof,
                "free_dof": sum(problem.fes.FreeDofs()),
            },
            "timings": {
                "mesh_generation": t_mesh,
                "assembly": t_assembly,
                "solve": t_solve,
                "total": t_total,
            },
        }

        with open(output_path, "w") as f:
            json.dump(results, f, indent=2)

        if not args.quiet:
            print(f"\nResults saved to: {output_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
