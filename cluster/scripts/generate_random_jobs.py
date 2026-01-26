#!/usr/bin/env python3
"""
Generate random job parameter file from YAML configuration.

This script reads a YAML configuration file defining a parametric sweep
and generates a JSON file containing all parameter combinations for each job.

Usage:
    python generate_jobs.py configs/spherical_sweep.yaml
    python generate_jobs.py configs/my_sweep.yaml --validate-only
"""

import argparse
import sys
import json
import yaml
import numpy as np
from pathlib import Path
from itertools import product
from typing import Dict, List, Any

def generate_wave_samples(num_waves: int = 100, wavelength_range: tuple = (0.33, 2.0)) -> List[Dict[str, Any]]:
    # Wavelength sweep
    _wavelengths = np.random.random(num_waves) * (wavelength_range[1] - wavelength_range[0]) + wavelength_range[0]

    # Wave configurations (incident direction, polarization)
    _directions = np.random.random(size=(num_waves, 3))
    _directions *= 1.0 / np.linalg.norm(_directions, axis=1, keepdims=True)

    _polarizations = np.random.random(size=(num_waves, 3))
    _polarizations *= 1.0 / np.linalg.norm(_polarizations, axis=1, keepdims=True)

    wave_configs = [ {
        "wavelength": float(wavelength),
        "propagation_dir": [i.item() for i in id], 
        "polarization": [i.item() for i in ip], 
        "label": f"wave_{i}"
        } for i, (wavelength,id,ip) in enumerate(zip(_wavelengths, _directions, _polarizations)) ]
    
    return wave_configs

def generate_scattering_random_ellipsoid(
        num_geometries: int = 100, # number of geometry configurations
        radius_range: tuple = (0.05, 0.5), # scatterer radius range in meters
        wave_configs: List[Dict[str, Any]] = None, # list of wave configurations
        R: float = 1.0, # domain outer radius
        PMLw: float = 0.25, # PML-Layer width
        h_max: float = 0.08, # max mesh size
        curve_order: int = 5 # curve order of geometry
) -> List[Dict[str, Any]]:
    """
    Generate 1000 scattering problem configurations.

    Strategy: 5 wavelengths × 40 geometries × 5 wave configs = 1000 cases

    Returns:
        List of parameter dictionaries
    """

    configs = []

    if wave_configs is None:
        print("Warning: No wave configurations provided, generating default 100 wave samples.")
        wave_configs = generate_wave_samples()
    num_waves = len(wave_configs)

    # Geometry configurations
    geometries = []

    # Tri-axial ellipsoids
    _axis = np.random.random(size=(num_geometries, 3)) * (radius_range[1] - radius_range[0]) + radius_range[0]

    for a, b, c in _axis:
        geometries.append({
            "type": "ellipsoid",
            "semi_axis_a_factor": float(a),
            "semi_axis_b_factor": float(b),
            "semi_axis_c_factor": float(c),
            "label": f"ellipsoid_{a:.2f}_{b:.2f}_{c:.2f}λ"
        })

    assert len(geometries) == num_geometries, f"Expected {num_geometries} geometries, got {len(geometries)}"

    # Generate all combinations
    job_id = 0
    for wave in wave_configs:
        for geom in geometries:
            config = {
                "job_id": job_id,
                "problem_type": "scattering",

                # Physical parameters
                "parameters": {
                    # "wavelength": float(wavelength),
                    
                    "ellipsoid_semi_axis_a": geom["semi_axis_a_factor"],
                    "ellipsoid_semi_axis_b": geom["semi_axis_b_factor"],
                    "ellipsoid_semi_axis_c": geom["semi_axis_c_factor"],

                    # # Incident wave
                    # "propagation_dir": wave["direction"],
                    # "polarization": wave["polarization"]
                },

                # Geometry
                "geometry": {
                    "type": geom["type"],

                    "R": max(R, 1.5*wave["wavelength"]),
                    "PMLw": max(PMLw, 0.375*wave["wavelength"]),
                    "h_max": max(h_max, wave["wavelength"] / 8.0),

                    "curve_order": curve_order
                },

                # Solver
                "solver": {
                    "method": "gmres",
                    "preconditioner": "block_jacobi",
                    "fes_order": 3,
                    "maxiter": 2000,
                    "tol": 1e-6
                },

                # Output
                "output": {
                    "save_solution": True,
                },
            }

            config["parameters"].update(wave)

            configs.append(config)
            job_id += 1

    assert len(configs) == num_geometries*num_waves, f"Expected {num_geometries*num_waves} configs, got {len(configs)}"
    return configs


def main():
    """Generate sweep job list."""

    # Create the argument parser
    parser = argparse.ArgumentParser(description="Create the jobs list (the configurations for each job).")

    parser.add_argument("--filename",  type=str,  default=None, help="Name of the output file (default: random_'problem'_'object'_jobs).")
    parser.add_argument("--problem",  type=str,  default="scatterer",        help="Type of problem: antenna or scatterer (default: scatterer).")
    parser.add_argument("--object",    type=str,  default="ellipsoid", help="Type of geometry object for the simulation: ellipsoid (default: ellipsoid).")
    parser.add_argument("--num-geometry",    type=int,  default=10,             help="Number of geometry configurations (default: 10).")
    parser.add_argument("--num-sources",    type=int,  default=10,             help="Number of incident waves/sources (default: 10).")
    # parser.add_argument("--idx_end",      type=int,  default=None,          help="Ending Job-Index (default: None).")
    # parser.add_argument("--verbose",      type=bool, default=True,          help="Create more output (default: True)")

    args = parser.parse_args()

    print("=" * 70)
    print("Generating Random Sweep Jobs Configurations")
    print(f"  for {args.problem} problem with {args.object} geometries.")
    print("=" * 70)
    print()

    # Create configs directory if it doesn't exist
    configs_dir = Path(__file__).parent.parent / "configs"
    configs_dir.mkdir(exist_ok=True)

    # Create wave configurations
    print(f"(1) Generating random incident wave configurations ({args.num_sources} waves)...", end="", flush=True)
    wave_configs = generate_wave_samples(num_waves=args.num_sources)
    print("  Done.", flush=True)
    print(f"   Wavelengths: 0.33m to 2.0m")
    print(f"   Wave configs: {args.num_sources} unique incident waves")

    # Generate full job configurations
    job_configs = None

    if args.problem.lower() == "scatterer":
        # (1) SCATTERER: Random Ellipsoids
        print(f"(2) Generating random scattering jobs with ellipsoids ({args.num_geometry * args.num_sources} jobs)...", end="", flush=True)
        job_configs = generate_scattering_random_ellipsoid(wave_configs=wave_configs, num_geometries=args.num_geometry)
    else:
        print(f"Error: Unknown problem type '{args.problem}'. Supported: 'scatterer'.")
        sys.exit(1)
    
    print("  Done.", flush=True)
    print(f"   Geometries: {args.num_geometry} geometreries")
    print(f"   Type: {args.object}")

    # Save to file
    print(f"(3) Saving job configurations to file...", end="", flush=True)
    if args.filename is None:
        args.filename = f"random_{args.problem}_{args.object}_jobs"
    
    output_file = (configs_dir / (args.filename + ".json")).resolve()
    with open(output_file, 'w') as f:
        json.dump(job_configs, f, indent=2)

    print("  Done.", flush=True)
    print(f"   Output file: {output_file}")

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
✓ Created 1 configuration files:

  {args.filename}.json  - {len(job_configs)} {args.problem} problems

Next steps:
  1. Run pilot studies to validate convergence
  2. Review results and adjust parameters if needed
  3. Deploy full 10'000-case sweeps on cluster

Usage:
  # Local test (single job)
  python3 cluster/scripts/run_simulation.py pilot_scattering 0

  # Cluster submission (example for pilot)
  sbatch --array=0-8 cluster/slurm/array_job.slurm pilot_scattering

  # Full sweep (do NOT run yet - validation first!)
  sbatch --array=0-999 cluster/slurm/array_job.slurm scattering_1000_sweep
""")


if __name__ == "__main__":
    main()
