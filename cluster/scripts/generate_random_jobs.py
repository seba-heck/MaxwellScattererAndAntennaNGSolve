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

    # Wave configurations - incident direction
    _directions = 2.0* np.random.random(size=(num_waves, 3)) - 1.0
    _directions *= 1.0 / np.linalg.norm(_directions, axis=1, keepdims=True)

    # polarization vectors
    _polarizations = np.random.random(size=(num_waves, 3))

    # replace parallel polarizations
    for d,p in zip(_directions, _polarizations):
        while abs(np.dot(d,p)) > 0.99:
            p = np.random.random(size=(3,))
        
    # orthogonalize polarization
    _polarizations = _polarizations - np.sum(_polarizations * _directions, axis=1, keepdims=True) * _directions

    if np.sum( np.abs( np.sum(_polarizations * _directions, axis=1, keepdims=True)) ) > 1e-6:
        raise ValueError("Polarization vectors are not orthogonal to propagation directions!")

    _polarizations *= 1.0 / np.linalg.norm(_polarizations, axis=1, keepdims=True)

    wave_configs = [ {
        "wavelength": float(wavelength),
        "propagation_dir": [i.item() for i in id], 
        "polarization": [i.item() for i in ip], 
        "label": f"wave_{i}"
        } for i, (wavelength,id,ip) in enumerate(zip(_wavelengths, _directions, _polarizations)) ]
    
    return wave_configs

def assemble_parameter_dict(wave_configs: List[Dict[str, Any]], geometries: List[Dict[str, Any]], geom_type: str, R: float, PMLw: float, h_max: float, curve_order: int, start_job_id: int = 0) -> List[Dict[str, Any]]:
    """Assemble parameter dictionary from names and values."""
    configs = []
    job_id = start_job_id

    # Generate all combinations
    for wave in wave_configs:
        for geom in geometries:
            config = {
                "job_id": job_id,
                "problem_type": "scattering",

                # Physical parameters
                "parameters": {
                    # "wavelength": float(wavelength),
                    
                    # "axis_a": geom["semi_axis_a_factor"],
                    # "axis_b": geom["semi_axis_b_factor"],
                    # "axis_c": geom["semi_axis_c_factor"],

                    # "edge_radius": geom["edge_radius"],

                    # # Incident wave
                    # "propagation_dir": wave["direction"],
                    # "polarization": wave["polarization"]
                },

                # Geometry
                "geometry": {
                    "type": geom_type,

                    "R": max(R, 1.5*wave["wavelength"]),
                    "PMLw": max(PMLw, 0.375*wave["wavelength"]),
                    "h_max": max(h_max, wave["wavelength"] / 5.0),

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
            config["parameters"].update(geom)

            configs.append(config)
            job_id += 1

    return configs

def generate_scattering_random_ellipsoid(
        num_geometries: int = 100, # number of geometry configurations
        radius_range: tuple = (0.05, 0.5), # scatterer radius range in meters
        wave_configs: List[Dict[str, Any]] = None, # list of wave configurations
        R: float = 1.0, # domain outer radius
        PMLw: float = 0.25, # PML-Layer width
        h_max: float = 0.08, # max mesh size
        curve_order: int = 5, # curve order of geometry
        start_job_id: int = 0
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
            "ellipsoid_semi_axis_a": float(a),
            "ellipsoid_semi_axis_b": float(b),
            "ellipsoid_semi_axis_c": float(c),
            "label": f"ellipsoid_{a:.2f}_{b:.2f}_{c:.2f}λ"
        })

    assert len(geometries) == num_geometries, f"Expected {num_geometries} geometries, got {len(geometries)}"

    # # Generate all combinations
    # job_id = 0
    # for wave in wave_configs:
    #     for geom in geometries:
    #         config = {
    #             "job_id": job_id,
    #             "problem_type": "scattering",

    #             # Physical parameters
    #             "parameters": {
    #                 # "wavelength": float(wavelength),
                    
    #                 "ellipsoid_semi_axis_a": geom["semi_axis_a_factor"],
    #                 "ellipsoid_semi_axis_b": geom["semi_axis_b_factor"],
    #                 "ellipsoid_semi_axis_c": geom["semi_axis_c_factor"],

    #                 # # Incident wave
    #                 # "propagation_dir": wave["direction"],
    #                 # "polarization": wave["polarization"]
    #             },

    #             # Geometry
    #             "geometry": {
    #                 "type": geom["type"],

    #                 "R": max(R, 1.5*wave["wavelength"]),
    #                 "PMLw": max(PMLw, 0.375*wave["wavelength"]),
    #                 "h_max": max(h_max, wave["wavelength"] / 8.0),

    #                 "curve_order": curve_order
    #             },

    #             # Solver
    #             "solver": {
    #                 "method": "gmres",
    #                 "preconditioner": "block_jacobi",
    #                 "fes_order": 3,
    #                 "maxiter": 2000,
    #                 "tol": 1e-6
    #             },

    #             # Output
    #             "output": {
    #                 "save_solution": True,
    #             },
    #         }

    #         config["parameters"].update(wave)

    #         configs.append(config)
    #         job_id += 1

    configs = assemble_parameter_dict(
        wave_configs=wave_configs,
        geometries=geometries,
        geom_type="ellipsoid",
        R=R, PMLw=PMLw, h_max=h_max, curve_order=curve_order
    )

    assert len(configs) == num_geometries*num_waves, f"Expected {num_geometries*num_waves} configs, got {len(configs)}"
    return configs

def generate_scattering_random_box(
        num_geometries: int = 100, # number of geometry configurations
        axis_range: tuple = (0.05, 0.5), # scatterer axis range in meters
        radius_range: tuple = (0.05, 0.2), # round edge radius range in meters
        wave_configs: List[Dict[str, Any]] = None, # list of wave configurations
        R: float = 1.0, # domain outer radius
        PMLw: float = 0.25, # PML-Layer width
        h_max: float = 0.08, # max mesh size
        curve_order: int = 5, # curve order of geometry
        start_job_id: int = 0
) -> List[Dict[str, Any]]:
    """
    Generate random scattering problem configurations.

    Strategy: 'num_geometries' geometries x 'num_waves' wave configs = total cases

    Returns:
        List of parameter dictionaries
    """

    if wave_configs is None:
        print("Warning: No wave configurations provided, generating default 100 wave samples.")
        wave_configs = generate_wave_samples()
    num_waves = len(wave_configs)

    # Geometry configurations
    geometries = []

    # (Tri-axial) box
    _axis = np.random.random(size=(num_geometries, 3)) * (axis_range[1] - axis_range[0]) + axis_range[0]
    _rad = np.random.random(size=(num_geometries, )) * (radius_range[1] - radius_range[0]) + radius_range[0]

    for a, b, c, r in zip(_axis[:, 0], _axis[:, 1], _axis[:, 2], _rad):
        geometries.append({
            "axis_a": float(a),
            "axis_b": float(b),
            "axis_c": float(c),
            "edge_radius": float( min(r, 0.49*min(a,b,c)) ),
            "label": f"box_{a:.2f}_{b:.2f}_{c:.2f}_r{r:.2f}",
        })

    assert len(geometries) == num_geometries, f"Expected {num_geometries} geometries, got {len(geometries)}"

    configs = assemble_parameter_dict(
        wave_configs=wave_configs,
        geometries=geometries,
        geom_type="box",
        R=R, PMLw=PMLw, h_max=h_max, curve_order=curve_order,
        start_job_id=start_job_id
    )

    assert len(configs) == num_geometries*num_waves, f"Expected {num_geometries*num_waves} configs, got {len(configs)}"
    return configs

def generate_scattering_random_cylinder(
        num_geometries: int = 100, # number of geometry configurations
        axis_range: tuple = (0.05, 0.5), # scatterer axis range in meters
        radius_range: tuple = (0.05, 0.2), # round edge radius range in meters
        wave_configs: List[Dict[str, Any]] = None, # list of wave configurations
        R: float = 1.0, # domain outer radius
        PMLw: float = 0.25, # PML-Layer width
        h_max: float = 0.08, # max mesh size
        curve_order: int = 5, # curve order of geometry
        start_job_id: int = 0
) -> List[Dict[str, Any]]:
    """
    Generate random scattering problem configurations.

    Strategy: 'num_geometries' geometries x 'num_waves' wave configs = total cases

    Returns:
        List of parameter dictionaries
    """

    if wave_configs is None:
        print("Warning: No wave configurations provided, generating default 100 wave samples.")
        wave_configs = generate_wave_samples()
    num_waves = len(wave_configs)

    # Geometry configurations
    geometries = []

    # cylinder
    _axis = np.random.random(size=(num_geometries, 2)) * (axis_range[1] - axis_range[0]) + axis_range[0]
    _height = (np.random.random(size=(num_geometries, )) * (axis_range[1] - axis_range[0]) + axis_range[0]) * 2.0
    _rad = np.random.random(size=(num_geometries, )) * (radius_range[1] - radius_range[0]) + radius_range[0]

    for a, b, c, r in zip(_axis[:, 0], _axis[:, 1], _height[:], _rad):
        geometries.append({
            "radius_major": float(a),
            "radius_minor": float(b),
            "height": float(c),
            "edge_radius": float( min(r, 0.49*min(a,b,c)) ),
            "label": f"cylinder_{a:.2f}_{b:.2f}_{c:.2f}_r{r:.2f}",
        })

    assert len(geometries) == num_geometries, f"Expected {num_geometries} geometries, got {len(geometries)}"

    configs = assemble_parameter_dict(
        wave_configs=wave_configs,
        geometries=geometries,
        geom_type="cylinder",
        R=R, PMLw=PMLw, h_max=h_max, curve_order=curve_order,
        start_job_id=start_job_id
    )

    assert len(configs) == num_geometries*num_waves, f"Expected {num_geometries*num_waves} configs, got {len(configs)}"
    return configs


def main():
    """Generate sweep job list."""

    # Create the argument parser
    parser = argparse.ArgumentParser(description="Create the jobs list (the configurations for each job).")

    parser.add_argument("--filename",  type=str,  default=None, help="Name of the output file (default: random_'problem'_'object'_jobs).")
    parser.add_argument("--problem",  type=str,  default="scatterer",        help="Type of problem: antenna or scatterer (default: scatterer).")
    # parser.add_argument("--object",    type=str,  default="ellipsoid", help="Type of geometry object for the simulation: ellipsoid (default: ellipsoid).")
    parser.add_argument("--ellipsoid", action='store_true', help='Make jobs with ellipsoid geometries.')
    parser.add_argument("--box", action='store_true', help='Make jobs with box geometries.')
    parser.add_argument("--cylinder", action='store_true', help='Make jobs with cylinder geometries.')
    parser.add_argument("--num-geometry",    type=str,  default='10',             help="Comma separated number of geometry configurations (default: 10).")
    parser.add_argument("--num-sources",    type=int,  default=10,             help="Number of incident waves/sources (default: 10).")
    # parser.add_argument("--idx_end",      type=int,  default=None,          help="Ending Job-Index (default: None).")
    # parser.add_argument("--verbose",      type=bool, default=True,          help="Create more output (default: True)")

    args = parser.parse_args()

    print("=" * 70)
    print("Generating Random Sweep Jobs Configurations")
    print(f"  for {args.problem} problem.")
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
    job_configs = []
    num_geometries = [int(x.strip()) for x in args.num_geometry.split(',')]
    sum_geometries = sum(num_geometries)

    if args.problem.lower() == "scatterer":
        # (1) SCATTERER: Random Job Configurations
        print(f"(2) Generating random scattering jobs ({sum_geometries * args.num_sources} jobs)", flush=True)

        if args.ellipsoid:
            print(f"   - {num_geometries[0]} Ellipsoid geometries...", flush=True)
            job_configs += generate_scattering_random_ellipsoid(wave_configs=wave_configs, num_geometries=num_geometries.pop(0))
        if args.box:
            print(f"   - {num_geometries[0]} Box geometries...", flush=True)
            job_configs += generate_scattering_random_box(wave_configs=wave_configs, num_geometries=num_geometries.pop(0), start_job_id=len(job_configs))
        if args.cylinder:
            print(f"   - {num_geometries[0]} Cylinder geometries...", flush=True)
            job_configs += generate_scattering_random_cylinder(wave_configs=wave_configs, num_geometries=num_geometries.pop(0), start_job_id=len(job_configs))
    else:
        print(f"Error: Unknown problem type '{args.problem}'. Supported: 'scatterer'.")
        sys.exit(1)
    
    print(f" > Finished generating job configurations.")
    print(f"   Geometries: {sum_geometries} geometreries")
    print(f"   Type: {args.ellipsoid * 'Ellipsoid '} {args.box * 'Box '} {args.cylinder * 'Cylinder '}")

    # Save to file
    print(f"(3) Saving job configurations to file...", end="", flush=True)
    if args.filename is None:
        object_str = ("ellipsoid_" if args.ellipsoid else "") + ("box_" if args.box else "") + ("cylinder_" if args.cylinder else "")
        args.filename = f"random_{args.problem}_{object_str}jobs"
    
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
