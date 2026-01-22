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

def generate_scattering_random_ellipsoid(
        num_geometries: int = 100, radius_range: tuple = (0.05, 0.5),
        num_waves: int = 100, wave_depth_range: tuple = (1, 10)
) -> List[Dict[str, Any]]:
    """
    Generate 1000 scattering problem configurations.

    Strategy: 5 wavelengths × 40 geometries × 5 wave configs = 1000 cases

    Returns:
        List of parameter dictionaries
    """
    R = 1.0
    PMLw = 0.25
    h_max = 0.01
    curve_order = 5

    configs = []

    # Wavelength sweep
    _c = np.random.random(num_waves) * (wave_depth_range[1] - wave_depth_range[0]) + wave_depth_range[0]
    _d = 1.0
    wavelengths = _d / _c

    # Wave configurations (incident direction, polarization)
    _directions = np.random.random(size=(num_waves, 3))
    _directions *= 1.0 / np.linalg.norm(_directions, axis=1, keepdims=True)

    _polarizations = np.random.random(size=(num_waves, 3))
    _polarizations *= 1.0 / np.linalg.norm(_polarizations, axis=1, keepdims=True)

    wave_configs = [ {"direction": [i.item() for i in id], "polarization": [i.item() for i in ip], "label": f"wave_{i}"} for i, (id,ip) in enumerate(zip(_directions, _polarizations)) ]

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
    for wavelength,wave in zip(wavelengths, wave_configs):
        for geom in geometries:
            config = {
                "job_id": job_id,
                "problem_type": "scattering",

                # Physical parameters
                "parameters": {
                    "wavelength": float(wavelength),
                    
                    "ellipsoid_semi_axis_a": geom["semi_axis_a_factor"],
                    "ellipsoid_semi_axis_b": geom["semi_axis_b_factor"],
                    "ellipsoid_semi_axis_c": geom["semi_axis_c_factor"],

                    # Incident wave
                    "propagation_dir": wave["direction"],
                    "polarization": wave["polarization"]
                },

                # Geometry
                "geometry": {
                    "type": geom["type"],

                    "R": R,
                    "PMLw": PMLw,
                    "h_max": min(h_max, float(wavelength) / 6.0),

                    "curve_order": curve_order
                },

                # Solver
                "solver": {
                    "method": "direct",
                    "preconditioner": "block_jacobi",
                    "fes_order": 3,
                    "maxiter": 1000,
                    "tol": 1e-6
                },

                # Output
                "output": {
                    "save_solution": True,
                },
            }

            configs.append(config)
            job_id += 1

    assert len(configs) == num_geometries*num_waves, f"Expected {num_geometries*num_waves} configs, got {len(configs)}"
    return configs


def main():
    """Generate sweep configurations."""

    print("=" * 70)
    print("Generating Random Ellipsoid Sweep Configurations")
    print("=" * 70)
    print()

    # Create configs directory if it doesn't exist
    configs_dir = Path(__file__).parent.parent / "configs"
    configs_dir.mkdir(exist_ok=True)

    # Generate scattering sweep
    print("1. Generating scattering random ellipsoid (10'000 jobs)...")
    scattering_configs = generate_scattering_random_ellipsoid()
    with open(configs_dir / "scattering_random_ellipsoid_sweep_jobs.json", 'w') as f:
        json.dump(scattering_configs, f, indent=2)
    print(f"   Wavelengths: 0.1m to 1.0m")
    print(f"   Geometries: 100 types (ellipsoids)")
    print(f"   Wave configs: 100 incident directions")
    print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
✓ Created 1 configuration files:

  scattering_random_ellipsoid_sweep_jobs.json  - 10'000 scattering problems

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
