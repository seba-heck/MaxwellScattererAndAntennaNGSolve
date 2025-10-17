#!/usr/bin/env python3
"""
Generate 1000-case parameter sweep configurations.

Creates YAML configuration files for:
1. Scattering problems (1000 cases)
2. Antenna problems (1000 cases)

Each configuration explores diverse physical regimes through
systematic variation of dimensionless parameters.
"""

import yaml
import itertools
from pathlib import Path
from typing import List, Dict, Any
import numpy as np


def generate_scattering_sweep() -> List[Dict[str, Any]]:
    """
    Generate 1000 scattering problem configurations.

    Strategy: 5 wavelengths × 40 geometries × 5 wave configs = 1000 cases

    Returns:
        List of parameter dictionaries
    """
    configs = []

    # Wavelength sweep
    wavelengths = [0.5, 0.75, 1.0, 1.5, 2.0]

    # Wave configurations (incident direction, polarization)
    wave_configs = [
        {"direction": [0, 0, 1], "polarization": [1, 0, 0], "label": "z_prop"},
        {"direction": [1, 0, 0], "polarization": [0, 1, 0], "label": "x_prop"},
        {"direction": [0, 1, 0], "polarization": [0, 0, 1], "label": "y_prop"},
        {"direction": [0.7071, 0.7071, 0], "polarization": [0, 0, 1], "label": "xy_diag"},
        {"direction": [0.5774, 0.5774, 0.5774], "polarization": [0.7071, -0.7071, 0], "label": "xyz_diag"},
    ]

    # Geometry configurations (40 total)
    geometries = []

    # 1. Spheres (5 configurations)
    for r_factor in [0.1, 0.15, 0.2, 0.25, 0.3]:
        geometries.append({
            "type": "sphere",
            "radius_factor": r_factor,
            "label": f"sphere_r{r_factor:.2f}λ"
        })

    # 2. Prolate spheroids - "cigar" (10 configurations)
    prolate_configs = [
        (0.1, 0.3), (0.1, 0.4), (0.1, 0.5),
        (0.15, 0.3), (0.15, 0.4), (0.15, 0.5),
        (0.2, 0.3), (0.2, 0.4), (0.2, 0.5),
        (0.15, 0.45),  # Additional intermediate
    ]
    for a_eq, c_polar in prolate_configs:
        geometries.append({
            "type": "spheroid",
            "equatorial_radius_factor": a_eq,
            "polar_radius_factor": c_polar,
            "label": f"prolate_a{a_eq:.2f}_c{c_polar:.2f}λ"
        })

    # 3. Oblate spheroids - "pancake" (10 configurations)
    oblate_configs = [
        (0.3, 0.1), (0.3, 0.15), (0.3, 0.2),
        (0.4, 0.1), (0.4, 0.15), (0.4, 0.2),
        (0.5, 0.1), (0.5, 0.15), (0.5, 0.2),
        (0.35, 0.15),  # Additional intermediate
    ]
    for a_eq, c_polar in oblate_configs:
        geometries.append({
            "type": "spheroid",
            "equatorial_radius_factor": a_eq,
            "polar_radius_factor": c_polar,
            "label": f"oblate_a{a_eq:.2f}_c{c_polar:.2f}λ"
        })

    # 4. Tri-axial ellipsoids (15 configurations)
    # Explore different aspect ratios: a:b:c
    ellipsoid_configs = [
        # Moderate asymmetry (1:1.5:2 ratio)
        (0.1, 0.15, 0.2),
        (0.15, 0.225, 0.3),
        (0.2, 0.3, 0.4),
        # Strong asymmetry (1:2:3 ratio)
        (0.1, 0.2, 0.3),
        (0.12, 0.24, 0.36),
        (0.15, 0.3, 0.45),
        # Mixed ratios (1:1.5:3)
        (0.1, 0.15, 0.3),
        (0.12, 0.18, 0.36),
        (0.15, 0.225, 0.45),
        # Biaxial variations
        (0.1, 0.1, 0.3),   # Prolate-like
        (0.3, 0.3, 0.1),   # Oblate-like
        (0.15, 0.2, 0.25), # Nearly spherical
        (0.1, 0.2, 0.4),   # Strong variation
        (0.2, 0.25, 0.3),  # Mild variation
        (0.15, 0.25, 0.35), # Intermediate
    ]
    for a, b, c in ellipsoid_configs:
        geometries.append({
            "type": "ellipsoid",
            "semi_axis_a_factor": a,
            "semi_axis_b_factor": b,
            "semi_axis_c_factor": c,
            "label": f"ellipsoid_{a:.2f}_{b:.2f}_{c:.2f}λ"
        })

    assert len(geometries) == 40, f"Expected 40 geometries, got {len(geometries)}"

    # Generate all combinations
    job_id = 0
    for wavelength in wavelengths:
        for geom in geometries:
            for wave in wave_configs:
                config = {
                    "job_id": job_id,
                    "problem_type": "scattering",

                    # Physical parameters
                    "wavelength": wavelength,
                    "wavenumber": 2 * np.pi / wavelength,

                    # Geometry
                    "geometry_type": geom["type"],

                    # Incident wave
                    "incident_direction": wave["direction"],
                    "incident_polarization": wave["polarization"],

                    # Domain
                    "domain_radius": 2.0 * wavelength,
                    "pml_width": 0.25 * wavelength,

                    # Mesh
                    "max_mesh_size": wavelength / 15.0,

                    # Solver
                    "fes_order": 3,
                    "solver_method": "direct",

                    # Output
                    "save_solution": False,  # VTK off by default for large sweeps

                    # Labels for analysis
                    "geometry_label": geom["label"],
                    "wave_label": wave["label"],
                }

                # Add geometry-specific parameters
                if geom["type"] == "sphere":
                    config["sphere_radius"] = geom["radius_factor"] * wavelength

                elif geom["type"] == "spheroid":
                    config["spheroid_equatorial_radius"] = geom["equatorial_radius_factor"] * wavelength
                    config["spheroid_polar_radius"] = geom["polar_radius_factor"] * wavelength

                elif geom["type"] == "ellipsoid":
                    config["ellipsoid_semi_axis_a"] = geom["semi_axis_a_factor"] * wavelength
                    config["ellipsoid_semi_axis_b"] = geom["semi_axis_b_factor"] * wavelength
                    config["ellipsoid_semi_axis_c"] = geom["semi_axis_c_factor"] * wavelength
                    config["ellipsoid_orientation"] = "z"

                configs.append(config)
                job_id += 1

    assert len(configs) == 1000, f"Expected 1000 configs, got {len(configs)}"
    return configs


def generate_antenna_sweep() -> List[Dict[str, Any]]:
    """
    Generate 1000 antenna problem configurations.

    Strategy: 5 wavelengths × 10 lengths × 4 radii × 5 excitation/orientation = 1000 cases

    Returns:
        List of parameter dictionaries
    """
    configs = []

    # Wavelength sweep
    wavelengths = [0.5, 0.75, 1.0, 1.5, 2.0]

    # Dipole lengths (as fraction of wavelength)
    length_factors = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0]

    # Dipole radii (as fraction of wavelength)
    radius_factors = [0.005, 0.01, 0.015, 0.02]

    # Excitation and orientation configurations (5 total)
    excitation_configs = [
        {"orientation": "z", "amplitude": 1.0, "feed": "center", "label": "z_I1.0"},
        {"orientation": "z", "amplitude": 2.0, "feed": "center", "label": "z_I2.0"},
        {"orientation": "x", "amplitude": 1.0, "feed": "center", "label": "x_I1.0"},
        {"orientation": "y", "amplitude": 1.0, "feed": "center", "label": "y_I1.0"},
        {"orientation": "z", "amplitude": 1.0, "feed": "offset", "label": "z_offset"},
    ]

    # Generate all combinations
    job_id = 0
    for wavelength in wavelengths:
        for L_factor in length_factors:
            for r_factor in radius_factors:
                for excite in excitation_configs:

                    dipole_length = L_factor * wavelength
                    dipole_radius = r_factor * wavelength

                    config = {
                        "job_id": job_id,
                        "problem_type": "antenna",

                        # Physical parameters
                        "wavelength": wavelength,
                        "wavenumber": 2 * np.pi / wavelength,

                        # Geometry
                        "geometry_type": "dipole",
                        "dipole_length": dipole_length,
                        "dipole_radius": dipole_radius,
                        "dipole_orientation": excite["orientation"],

                        # Excitation
                        "source_amplitude": excite["amplitude"],
                        "feed_position": excite["feed"],

                        # Domain
                        "domain_radius": 2.0 * wavelength,
                        "pml_width": 0.25 * wavelength,

                        # Mesh
                        "max_mesh_size": wavelength / 15.0,

                        # Solver
                        "fes_order": 3,
                        "solver_method": "direct",

                        # Output
                        "save_solution": False,

                        # Labels for analysis
                        "length_label": f"L{L_factor:.1f}λ",
                        "radius_label": f"r{r_factor:.3f}λ",
                        "aspect_ratio": dipole_length / dipole_radius,
                        "excitation_label": excite["label"],
                    }

                    configs.append(config)
                    job_id += 1

    assert len(configs) == 1000, f"Expected 1000 configs, got {len(configs)}"
    return configs


def generate_pilot_study() -> Dict[str, List[Dict[str, Any]]]:
    """
    Generate pilot study configurations for validation.

    Selects representative cases from each physical regime:
    - Scattering: Small (ka~0.5), Medium (ka~1.5), Large (ka~2.5)
    - Antenna: Short (L~0.3λ), Resonant (L~0.5λ), Long (L~1.5λ)

    Returns:
        Dictionary with 'scattering' and 'antenna' pilot configs
    """
    pilot = {"scattering": [], "antenna": []}

    # Scattering pilots (9 cases)
    scattering_pilots = [
        # Small sphere (ka ~ 0.5)
        {"wavelength": 1.0, "type": "sphere", "radius": 0.1, "direction": [0, 0, 1]},
        {"wavelength": 1.0, "type": "sphere", "radius": 0.1, "direction": [1, 0, 0]},
        {"wavelength": 1.0, "type": "sphere", "radius": 0.1, "direction": [0.5774, 0.5774, 0.5774]},

        # Medium ellipsoid (ka ~ 1.5)
        {"wavelength": 1.0, "type": "ellipsoid", "semi_axes": (0.15, 0.225, 0.3), "direction": [0, 0, 1]},
        {"wavelength": 1.0, "type": "ellipsoid", "semi_axes": (0.15, 0.225, 0.3), "direction": [1, 0, 0]},

        # Large spheroid (ka ~ 2.5)
        {"wavelength": 1.0, "type": "spheroid", "equatorial": 0.3, "polar": 0.5, "direction": [0, 0, 1]},
        {"wavelength": 1.0, "type": "spheroid", "equatorial": 0.3, "polar": 0.5, "direction": [1, 0, 0]},
        {"wavelength": 1.0, "type": "spheroid", "equatorial": 0.5, "polar": 0.3, "direction": [0, 0, 1]},
        {"wavelength": 1.0, "type": "spheroid", "equatorial": 0.5, "polar": 0.3, "direction": [1, 0, 0]},
    ]

    job_id = 0
    for p in scattering_pilots:
        wavelength = p["wavelength"]
        config = {
            "job_id": job_id,
            "problem_type": "scattering",
            "wavelength": wavelength,
            "wavenumber": 2 * np.pi / wavelength,
            "geometry_type": p["type"],
            "incident_direction": p["direction"],
            "incident_polarization": [1, 0, 0] if p["direction"] != [1, 0, 0] else [0, 1, 0],
            "domain_radius": 2.0 * wavelength,
            "pml_width": 0.25 * wavelength,
            "max_mesh_size": wavelength / 15.0,
            "fes_order": 3,
            "solver_method": "direct",
            "save_solution": True,  # Enable VTK for pilots
        }

        if p["type"] == "sphere":
            config["sphere_radius"] = p["radius"] * wavelength
        elif p["type"] == "ellipsoid":
            a, b, c = p["semi_axes"]
            config["ellipsoid_semi_axis_a"] = a * wavelength
            config["ellipsoid_semi_axis_b"] = b * wavelength
            config["ellipsoid_semi_axis_c"] = c * wavelength
            config["ellipsoid_orientation"] = "z"
        elif p["type"] == "spheroid":
            config["spheroid_equatorial_radius"] = p["equatorial"] * wavelength
            config["spheroid_polar_radius"] = p["polar"] * wavelength

        pilot["scattering"].append(config)
        job_id += 1

    # Antenna pilots (12 cases)
    antenna_pilots = [
        # Short dipole (L = 0.3λ)
        {"wavelength": 1.0, "length_factor": 0.3, "radius_factor": 0.01, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 0.3, "radius_factor": 0.02, "orientation": "z"},

        # Resonant (L = 0.5λ) - most important
        {"wavelength": 1.0, "length_factor": 0.5, "radius_factor": 0.005, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 0.5, "radius_factor": 0.01, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 0.5, "radius_factor": 0.015, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 0.5, "radius_factor": 0.01, "orientation": "x"},

        # Full-wave (L = 1.0λ)
        {"wavelength": 1.0, "length_factor": 1.0, "radius_factor": 0.01, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 1.0, "radius_factor": 0.015, "orientation": "z"},

        # Long dipole (L = 1.5λ)
        {"wavelength": 1.0, "length_factor": 1.5, "radius_factor": 0.01, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 1.5, "radius_factor": 0.015, "orientation": "z"},

        # Very long (L = 2.0λ)
        {"wavelength": 1.0, "length_factor": 2.0, "radius_factor": 0.01, "orientation": "z"},
        {"wavelength": 1.0, "length_factor": 2.0, "radius_factor": 0.015, "orientation": "z"},
    ]

    job_id = 0
    for p in antenna_pilots:
        wavelength = p["wavelength"]
        config = {
            "job_id": job_id,
            "problem_type": "antenna",
            "wavelength": wavelength,
            "wavenumber": 2 * np.pi / wavelength,
            "geometry_type": "dipole",
            "dipole_length": p["length_factor"] * wavelength,
            "dipole_radius": p["radius_factor"] * wavelength,
            "dipole_orientation": p["orientation"],
            "source_amplitude": 1.0,
            "feed_position": "center",
            "domain_radius": 2.0 * wavelength,
            "pml_width": 0.25 * wavelength,
            "max_mesh_size": wavelength / 15.0,
            "fes_order": 3,
            "solver_method": "direct",
            "save_solution": True,  # Enable VTK for pilots
            "length_label": f"L{p['length_factor']:.1f}λ",
            "radius_label": f"r{p['radius_factor']:.3f}λ",
        }

        pilot["antenna"].append(config)
        job_id += 1

    return pilot


def save_sweep_config(configs: List[Dict[str, Any]], filename: str):
    """Save sweep configuration as YAML file."""
    output_path = Path(__file__).parent.parent / "configs" / filename

    # Create sweep metadata
    sweep_data = {
        "metadata": {
            "total_jobs": len(configs),
            "problem_type": configs[0]["problem_type"],
            "description": f"{len(configs)}-case parameter sweep",
        },
        "jobs": configs
    }

    with open(output_path, 'w') as f:
        yaml.dump(sweep_data, f, default_flow_style=False, sort_keys=False)

    print(f"✓ Created {filename} with {len(configs)} configurations")
    print(f"  Saved to: {output_path}")


def main():
    """Generate all sweep configurations."""

    print("=" * 70)
    print("Generating 1000-Case Parameter Sweep Configurations")
    print("=" * 70)
    print()

    # Create configs directory if it doesn't exist
    configs_dir = Path(__file__).parent.parent / "configs"
    configs_dir.mkdir(exist_ok=True)

    # Generate scattering sweep
    print("1. Generating scattering sweep (1000 cases)...")
    scattering_configs = generate_scattering_sweep()
    save_sweep_config(scattering_configs, "scattering_1000_sweep.yaml")
    print(f"   Wavelengths: 0.5m to 2.0m")
    print(f"   Geometries: 40 types (spheres, spheroids, ellipsoids)")
    print(f"   Wave configs: 5 incident directions")
    print()

    # Generate antenna sweep
    print("2. Generating antenna sweep (1000 cases)...")
    antenna_configs = generate_antenna_sweep()
    save_sweep_config(antenna_configs, "antenna_1000_sweep.yaml")
    print(f"   Wavelengths: 0.5m to 2.0m")
    print(f"   Dipole lengths: 0.3λ to 2.0λ (10 values)")
    print(f"   Dipole radii: 0.005λ to 0.02λ (4 values)")
    print(f"   Excitations: 5 orientation/amplitude combos")
    print()

    # Generate pilot study
    print("3. Generating pilot study configurations...")
    pilot_configs = generate_pilot_study()
    save_sweep_config(pilot_configs["scattering"], "pilot_scattering.yaml")
    save_sweep_config(pilot_configs["antenna"], "pilot_antenna.yaml")
    print(f"   Scattering pilots: {len(pilot_configs['scattering'])} representative cases")
    print(f"   Antenna pilots: {len(pilot_configs['antenna'])} representative cases")
    print()

    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"""
✓ Created 4 configuration files:

  scattering_1000_sweep.yaml  - 1000 scattering problems
  antenna_1000_sweep.yaml     - 1000 antenna problems
  pilot_scattering.yaml       - 9 validation cases
  pilot_antenna.yaml          - 12 validation cases

Next steps:
  1. Run pilot studies to validate convergence
  2. Review results and adjust parameters if needed
  3. Deploy full 1000-case sweeps on cluster

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
