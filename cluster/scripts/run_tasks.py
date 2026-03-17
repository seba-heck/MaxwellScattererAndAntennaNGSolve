#!/usr/bin/env python3
"""
Run a single electromagnetic simulation with given parameters.

This script supports both SLURM cluster execution and local testing.
It automatically detects the environment and adapts accordingly.

Usage:
    # On cluster (called by SLURM)
    python run_simulation.py --job-file jobs.json --job-id 0 --output-dir results/

    # Local testing
    python run_simulation.py --job-file jobs.json --job-id 0 --output-dir results/ --local
"""

import argparse
import sys
import json
import time
import os
from pathlib import Path
from typing import Dict, Any
import numpy as np
from scipy import constants
import pyvista as pv

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from ngsolve import TaskManager, SetNumThreads, pi, VTKOutput
from ngsolve.webgui import Draw
from src import (
    create_spherical_geometry,
    create_ellipsoid_scatterer_geometry,
    create_box_scatterer_geometry,
    create_cylinder_scatterer_geometry,
    create_spheroid_scatterer_geometry,
    create_dipole_antenna_geometry,
    create_incident_wave,
    create_antenna_source,
    MaxwellProblem,
    solve_gmres,
    solve_direct
)


def is_slurm_environment() -> bool:
    """Check if running in SLURM environment."""
    return 'SLURM_JOB_ID' in os.environ


def load_job_parameters(job_file: Path, job_id: int) -> Dict[str, Any]:
    """Load parameters for specific job from JSON file."""
    with open(job_file, 'r') as f:
        jobs = json.load(f)

    if job_id < 0 or job_id >= len(jobs):
        raise ValueError(f"Invalid job_id {job_id}. Valid range: 0-{len(jobs)-1}")

    return jobs[job_id]


def create_geometry(geometry_config: Dict[str, Any], params: Dict[str, Any], verbose: bool = True, surface_mesh: bool = True):
    """Create geometry and mesh based on configuration."""
    wavelength = params['wavelength']
    
    if geometry_config["type"] == "ellipsoid" or 'ellipsoid_semi_axis_a' in params:
        # Tri-axial ellipsoid scatterer
        if verbose: print(f"  → Creating tri-axial ellipsoid scatterer geometry")
        mesh = create_ellipsoid_scatterer_geometry(
            wavelength=wavelength,
            semi_axis_a=params['ellipsoid_semi_axis_a'],
            semi_axis_b=params['ellipsoid_semi_axis_b'],
            semi_axis_c=params['ellipsoid_semi_axis_c'],
            domain_radius=geometry_config.get('R', 1.0),
            pml_width=geometry_config.get('PMLw', 0.25),
            max_mesh_size=geometry_config.get('h_max', None),
            orientation=params.get('ellipsoid_orientation', 'z'),
            curve_order=geometry_config.get('curve_order', 5),
            surface_mesh=surface_mesh
        )
        geometry_type = 'ellipsoid'

    elif geometry_config["type"] == "box":
        # Box geometry
        if verbose: print(f"  → Creating box geometry")
        mesh = create_box_scatterer_geometry(
            wavelength=wavelength,
            axis_a=params['axis_a'],
            axis_b=params['axis_b'],
            axis_c=params['axis_c'],
            box_radius=params['edge_radius'],
            domain_radius=geometry_config.get('R', 1.0),
            pml_width=geometry_config.get('PMLw', 0.25),
            max_mesh_size=geometry_config.get('h_max', None),
            curve_order=geometry_config.get('curve_order', 5),
            surface_mesh=surface_mesh
        )
        geometry_type = 'box'

    elif geometry_config["type"] == "cylinder":
        # Cylinder geometry
        if verbose: print(f"  → Creating cylinder geometry")
        mesh = create_cylinder_scatterer_geometry(
            wavelength=wavelength,
            height=params['height'],
            radius=params['radius_major'],
            radius_2=params.get('radius_minor', None),
            box_radius=params['edge_radius'],
            domain_radius=geometry_config.get('R', 1.0),
            pml_width=geometry_config.get('PMLw', 0.25),
            max_mesh_size=geometry_config.get('h_max', None),
            curve_order=geometry_config.get('curve_order', 5),
            surface_mesh=surface_mesh
        )
        geometry_type = 'cylinder'

    elif 'dipole_length_factor' in params:
        # Dipole antenna geometry
        if verbose: print(f"  → Creating cylindrical dipole antenna geometry")
        mesh = create_dipole_antenna_geometry(
            wavelength=wavelength,
            length_factor=params['dipole_length_factor'],
            radius_factor=params.get('dipole_radius_factor', 0.01),
            domain_radius=geometry_config.get('R', 1.0),
            pml_width=geometry_config.get('PMLw', 0.25),
            max_mesh_size=geometry_config.get('h_max', None),
            orientation=params.get('dipole_orientation', 'z'),
            curve_order=geometry_config.get('curve_order', 5)
        )
        geometry_type = 'dipole'

    elif 'spheroid_equatorial_radius' in params:
        # Spheroid scatterer
        if verbose: print(f"  → Creating spheroid scatterer geometry")
        mesh = create_spheroid_scatterer_geometry(
            wavelength=wavelength,
            equatorial_radius=params['spheroid_equatorial_radius'],
            polar_radius=params['spheroid_polar_radius'],
            domain_radius=geometry_config.get('R', 1.0),
            pml_width=geometry_config.get('PMLw', 0.25),
            max_mesh_size=geometry_config.get('h_max', None),
            orientation=params.get('spheroid_orientation', 'z'),
            curve_order=geometry_config.get('curve_order', 5)
        )
        geometry_type = 'spheroid'

    else:
        # Default: spherical geometry (backward compatibility)
        if verbose: print(f"  → Creating spherical geometry (default)")
        mesh = create_spherical_geometry(
            R=geometry_config.get('R', 1.0),
            PMLw=geometry_config.get('PMLw', 0.25),
            r=geometry_config.get('r', 0.1),
            h_max=geometry_config.get('h_max', 0.5),
            curve_order=geometry_config.get('curve_order', 5)
        )
        geometry_type = 'sphere'
    
    return mesh, geometry_type

def run_task(job_params: Dict[str, Any], job_path: str) -> Dict[str, Any]:
    """
    Run electromagnetic simulation with given parameters.

    Returns:
        Dictionary with results including timing and solver info
    """

    # Extract parameters
    params = job_params['parameters']
    geometry_config = job_params['geometry']
    solver_config = job_params['solver']

    # 1. Create geometry and mesh
    print(f"\n[1/4] Generating surface mesh...")
    t0 = time.perf_counter()

    wavelength = params['wavelength']

    mesh,geometry_type = create_geometry(geometry_config, params, surface_mesh=True)

    # Redirect stdout to null
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')

    # Output
    filepath = job_path + f"mesh_scatterer.stl"
    mesh.Export(filepath, "STL Format")

    # Restore stdout
    sys.stdout = original_stdout
    sys.stderr = original_stderr

    mesh,geometry_type = create_geometry(geometry_config, params, surface_mesh=False)

    print(f"  ✓ Mesh generated ({geometry_type}): {mesh.ne} elements in {time.perf_counter() - t0:.2f}s", flush=True)

    # 2. Compute wavenumber from wavelength
    wavelength = params['wavelength']
    k = 2 * pi / wavelength

    # 3. Auto-detect problem type and create source field
    print(f"\n[2/4] Setting up source field...")
    t0 = time.perf_counter()

    # Auto-detect: if propagation_dir exists → scattering, if amplitude exists → antenna
    if 'propagation_dir' in params:
        # Scattering problem: incident wave
        source = create_incident_wave(
            k=k,
            propagation_dir=tuple(params['propagation_dir']),
            polarization=tuple(params['polarization'])
        )
        problem_type = 'scattering'
        print(f"  ✓ Incident wave (scattering) configured in {time.perf_counter() - t0:.3f}s", flush=True)
    elif 'amplitude' in params:
        # Antenna problem: antenna excitation
        source = create_antenna_source(
            polarization=tuple(params['polarization']),
            amplitude=params['amplitude']
        )
        problem_type = 'antenna'
        print(f"  ✓ Antenna source configured in {time.perf_counter() - t0:.3f}s", flush=True)
    else:
        raise ValueError(
            "Config must specify either 'propagation_dir' (scattering) or 'amplitude' (antenna)"
        )
    
    # Export incident wave (real and imaginary parts for complex field)
    if False:
        vtk_path = job_path + "data_E_inc"

        vtk = VTKOutput(
            ma=mesh,
            coefs=[source.real[0], source.real[1], source.real[2], source.imag[0], source.imag[1], source.imag[2]],
            names=["E_inc_real_x", "E_inc_real_y", "E_inc_real_z", "E_inc_imag_x", "E_inc_imag_y", "E_inc_imag_z"],
            filename=vtk_path,
            subdivision=2,
            legacy=True
        )
        vtk.Do()

        vtu = VTKOutput(
            ma=mesh,
            coefs=[source.real[0], source.real[1], source.real[2], source.imag[0], source.imag[1], source.imag[2]],
            names=["E_inc_real_x", "E_inc_real_y", "E_inc_real_z", "E_inc_imag_x", "E_inc_imag_y", "E_inc_imag_z"],
            filename=vtk_path,
            subdivision=2  # Smoother visualization in ParaView
        )
        vtu.Do()

    # 4. Setup unified Maxwell problem
    print(f"\n[3/4] Assembling FEM system ({problem_type} problem)...")
    t0 = time.perf_counter()
    problem = MaxwellProblem(
        mesh=mesh,
        k=k,
        source=source,
        fes_order=solver_config.get('fes_order', 5),
        use_type1=solver_config.get('use_type1', False)
    )
    problem.assemble_system()
    print(f"  ✓ System assembled in {time.perf_counter() - t0}s", flush=True)

    mu0 = 4*pi*1e-7
    sigma = mesh.MaterialCF({"inner":1e16}, default=0)
    mu = mesh.MaterialCF({"inner":0}, default=1)
    Js = mesh.MaterialCF({"inner":0}, default=0)
    ak = mesh.MaterialCF({"inner":0}, default=0)

    # Export material parameter (real and imaginary parts for complex field)
    vtk_path = job_path + "data_params"

    vtk = VTKOutput(
        ma=problem.fes.mesh,
        coefs=[sigma,mu,Js,ak],
        names=["conductivity", "relative_permeability", "saturation_level", "knee_shape"],
        filename=vtk_path,
        subdivision=2,
        legacy=True
    )
    vtk.Do()

    vtu = VTKOutput(
        ma=problem.fes.mesh,
        coefs=[sigma,mu,Js,ak],
        names=["conductivity", "relative_permeability", "saturation_level", "knee_shape"],
        filename=vtk_path,
        subdivision=2  # Smoother visualization in ParaView
    )
    vtu.Do()

    # rename E_field files to 'data-E_field' for consistency
    if os.path.exists(job_path + "E_field.vtk"):
        os.rename(job_path + "E_field.vtk", job_path + "data-E_field.vtk")
    if os.path.exists(job_path + "E_field.vtu"):
        os.rename(job_path + "E_field.vtu", job_path + "data-E_field.vtu")

    print(f"  ✓ Saved material parameters in {time.perf_counter() - t0}s", flush=True)
    print(f"\n[4/4] Making inputs / natural seeds list...")

    # Check results and meshes
    # mesh_E_inc = pv.read(job_path+"data-E_inc.vtu")
    mesh_E_fld = pv.read(job_path+"data_E_field.vtu")
    mesh_param = pv.read(job_path+"data_params.vtu")

    print(f"  Array names in params mesh: {mesh_param.array_names}")
    print(f"  Array names in E_field mesh: {mesh_E_fld.array_names}")

    assert mesh_E_fld.n_points == mesh_param.n_points, "✗ Mesh point counts do not match between E_inc, E_field, and parameters"
    assert mesh_E_fld.n_cells == mesh_param.n_cells, "✗ Mesh cell counts do not match between E_inc, E_field, and parameters"

    print(f"  ✓ Nbr. of points: {mesh_E_fld.n_points}, {mesh_param.n_points}")
    print(f"  ✓ Nbr. of cells: {mesh_E_fld.n_cells}, {mesh_param.n_cells}")

    # Collect inputs
    inputs = {
        'wavelength': wavelength,
        'wavenumber': float(k),
        'frequency': float(constants.c / wavelength),
        'amplitude': 1.0,
        'polarization': params['polarization'],
        'polarization_x': params['polarization'][0],
        'polarization_y': params['polarization'][1],
        'polarization_z': params['polarization'][2],
        'propagation_dir': params['propagation_dir'],
        'propagation_dir_x': params['propagation_dir'][0],
        'propagation_dir_y': params['propagation_dir'][1],
        'propagation_dir_z': params['propagation_dir'][2],
        'permeability_vacuum': 1.256637061e-6,
        'turns_coil': 0.0
    }

    metadata_file = job_path + "inputs.json"
    with open(metadata_file, 'w') as f:
        json.dump(inputs, f, indent=2)

    print(f"  ✓ Saved inputs list.", flush=True)

    return inputs, problem


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(description="Run electromagnetic simulation with specified parameters")

    parser.add_argument('--job-file',      type=str, required=True, help='JSON file containing job parameters')
    parser.add_argument('--job-id',        type=int, required=True, help='Job ID (index into job list)')
    parser.add_argument('--output-dir',    type=str, required=True, help='Output directory for results')
    parser.add_argument('--num-threads',   type=int, default=None,  help='Number of threads (default: from SLURM_CPUS_PER_TASK or 1)')
    parser.add_argument('--save-solution', action='store_true',     help='Save field solution (warning: can be large)')
    parser.add_argument('--local',         action='store_true',     help='Force local mode (skip SLURM detection)')

    args = parser.parse_args()

    # Determine execution mode
    is_cluster = is_slurm_environment() and not args.local

    # Determine number of threads
    if args.num_threads is not None:
        num_threads = args.num_threads
    else:
        num_threads = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    
    print("=" * 70)
    print(f"MaxwellScattererAndAntennaNGSolve Prepare Sample")
    print("=" * 70)
    print(f"Mode:         {'CLUSTER (SLURM)' if is_cluster else 'LOCAL'}")
    print(f"Job ID:       {args.job_id}")
    print(f"Job file:     {args.job_file}")
    print(f"Output dir:   {args.output_dir}")
    print(f"Threads:      {num_threads}")
    if is_cluster:
        print(f"SLURM Job:    {os.environ.get('SLURM_JOB_ID', 'N/A')}")
        print(f"Node:         {os.environ.get('SLURM_NODELIST', 'N/A')}")
    print("=" * 70)

    # Load job parameters
    try:
        job_file = Path(args.job_file)
        job_params = load_job_parameters(job_file, args.job_id)
        params = job_params['parameters']
        print(f"\n✓ Loaded parameters for job {args.job_id}")
        print(f"  Wavelength:   {params['wavelength']}")
        if 'propagation_dir' in params:
            print(f"  Type:         Scattering")
            print(f"  Propagation:  {params['propagation_dir']}")
        elif 'amplitude' in params:
            print(f"  Type:         Antenna")
            print(f"  Amplitude:    {params['amplitude']}")
        print(f"  Polarization: {params['polarization']}", flush=True)
    except Exception as e:
        print(f"✗ Error loading job parameters: {e}", file=sys.stderr)
        return 1

    # Run simulation
    try:
        # Set number of threads
        SetNumThreads(num_threads)

        job_dir = args.output_dir + f"job_{args.job_id:04d}/"

        print("\n" + "=" * 70)
        print("Starting tasks...")
        print("=" * 70)
        results, problem = run_task(job_params, job_dir)
        print("\n" + "=" * 70)
        print("✓ Tasks completed successfully!")
        print("=" * 70)
    except Exception as e:
        print(f"\n✗ Error during tasks: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    print("\n✓ Job completed successfully\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
