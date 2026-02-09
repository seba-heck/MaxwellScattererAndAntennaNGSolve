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

# Add parent directory to path to import src modules
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from ngsolve import TaskManager, SetNumThreads, pi
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


def create_geometry(geometry_config: Dict[str, Any], params: Dict[str, Any], verbose: bool = True):
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
            curve_order=geometry_config.get('curve_order', 5)
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
            curve_order=geometry_config.get('curve_order', 5)
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
            curve_order=geometry_config.get('curve_order', 5)
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

def run_simulation(job_params: Dict[str, Any], num_threads: int) -> Dict[str, Any]:
    """
    Run electromagnetic simulation with given parameters.

    Returns:
        Dictionary with results including timing and solver info
    """
    # Set number of threads
    SetNumThreads(num_threads)

    # Extract parameters
    params = job_params['parameters']
    geometry_config = job_params['geometry']
    solver_config = job_params['solver']

    # Timing breakdown
    timings = {}

    # 1. Create geometry and mesh (auto-detect geometry type)
    print(f"\n[1/4] Generating mesh...")
    t0 = time.perf_counter()

    # Auto-detect geometry type from parameters
    wavelength = params['wavelength']

    mesh,geometry_type = create_geometry(geometry_config, params)

    # if geometry_config["type"] == "ellipsoid" or 'ellipsoid_semi_axis_a' in params:
    #     # Tri-axial ellipsoid scatterer
    #     print(f"  → Creating tri-axial ellipsoid scatterer geometry")
    #     mesh = create_ellipsoid_scatterer_geometry(
    #         wavelength=wavelength,
    #         semi_axis_a=params['ellipsoid_semi_axis_a'],
    #         semi_axis_b=params['ellipsoid_semi_axis_b'],
    #         semi_axis_c=params['ellipsoid_semi_axis_c'],
    #         domain_radius=geometry_config.get('R', 1.0),
    #         pml_width=geometry_config.get('PMLw', 0.25),
    #         max_mesh_size=geometry_config.get('h_max', None),
    #         orientation=params.get('ellipsoid_orientation', 'z'),
    #         curve_order=geometry_config.get('curve_order', 5)
    #     )
    #     geometry_type = 'ellipsoid'

    # elif geometry_config["type"] == "box":
    #     # Box geometry
    #     print(f"  → Creating box geometry")
    #     mesh = create_box_scatterer_geometry(
    #         wavelength=wavelength,
    #         axis_a=params['axis_a'],
    #         axis_b=params['axis_b'],
    #         axis_c=params['axis_c'],
    #         box_radius=params['edge_radius'],
    #         domain_radius=geometry_config.get('R', 1.0),
    #         pml_width=geometry_config.get('PMLw', 0.25),
    #         max_mesh_size=geometry_config.get('h_max', None),
    #         curve_order=geometry_config.get('curve_order', 5)
    #     )
    #     geometry_type = 'box'

    # elif geometry_config["type"] == "cylinder":
    #     # Cylinder geometry
    #     print(f"  → Creating cylinder geometry")
    #     mesh = create_cylinder_scatterer_geometry(
    #         wavelength=wavelength,
    #         height=params['height'],
    #         radius=params['radius_major'],
    #         radius_2=params.get('radius_minor', None),
    #         box_radius=params['edge_radius'],
    #         domain_radius=geometry_config.get('R', 1.0),
    #         pml_width=geometry_config.get('PMLw', 0.25),
    #         max_mesh_size=geometry_config.get('h_max', None),
    #         curve_order=geometry_config.get('curve_order', 5)
    #     )
    #     geometry_type = 'cylinder'

    # elif 'dipole_length_factor' in params:
    #     # Dipole antenna geometry
    #     print(f"  → Creating cylindrical dipole antenna geometry")
    #     mesh = create_dipole_antenna_geometry(
    #         wavelength=wavelength,
    #         length_factor=params['dipole_length_factor'],
    #         radius_factor=params.get('dipole_radius_factor', 0.01),
    #         domain_radius=geometry_config.get('R', 1.0),
    #         pml_width=geometry_config.get('PMLw', 0.25),
    #         max_mesh_size=geometry_config.get('h_max', None),
    #         orientation=params.get('dipole_orientation', 'z'),
    #         curve_order=geometry_config.get('curve_order', 5)
    #     )
    #     geometry_type = 'dipole'

    # elif 'spheroid_equatorial_radius' in params:
    #     # Spheroid scatterer
    #     print(f"  → Creating spheroid scatterer geometry")
    #     mesh = create_spheroid_scatterer_geometry(
    #         wavelength=wavelength,
    #         equatorial_radius=params['spheroid_equatorial_radius'],
    #         polar_radius=params['spheroid_polar_radius'],
    #         domain_radius=geometry_config.get('R', 1.0),
    #         pml_width=geometry_config.get('PMLw', 0.25),
    #         max_mesh_size=geometry_config.get('h_max', None),
    #         orientation=params.get('spheroid_orientation', 'z'),
    #         curve_order=geometry_config.get('curve_order', 5)
    #     )
    #     geometry_type = 'spheroid'

    # else:
    #     # Default: spherical geometry (backward compatibility)
    #     print(f"  → Creating spherical geometry (default)")
    #     mesh = create_spherical_geometry(
    #         R=geometry_config.get('R', 1.0),
    #         PMLw=geometry_config.get('PMLw', 0.25),
    #         r=geometry_config.get('r', 0.1),
    #         h_max=geometry_config.get('h_max', 0.5),
    #         curve_order=geometry_config.get('curve_order', 5)
    #     )
    #     geometry_type = 'sphere'

    timings['mesh_gen'] = time.perf_counter() - t0
    print(f"  ✓ Mesh generated ({geometry_type}): {mesh.ne} elements in {timings['mesh_gen']:.2f}s")

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
        print(f"  ✓ Incident wave (scattering) configured in {time.perf_counter() - t0:.3f}s")
    elif 'amplitude' in params:
        # Antenna problem: antenna excitation
        source = create_antenna_source(
            polarization=tuple(params['polarization']),
            amplitude=params['amplitude']
        )
        problem_type = 'antenna'
        print(f"  ✓ Antenna source configured in {time.perf_counter() - t0:.3f}s")
    else:
        raise ValueError(
            "Config must specify either 'propagation_dir' (scattering) or 'amplitude' (antenna)"
        )

    timings['source_setup'] = time.perf_counter() - t0

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
    timings['assembly'] = time.perf_counter() - t0
    print(f"  ✓ System assembled in {timings['assembly']:.2f}s")

    # 5. Solve
    print(f"\n[4/4] Solving linear system...")
    t0 = time.perf_counter()
    solver_method = solver_config.get('method', 'gmres')

    with TaskManager():
        if solver_method == 'gmres':
            solution = solve_gmres(
                problem.a, problem.l, problem.fes,
                preconditioner=solver_config.get('preconditioner', 'block_jacobi'),
                maxsteps=solver_config.get('maxiter', 1000),
                tol=solver_config.get('tol', 1e-6),
                printrates=True
            )
        elif solver_method == 'direct':
            solution = solve_direct(
                problem.a, problem.l, problem.fes,
                print_info=True
            )
        else:
            raise ValueError(f"Unknown solver method: {solver_method}")

    timings['solve'] = time.perf_counter() - t0
    timings['total'] = sum(timings.values())
    print(f"  ✓ Solved in {timings['solve']:.2f}s")

    # Store solution in problem
    problem.set_solution(solution)

    # Collect results
    results = {
        'job_id': job_params['job_id'],
        'problem_type': problem_type,
        'geometry_type': geometry_type,
        'parameters': {
            'wavelength': wavelength,
            'wavenumber': float(k),
            'polarization': params['polarization'],
            'geometry': geometry_config,
        },
        'solver': {
            'method': solver_method,
            'preconditioner': solver_config.get('preconditioner', 'N/A'),
            'fes_order': solver_config.get('fes_order', 5),
        },
        'mesh': {
            'elements': mesh.ne,
            'vertices': mesh.nv,
            'ndof': problem.fes.ndof,
            'free_dofs': sum(problem.fes.FreeDofs()),
        },
        'timings': timings,
        'resources': {
            'cpus': num_threads,
        },
        'environment': {
            'is_slurm': is_slurm_environment(),
            'slurm_job_id': os.environ.get('SLURM_JOB_ID', 'N/A'),
        }
    }

    return results, problem


def save_results(results: Dict[str, Any], output_dir: Path, job_id: int,
                 params: Dict[str, Any],
                save_solution: bool = False, problem=None):
    """Save simulation results to output directory."""
    # Create job-specific directory
    job_dir = output_dir / f"job_{job_id:04d}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # Optionally save VTK solution (can be large!)
    vtk_files = None
    if save_solution and problem is not None:
        try:
            from ngsolve import VTKOutput

            # Get the solution GridFunction
            gfu = problem.get_solution()

            # Export E-field (real and imaginary parts for complex field)
            vtk_path = str(job_dir / "E_field")

            vtk = VTKOutput(
                ma=problem.fes.mesh,
                coefs=[gfu.real, gfu.imag],
                names=["E_real", "E_imag"],
                filename=vtk_path,
                subdivision=2  # Smoother visualization in ParaView
            )
            vtk.Do()

            # NGSolve creates files like: E_field.vtu or E_field_000000.vtu
            # Record which files were created
            vtk_files = ["E_field.vtu"]
            print(f"  ✓ VTK output saved: E_field.vtu")

        except Exception as e:
            print(f"  ⚠ Warning: Could not save VTK output: {e}")
            vtk_files = None

    # Add VTK info to results
    results['vtk_output'] = vtk_files

    # Save metadata
    metadata_file = job_dir / "metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(results, f, indent=2)

    # Save timings separately for easy aggregation
    timings_file = job_dir / "timings.json"
    with open(timings_file, 'w') as f:
        json.dump(results['timings'], f, indent=2)

    # Draw simulation results
    draw_file = job_dir / "field_visualization.html"
    clipping = { "function" : True,  "pnt" : (0,0.0,0), "vec" : (0,0,-1) }

    line_1 = { "type": "lines", "position": [params['wavelength'],params['wavelength'],0, params['wavelength']+0.5*params['propagation_dir'][0], params['wavelength']+0.5*params['propagation_dir'][1], 0.5*params['propagation_dir'][2]], "name": "propagation direction", "color": "red",}
    line_2 = { "type": "lines", "position": [params['wavelength'],params['wavelength'],0, params['wavelength']+0.5*params['polarization'][0], params['wavelength']+0.5*params['polarization'][1], 0.5*params['polarization'][2]], "name": "polarization direction", "color": "blue"}
    points = { "type": "points", "position": [params['wavelength'],params['wavelength'],0], "size":20, "color": "black", "name": "origin"}
    text_1 = { "type": "text", "name": "info1", "text": f" wavelength = {params['wavelength']}, outer radius = {params['R']}, PML width = {params['PMLw']}, mesh size = {params['h_max']}", "position": [-params['wavelength'],-params['wavelength']-0.2,0]}
    text_2 = { "type": "text", "name": "info2", "text": f" elements = {problem.fes.mesh.ne}, vertices = {problem.fes.mesh.nv}, free DOFs = {sum(problem.fes.FreeDofs())}", "position": [-params['wavelength'],-params['wavelength']-0.3,0]}
    
    Draw(problem.solution, problem.fes.mesh, "B", objects=[line_1,line_2,points,text_1,text_2], clipping=clipping, max = 10e-3, min = 0, draw_surf=False, filename=draw_file)
    
    print(f"\n✓ Results saved to: {job_dir}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Run electromagnetic simulation with specified parameters"
    )
    parser.add_argument(
        '--job-file',
        type=str,
        required=True,
        help='JSON file containing job parameters'
    )
    parser.add_argument(
        '--job-id',
        type=int,
        required=True,
        help='Job ID (index into job list)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory for results'
    )
    parser.add_argument(
        '--num-threads',
        type=int,
        default=None,
        help='Number of threads (default: from SLURM_CPUS_PER_TASK or 1)'
    )
    parser.add_argument(
        '--save-solution',
        action='store_true',
        help='Save field solution (warning: can be large)'
    )
    parser.add_argument(
        '--local',
        action='store_true',
        help='Force local mode (skip SLURM detection)'
    )

    args = parser.parse_args()

    # Determine execution mode
    is_cluster = is_slurm_environment() and not args.local

    # Determine number of threads
    if args.num_threads is not None:
        num_threads = args.num_threads
    else:
        num_threads = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))

    print("=" * 70)
    print(f"MaxwellScattererAndAntennaNGSolve Simulation")
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
        print(f"  Polarization: {params['polarization']}")
    except Exception as e:
        print(f"Error loading job parameters: {e}", file=sys.stderr)
        return 1

    # Run simulation
    try:
        print("\n" + "=" * 70)
        print("Starting simulation...")
        print("=" * 70)
        results, problem = run_simulation(job_params, num_threads)
        print("\n" + "=" * 70)
        print("✓ Simulation completed successfully!")
        print("=" * 70)
        print(f"Total time: {results['timings']['total']:.2f} s")
        print(f"  Mesh:     {results['timings']['mesh_gen']:.2f} s")
        print(f"  Assembly: {results['timings']['assembly']:.2f} s")
        print(f"  Solve:    {results['timings']['solve']:.2f} s")
    except Exception as e:
        print(f"\n✗ Error during simulation: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    # Save results
    try:
        output_dir = Path(args.output_dir)
        save_solution = args.save_solution or job_params['output'].get('save_solution', False)
        save_results(results, output_dir, args.job_id, params, save_solution, problem)
    except Exception as e:
        print(f"Error saving results: {e}", file=sys.stderr)
        return 1

    print("\n✓ Job completed successfully\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
