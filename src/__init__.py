"""
MaxwellScattererAndAntennaNGSolve - Electromagnetic Antenna and Scattering Simulations

This package provides tools for solving electromagnetic scattering and antenna
problems using the NGSolve finite element library.

Main Components:
    - utils: Geometry generation and helper functions
    - scatterer: Scattering problem formulation
    - antenna: Antenna radiation problem formulation
    - solvers: Direct and iterative solver strategies

Example usage (Scattering):
    >>> from src import create_spherical_geometry, create_incident_wave
    >>> from src import ScattererProblem, solve_gmres
    >>> from ngsolve import TaskManager, pi
    >>>
    >>> mesh = create_spherical_geometry(R=1.0, PMLw=0.25, r=0.1, h_max=0.5)
    >>> k = 2*pi/0.65
    >>> E_inc = create_incident_wave(k, (0,0,1), (1,0,0))
    >>> problem = ScattererProblem(mesh, k, E_inc)
    >>> problem.assemble_system()
    >>> with TaskManager():
    >>>     solution = solve_gmres(problem.a, problem.l, problem.fes)

Example usage (Antenna):
    >>> from src import create_spherical_geometry, create_antenna_source
    >>> from src import AntennaProblem, solve_gmres
    >>> from ngsolve import TaskManager, pi
    >>>
    >>> mesh = create_spherical_geometry(R=1.0, PMLw=0.25, r=0.1, h_max=0.5)
    >>> k = 2*pi/0.65
    >>> E_source = create_antenna_source(polarization=(1,0,0), amplitude=1.0)
    >>> problem = AntennaProblem(mesh, k, E_source)
    >>> problem.assemble_system()
    >>> with TaskManager():
    >>>     solution = solve_gmres(problem.a, problem.l, problem.fes)
"""

# Version information
__version__ = "0.1.0"
__author__ = "Camilo & Maximilian"

# Import main components from submodules
from .utils import (
    create_incident_wave,
    create_antenna_source,
    create_zero_source,
    normalize_vector,
    compute_wavelength,
    compute_frequency,
    PhysicalParameters,
)

from .geometry import (
    create_ellipsoid_scatterer_geometry,
    create_spheroid_scatterer_geometry,
    create_dipole_antenna_geometry,
    create_spherical_geometry,

    create_box_scatterer_geometry,
    create_two_box_scatterer_geometry,
    create_two_ellipsoid_scatterer_geometry,
    create_cylinder_scatterer_geometry,
    create_empty_geometry,
)

from .scatterer import ScattererProblem

from .antenna import AntennaProblem

from .maxwell import MaxwellProblem

from .solvers import (
    solve_direct,
    solve_gmres,
    solve_cg,
    solve_bvp,
    solve_with_taskmanager,
    create_preconditioner,
    IterativeSolver,
)

# Define public API
__all__ = [
    # Version
    "__version__",
    "__author__",

    # Geometry generation
    "create_ellipsoid_scatterer_geometry",
    "create_spheroid_scatterer_geometry",
    "create_dipole_antenna_geometry",
    "create_spherical_geometry",

    "create_box_scatterer_geometry",
    "create_two_box_scatterer_geometry",
    "create_two_ellipsoid_scatterer_geometry",
    "create_cylinder_scatterer_geometry",
    "create_empty_geometry",

    # Utilities
    "create_incident_wave",
    "create_antenna_source",
    "create_zero_source",
    "normalize_vector",
    "compute_wavelength",
    "compute_frequency",
    "PhysicalParameters",

    # Problem formulations
    "ScattererProblem",
    "AntennaProblem",
    "MaxwellProblem",

    # Solvers
    "solve_direct",
    "solve_gmres",
    "solve_cg",
    "solve_bvp",
    "solve_with_taskmanager",
    "create_preconditioner",
    "IterativeSolver",
]
