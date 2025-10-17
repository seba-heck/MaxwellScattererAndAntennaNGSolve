"""
Utility functions for electromagnetic scattering simulations.

This module provides geometry generation, incident wave definitions,
and helper functions for NGSolve-based electromagnetic simulations.
"""

from dataclasses import dataclass
from typing import Tuple, Optional
from ngsolve import CF, exp, pi, x, y, z
from netgen.occ import Sphere, Pnt, Glue, OCCGeometry
from ngsolve import Mesh, pml


@dataclass
class PhysicalParameters:
    """Container for physical parameters of electromagnetic problems."""

    wavelength: float  # Wavelength λ in meters
    outer_radius: float = 1.0  # Outer computational domain radius
    scatterer_radius: float = 0.1  # Inner scatterer radius
    pml_width: float = 0.25  # PML layer thickness
    max_mesh_size: float = 0.5  # Maximum mesh element size
    curve_order: int = 5  # Mesh curvature order

    @property
    def wavenumber(self) -> float:
        """Compute wavenumber k = 2π/λ."""
        return 2 * pi / self.wavelength


def create_spherical_geometry(
    R: float,
    PMLw: float,
    r: float,
    h_max: float,
    curve_order: int = 5
) -> Mesh:
    """
    Create a spherical scatterer geometry with PML boundary layer.

    The geometry consists of three concentric spheres:
    - Inner sphere (radius r): Scatterer surface
    - Middle sphere (radius R-PMLw): PML inner boundary
    - Outer sphere (radius R): Computational domain boundary

    Args:
        R: Outer radius of computational domain
        PMLw: Width of PML layer
        r: Radius of inner scatterer
        h_max: Maximum mesh element size
        curve_order: Order of mesh curvature (default: 5)

    Returns:
        NGSolve Mesh object with PML region configured

    Example:
        >>> mesh = create_spherical_geometry(R=1.0, PMLw=0.25, r=0.1, h_max=0.5)
    """
    # Create three concentric spheres
    sphere_R = Sphere(Pnt(0, 0, 0), R)          # Largest sphere (outer boundary)
    sphere_Rr = Sphere(Pnt(0, 0, 0), R - PMLw)  # Middle sphere (PML inner boundary)
    sphere_rr = Sphere(Pnt(0, 0, 0), r)         # Smallest sphere (scatterer)

    # Name the boundaries
    sphere_R.faces.name = "outer"
    sphere_rr.faces.name = "inner"

    # Create PML and vacuum regions
    PML = sphere_R - sphere_Rr
    vacuum = sphere_Rr - sphere_rr

    PML.name = "PML"
    vacuum.name = "vacuum"

    # Glue regions together
    domain = Glue([PML, vacuum])

    # Convert to OpenCascade geometry and generate mesh
    geo = OCCGeometry(domain, dim=3)
    mesh = Mesh(geo.GenerateMesh(maxh=h_max))
    mesh.Curve(curve_order)  # Curved elements to capture sphere curvature

    # Configure PML with radial attenuation
    mesh.SetPML(pml.Radial(rad=R - PMLw, alpha=1j, origin=(0, 0, 0)), "PML")

    return mesh


def create_incident_wave(
    k: float,
    propagation_dir: Tuple[float, float, float] = (0, 0, 1),
    polarization: Tuple[float, float, float] = (1, 0, 0)
) -> CF:
    """
    Create an incident plane wave field.

    The incident field is given by:
        E_inc = pol * exp(i * k * p · r)

    where:
        - pol is the polarization vector
        - p is the propagation direction (unit vector)
        - r is the position vector
        - k is the wavenumber

    Args:
        k: Wavenumber (2π/λ)
        propagation_dir: Propagation direction vector (px, py, pz)
        polarization: Polarization vector (Ex, Ey, Ez)

    Returns:
        NGSolve CoefficientFunction representing the incident field

    Example:
        >>> k = 2*pi/0.65
        >>> E_inc = create_incident_wave(k, prop_dir=(0,0,1), pol=(1,0,0))
    """
    p = CF(propagation_dir)      # Propagation direction (Poynting vector direction)
    pol = CF(polarization)       # Polarization vector
    coord_vec = CF((x, y, z))    # Position vector

    # Phase factor: exp(i k p·r)
    phase = exp(1j * k * p * coord_vec)

    # Incident field: E_inc = pol * phase
    E_inc = pol * phase

    return E_inc


def create_antenna_source(
    polarization: Tuple[float, float, float] = (1, 0, 0),
    amplitude: float = 1.0
) -> CF:
    """
    Create a simple uniform antenna excitation field.

    This creates a constant field at the antenna surface, suitable for
    basic antenna simulations. The field is uniform across the antenna.

    Args:
        polarization: Direction of excitation field (Ex, Ey, Ez)
        amplitude: Amplitude of the excitation field

    Returns:
        NGSolve CoefficientFunction representing the antenna source

    Example:
        >>> E_source = create_antenna_source(polarization=(1,0,0), amplitude=1.0)
    """
    pol = CF(polarization)
    return amplitude * pol


def create_zero_source() -> CF:
    """
    Create a zero volumetric source term.

    Returns:
        Zero vector coefficient function
    """
    return CF((0, 0, 0))


def normalize_vector(vec: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Normalize a 3D vector to unit length.

    Args:
        vec: Input vector (x, y, z)

    Returns:
        Normalized unit vector

    Example:
        >>> normalize_vector((3, 4, 0))
        (0.6, 0.8, 0.0)
    """
    magnitude = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
    if magnitude < 1e-12:
        raise ValueError("Cannot normalize zero vector")
    return (vec[0]/magnitude, vec[1]/magnitude, vec[2]/magnitude)


def compute_wavelength(frequency: float, c: float = 299792458.0) -> float:
    """
    Compute wavelength from frequency.

    Args:
        frequency: Frequency in Hz
        c: Speed of light in m/s (default: exact value)

    Returns:
        Wavelength in meters

    Example:
        >>> compute_wavelength(1e9)  # 1 GHz
        0.299792458
    """
    return c / frequency


def compute_frequency(wavelength: float, c: float = 299792458.0) -> float:
    """
    Compute frequency from wavelength.

    Args:
        wavelength: Wavelength in meters
        c: Speed of light in m/s (default: exact value)

    Returns:
        Frequency in Hz

    Example:
        >>> compute_frequency(0.65)
        461219166.15384614
    """
    return c / wavelength
