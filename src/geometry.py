"""
Parametrizable geometry generation for antenna and scatterer problems.

This module provides functions to create various electromagnetic geometries:
- Ellipsoidal scatterers (tri-axial and spheroid)
- Cylindrical dipole antennas

All geometries include:
- Computational domain (sphere)
- PML (Perfectly Matched Layer) boundary region
- Proper material and boundary labeling
- Mesh size control

Physical conventions:
- Units: meters
- Coordinate system: Right-handed Cartesian (x, y, z)
- Default orientation: Major axis along z-axis
- Origin: Geometry center at (0, 0, 0)
"""

from netgen.occ import (
    Pnt, gp_Dir, Axes, Vec,
    Ellipsoid, Cylinder, Sphere, Box,
    Glue, OCCGeometry,Rectangle,Circle,Ellipse,WorkPlane
)
from netgen.meshing import MeshingParameters
from ngsolve import Mesh
import sys
from typing import Optional, Tuple
import numpy as np

def setup_basic_geometry(R: float, R_pml: float, scatterer, max_mesh_size: float):
    # Create spherical computational domain
    outer_sphere = Sphere(Pnt(0, 0, 0), R)
    pml_sphere = Sphere(Pnt(0, 0, 0), R_pml)

    # Create vacuum region (between scatterer and PML)
    vacuum_region = pml_sphere - scatterer
    vacuum_region.mat("vacuum")
    vacuum_region.maxh = max_mesh_size
    vacuum_region.faces.col = (0.7, 0.7, 1.0)  # Light blue

    # Create PML region (absorbing boundary layer)
    pml_region = outer_sphere - pml_sphere
    pml_region.mat("pml")
    pml_region.maxh = max_mesh_size
    pml_region.faces.col = (0.5, 0.5, 0.5)  # Gray

    # Name the outer boundary
    for face in outer_sphere.faces:
        face.name = "outer"
        face.maxh = max_mesh_size * 1.5  # Slightly coarser on outer boundary

    return vacuum_region, pml_region

def make_mesh_from_geometry(geometry, max_mesh_size: float, curve_order: int) -> Mesh:
    # Generate mesh
    occ_geo = OCCGeometry(geometry)
    ngmesh = occ_geo.GenerateMesh(maxh=max_mesh_size)
    ngmesh.SetGeometry(occ_geo)

    # Convert to NGSolve mesh
    mesh = Mesh(ngmesh)
    mesh.Curve(curve_order)  # Higher-order geometry approximation

    return mesh


def create_ellipsoid_scatterer_geometry(
    wavelength: float,
    semi_axis_a: float,
    semi_axis_b: float,
    semi_axis_c: float,
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    scatterer_mesh_size: Optional[float] = None,
    orientation: str = 'z',
    curve_order: int = 5
) -> Mesh:
    """
    Create ellipsoidal scatterer geometry with PML boundary conditions.

    The ellipsoid is defined by three semi-axes (a, b, c) corresponding to
    the x, y, z directions respectively (after applying orientation).

    Geometry structure:
    - Inner region: Ellipsoidal scatterer (perfect conductor)
    - Middle region: Vacuum (field propagation)
    - Outer region: PML layer (absorbing boundary)

    Args:
        wavelength: Wavelength in meters (for mesh sizing reference)
        semi_axis_a: Semi-axis length along first direction (meters)
        semi_axis_b: Semi-axis length along second direction (meters)
        semi_axis_c: Semi-axis length along third direction (meters)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Major axis orientation ('x', 'y', or 'z')
        curve_order: Mesh curve order (higher = better geometry approximation)

    Returns:
        NGSolve Mesh object with labeled regions and boundaries

    Example:
        >>> # Prolate spheroid (cigar shape, c > a = b)
        >>> mesh = create_ellipsoid_scatterer_geometry(
        ...     wavelength=0.65,
        ...     semi_axis_a=0.05,
        ...     semi_axis_b=0.05,
        ...     semi_axis_c=0.15
        ... )
        >>> # Oblate spheroid (pancake shape, c < a = b)
        >>> mesh = create_ellipsoid_scatterer_geometry(
        ...     wavelength=0.65,
        ...     semi_axis_a=0.15,
        ...     semi_axis_b=0.15,
        ...     semi_axis_c=0.05
        ... )
        >>> # Tri-axial ellipsoid (all different)
        >>> mesh = create_ellipsoid_scatterer_geometry(
        ...     wavelength=0.65,
        ...     semi_axis_a=0.05,
        ...     semi_axis_b=0.10,
        ...     semi_axis_c=0.15
        ... )
    """
    # Set default mesh size (10-15 elements per wavelength)
    if max_mesh_size is None:
        max_mesh_size = wavelength / 15.0

    # Validate parameters
    if pml_width >= domain_radius:
        raise ValueError(f"PML width ({pml_width}) must be less than domain radius ({domain_radius})")

    max_scatterer_dim = max(semi_axis_a, semi_axis_b, semi_axis_c)
    pml_inner_radius = domain_radius - pml_width

    if max_scatterer_dim >= pml_inner_radius:
        raise ValueError(
            f"Scatterer max dimension ({max_scatterer_dim:.3f}) must be less than "
            f"PML inner radius ({pml_inner_radius:.3f})"
        )

    # Create coordinate system based on orientation
    center = Pnt(0, 0, 0)
    if orientation == 'z':
        # Standard: z-axis is first semi-axis direction
        axes = Axes(center, n=gp_Dir(0, 0, 1), h=gp_Dir(1, 0, 0))
        # Semi-axes map: a→z, b→x, c→y
        r1, r2, r3 = semi_axis_a, semi_axis_b, semi_axis_c
    elif orientation == 'x':
        # x-axis is first semi-axis direction
        axes = Axes(center, n=gp_Dir(1, 0, 0), h=gp_Dir(0, 1, 0))
        # Semi-axes map: a→x, b→y, c→z
        r1, r2, r3 = semi_axis_a, semi_axis_b, semi_axis_c
    elif orientation == 'y':
        # y-axis is first semi-axis direction
        axes = Axes(center, n=gp_Dir(0, 1, 0), h=gp_Dir(0, 0, 1))
        # Semi-axes map: a→y, b→z, c→x
        r1, r2, r3 = semi_axis_a, semi_axis_b, semi_axis_c
    else:
        raise ValueError(f"orientation must be 'x', 'y', or 'z', got '{orientation}'")

    # Create ellipsoid scatterer
    scatterer = Ellipsoid(axes, r1, r2, r3)
    scatterer.faces.name = "inner"
    scatterer.faces.maxh = max_mesh_size if scatterer_mesh_size == None else scatterer_mesh_size
    scatterer.faces.col = (1, 0, 0)  # Red for scatterer

    # Setup basic geometry regions
    vacuum_region, pml_region = setup_basic_geometry(domain_radius, pml_inner_radius, scatterer, max_mesh_size)

    # Combine regions
    geometry = Glue([vacuum_region, pml_region])

    mesh = make_mesh_from_geometry(geometry, max_mesh_size, curve_order)

    return mesh


def create_box_scatterer_geometry(
    wavelength: float,
    axis_a: float,
    axis_b: float,
    axis_c: float,
    box_radius: float = None,
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    scatterer_mesh_size: Optional[float] = None,
    orientation: str = 'z',
    curve_order: int = 5
) -> Mesh:
    """
    Create ellipsoidal scatterer geometry with PML boundary conditions.

    The ellipsoid is defined by three semi-axes (a, b, c) corresponding to
    the x, y, z directions respectively (after applying orientation).

    Geometry structure:
    - Inner region: Ellipsoidal scatterer (perfect conductor)
    - Middle region: Vacuum (field propagation)
    - Outer region: PML layer (absorbing boundary)

    Args:
        wavelength: Wavelength in meters (for mesh sizing reference)
        semi_axis_a: Semi-axis length along first direction (meters)
        semi_axis_b: Semi-axis length along second direction (meters)
        semi_axis_c: Semi-axis length along third direction (meters)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Major axis orientation ('x', 'y', or 'z')
        curve_order: Mesh curve order (higher = better geometry approximation)

    Returns:
        NGSolve Mesh object with labeled regions and boundaries
    """
    # Set default mesh size (10-15 elements per wavelength)
    if max_mesh_size is None:
        max_mesh_size = wavelength / 15.0

    # Validate parameters
    if pml_width >= domain_radius:
        raise ValueError(f"PML width ({pml_width}) must be less than domain radius ({domain_radius})")

    max_scatterer_dim = max(axis_a, axis_b, axis_c)
    pml_inner_radius = domain_radius - pml_width

    if box_radius is None:
        box_radius = max_scatterer_dim / 4
    
    if max_scatterer_dim >= pml_inner_radius:
        raise ValueError(
            f"Scatterer max dimension ({max_scatterer_dim:.3f}) must be less than "
            f"PML inner radius ({pml_inner_radius:.3f})"
        )

    # Create coordinate system based on orientation
    # center = Pnt(0, 0, 0)
    # axes = Axes(center, n=gp_Dir(0, 0, 1), h=gp_Dir(1, 0, 0))

    # Create box scatterer
    if box_radius == 0.0:
        scatterer = Box(Pnt(-axis_a/2,-axis_b/2,-axis_c/2),
                        Pnt( axis_a/2, axis_b/2, axis_c/2))
    else:
        rect = Rectangle(axis_a,axis_b).Face()
        body = rect.Extrude(axis_c)
        scatterer = body.MakeFillet(body.edges, box_radius).Move((-axis_a/2,-axis_b/2,-axis_c/2))
    scatterer.faces.name = "inner"
    scatterer.faces.maxh = max_mesh_size if scatterer_mesh_size == None else scatterer_mesh_size
    scatterer.faces.col = (1, 0, 0)  # Red for scatterer

    # Setup basic geometry regions
    vacuum_region, pml_region = setup_basic_geometry(domain_radius, pml_inner_radius, scatterer, max_mesh_size)

    # Combine regions
    geometry = Glue([vacuum_region, pml_region])

    if box_radius == 0.0:
        edges = set([tuple(sorted([(e.start.x, e.start.y, e.start.z), (e.end.x, e.end.y, e.end.z)])) for e in scatterer.edges])
        n_points = 20
        mp = MeshingParameters(maxh=max_mesh_size)

        for p_start,p_end in edges:
            
            v = [_e - _s for _e,_s in zip(p_end,p_start)]
            steps = [_i/n_points for _i in range(n_points+1)]

            for _m in steps:
                x,y,z = [_p + _m*_v for _p,_v in zip(p_start,v)]
                mp.RestrictH(x, y, z, h=max_mesh_size/10)    

        # Generate mesh
        occ_geo = OCCGeometry(geometry)
        # ngmesh = occ_geo.GenerateMesh(mp=mp)
        ngmesh = occ_geo.GenerateMesh(maxh=max_mesh_size)
        ngmesh.SetGeometry(occ_geo)

        # Convert to NGSolve mesh
        mesh = Mesh(ngmesh)
        mesh.Curve(curve_order)  # Higher-order geometry approximation
    else:
        mesh = make_mesh_from_geometry(geometry, max_mesh_size, curve_order)

    return mesh

def create_cylinder_scatterer_geometry(
    wavelength: float,
    height: float,
    radius: float,
    box_radius: float,
    radius_2: float = None,
    origin: tuple = (0,0,0),
    direction: tuple = (0,0,1),
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    curve_order: int = 5
) -> Mesh:
    """
    Create ellipsoidal scatterer geometry with PML boundary conditions.

    The ellipsoid is defined by three semi-axes (a, b, c) corresponding to
    the x, y, z directions respectively (after applying orientation).

    Geometry structure:
    - Inner region: Ellipsoidal scatterer (perfect conductor)
    - Middle region: Vacuum (field propagation)
    - Outer region: PML layer (absorbing boundary)

    Args:
        wavelength: Wavelength in meters (for mesh sizing reference)
        semi_axis_a: Semi-axis length along first direction (meters)
        semi_axis_b: Semi-axis length along second direction (meters)
        semi_axis_c: Semi-axis length along third direction (meters)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Major axis orientation ('x', 'y', or 'z')
        curve_order: Mesh curve order (higher = better geometry approximation)

    Returns:
        NGSolve Mesh object with labeled regions and boundaries
    """
    # Set default mesh size (10-15 elements per wavelength)
    if max_mesh_size is None:
        max_mesh_size = wavelength / 15.0

    # Validate parameters
    if pml_width >= domain_radius:
        raise ValueError(f"PML width ({pml_width}) must be less than domain radius ({domain_radius})")

    # max_scatterer_dim = max(axis_a, axis_b, axis_c)
    pml_inner_radius = domain_radius - pml_width
    
    # if max_scatterer_dim >= pml_inner_radius:
    #     raise ValueError(
    #         f"Scatterer max dimension ({max_scatterer_dim:.3f}) must be less than "
    #         f"PML inner radius ({pml_inner_radius:.3f})"
    #     )

    # Create coordinate system based on orientation

    # Create box scatterer
    if box_radius == 0.0:
        origin_ = Pnt(origin[0],origin[1],origin[2]-height/2)
        scatterer = Cylinder(origin_, direction, r=radius, h=height)
    else:
        origin_ = Pnt(origin[0],origin[1])
        if radius_2 is None:
            rect = Circle(origin_, radius).Face()
        else:
            wpplate = WorkPlane(Axes())
            rect = wpplate.Ellipse(radius, radius_2).Face()
        body = rect.Extrude(height)
        # Try fillet with decreasing radii; if it fails, fall back to no fillet.
        scatterer_body = None
        try:
            if len(body.edges) == 0:
                raise RuntimeError("No edges to fillet")
            scatterer_body = body.MakeFillet(body.edges, box_radius)
        except Exception as e:
            # Retry with smaller radii
            fillet_radius = float(box_radius)
            for attempt in range(5):
                fillet_radius *= 0.5
                if fillet_radius <= 1e-6:
                    break
                try:
                    scatterer_body = body.MakeFillet(body.edges, fillet_radius)
                    break
                except Exception:
                    continue
            if scatterer_body is None:
                print(f"Warning: fillet failed ({e}); proceeding without fillet.", file=sys.stderr)
                scatterer_body = body
        scatterer = scatterer_body.Move((0,0,-height/2))
    scatterer.faces.name = "inner"

    scatterer.faces.maxh = max_mesh_size / 4
    scatterer.faces.col = (1, 0, 0)  # Red for scatterer

    # Setup basic geometry regions
    vacuum_region, pml_region = setup_basic_geometry(domain_radius, pml_inner_radius, scatterer, max_mesh_size)

    # Combine regions
    geometry = Glue([vacuum_region, pml_region])

    mesh = make_mesh_from_geometry(geometry, max_mesh_size, curve_order)

    return mesh

def create_two_box_scatterer_geometry(
    wavelength: float,
    dist: float,
    b1_axis_a: float,
    b1_axis_b: float,
    b1_axis_c: float,
    b2_axis_a: float,
    b2_axis_b: float,
    b2_axis_c: float,
    box_radius: float = None,
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    orientation: str = 'z',
    curve_order: int = 5
) -> Mesh:
    """
    Create ellipsoidal scatterer geometry with PML boundary conditions.

    The ellipsoid is defined by three semi-axes (a, b, c) corresponding to
    the x, y, z directions respectively (after applying orientation).

    Geometry structure:
    - Inner region: Ellipsoidal scatterer (perfect conductor)
    - Middle region: Vacuum (field propagation)
    - Outer region: PML layer (absorbing boundary)

    Args:
        wavelength: Wavelength in meters (for mesh sizing reference)
        semi_axis_a: Semi-axis length along first direction (meters)
        semi_axis_b: Semi-axis length along second direction (meters)
        semi_axis_c: Semi-axis length along third direction (meters)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Major axis orientation ('x', 'y', or 'z')
        curve_order: Mesh curve order (higher = better geometry approximation)

    Returns:
        NGSolve Mesh object with labeled regions and boundaries
    """
    # Set default mesh size (10-15 elements per wavelength)
    if max_mesh_size is None:
        max_mesh_size = wavelength / 15.0

    # Validate parameters
    if pml_width >= domain_radius:
        raise ValueError(f"PML width ({pml_width}) must be less than domain radius ({domain_radius})")

    max_scatterer_dim = max(b1_axis_a, b1_axis_b, b1_axis_c)
    pml_inner_radius = domain_radius - pml_width

    if box_radius is None:
        box_radius = max_scatterer_dim / 4
    
    if max_scatterer_dim >= pml_inner_radius:
        raise ValueError(
            f"Scatterer max dimension ({max_scatterer_dim:.3f}) must be less than "
            f"PML inner radius ({pml_inner_radius:.3f})"
        )

    # Create coordinate system based on orientation
    center = Pnt(0, 0, 0)
    axes = Axes(center, n=gp_Dir(0, 0, 1), h=gp_Dir(1, 0, 0))

    # Create box scatterer
    rect = Rectangle(b1_axis_a,b1_axis_b).Face()
    body = rect.Extrude(b1_axis_c)
    box1 = body.MakeFillet(body.edges, box_radius).Move((-b1_axis_a/2-dist/2,-b1_axis_b/2,-b1_axis_c/2))

    # Create second box
    rect2 = Rectangle(b2_axis_a,b2_axis_b).Face()
    body2 = rect2.Extrude(b2_axis_c)
    box2 = body2.MakeFillet(body2.edges, box_radius).Move((-b2_axis_a/2+dist/2,-b2_axis_b/2,-b2_axis_c/2))

    scatterer = box1 + box2
    scatterer.faces.name = "inner"
    scatterer.faces.maxh = max_mesh_size / 4
    scatterer.faces.col = (1, 0, 0)  # Red for scatterer

    # Setup basic geometry regions
    vacuum_region, pml_region = setup_basic_geometry(domain_radius, pml_inner_radius, scatterer, max_mesh_size)

    # Combine regions
    geometry = Glue([vacuum_region, pml_region])

    mesh = make_mesh_from_geometry(geometry, max_mesh_size, curve_order)

    return mesh

def create_two_ellipsoid_scatterer_geometry(
    wavelength: float,
    dist: float,
    b1_axis_a: float,
    b1_axis_b: float,
    b1_axis_c: float,
    b2_axis_a: float,
    b2_axis_b: float,
    b2_axis_c: float,
    box_radius: float = None,
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    orientation: str = 'z',
    curve_order: int = 5
) -> Mesh:
    """
    Create ellipsoidal scatterer geometry with PML boundary conditions.

    The ellipsoid is defined by three semi-axes (a, b, c) corresponding to
    the x, y, z directions respectively (after applying orientation).

    Geometry structure:
    - Inner region: Ellipsoidal scatterer (perfect conductor)
    - Middle region: Vacuum (field propagation)
    - Outer region: PML layer (absorbing boundary)

    Args:
        wavelength: Wavelength in meters (for mesh sizing reference)
        semi_axis_a: Semi-axis length along first direction (meters)
        semi_axis_b: Semi-axis length along second direction (meters)
        semi_axis_c: Semi-axis length along third direction (meters)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Major axis orientation ('x', 'y', or 'z')
        curve_order: Mesh curve order (higher = better geometry approximation)

    Returns:
        NGSolve Mesh object with labeled regions and boundaries
    """
    # Set default mesh size (10-15 elements per wavelength)
    if max_mesh_size is None:
        max_mesh_size = wavelength / 15.0

    # Validate parameters
    if pml_width >= domain_radius:
        raise ValueError(f"PML width ({pml_width}) must be less than domain radius ({domain_radius})")

    max_scatterer_dim = max(b1_axis_a, b1_axis_b, b1_axis_c)
    pml_inner_radius = domain_radius - pml_width

    if box_radius is None:
        box_radius = max_scatterer_dim / 4
    
    if max_scatterer_dim >= pml_inner_radius:
        raise ValueError(
            f"Scatterer max dimension ({max_scatterer_dim:.3f}) must be less than "
            f"PML inner radius ({pml_inner_radius:.3f})"
        )

    # Create coordinate system based on orientation
    center = Pnt(0, 0, 0)
    axes = Axes(center)
    # axes = Axes(center, n=gp_Dir(1, 0, 0), h=gp_Dir(0, 1, 0))

    # Create ellipsoid scatterer
    box1 = Ellipsoid(axes, b1_axis_a, b1_axis_b, b1_axis_c).Move((-b1_axis_b-dist/2,0,0))
    box2 = Ellipsoid(axes, b2_axis_a, b2_axis_b, b2_axis_c).Move((b2_axis_b+dist/2,0,0))
    scatterer = box1 + box2
    scatterer.faces.name = "inner"
    scatterer.faces.maxh = max_mesh_size
    scatterer.faces.col = (1, 0, 0)  # Red for scatterer

    # Setup basic geometry regions
    vacuum_region, pml_region = setup_basic_geometry(domain_radius, pml_inner_radius, scatterer, max_mesh_size)

    # Combine regions
    geometry = Glue([vacuum_region, pml_region])

    mesh = make_mesh_from_geometry(geometry, max_mesh_size, curve_order)

    return mesh

def create_spheroid_scatterer_geometry(
    wavelength: float,
    equatorial_radius: float,
    polar_radius: float,
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    orientation: str = 'z',
    curve_order: int = 5
) -> Mesh:
    """
    Create spheroidal (ellipsoid of revolution) scatterer geometry.

    A spheroid is an ellipsoid with two equal semi-axes (rotational symmetry).
    - Prolate: equatorial_radius < polar_radius (cigar/football shape)
    - Oblate: equatorial_radius > polar_radius (pancake/lens shape)

    Args:
        wavelength: Wavelength in meters
        equatorial_radius: Radius in the plane perpendicular to symmetry axis (a = b)
        polar_radius: Radius along symmetry axis (c)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Symmetry axis direction ('x', 'y', or 'z')
        curve_order: Mesh curve order

    Returns:
        NGSolve Mesh object

    Example:
        >>> # Prolate spheroid (elongated along z)
        >>> mesh = create_spheroid_scatterer_geometry(
        ...     wavelength=0.65,
        ...     equatorial_radius=0.05,
        ...     polar_radius=0.15
        ... )
    """
    return create_ellipsoid_scatterer_geometry(
        wavelength=wavelength,
        semi_axis_a=equatorial_radius,
        semi_axis_b=equatorial_radius,
        semi_axis_c=polar_radius,
        domain_radius=domain_radius,
        pml_width=pml_width,
        max_mesh_size=max_mesh_size,
        orientation=orientation,
        curve_order=curve_order
    )


def create_dipole_antenna_geometry(
    wavelength: float,
    length_factor: float = 0.5,
    radius_factor: float = 0.01,
    domain_radius: float = 1.0,
    pml_width: float = 0.25,
    max_mesh_size: Optional[float] = None,
    orientation: str = 'z',
    curve_order: int = 5
) -> Mesh:
    """
    Create cylindrical dipole antenna geometry with PML boundary conditions.

    The dipole is a thin cylindrical conductor, typically resonant at half-wavelength
    (length_factor = 0.5 gives a half-wave dipole).

    Geometry structure:
    - Inner region: Cylindrical dipole (antenna surface for excitation)
    - Middle region: Vacuum (field propagation)
    - Outer region: PML layer (absorbing boundary)

    Args:
        wavelength: Wavelength in meters
        length_factor: Dipole length as fraction of wavelength (L = factor × λ)
                      Common values:
                      - 0.5: Half-wave dipole (fundamental resonance, ~73Ω impedance)
                      - 0.25: Quarter-wave (with ground plane)
                      - 0.3-0.7: Typical parametric study range
        radius_factor: Dipole radius as fraction of wavelength (r = factor × λ)
                      Common values:
                      - 0.01: Thin wire approximation (L/r ~ 50)
                      - 0.005: Very thin wire (L/r ~ 100)
                      - 0.02: Fat dipole (increased bandwidth)
        domain_radius: Outer sphere radius (meters)
        pml_width: PML layer thickness (meters)
        max_mesh_size: Maximum element size (meters). Default: wavelength/15
        orientation: Dipole axis direction ('x', 'y', or 'z')
        curve_order: Mesh curve order

    Returns:
        NGSolve Mesh object with labeled regions and boundaries

    Example:
        >>> # Half-wave dipole at 461 MHz (λ = 0.65 m)
        >>> mesh = create_dipole_antenna_geometry(
        ...     wavelength=0.65,
        ...     length_factor=0.5,    # L = 0.325 m
        ...     radius_factor=0.01    # r = 0.0065 m (L/r = 50)
        ... )
    """
    # Set default mesh size
    if max_mesh_size is None:
        max_mesh_size = wavelength / 15.0

    # Calculate dipole dimensions
    length = length_factor * wavelength
    radius = radius_factor * wavelength

    # Validate parameters
    if radius >= length / 10:
        import warnings
        warnings.warn(
            f"Dipole radius ({radius:.4f}) is large compared to length ({length:.4f}). "
            f"L/r = {length/radius:.1f} < 10. Thin wire approximation may not be valid."
        )

    if pml_width >= domain_radius:
        raise ValueError(f"PML width ({pml_width}) must be less than domain radius ({domain_radius})")

    pml_inner_radius = domain_radius - pml_width
    if length / 2 >= pml_inner_radius:
        raise ValueError(
            f"Dipole half-length ({length/2:.3f}) must be less than "
            f"PML inner radius ({pml_inner_radius:.3f})"
        )

    # Create dipole cylinder centered at origin
    if orientation == 'z':
        # Cylinder along z-axis
        p_start = Pnt(0, 0, -length/2)
        direction = gp_Dir(0, 0, 1)
    elif orientation == 'x':
        # Cylinder along x-axis
        p_start = Pnt(-length/2, 0, 0)
        direction = gp_Dir(1, 0, 0)
    elif orientation == 'y':
        # Cylinder along y-axis
        p_start = Pnt(0, -length/2, 0)
        direction = gp_Dir(0, 1, 0)
    else:
        raise ValueError(f"orientation must be 'x', 'y', or 'z', got '{orientation}'")

    dipole = Cylinder(p_start, direction, r=radius, h=length)
    dipole.faces.name = "antenna"
    dipole.faces.maxh = max_mesh_size * 0.5  # Finer mesh on antenna surface
    dipole.faces.col = (0, 0, 1)  # Blue for antenna

    # Setup basic geometry regions
    vacuum_region, pml_region = setup_basic_geometry(domain_radius, pml_inner_radius, scatterer, max_mesh_size)

    # Combine regions
    geometry = Glue([vacuum_region, pml_region])

    mesh = make_mesh_from_geometry(geometry, max_mesh_size, curve_order)

    return mesh


# Convenience function for backward compatibility with existing code
def create_spherical_geometry(
    R: float = 1.0,
    PMLw: float = 0.25,
    r: float = 0.1,
    h_max: float = 0.5,
    curve_order: int = 5
) -> Mesh:
    """
    Create simple spherical scatterer geometry (backward compatibility).

    This maintains compatibility with the existing codebase that uses
    spherical geometries.

    Args:
        R: Outer domain radius
        PMLw: PML width
        r: Scatterer sphere radius
        h_max: Maximum mesh element size
        curve_order: Mesh curve order

    Returns:
        NGSolve Mesh object

    Note:
        This is a simplified interface. For new code, consider using
        create_ellipsoid_scatterer_geometry or create_dipole_antenna_geometry
        with more explicit parameters.
    """
    # Create as spheroid with equal radii (sphere)
    # Use a dummy wavelength for mesh sizing
    wavelength = h_max * 15.0  # Reverse of default mesh sizing

    return create_spheroid_scatterer_geometry(
        wavelength=wavelength,
        equatorial_radius=r,
        polar_radius=r,
        domain_radius=R,
        pml_width=PMLw,
        max_mesh_size=h_max,
        orientation='z',
        curve_order=curve_order
    )
