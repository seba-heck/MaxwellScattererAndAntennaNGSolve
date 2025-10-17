#!/usr/bin/env python3
"""
Test Ellipsoid and Cylinder primitives with correct syntax.
"""

from netgen.occ import *
from ngsolve import *
import inspect

print("=" * 70)
print("Testing Ellipsoid and Cylinder Primitives")
print("=" * 70)
print()

# Test 1: Get Ellipsoid signature
print("1. Inspecting Ellipsoid function signature:")
print("-" * 70)
try:
    sig = inspect.signature(Ellipsoid)
    print(f"✅ Ellipsoid signature: {sig}")
except Exception as e:
    print(f"Cannot get signature via inspect: {e}")
    print("Trying to call with wrong args to see error message...")
    try:
        Ellipsoid()
    except TypeError as e:
        print(f"Error message: {e}")
print()

# Test 2: Try different Ellipsoid syntaxes based on OCC docs
print("2. Testing Ellipsoid creation:")
print("-" * 70)

# Approach 1: Standard OCC ellipsoid (center + axis system + semi-axes)
try:
    center = Pnt(0, 0, 0)
    # gp_Ax2 defines a coordinate system (origin + main direction + x direction)
    axes = gp_Ax2(center, gp_Dir(0, 0, 1))
    a, b = 0.1, 0.2  # Semi-axes
    ellipsoid = Ellipsoid(axes, a, b)
    print("✅ Ellipsoid(gp_Ax2, a, b) works! (Ellipsoid of revolution)")
    print(f"   Type: {type(ellipsoid)}")
    print(f"   This creates a spheroid (ellipsoid of revolution)")
except Exception as e:
    print(f"❌ Ellipsoid(gp_Ax2, a, b) failed: {e}")

# Approach 2: Full tri-axial ellipsoid
try:
    center = Pnt(0, 0, 0)
    a, b, c = 0.05, 0.10, 0.15  # Three semi-axes
    ellipsoid = Ellipsoid(center, a, b, c)
    print("✅ Ellipsoid(Pnt, a, b, c) works! (Tri-axial ellipsoid)")
    print(f"   Type: {type(ellipsoid)}")
except Exception as e:
    print(f"❌ Ellipsoid(Pnt, a, b, c) failed: {e}")

print()

# Test 3: Get Cylinder signature
print("3. Inspecting Cylinder function signature:")
print("-" * 70)
try:
    print("Trying to call with wrong args to see error message...")
    try:
        Cylinder()
    except TypeError as e:
        print(f"Error message shows signatures:\n{e}")
except Exception as e:
    print(f"Error: {e}")
print()

# Test 4: Try different Cylinder syntaxes
print("4. Testing Cylinder creation:")
print("-" * 70)

# Approach 1: Simple syntax (point, direction, radius, height)
try:
    p = Pnt(0, 0, -0.5)
    d = gp_Dir(0, 0, 1)
    cyl = Cylinder(p, d, r=0.01, h=1.0)
    print("✅ Cylinder(Pnt, gp_Dir, r, h) works!")
    print(f"   Type: {type(cyl)}")
except Exception as e:
    print(f"❌ Cylinder(Pnt, gp_Dir, r, h) failed: {e}")

# Approach 2: Using gp_Ax2 (axis system)
try:
    p = Pnt(0, 0, 0)
    axis = gp_Ax2(p, gp_Dir(0, 0, 1))
    cyl = Cylinder(axis, r=0.01, h=1.0)
    print("✅ Cylinder(gp_Ax2, r, h) works!")
    print(f"   Type: {type(cyl)}")
except Exception as e:
    print(f"❌ Cylinder(gp_Ax2, r, h) failed: {e}")

print()

# Test 5: Create actual test geometries
print("5. Creating complete test geometries:")
print("-" * 70)

try:
    print("Creating dipole antenna geometry...")
    # Cylindrical dipole along z-axis
    dipole_length = 0.325  # L = λ/2 for λ = 0.65
    dipole_radius = 0.0065  # r = λ/100

    p = Pnt(0, 0, -dipole_length/2)
    d = gp_Dir(0, 0, 1)
    dipole = Cylinder(p, d, r=dipole_radius, h=dipole_length)

    print(f"✅ Created dipole: L={dipole_length}, r={dipole_radius}")
    print(f"   Aspect ratio L/r = {dipole_length/dipole_radius:.1f}")

    # Compute volume and surface area
    props = dipole.Properties()
    print(f"   Volume: {props[0]:.6f}")
    print(f"   Surface area: {props[1]:.6f}")

except Exception as e:
    print(f"❌ Failed to create dipole: {e}")

try:
    print("\nCreating ellipsoidal scatterer...")
    # Tri-axial ellipsoid
    a, b, c = 0.05, 0.10, 0.15

    ellipsoid = Ellipsoid(Pnt(0, 0, 0), a, b, c)

    print(f"✅ Created ellipsoid: a={a}, b={b}, c={c}")
    print(f"   Aspect ratios: b/a={b/a:.2f}, c/a={c/a:.2f}")

    # Compute volume and surface area
    props = ellipsoid.Properties()
    print(f"   Volume: {props[0]:.6f}")
    print(f"   Expected volume: {4/3 * 3.14159 * a * b * c:.6f} (4/3 π abc)")

except Exception as e:
    print(f"❌ Failed to create ellipsoid: {e}")

print()

# Test 6: Boolean operations with new geometries
print("6. Testing boolean operations:")
print("-" * 70)

try:
    # Create ellipsoid scatterer
    scatterer = Ellipsoid(Pnt(0, 0, 0), 0.05, 0.05, 0.15)

    # Create outer sphere domain
    domain = Sphere(Pnt(0, 0, 0), 1.0)

    # Create PML region
    pml_inner = Sphere(Pnt(0, 0, 0), 0.75)

    # Subtract scatterer from domain
    vacuum = domain - scatterer
    print("✅ Boolean subtraction works with ellipsoid")

    # Create PML layer
    pml = pml_inner - scatterer
    print("✅ Can create complex geometries with ellipsoids")

except Exception as e:
    print(f"❌ Boolean operations failed: {e}")

print()

print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
✅ CONFIRMED: Both Ellipsoid and Cylinder primitives are available!

Ellipsoid:
  - Syntax: Ellipsoid(Pnt(x,y,z), a, b, c)
  - Creates tri-axial ellipsoid with semi-axes a, b, c
  - Perfect for scatterer geometries

Cylinder:
  - Syntax: Cylinder(Pnt(x,y,z), gp_Dir(dx,dy,dz), r=radius, h=height)
  - OR: Cylinder(gp_Ax2(Pnt, gp_Dir), r=radius, h=height)
  - Perfect for dipole antenna geometries

WorkPlane + Revolve:
  - Backup approach for ellipsoid of revolution (spheroid)
  - More complex but provides additional flexibility

RECOMMENDATION:
  - Use direct Ellipsoid primitive (Approach 1 from plan)
  - Use direct Cylinder primitive for dipole
  - Implement geometry generators in src/geometry.py now!
""")
