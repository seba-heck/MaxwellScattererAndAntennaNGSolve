#!/usr/bin/env python3
"""
Test correct OCC syntax for geometry creation.
"""

from netgen.occ import *
from ngsolve import *

print("=" * 70)
print("Testing Correct OCC Syntax")
print("=" * 70)
print()

# Test 1: Find available geometry constructors
print("1. Looking for geometry primitives...")
print("-" * 70)

# Check what's actually available in the module
import netgen.occ as occ
all_names = dir(occ)

# Look for primitive shapes
primitives = []
for name in all_names:
    if name[0].isupper() and not name.startswith('_'):
        obj = getattr(occ, name)
        if callable(obj):
            primitives.append(name)

print(f"Found {len(primitives)} callable objects (potential primitives):")
for p in sorted(primitives):
    print(f"  - {p}")
print()

# Test 2: Try Cylinder with correct syntax
print("2. Testing Cylinder creation with different approaches:")
print("-" * 70)

# Approach 1: Using Axis constructor properly
try:
    p = Pnt(0, 0, -0.5)
    d = gp_Dir(0, 0, 1)
    axis = Axis(p, d)
    cyl = Cylinder(axis, r=0.01, h=1.0)
    print("✅ Cylinder(Axis(Pnt, gp_Dir), r, h) works!")
    print(f"   Type: {type(cyl)}")
except Exception as e:
    print(f"❌ Failed: {e}")

# Approach 2: Using WorkPlane + Extrude
try:
    wp = WorkPlane()
    wp.MoveTo(0, 0).Circle(0.01).Close()
    face = wp.Face()
    cyl2 = face.Extrude(1.0 * gp_Vec(0, 0, 1))
    print("✅ WorkPlane().Circle().Face().Extrude() works!")
    print(f"   Type: {type(cyl2)}")
except Exception as e:
    print(f"❌ Failed: {e}")

print()

# Test 3: Try creating ellipsoid using WorkPlane + Revolve
print("3. Testing ellipsoid via WorkPlane + Revolve:")
print("-" * 70)

try:
    # Create an ellipse profile in 2D
    wp = WorkPlane()
    # Create ellipse by parametric approach
    # Not all WorkPlane methods support ellipses directly, so let's try

    # Try to see if Ellipse method exists
    if hasattr(wp, 'Ellipse'):
        print("✅ WorkPlane has Ellipse method!")
    else:
        print("❌ WorkPlane doesn't have Ellipse method")

    # Try to create using spline approximation
    import numpy as np
    theta = np.linspace(0, 2*np.pi, 50)
    a, c = 0.1, 0.2  # semi-axes for ellipse in x-z plane

    # Start at first point
    wp.MoveTo(a, 0)
    for i in range(1, len(theta)):
        x = a * np.cos(theta[i])
        z = c * np.sin(theta[i])
        wp.LineTo(x, z)
    wp.Close()

    profile = wp.Face()
    print("✅ Created ellipse profile via spline approximation")

    # Now revolve around z-axis
    ellipsoid = profile.Revolve(Axis(Pnt(0,0,0), gp_Dir(0,0,1)), 360)
    print("✅ Revolved ellipse profile to create ellipsoid of revolution!")
    print(f"   Type: {type(ellipsoid)}")

except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()

print()

# Test 4: Test Box and Scale
print("4. Testing Box creation and transformation:")
print("-" * 70)

try:
    # Create unit box
    unit_box = Box(Pnt(-1, -1, -1), Pnt(1, 1, 1))
    print("✅ Created unit box")

    # Scale uniformly
    scaled_box = unit_box.Scale(Pnt(0, 0, 0), 2.0)
    print("✅ Scaled box uniformly by 2.0")

    # Try to create ellipsoid-like shape by scaling box then... no that won't work

    # What about creating a sphere and trying transformations?
    unit_sphere = Sphere(Pnt(0, 0, 0), 1.0)

    # Check if there's a transformation matrix we can apply
    print(f"\nChecking transformation methods on sphere:")
    transform_methods = [m for m in dir(unit_sphere) if 'trsf' in m.lower() or 'trans' in m.lower()]
    print(f"  Found: {transform_methods}")

except Exception as e:
    print(f"❌ Failed: {e}")

print()

# Test 5: Check the Located method (might allow transformations)
print("5. Testing Located method (for transformations):")
print("-" * 70)

try:
    sphere = Sphere(Pnt(0, 0, 0), 1.0)

    # Try creating a transformation
    trsf = gp_Trsf()
    print(f"✅ Created gp_Trsf object")
    print(f"   Available methods: {[m for m in dir(trsf) if not m.startswith('_')]}")

    # Try to set translation
    trsf.SetTranslation(gp_Vec(1, 0, 0))
    transformed = sphere.Located(trsf)
    print(f"✅ sphere.Located(gp_Trsf) works for translation!")

except Exception as e:
    print(f"❌ Failed: {e}")

print()

# Test 6: Check for affine transformations with gp_GTrsf
print("6. Testing gp_GTrsf (general transformation - non-uniform scale):")
print("-" * 70)

try:
    sphere = Sphere(Pnt(0, 0, 0), 1.0)

    # Create general transformation (allows non-uniform scaling)
    gtrsf = gp_GTrsf()
    print(f"✅ Created gp_GTrsf object")
    print(f"   Methods: {[m for m in dir(gtrsf) if 'Scale' in m or 'Matrix' in m]}")

    # Try to set non-uniform scaling
    # gp_GTrsf can represent affine transformations including non-uniform scaling
    # But need to check if we can apply it to shapes

    # Check if shape has method to apply gp_GTrsf
    if hasattr(sphere, 'Located'):
        print("   Shape has Located method")
    if hasattr(sphere, 'Transformed'):
        print("   Shape has Transformed method")

except Exception as e:
    print(f"❌ Failed: {e}")

print()

print("=" * 70)
print("FINDINGS:")
print("=" * 70)
print("""
Based on tests:
1. Cylinder: Requires Axis(Pnt, gp_Dir) syntax
2. No direct Ellipsoid primitive
3. Scale only works uniformly (one scale factor)
4. WorkPlane + Revolve can create ellipsoid of revolution (spheroid)
5. gp_Trsf exists but only for rigid transformations
6. gp_GTrsf exists (general transformation) - might support non-uniform scale

Next: Test if gp_GTrsf can be applied to create general ellipsoids.
""")
