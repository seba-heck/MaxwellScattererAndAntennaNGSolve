#!/usr/bin/env python3
"""
Test script to verify NGSolve OCC capabilities for parametric geometries.
Run this BEFORE implementing geometry generators.

This script tests:
1. Available classes in netgen.occ
2. Ellipsoid primitive (if available)
3. Scale transformation (for creating ellipsoids from spheres)
4. Cylinder primitive (for dipole antennas)
5. WorkPlane and revolution capabilities
6. Available methods on Shape objects
"""

import sys
import inspect
from ngsolve import *

print("=" * 70)
print("NGSolve OCC Capability Verification")
print("=" * 70)
print()

# Import netgen.occ
try:
    from netgen.occ import *
    import netgen.occ as occ
    print("✅ netgen.occ module imported successfully")
    print()
except ImportError as e:
    print(f"❌ Failed to import netgen.occ: {e}")
    sys.exit(1)

# Test 1: List all available classes in netgen.occ
print("-" * 70)
print("1. Available netgen.occ classes:")
print("-" * 70)
occ_classes = [name for name, obj in inspect.getmembers(occ)
               if inspect.isclass(obj)]
print(f"Found {len(occ_classes)} classes:")
for cls in sorted(occ_classes):
    print(f"  - {cls}")
print()

# Test 2: Check for Ellipsoid primitive
print("-" * 70)
print("2. Testing Ellipsoid primitive:")
print("-" * 70)
try:
    # Try different possible signatures
    ellipsoid = None

    # Attempt 1: Direct constructor
    try:
        ellipsoid = Ellipsoid(Pnt(0, 0, 0), 1.0, 0.5, 0.3)
        print("✅ Ellipsoid(Pnt, a, b, c) constructor works!")
    except (NameError, AttributeError, TypeError) as e1:
        print(f"❌ Ellipsoid(Pnt, a, b, c) failed: {type(e1).__name__}")

        # Attempt 2: Different signature
        try:
            ellipsoid = Ellipsoid(Pnt(0, 0, 0), Pnt(1, 0, 0), Pnt(0, 0.5, 0), Pnt(0, 0, 0.3))
            print("✅ Ellipsoid with point-based constructor works!")
        except (NameError, AttributeError, TypeError) as e2:
            print(f"❌ Ellipsoid point-based constructor failed: {type(e2).__name__}")

    if ellipsoid is not None:
        print(f"   Ellipsoid type: {type(ellipsoid)}")
except Exception as e:
    print(f"❌ Ellipsoid primitive NOT AVAILABLE: {e}")
print()

# Test 3: Check for Scale method on Shape objects
print("-" * 70)
print("3. Testing Scale transformation:")
print("-" * 70)
try:
    sphere = Sphere(Pnt(0, 0, 0), 1.0)
    print("✅ Created test sphere")

    # Try different scale signatures
    scaled = None

    # Attempt 1: Non-uniform scaling with tuple
    try:
        scaled = sphere.Scale(Pnt(0, 0, 0), (2.0, 1.0, 0.5))
        print("✅ sphere.Scale(Pnt, (sx, sy, sz)) works!")
        print(f"   Scaled object type: {type(scaled)}")
    except (AttributeError, TypeError) as e1:
        print(f"❌ sphere.Scale(Pnt, tuple) failed: {type(e1).__name__}: {e1}")

        # Attempt 2: Uniform scaling
        try:
            scaled = sphere.Scale(Pnt(0, 0, 0), 2.0)
            print("⚠️  sphere.Scale(Pnt, factor) works (uniform only)")
            print(f"   Scaled object type: {type(scaled)}")
        except (AttributeError, TypeError) as e2:
            print(f"❌ sphere.Scale(Pnt, factor) also failed: {type(e2).__name__}")

except Exception as e:
    print(f"❌ Scale transformation NOT AVAILABLE: {e}")
print()

# Test 4: Cylinder creation (should work based on docs)
print("-" * 70)
print("4. Testing Cylinder primitive (for dipole antenna):")
print("-" * 70)
try:
    # Create cylinder along z-axis
    axis_start = Pnt(0, 0, -0.5)
    axis_direction = Axis(Pnt(0, 0, 0), direction=(0, 0, 1))
    cyl = Cylinder(axis_start, axis_direction, r=0.01, h=1.0)
    print("✅ Cylinder creation successful!")
    print(f"   Cylinder type: {type(cyl)}")
    print("   Parameters: r=0.01, h=1.0, along z-axis")
except Exception as e:
    print(f"❌ Cylinder creation failed: {e}")
print()

# Test 5: List methods on Sphere object
print("-" * 70)
print("5. Available methods on Sphere/Shape object:")
print("-" * 70)
try:
    sphere = Sphere(Pnt(0, 0, 0), 1.0)
    sphere_methods = [m for m in dir(sphere) if not m.startswith('_')]
    print(f"Found {len(sphere_methods)} public methods:")
    for method in sorted(sphere_methods):
        print(f"  - {method}")
except Exception as e:
    print(f"❌ Could not inspect Sphere methods: {e}")
print()

# Test 6: Check for WorkPlane and revolution
print("-" * 70)
print("6. Testing WorkPlane for profile revolution:")
print("-" * 70)
try:
    wp = WorkPlane()
    print("✅ WorkPlane exists")

    # Try to inspect WorkPlane methods
    wp_methods = [m for m in dir(wp) if not m.startswith('_')]
    print(f"   WorkPlane has {len(wp_methods)} public methods:")
    for method in sorted(wp_methods):
        print(f"     - {method}")

except (NameError, AttributeError) as e:
    print(f"❌ WorkPlane NOT AVAILABLE: {type(e).__name__}")
print()

# Test 7: Test transformation methods (Move, Rotate, Mirror)
print("-" * 70)
print("7. Testing transformation methods:")
print("-" * 70)
sphere = Sphere(Pnt(0, 0, 0), 0.1)

# Test Move
try:
    moved = sphere.Move(Vec(1, 0, 0))
    print("✅ sphere.Move(Vec) works")
except Exception as e:
    print(f"❌ sphere.Move failed: {type(e).__name__}")

# Test Rotate
try:
    rotated = sphere.Rotate(Axis(Pnt(0, 0, 0), direction=(0, 0, 1)), 45)
    print("✅ sphere.Rotate(Axis, angle) works")
except Exception as e:
    print(f"❌ sphere.Rotate failed: {type(e).__name__}")

# Test Mirror
try:
    mirrored = sphere.Mirror(Axis(Pnt(0, 0, 0), direction=(0, 0, 1)))
    print("✅ sphere.Mirror(Axis) works")
except Exception as e:
    print(f"❌ sphere.Mirror failed: {type(e).__name__}")
print()

# Test 8: Boolean operations (should work)
print("-" * 70)
print("8. Testing Boolean operations:")
print("-" * 70)
try:
    s1 = Sphere(Pnt(0, 0, 0), 0.1)
    s2 = Sphere(Pnt(0.05, 0, 0), 0.08)

    union = s1 + s2
    print("✅ Boolean union (s1 + s2) works")

    intersection = s1 * s2
    print("✅ Boolean intersection (s1 * s2) works")

    difference = s1 - s2
    print("✅ Boolean subtraction (s1 - s2) works")

except Exception as e:
    print(f"❌ Boolean operations failed: {e}")
print()

# Test 9: Check for Box primitive
print("-" * 70)
print("9. Testing Box primitive:")
print("-" * 70)
try:
    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    print("✅ Box(Pnt1, Pnt2) works")
except Exception as e:
    print(f"❌ Box creation failed: {e}")
print()

# Test 10: Try to access underlying OCC Core
print("-" * 70)
print("10. Testing access to OCC.Core (low-level):")
print("-" * 70)
try:
    from OCC.Core.gp import gp_GTrsf, gp_XYZ
    print("✅ OCC.Core.gp module accessible")

    from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
    print("✅ OCC.Core.BRepBuilderAPI accessible")

    print("   (Low-level OCC transformations may be possible)")
except ImportError as e:
    print(f"❌ OCC.Core not accessible: {e}")
    print("   (Direct OCC access not available)")
print()

# Summary and recommendations
print("=" * 70)
print("SUMMARY & RECOMMENDATIONS")
print("=" * 70)
print()

recommendations = []

# Check Ellipsoid availability
if 'ellipsoid' in locals() and ellipsoid is not None:
    recommendations.append("✅ Use direct Ellipsoid primitive (Approach 1)")
else:
    recommendations.append("❌ No direct Ellipsoid primitive found")

# Check Scale availability
if 'scaled' in locals() and scaled is not None:
    recommendations.append("✅ Use Scale transformation for ellipsoids (Approach 2)")
else:
    recommendations.append("❌ Scale transformation not available or limited")

# Check Cylinder
if 'cyl' in locals():
    recommendations.append("✅ Use Cylinder for dipole antenna")

# Check WorkPlane
if 'wp' in locals():
    recommendations.append("✅ WorkPlane available for complex profiles")

# Fallback option
recommendations.append("🔄 Fallback: Use sphere-with-dent geometry (boolean ops work)")

print("Implementation recommendations:")
for i, rec in enumerate(recommendations, 1):
    print(f"{i}. {rec}")
print()

print("Next steps:")
print("- Review the output above")
print("- Choose implementation approach based on available features")
print("- Implement geometry generators in src/geometry.py")
print()
print("=" * 70)
