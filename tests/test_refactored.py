"""
Test script for refactored MaxwellScattererAndAntennaNGSolve package.

This script tests the basic functionality of the refactored code
without running the full simulation (which requires NGSolve).
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")

    try:
        from src import (
            create_spherical_geometry,
            create_incident_wave,
            ScattererProblem,
            AntennaProblem,
            solve_gmres,
            solve_direct,
            PhysicalParameters
        )
        print("✓ All imports successful")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        return False


def test_physical_parameters():
    """Test PhysicalParameters dataclass."""
    print("\nTesting PhysicalParameters...")

    try:
        from src import PhysicalParameters

        params = PhysicalParameters(wavelength=0.65)
        assert abs(params.wavenumber - 2*3.14159/0.65) < 0.01
        print(f"✓ PhysicalParameters works")
        print(f"  Wavelength: {params.wavelength}")
        print(f"  Wavenumber: {params.wavenumber:.4f}")
        return True
    except Exception as e:
        print(f"✗ PhysicalParameters test failed: {e}")
        return False


def test_utility_functions():
    """Test utility functions."""
    print("\nTesting utility functions...")

    try:
        from src import compute_wavelength, compute_frequency, normalize_vector

        # Test wavelength/frequency conversion
        freq = 1e9  # 1 GHz
        wl = compute_wavelength(freq)
        freq_back = compute_frequency(wl)
        assert abs(freq - freq_back) < 1.0
        print(f"✓ Wavelength/frequency conversion works")
        print(f"  1 GHz → λ = {wl:.4f} m")

        # Test vector normalization
        vec = (3.0, 4.0, 0.0)
        norm_vec = normalize_vector(vec)
        magnitude = (norm_vec[0]**2 + norm_vec[1]**2 + norm_vec[2]**2)**0.5
        assert abs(magnitude - 1.0) < 1e-10
        print(f"✓ Vector normalization works")
        print(f"  {vec} → {norm_vec}")

        return True
    except Exception as e:
        print(f"✗ Utility functions test failed: {e}")
        return False


def test_module_structure():
    """Test that all expected modules and classes exist."""
    print("\nTesting module structure...")

    try:
        import src

        # Check version info
        assert hasattr(src, '__version__')
        print(f"✓ Package version: {src.__version__}")

        # Check all expected exports
        expected = [
            'create_spherical_geometry',
            'create_incident_wave',
            'ScattererProblem',
            'AntennaProblem',
            'solve_gmres',
            'solve_direct',
            'PhysicalParameters',
        ]

        for name in expected:
            assert hasattr(src, name), f"Missing export: {name}"

        print(f"✓ All expected exports present ({len(expected)} items)")
        return True
    except Exception as e:
        print(f"✗ Module structure test failed: {e}")
        return False


def main():
    """Run all tests."""
    print("=" * 60)
    print("MaxwellScattererAndAntennaNGSolve Refactored Code Test Suite")
    print("=" * 60)

    tests = [
        test_imports,
        test_physical_parameters,
        test_utility_functions,
        test_module_structure,
    ]

    results = []
    for test in tests:
        results.append(test())

    print("\n" + "=" * 60)
    print(f"Test Results: {sum(results)}/{len(results)} passed")
    print("=" * 60)

    if all(results):
        print("\n✓ All tests passed! The refactored code structure is working.")
        print("\nNote: Full simulation test requires NGSolve to be installed.")
        print("      Use the Jupyter notebook to run a complete simulation.")
        return 0
    else:
        print("\n✗ Some tests failed. Check the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
