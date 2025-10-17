#!/usr/bin/env python3
"""
Unified test runner for MaxwellScattererAndAntennaNGSolve.

Runs all test suites in the tests/ directory in sequence.
Reports overall success/failure status.

Usage:
    python3 run_all_tests.py [--verbose]
"""

import sys
import subprocess
from pathlib import Path
import argparse
from typing import Tuple

# Test files in order of execution
TEST_FILES = [
    "tests/test_refactored.py",
    "tests/test_occ_syntax.py",
    "tests/test_ellipsoid.py",
    "tests/test_final_syntax.py",
    "tests/test_geometry_generators.py",
    "tests/test_environment_setup.py",
]

def run_test(test_file: Path, verbose: bool = False) -> Tuple[bool, str]:
    """
    Run a single test file.

    Returns:
        (success, output): Boolean success status and output string
    """
    print(f"\n{'=' * 70}")
    print(f"Running: {test_file.name}")
    print(f"{'=' * 70}")

    try:
        result = subprocess.run(
            ["python3", str(test_file)],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout per test
        )

        output = result.stdout + result.stderr

        if verbose or result.returncode != 0:
            print(output)

        if result.returncode == 0:
            print(f"✅ {test_file.name} PASSED")
            return True, output
        else:
            print(f"❌ {test_file.name} FAILED (exit code: {result.returncode})")
            return False, output

    except subprocess.TimeoutExpired:
        print(f"❌ {test_file.name} TIMEOUT (>5 minutes)")
        return False, "Test timeout"
    except Exception as e:
        print(f"❌ {test_file.name} ERROR: {e}")
        return False, str(e)


def main():
    """Run all tests and report results."""
    parser = argparse.ArgumentParser(description="Run all MaxwellScattererAndAntennaNGSolve tests")
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show output for all tests (not just failures)'
    )
    args = parser.parse_args()

    print("=" * 70)
    print("MaxwellScattererAndAntennaNGSolve Test Suite")
    print("=" * 70)
    print(f"Running {len(TEST_FILES)} test suites...")

    project_root = Path(__file__).parent
    results = []

    for test_file in TEST_FILES:
        test_path = project_root / test_file
        if not test_path.exists():
            print(f"\n⚠️  Warning: Test file not found: {test_file}")
            results.append((test_file, False, "File not found"))
            continue

        success, output = run_test(test_path, verbose=args.verbose)
        results.append((test_file, success, output))

    # Summary
    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)

    passed = sum(1 for _, success, _ in results if success)
    failed = len(results) - passed

    for test_file, success, _ in results:
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status}  {test_file}")

    print("=" * 70)
    print(f"Total: {len(results)} tests")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print("=" * 70)

    if failed > 0:
        print("\n❌ Some tests failed!")
        return 1
    else:
        print("\n✅ All tests passed!")
        return 0


if __name__ == '__main__':
    sys.exit(main())
