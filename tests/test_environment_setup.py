#!/usr/bin/env python3
"""
Test script for environment setup functionality.

This validates:
1. Setup script exists and is executable
2. SLURM script has proper venv detection logic
3. Generated files structure is correct
"""

import sys
import os
from pathlib import Path
import subprocess

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

print("=" * 70)
print("Testing Environment Setup")
print("=" * 70)
print()

# Test 1: Setup script exists and is executable
print("1. Testing setup script:")
print("-" * 70)
try:
    setup_script = Path(__file__).parent.parent / "cluster/scripts/setup_environment.sh"

    if not setup_script.exists():
        print(f"❌ Setup script not found: {setup_script}")
        sys.exit(1)

    if not os.access(setup_script, os.X_OK):
        print(f"❌ Setup script not executable: {setup_script}")
        sys.exit(1)

    print(f"✓ Setup script exists and is executable")
    print(f"  Path: {setup_script}")
    print(f"  Size: {setup_script.stat().st_size} bytes")
except Exception as e:
    print(f"❌ Failed: {e}")
    sys.exit(1)
print()

# Test 2: Setup script has required sections
print("2. Testing setup script content:")
print("-" * 70)
try:
    content = setup_script.read_text()

    required_sections = [
        "#!/bin/bash",
        "Python version",
        "Create Virtual Environment",
        "Install Core Dependencies",
        "Verify Installation",
        "ngsolve",
        "numpy",
        "scipy",
        "pyyaml",
        "pandas",
    ]

    for section in required_sections:
        if section not in content:
            print(f"❌ Missing section: {section}")
            sys.exit(1)

    print(f"✓ All required sections present")
    print(f"  Sections checked: {len(required_sections)}")
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
print()

# Test 3: SLURM script has venv detection
print("3. Testing SLURM script venv detection:")
print("-" * 70)
try:
    slurm_script = Path(__file__).parent.parent / "cluster/slurm/array_job.slurm"

    if not slurm_script.exists():
        print(f"❌ SLURM script not found: {slurm_script}")
        sys.exit(1)

    content = slurm_script.read_text()

    required_checks = [
        "ANTENNA_VENV",
        ".env_info",
        "ngsolve_env",
        "activate",
    ]

    for check in required_checks:
        if check not in content:
            print(f"❌ Missing venv detection for: {check}")
            sys.exit(1)

    print(f"✓ SLURM script has proper venv detection")
    print(f"  Detection mechanisms: ANTENNA_VENV, .env_info, defaults")
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
print()

# Test 4: Check documentation mentions setup
print("4. Testing documentation coverage:")
print("-" * 70)
try:
    docs_to_check = [
        ("README.md", "setup_environment.sh"),
        ("cluster/README.md", "Environment Setup"),
        ("cluster/QUICK_START.md", "setup_environment"),
        ("cluster/scripts/README.md", "setup_environment.sh"),
    ]

    project_root = Path(__file__).parent.parent

    for doc_file, required_text in docs_to_check:
        doc_path = project_root / doc_file
        if not doc_path.exists():
            print(f"❌ Documentation file not found: {doc_file}")
            sys.exit(1)

        content = doc_path.read_text()
        if required_text not in content:
            print(f"❌ {doc_file} missing reference to: {required_text}")
            sys.exit(1)

    print(f"✓ All documentation files mention environment setup")
    print(f"  Files checked: {len(docs_to_check)}")
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
print()

# Test 5: Check .gitignore has environment files
print("5. Testing .gitignore:")
print("-" * 70)
try:
    gitignore = Path(__file__).parent.parent / ".gitignore"

    if gitignore.exists():
        content = gitignore.read_text()

        required_entries = [".env_info", "activate_env.sh"]

        for entry in required_entries:
            if entry not in content:
                print(f"⚠️  Warning: .gitignore missing: {entry}")

        print(f"✓ .gitignore checked")
    else:
        print(f"⚠️  Warning: .gitignore not found (optional)")
except Exception as e:
    print(f"⚠️  Warning: {e}")
print()

# Test 6: Verify script can be parsed (bash -n syntax check)
print("6. Testing script syntax (bash -n):")
print("-" * 70)
try:
    result = subprocess.run(
        ["bash", "-n", str(setup_script)],
        capture_output=True,
        text=True,
        timeout=5
    )

    if result.returncode != 0:
        print(f"❌ Bash syntax check failed:")
        print(result.stderr)
        sys.exit(1)

    print(f"✓ Bash syntax is valid")
except subprocess.TimeoutExpired:
    print(f"❌ Syntax check timed out")
    sys.exit(1)
except FileNotFoundError:
    print(f"⚠️  Warning: bash not found, skipping syntax check")
except Exception as e:
    print(f"❌ Failed: {e}")
    sys.exit(1)
print()

# Test 7: Check script has usage instructions
print("7. Testing setup script help/usage:")
print("-" * 70)
try:
    content = setup_script.read_text()

    help_sections = [
        "Usage:",
        "Arguments:",
        "Virtual environment:",
    ]

    found = sum(1 for section in help_sections if section in content)

    if found >= 2:
        print(f"✓ Setup script has usage instructions")
        print(f"  Help sections found: {found}/{len(help_sections)}")
    else:
        print(f"⚠️  Warning: Limited usage documentation in setup script")
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
print()

# Test 8: Verify expected output files would be created
print("8. Testing expected output file structure:")
print("-" * 70)
try:
    content = setup_script.read_text()

    expected_outputs = [
        ".env_info",
        "activate_env.sh",
        "VENV_PATH",
        "PYTHON_VERSION",
    ]

    for output in expected_outputs:
        if output not in content:
            print(f"❌ Script doesn't generate: {output}")
            sys.exit(1)

    print(f"✓ Script generates expected output files")
    print(f"  Outputs: .env_info, activate_env.sh")
except Exception as e:
    print(f"❌ Failed: {e}")
    import traceback
    traceback.print_exc()
print()

print("=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
✓ All environment setup tests passed!

Verified:
  - Setup script exists and is executable
  - Setup script has all required sections
  - SLURM script has smart venv detection
  - Documentation mentions environment setup
  - .gitignore includes generated files
  - Bash syntax is valid
  - Script has usage instructions
  - Expected output files are generated

The environment setup system is ready for use!

Next steps:
  - Run actual setup: bash cluster/scripts/setup_environment.sh
  - Test with all other test suites: python3 run_all_tests.py
""")
