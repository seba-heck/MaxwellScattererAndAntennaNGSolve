#!/bin/bash
#
# Environment Setup Script for MaxwellScattererAndAntennaNGSolve
#
# This script creates a Python virtual environment and installs all dependencies
# required for running MaxwellScattererAndAntennaNGSolve simulations on an HPC cluster.
#
# Usage:
#   bash cluster/scripts/setup_environment.sh [venv_path]
#
# Arguments:
#   venv_path  - Optional path for virtual environment (default: $HOME/.venvs/ngsolve_env)
#

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

VENV_PATH="${1:-$HOME/.venvs/ngsolve_env}"
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

echo "========================================================================"
echo "MaxwellScattererAndAntennaNGSolve Environment Setup"
echo "========================================================================"
echo "Virtual environment: ${VENV_PATH}"
echo "Project root:        ${PROJECT_ROOT}"
echo ""

# ============================================================================
# Check Python Version
# ============================================================================

echo "Checking Python version..."
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
PYTHON_MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
PYTHON_MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)

echo "  Found Python ${PYTHON_VERSION}"

if [ "$PYTHON_MAJOR" -lt 3 ] || { [ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 8 ]; }; then
    echo "  ❌ ERROR: Python 3.8 or higher required"
    echo ""
    echo "  Please load a compatible Python module:"
    echo "    module load python/3.11"
    echo "    module load python/3.10"
    exit 1
fi

echo "  ✓ Python version is compatible"
echo ""

# ============================================================================
# Create Virtual Environment
# ============================================================================

if [ -d "${VENV_PATH}" ]; then
    echo "⚠️  Virtual environment already exists at: ${VENV_PATH}"
    read -p "   Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "   Removing existing environment..."
        rm -rf "${VENV_PATH}"
    else
        echo "   Keeping existing environment, will try to update packages..."
    fi
fi

if [ ! -d "${VENV_PATH}" ]; then
    echo "Creating virtual environment..."
    python3 -m venv "${VENV_PATH}"
    echo "  ✓ Virtual environment created"
else
    echo "Using existing virtual environment..."
fi

echo ""

# ============================================================================
# Activate Virtual Environment
# ============================================================================

echo "Activating virtual environment..."
source "${VENV_PATH}/bin/activate"

if [ $? -ne 0 ]; then
    echo "  ❌ ERROR: Failed to activate virtual environment"
    exit 1
fi

echo "  ✓ Virtual environment activated"
echo ""

# ============================================================================
# Upgrade pip
# ============================================================================

echo "Upgrading pip, setuptools, and wheel..."
pip install --upgrade pip setuptools wheel

if [ $? -ne 0 ]; then
    echo "  ❌ ERROR: Failed to upgrade pip"
    exit 1
fi

echo "  ✓ pip upgraded"
echo ""

# ============================================================================
# Install Core Dependencies
# ============================================================================

echo "Installing core dependencies..."
echo "  This may take several minutes, especially for NGSolve..."
echo ""

# Core scientific computing
echo "1/7 Installing numpy..."
pip install numpy

echo "2/7 Installing scipy..."
pip install scipy

echo "3/7 Installing matplotlib..."
pip install matplotlib

# NGSolve (this is the big one)
echo "4/7 Installing ngsolve (this may take a while)..."
pip install ngsolve

# Cluster workflow dependencies
echo "5/7 Installing pyyaml..."
pip install pyyaml

echo "6/7 Installing pandas..."
pip install pandas

# Optional but useful
echo "7/7 Installing jupyter (optional, for notebooks)..."
pip install jupyter ipython

echo ""
echo "  ✓ All dependencies installed"
echo ""

# ============================================================================
# Verify Installation
# ============================================================================

echo "Verifying installation..."
echo ""

python3 << 'EOF'
import sys
import importlib

packages = {
    'ngsolve': 'NGSolve (FEM solver)',
    'netgen': 'Netgen (mesh generation)',
    'numpy': 'NumPy (arrays)',
    'scipy': 'SciPy (scientific computing)',
    'yaml': 'PyYAML (config files)',
    'pandas': 'Pandas (data processing)',
    'matplotlib': 'Matplotlib (plotting)',
}

print("Package Verification:")
print("-" * 70)

all_ok = True
for module, description in packages.items():
    try:
        mod = importlib.import_module(module)
        version = getattr(mod, '__version__', 'unknown')
        print(f"  ✓ {description:30} v{version}")
    except ImportError:
        print(f"  ❌ {description:30} NOT FOUND")
        all_ok = False

print("-" * 70)

if all_ok:
    print("  ✓ All packages verified successfully!")
    sys.exit(0)
else:
    print("  ❌ Some packages failed to import")
    sys.exit(1)
EOF

VERIFY_EXIT=$?

echo ""

if [ $VERIFY_EXIT -ne 0 ]; then
    echo "❌ Verification failed!"
    echo ""
    echo "Some packages could not be imported. Please check the error messages above."
    exit 1
fi

# ============================================================================
# Test MaxwellScattererAndAntennaNGSolve Import
# ============================================================================

echo "Testing MaxwellScattererAndAntennaNGSolve module imports..."
echo ""

cd "${PROJECT_ROOT}"

python3 << 'EOF'
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path.cwd()))

try:
    from src import (
        create_spherical_geometry,
        create_ellipsoid_scatterer_geometry,
        create_spheroid_scatterer_geometry,
        create_dipole_antenna_geometry,
        create_incident_wave,
        create_antenna_source,
        MaxwellProblem,
        solve_gmres,
        solve_direct,
    )
    print("  ✓ All MaxwellScattererAndAntennaNGSolve modules imported successfully")
    print("")
    print("Available geometry generators:")
    print("  - create_spherical_geometry()")
    print("  - create_ellipsoid_scatterer_geometry()")
    print("  - create_spheroid_scatterer_geometry()")
    print("  - create_dipole_antenna_geometry()")
    print("")
    print("Available problem formulations:")
    print("  - MaxwellProblem (unified interface)")
    print("  - ScattererProblem (legacy)")
    print("  - AntennaProblem (legacy)")

except ImportError as e:
    print(f"  ❌ Failed to import MaxwellScattererAndAntennaNGSolve modules: {e}")
    sys.exit(1)
EOF

if [ $? -ne 0 ]; then
    echo ""
    echo "❌ MaxwellScattererAndAntennaNGSolve module test failed!"
    echo ""
    echo "This is not critical if you're only setting up the environment."
    echo "Make sure the src/ directory exists in the project root."
fi

# ============================================================================
# Create Activation Helper
# ============================================================================

ACTIVATE_SCRIPT="${PROJECT_ROOT}/activate_env.sh"

cat > "${ACTIVATE_SCRIPT}" << EOF
#!/bin/bash
# Quick activation script for MaxwellScattererAndAntennaNGSolve environment
# Usage: source activate_env.sh

source "${VENV_PATH}/bin/activate"

if [ \$? -eq 0 ]; then
    echo "✓ MaxwellScattererAndAntennaNGSolve environment activated"
    echo "  Python: \$(which python3)"
    echo "  Virtual env: ${VENV_PATH}"
else
    echo "❌ Failed to activate environment"
fi
EOF

chmod +x "${ACTIVATE_SCRIPT}"

echo ""
echo "  ✓ Created activation helper: activate_env.sh"

# ============================================================================
# Create Environment Info File
# ============================================================================

ENV_INFO="${PROJECT_ROOT}/.env_info"

cat > "${ENV_INFO}" << EOF
# MaxwellScattererAndAntennaNGSolve Environment Configuration
# Auto-generated by setup_environment.sh on $(date)

VENV_PATH=${VENV_PATH}
PYTHON_VERSION=${PYTHON_VERSION}
SETUP_DATE=$(date -Iseconds)
SETUP_HOST=$(hostname)
EOF

echo "  ✓ Created environment info: .env_info"

# ============================================================================
# Update SLURM Script Reference
# ============================================================================

SLURM_SCRIPT="${PROJECT_ROOT}/cluster/slurm/array_job.slurm"

if [ -f "${SLURM_SCRIPT}" ]; then
    # Check if SLURM script references the default venv path
    if grep -q "\.venvs/ngs311" "${SLURM_SCRIPT}"; then
        echo ""
        echo "⚠️  Note: cluster/slurm/array_job.slurm references a hardcoded venv path"
        echo "   Current venv: ${VENV_PATH}"
        echo ""
        if [ "${VENV_PATH}" != "$HOME/.venvs/ngs311" ]; then
            echo "   You may need to update the SLURM script to use your venv path,"
            echo "   or create a symlink:"
            echo "     mkdir -p $HOME/.venvs"
            echo "     ln -s ${VENV_PATH} $HOME/.venvs/ngs311"
        fi
    fi
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "========================================================================"
echo "✓ Environment Setup Complete!"
echo "========================================================================"
echo ""
echo "Virtual environment: ${VENV_PATH}"
echo ""
echo "To activate the environment manually:"
echo "  source ${VENV_PATH}/bin/activate"
echo ""
echo "Or use the convenience script:"
echo "  source activate_env.sh"
echo ""
echo "The SLURM job scripts will automatically activate this environment."
echo ""
echo "Next steps:"
echo "  1. Test locally:"
echo "       source activate_env.sh"
echo "       python3 run_all_tests.py"
echo ""
echo "  2. Run pilot study:"
echo "       cd cluster"
echo "       python3 scripts/run_simulation.py pilot_scattering 0"
echo ""
echo "  3. Submit to cluster (after pilot validation):"
echo "       sbatch --array=0-8 slurm/array_job.slurm pilot_scattering"
echo ""
echo "See cluster/QUICK_START.md for more information."
echo "========================================================================"
