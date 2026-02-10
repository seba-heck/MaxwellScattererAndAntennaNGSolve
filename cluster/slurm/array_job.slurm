#!/bin/bash
#SBATCH --job-name=scatterer_sweep
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=10G
#SBATCH --time=06:00:00
#SBATCH --hint=nomultithread
#SBATCH --output=logs/%x/%a.out
#SBATCH --error=logs/%x/%a.err

# MaxwellScattererAndAntennaNGSolve SLURM Job Array Script
#
# This script runs parametric sweeps using SLURM job arrays.
# Each array task runs a single simulation with parameters from the job list.
#
# Usage:
#   sbatch --array=0-N slurm/array_job.slurm <sweep_name>

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

SWEEP_NAME="${1:-default_sweep}"
PROJECT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
CLUSTER_DIR="${PROJECT_DIR}/cluster"
SCRIPTS_DIR="${CLUSTER_DIR}/scripts"
CONFIGS_DIR="${CLUSTER_DIR}/configs"

# Allow custom results base directory via environment variable
# Default: cluster/results (backward compatible)
# Usage: export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna_results
RESULTS_BASE="${ANTENNA_RESULTS_BASE:-${CLUSTER_DIR}/results}"
RESULTS_DIR="${RESULTS_BASE}/${SWEEP_NAME}"
LOGS_DIR="${ANTENNA_RESULTS_BASE:-${CLUSTER_DIR}/logs/${SWEEP_NAME}}"

JOB_FILE="${CONFIGS_DIR}/${SWEEP_NAME}_jobs.json"
TASK_ID="${SLURM_ARRAY_TASK_ID}"
NTHREADS="${SLURM_CPUS_PER_TASK:-1}"

# ============================================================================
# Environment Setup
# ============================================================================

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export SLURM_CPU_BIND=verbose

# Module loads (adjust for your cluster)
module purge
module load stack/2024-06 gcc/12.2.0
module load python/3.11.6

# Activate Python virtual environment
# Priority order:
#   1. ANTENNA_VENV environment variable (most flexible)
#   2. .env_info file in project root (created by setup_environment.sh)
#   3. Default locations (backward compatible)

if [ -n "${ANTENNA_VENV:-}" ] && [ -d "${ANTENNA_VENV}" ]; then
    echo "Using virtual environment from ANTENNA_VENV: ${ANTENNA_VENV}"
    source "${ANTENNA_VENV}/bin/activate"
elif [ -f "${PROJECT_DIR}/.env_info" ]; then
    VENV_PATH=$(grep "^VENV_PATH=" "${PROJECT_DIR}/.env_info" | cut -d= -f2)
    if [ -n "${VENV_PATH}" ] && [ -d "${VENV_PATH}" ]; then
        echo "Using virtual environment from .env_info: ${VENV_PATH}"
        source "${VENV_PATH}/bin/activate"
    fi
elif [ -d "$HOME/.venvs/ngsolve_env" ]; then
    echo "Using default virtual environment: $HOME/.venvs/ngsolve_env"
    source "$HOME/.venvs/ngsolve_env/bin/activate"
elif [ -d "$HOME/.venvs/ngs311" ]; then
    echo "Using legacy virtual environment: $HOME/.venvs/ngs311"
    source "$HOME/.venvs/ngs311/bin/activate"
elif [ -d "$HOME/VENVS/ma" ]; then
    echo "Using custom virtual environment: $HOME/VENVS/ma"
    source "$HOME/VENVS/ma/bin/activate"
else
    echo "⚠️  Warning: No virtual environment found, using system Python"
    echo "   Run: bash cluster/scripts/setup_environment.sh"
fi

# ============================================================================
# Pre-flight Checks
# ============================================================================

echo "======================================================================"
echo "MaxwellScattererAndAntennaNGSolve Cluster Job Array"
echo "======================================================================"
echo "Job ID:          ${SLURM_JOB_ID}"
echo "Array Task ID:   ${TASK_ID}"
echo "Node:            ${SLURM_NODELIST}"
echo "CPUs per task:   ${NTHREADS}"
echo "Sweep:           ${SWEEP_NAME}"
echo "Job file:        ${JOB_FILE}"
echo "Results base:    ${RESULTS_BASE}"
echo "Results dir:     ${RESULTS_DIR}"
echo "======================================================================"

mkdir -p "${RESULTS_DIR}"
mkdir -p "${LOGS_DIR}"

if [ ! -f "${JOB_FILE}" ]; then
    echo "ERROR: Job file not found: ${JOB_FILE}"
    echo "Run: python scripts/generate_jobs.py configs/${SWEEP_NAME}.yaml"
    exit 1
fi

echo ""
echo "Verifying Python environment..."
python3 -V
python3 -c "import ngsolve, numpy; print(f'✓ ngsolve {ngsolve.__version__}, numpy {numpy.__version__}')"

if [ $? -ne 0 ]; then
    echo "ERROR: Environment check failed"
    exit 1
fi

# ============================================================================
# Run Simulation
# ============================================================================

echo ""
echo "======================================================================"
echo "Starting simulation for task ${TASK_ID}..."
echo "======================================================================"
echo ""

/usr/bin/time -v -o "${LOGS_DIR}/time_${TASK_ID}.log" \
    srun --cpu-bind=cores \
    python3 "${SCRIPTS_DIR}/run_simulation.py" \
        --job-file "${JOB_FILE}" \
        --job-id "${TASK_ID}" \
        --output-dir "${RESULTS_DIR}" \
        --num-threads "${NTHREADS}"

EXIT_CODE=$?

# ============================================================================
# Cleanup and Exit
# ============================================================================

echo ""
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "======================================================================"
    echo "✓ Task ${TASK_ID} completed successfully"
    echo "======================================================================"
else
    echo "======================================================================"
    echo "✗ Task ${TASK_ID} failed with exit code ${EXIT_CODE}"
    echo "======================================================================"
fi

exit ${EXIT_CODE}
