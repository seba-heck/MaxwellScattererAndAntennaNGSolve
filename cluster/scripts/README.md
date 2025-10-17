# Cluster Scripts

Automated workflow scripts for HPC parametric sweeps. All scripts support automatic detection of problem types (scattering/antenna) and geometry types (ellipsoid/spheroid/dipole/sphere).

## First-Time Setup

**Before using any scripts, set up the Python environment:**

```bash
bash scripts/setup_environment.sh
```

This creates a virtual environment with all dependencies and configures auto-detection for SLURM jobs. See `cluster/README.md` for details.

## Scripts

### 1. setup_environment.sh

Creates Python virtual environment and installs all dependencies.

```bash
bash scripts/setup_environment.sh [venv_path]
```

**Default venv location:** `$HOME/.venvs/ngsolve_env`

**What it installs:**
- ngsolve (FEM solver)
- numpy, scipy (scientific computing)
- pyyaml (config parsing)
- pandas (result processing)
- matplotlib (plotting)
- jupyter (notebooks, optional)

**Output:**
- Virtual environment at specified path
- `.env_info` file in project root (for SLURM auto-detection)
- `activate_env.sh` convenience script

**Run this once per cluster/machine.**

### 2. generate_sweep_configs.py

Generates large-scale parameter sweep configurations (1000+ cases).

```bash
python3 scripts/generate_sweep_configs.py
```

**Output:**
- `configs/scattering_1000_sweep.yaml` - 1000 scattering configurations
- `configs/antenna_1000_sweep.yaml` - 1000 antenna configurations
- `configs/pilot_scattering.yaml` - 9 validation cases
- `configs/pilot_antenna.yaml` - 12 validation cases

**Features:**
- Systematic parameter space exploration
- Dimensionless parameter coverage (ka, L/λ, aspect ratios)
- Pilot studies for validation
- Computational feasibility (20-60s per job)

See `PARAMETER_SWEEPS.md` for detailed documentation.

### 3. generate_jobs.py

Generates job list JSON from YAML parameter sweep configuration.

```bash
python3 scripts/generate_jobs.py configs/ellipsoid_geometry_sweep.yaml
```

**Output:** `configs/ellipsoid_geometry_sweep_jobs.json` containing all parameter combinations.

**Features:**
- Cartesian product of all parameter values
- Validates required parameters
- Supports all geometry types (auto-detected from parameters)
- Supports both scattering and antenna problems

### 4. run_simulation.py

Runs a single simulation job (called by SLURM array or locally).

```bash
python3 scripts/run_simulation.py \
    --job-file configs/sweep_jobs.json \
    --job-id 0 \
    --output-dir results/sweep \
    --num-threads 4 \
    --local
```

**Key Features:**
- **Auto-detects problem type**: Scattering (propagation_dir) vs Antenna (amplitude)
- **Auto-detects geometry type**: Ellipsoid, spheroid, dipole, or sphere
- **Environment detection**: SLURM vs local execution
- **Unified interface**: Same script for all problem/geometry types

**Geometry Detection Logic:**
```python
if 'dipole_length_factor' in params:
    → create_dipole_antenna_geometry()
elif 'ellipsoid_semi_axis_a' in params:
    → create_ellipsoid_scatterer_geometry()
elif 'spheroid_equatorial_radius' in params:
    → create_spheroid_scatterer_geometry()
else:
    → create_spherical_geometry()  # Backward compatible
```

### 5. post_process.py

Aggregates results from all jobs into summary CSV and statistics.

```bash
python3 scripts/post_process.py results/sweep_name
```

**Output:**
- `summary.csv` - All jobs in tabular format
- `summary_statistics.json` - Aggregate statistics

**Fields include:**
- All input parameters (wavelength, geometry parameters, etc.)
- `problem_type` field (scattering/antenna)
- `geometry_type` field (ellipsoid/spheroid/dipole/sphere)
- Mesh statistics (elements, DOFs)
- Timings (mesh, assembly, solve, total)
- Resource usage

### 6. run_local_test.py

Helper script for running full sweeps locally (for testing).

```bash
python3 scripts/run_local_test.py sweep_name
```

Runs all jobs sequentially on local machine.

## Supported Geometry Types

The scripts automatically detect geometry type from config parameters:

| Geometry | Detection Parameter | Generator Function |
|----------|---------------------|-------------------|
| Ellipsoid | `ellipsoid_semi_axis_a` | `create_ellipsoid_scatterer_geometry()` |
| Spheroid | `spheroid_equatorial_radius` | `create_spheroid_scatterer_geometry()` |
| Dipole | `dipole_length_factor` | `create_dipole_antenna_geometry()` |
| Sphere | (none of above) | `create_spherical_geometry()` |

## Example Workflow

```bash
# 1. Generate jobs
python3 scripts/generate_jobs.py configs/ellipsoid_geometry_sweep.yaml

# 2. Test one job locally
python3 scripts/run_simulation.py \
    --job-file configs/ellipsoid_geometry_sweep_jobs.json \
    --job-id 0 \
    --output-dir results/ellipsoid_geometry_sweep \
    --num-threads 4 \
    --local

# 3. Submit to cluster (if test passes)
cd ..
sbatch --array=0-47 cluster/slurm/array_job.slurm ellipsoid_geometry_sweep

# 4. Post-process results (after completion)
cd cluster
python3 scripts/post_process.py results/ellipsoid_geometry_sweep
```

## Integration

All scripts integrate with the `src/` module:
- `src/geometry.py` - Parametric geometry generators
- `src/maxwell.py` - Unified Maxwell problem formulation
- `src/utils.py` - Wave/source functions
- `src/solvers.py` - Iterative solvers
