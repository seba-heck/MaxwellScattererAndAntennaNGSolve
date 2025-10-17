# MaxwellScattererAndAntennaNGSolve Cluster Computing

HPC cluster workflow for running large-scale parametric sweeps of electromagnetic simulations.

**Supports both scattering and antenna problems with multiple geometry types** through a unified, auto-detecting workflow!

## First-Time Setup

**Before running any simulations, set up the Python environment:**

```bash
# Run setup script (creates virtual environment and installs dependencies)
bash cluster/scripts/setup_environment.sh

# Optional: specify custom venv location
bash cluster/scripts/setup_environment.sh /path/to/my/venv
```

This will:
- Create a Python virtual environment (default: `$HOME/.venvs/ngsolve_env`)
- Install all required dependencies (ngsolve, numpy, scipy, pyyaml, pandas, etc.)
- Verify the installation
- Create activation helpers

**The SLURM scripts will automatically detect and use this environment.**

See "Environment Setup" section below for more details.

## Quick Start

**Ellipsoid Scattering:**
1. **Define sweep**: Edit `configs/ellipsoid_geometry_sweep.yaml` (use `ellipsoid_semi_axis_*` + `propagation_dir`)
2. **Generate jobs**: `python scripts/generate_jobs.py configs/ellipsoid_geometry_sweep.yaml`
3. **Submit**: `sbatch --array=0-47 slurm/array_job.slurm ellipsoid_geometry_sweep`
4. **Post-process**: `python scripts/post_process.py results/ellipsoid_geometry_sweep`

**Dipole Antenna:**
1. **Define sweep**: Edit `configs/dipole_geometry_sweep.yaml` (use `dipole_length_factor` + `amplitude`)
2. **Generate jobs**: `python scripts/generate_jobs.py configs/dipole_geometry_sweep.yaml`
3. **Submit**: `sbatch --array=0-9 slurm/array_job.slurm dipole_geometry_sweep`
4. **Post-process**: `python scripts/post_process.py results/dipole_geometry_sweep`

**Spheroid Scattering:**
1. **Define sweep**: Edit `configs/spheroid_geometry_sweep.yaml` (use `spheroid_*_radius` + `propagation_dir`)
2. **Generate jobs**: `python scripts/generate_jobs.py configs/spheroid_geometry_sweep.yaml`
3. **Submit**: `sbatch --array=0-9 slurm/array_job.slurm spheroid_geometry_sweep`
4. **Post-process**: `python scripts/post_process.py results/spheroid_geometry_sweep`

> **Auto-Detection**: Problem type and geometry type are detected from config parameters:
> - **Problem**: `propagation_dir` → Scattering, `amplitude` → Antenna
> - **Geometry**: `ellipsoid_semi_axis_a` → Ellipsoid, `dipole_length_factor` → Dipole, etc.

## Large-Scale Parameter Sweeps (1000+ Cases)

For comprehensive parameter exploration:

```bash
# Generate sweep configurations
python3 scripts/generate_sweep_configs.py

# Run pilot studies first (recommended)
sbatch --array=0-8 slurm/array_job.slurm pilot_scattering
sbatch --array=0-11 slurm/array_job.slurm pilot_antenna

# After validation, run full sweeps
sbatch --array=0-999 slurm/array_job.slurm scattering_1000_sweep
sbatch --array=0-999 slurm/array_job.slurm antenna_1000_sweep
```

See `PARAMETER_SWEEPS.md` for complete documentation.

## Structure

```
cluster/
├── configs/                    # YAML parameter sweeps
│   ├── scattering_1000_sweep.yaml
│   ├── antenna_1000_sweep.yaml
│   ├── pilot_*.yaml           # Validation configs
│   └── *.yaml                 # Custom sweeps
├── scripts/                    # Python automation
│   ├── generate_sweep_configs.py
│   ├── generate_jobs.py
│   ├── run_simulation.py
│   └── post_process.py
├── slurm/                      # SLURM batch scripts
│   └── array_job.slurm
├── results/                    # Output (by sweep name)
└── logs/                       # Job logs
```

## Environment Setup

### Automatic Setup (Recommended)

Use the provided setup script:

```bash
bash cluster/scripts/setup_environment.sh
```

**What it does:**
1. Checks Python version (requires Python 3.8+)
2. Creates virtual environment at `$HOME/.venvs/ngsolve_env`
3. Installs all dependencies:
   - ngsolve (FEM solver)
   - numpy, scipy (scientific computing)
   - pyyaml (config files)
   - pandas (result processing)
   - matplotlib (plotting)
   - jupyter (optional, for notebooks)
4. Verifies installation
5. Tests MaxwellScattererAndAntennaNGSolve module imports
6. Creates `.env_info` file for SLURM auto-detection
7. Creates `activate_env.sh` convenience script

**Custom virtual environment location:**
```bash
bash cluster/scripts/setup_environment.sh /scratch/$USER/venv
```

### Manual Activation

After setup, activate the environment for local work:

```bash
# Option 1: Use convenience script
source activate_env.sh

# Option 2: Direct activation
source $HOME/.venvs/ngsolve_env/bin/activate
```

### SLURM Auto-Detection

The SLURM job scripts automatically detect the virtual environment in this priority order:

1. **`ANTENNA_VENV` environment variable** (most flexible):
   ```bash
   export ANTENNA_VENV=/custom/path/to/venv
   sbatch --array=0-8 slurm/array_job.slurm pilot_scattering
   ```

2. **`.env_info` file** (created by `setup_environment.sh`):
   - Automatically used if it exists in project root
   - No configuration needed

3. **Default locations** (backward compatible):
   - `$HOME/.venvs/ngsolve_env` (new default)
   - `$HOME/.venvs/ngs311` (legacy)

4. **System Python** (fallback with warning):
   - Used if no venv found
   - May fail if dependencies not installed system-wide

### Cluster-Specific Module Loads

Edit `cluster/slurm/array_job.slurm` lines 54-57 for your cluster:

```bash
# Module loads (adjust for your cluster)
module purge
module load stack/2024-06 gcc/12.2.0
module load python/3.11.6
```

Common alternatives:
```bash
# Example: Different cluster
module load python/3.10
module load gcc/11.2

# Example: Minimal setup
module load python
```

### Verifying Environment

Test the environment locally:

```bash
source activate_env.sh
python3 run_all_tests.py
```

Test on cluster (single job):

```bash
cd cluster
python3 scripts/run_simulation.py pilot_scattering 0
```

### Troubleshooting

**"No module named ngsolve":**
- Run `bash cluster/scripts/setup_environment.sh`
- Or install manually: `pip install ngsolve numpy scipy pyyaml pandas`

**"Virtual environment not found" in SLURM logs:**
- Check `.env_info` exists: `cat .env_info`
- Or set explicitly: `export ANTENNA_VENV=/path/to/venv`

**"Permission denied" when running setup:**
- Make executable: `chmod +x cluster/scripts/setup_environment.sh`

**NGSolve installation fails:**
- Check Python version: `python3 --version` (must be 3.8+)
- Try upgrading pip: `pip install --upgrade pip`
- Check available disk space

## Documentation

- `QUICK_START.md` - Quick reference guide
- `PARAMETER_SWEEPS.md` - Large-scale sweep documentation
- `scripts/README.md` - Script documentation
- `slurm/README.md` - SLURM documentation
