# MaxwellScattererAndAntennaNGSolve

⚠️ This is a fork of tellocam/MaxwellScattererAndAntennaNGSolve to modify and add some features.

Electromagnetic scattering and antenna simulations using NGSolve finite element method library with HPC cluster support.

## Overview

MaxwellScattererAndAntennaNGSolve is a modular framework for solving electromagnetic scattering and antenna radiation problems using the NGSolve/Netgen finite element library. It provides both interactive development (Jupyter notebooks) and scalable HPC cluster execution for large parametric studies.

**Key Features:**
- Modular Python package for electromagnetic simulations
- Scattering and antenna radiation problem formulations
- Maxwell-Helmholtz equation solver with PML boundary conditions
- Standalone runner scripts for easy simulation execution
- HPC cluster workflow with SLURM job arrays
- Local testing capability (no cluster required)
- Parametric sweep support (1000s of parameter combinations)
- Automated result aggregation and post-processing

## Project Structure

```
MaxwellScattererAndAntennaNGSolve/
├── src/                      # Core simulation modules
│   ├── geometry.py           # Parametric geometry generators (ellipsoid, spheroid, dipole, sphere)
│   ├── utils.py              # Wave/source functions
│   ├── maxwell.py            # Unified Maxwell problem (recommended)
│   ├── scatterer.py          # Scattering problem formulation (legacy)
│   ├── antenna.py            # Antenna radiation problem formulation (legacy)
│   └── solvers.py            # GMRes, direct solvers
├── tests/                    # Test suite
│   ├── test_refactored.py    # Core module tests
│   ├── test_occ_syntax.py    # OCC primitive verification
│   ├── test_ellipsoid.py     # Ellipsoid geometry tests
│   ├── test_final_syntax.py  # Complete workflow tests
│   └── test_geometry_generators.py  # Parametric generator tests (10 tests)
├── run_all_tests.py          # Unified test runner
├── run_scatterer.py          # Standalone scattering simulation script
├── run_antenna.py            # Standalone antenna simulation script
├── cluster/                  # HPC cluster workflow
│   ├── configs/              # YAML parameter sweeps
│   ├── scripts/              # Job generation, execution, post-processing
│   ├── slurm/                # SLURM batch scripts
│   ├── results/              # Simulation outputs
│   └── QUICK_START.md        # Cluster usage guide
├── helmholtzTest.ipynb       # Original interactive notebook
└── README.md                 # This file
```

## Quick Start

### Installation

**Option 1: Automated Setup (Recommended for Cluster)**

```bash
# Navigate to project
cd MaxwellScattererAndAntennaNGSolve

# Run setup script (creates venv and installs all dependencies)
bash cluster/scripts/setup_environment.sh

# Activate environment
source activate_env.sh
```

**Option 2: Manual Installation (Local Development)**

```bash
# Install NGSolve and dependencies
pip install ngsolve jupyter numpy scipy matplotlib pandas pyyaml

# Or in a virtual environment
python3 -m venv venv
source venv/bin/activate
pip install ngsolve jupyter numpy scipy matplotlib pandas pyyaml
```

The automated setup script is recommended for cluster deployment as it handles dependency installation and creates auto-detection files for SLURM jobs.

### Running Simulations

**Option 1: Standalone Scripts (Easiest)**

Run a scattering simulation:
```bash
./run_scatterer.py --wavelength 0.65 --mesh-size 0.3 --output results_scatterer.json
```

Run an antenna simulation:
```bash
./run_antenna.py --wavelength 0.65 --amplitude 1.0 --output results_antenna.json
```

Use `--help` to see all available options.

**Option 2: Interactive Python (Unified Approach - Recommended)**

For scattering:
```python
from src import create_ellipsoid_scatterer_geometry, create_incident_wave
from src import MaxwellProblem, solve_gmres
from ngsolve import TaskManager, pi

# Create ellipsoid scatterer mesh (or use other geometry generators)
mesh = create_ellipsoid_scatterer_geometry(
    wavelength=0.65,
    semi_axis_a=0.05,
    semi_axis_b=0.08,
    semi_axis_c=0.12,
    domain_radius=1.0,
    pml_width=0.25
)
k = 2*pi/0.65
source = create_incident_wave(k, propagation_dir=(0,0,1), polarization=(1,0,0))
problem = MaxwellProblem(mesh, k, source)
problem.assemble_system()

with TaskManager():
    solution = solve_gmres(problem.a, problem.l, problem.fes)
```

For antenna:
```python
from src import create_dipole_antenna_geometry, create_antenna_source
from src import MaxwellProblem, solve_gmres
from ngsolve import TaskManager, pi

# Create cylindrical dipole antenna mesh
mesh = create_dipole_antenna_geometry(
    wavelength=0.65,
    length_factor=0.5,      # L = 0.5λ (half-wave dipole)
    radius_factor=0.01,     # r = 0.01λ
    domain_radius=1.0,
    pml_width=0.25
)
k = 2*pi/0.65
source = create_antenna_source(polarization=(0,0,1), amplitude=1.0)
problem = MaxwellProblem(mesh, k, source)
problem.assemble_system()

with TaskManager():
    solution = solve_gmres(problem.a, problem.l, problem.fes)
```

> **Note**: `ScattererProblem` and `AntennaProblem` are still available for backwards compatibility, but `MaxwellProblem` is recommended as it provides a unified interface.

**Option 3: Jupyter Notebook**

```bash
jupyter notebook helmholtzTest.ipynb
```

### Local Testing (No Cluster)

Test the cluster workflow locally:

```bash
cd cluster

# Generate jobs from config
python3 scripts/generate_jobs.py configs/test_quick.yaml

# Run locally (single job, ~3 seconds)
python3 scripts/run_simulation.py \
    --job-file configs/test_quick_jobs.json \
    --job-id 0 \
    --output-dir results/test_quick \
    --local

# Check results
cat results/test_quick/job_0000/metadata.json
```

### Custom Results Directory

By default, results are stored in `cluster/results/`. For cluster deployments, you can specify a custom location:

```bash
# Set custom results base directory (e.g., cluster scratch)
export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna_results

# Submit jobs (automatically uses custom path)
cd cluster
sbatch --array=0-47 slurm/array_job.slurm ellipsoid_geometry_sweep

# Results will be in: /scratch/$USER/antenna_results/ellipsoid_geometry_sweep/
```

**Use cases:**
- Cluster scratch filesystems (faster I/O)
- Shared project directories
- Better quota management for large sweeps

### Cluster Deployment

For large parametric sweeps on HPC clusters:

**Scattering Study:**
```bash
cd cluster

# 1. Create parameter sweep config (or use template)
cp configs/template.yaml configs/my_scattering_sweep.yaml
# Edit to set wavelengths, propagation_dir, polarization

# 2. Generate job list
python3 scripts/generate_jobs.py configs/my_scattering_sweep.yaml

# 3. Submit to SLURM
sbatch --array=0-99 slurm/array_job.slurm my_scattering_sweep

# 4. Post-process results
python3 scripts/post_process.py results/my_scattering_sweep
```

**Antenna Study:**
```bash
cd cluster

# 1. Create antenna sweep config (or use template)
cp configs/antenna_template.yaml configs/my_antenna_sweep.yaml
# Edit to set wavelengths, amplitude, polarization

# 2. Generate and submit (same workflow!)
python3 scripts/generate_jobs.py configs/my_antenna_sweep.yaml
sbatch --array=0-23 slurm/array_job.slurm my_antenna_sweep

# 3. Post-process results
python3 scripts/post_process.py results/my_antenna_sweep
```

> **Auto-Detection**: The cluster workflow automatically detects problem type:
> - Configs with `propagation_dir` → Scattering problem
> - Configs with `amplitude` → Antenna problem

See `cluster/QUICK_START.md` for detailed usage examples.

## Architecture

### Unified Maxwell Problem

MaxwellScattererAndAntennaNGSolve uses a **unified architecture** where both scattering and antenna problems are solved with the same code:

- **Single problem class**: `MaxwellProblem` handles both problem types
- **Auto-detection**: System detects problem type from config parameters
  - `propagation_dir` → Scattering problem (incident wave)
  - `amplitude` → Antenna problem (excitation source)
- **Same equations**: Both problems use identical weak formulation, only the source term differs
- **Dual auto-detection**: Problem type AND geometry type are both detected automatically

**Why unified?**
- Eliminates code duplication (previously 95% identical code)
- Reflects underlying physics (same Maxwell equations)
- Simpler maintenance and testing
- One workflow for all simulations

**Legacy support**: `ScattererProblem` and `AntennaProblem` classes still available for backward compatibility.

## Physics & Formulation

### Maxwell-Helmholtz Equations

The code solves the time-harmonic Maxwell equations in frequency domain:

```
∇ × (∇ × E) - k²E = 0
```

with:
- **k = 2π/λ**: Wavenumber
- **E**: Electric field (complex-valued)
- **PML**: Perfectly Matched Layer for radiation boundary conditions

### Problem Types

**1. Scattering Problem**
- External plane wave illuminates a passive object
- Solves for the scattered field
- Source: Incident wave `E_inc` from infinity

**2. Antenna Problem**
- Active object radiates electromagnetic waves
- Solves for the radiated field
- Source: Excitation field `E_source` at antenna surface

Both use the same weak formulation with different source terms.

### Weak Formulation

Using HCurl finite element spaces (Nedelec elements):

```
Bilinear form:  a(E, E') = ∫_Ω (∇×E)·(∇×E') - k²E·E' dΩ - ik ∫_Γ (n×E)·(n×E') dΓ

Linear form:    l(E') = -ik ∫_Γ (n×E_input)·(n×E') dΓ
```

where `E_input` is either:
- `E_inc` (incident wave) for scattering
- `E_source` (antenna excitation) for radiation

### Geometry

**Multiple parametric geometry types available:**

1. **Tri-axial Ellipsoid Scatterers** (`create_ellipsoid_scatterer_geometry`)
   - Three independent semi-axes (a, b, c)
   - Suitable for studying aspect ratio effects
   - Axis-aligned orientations (x, y, or z)

2. **Spheroid Scatterers** (`create_spheroid_scatterer_geometry`)
   - Ellipsoid of revolution (two equal axes)
   - Models prolate (cigar) and oblate (pancake) shapes
   - Useful for comparing with analytical solutions

3. **Cylindrical Dipole Antennas** (`create_dipole_antenna_geometry`)
   - Parametric length and radius (relative to wavelength)
   - Typical: half-wave dipole (L=0.5λ)
   - Models thin-wire antennas

4. **Spherical Geometry** (`create_spherical_geometry`)
   - Simple sphere (backward compatible)
   - Fast meshing for basic tests

**All geometries include:**
- Outer domain sphere: radius `R` (computational boundary)
- PML layer: thickness `PMLw` at outer boundary
- Automatic mesh sizing: h_max ≈ wavelength/15 (configurable)

## Solver Capabilities

### Direct Solvers
- Sparse Cholesky factorization
- Suitable for problems up to ~100k DOFs

### Iterative Solvers
- GMRes with multiple preconditioners:
  - Block Jacobi (default, robust)
  - BDDC (domain decomposition)
  - Hiptmair-Xu AMG (HCurl-specific)
- Suitable for large problems (>100k DOFs)

## Parametric Studies

Define parameter sweeps in YAML. The system automatically detects both problem type and geometry type based on parameters.

**Ellipsoid Scattering Study:**
```yaml
sweep_name: "ellipsoid_aspect_ratio"

parameters:
  wavelength:
    values: [0.65]

  ellipsoid_semi_axis_a:  # ← Auto-detects ellipsoid geometry
    values: [0.05, 0.08, 0.10]

  ellipsoid_semi_axis_b:
    values: [0.05]

  ellipsoid_semi_axis_c:
    values: [0.05, 0.10, 0.15, 0.20]

  propagation_dir:  # ← Auto-detects scattering problem
    values:
      - [0, 0, 1]
      - [1, 0, 0]

  polarization:
    values:
      - [1, 0, 0]

# Total: 1 × 3 × 1 × 4 × 2 × 1 = 24 jobs
```

**Dipole Antenna Study:**
```yaml
sweep_name: "dipole_length_study"

parameters:
  wavelength:
    values: [0.65]

  dipole_length_factor:  # ← Auto-detects dipole geometry
    values: [0.3, 0.4, 0.5, 0.6, 0.7]  # L/λ ratios

  dipole_radius_factor:
    values: [0.01]

  amplitude:  # ← Auto-detects antenna problem
    values: [1.0]

  polarization:
    values:
      - [0, 0, 1]  # z-polarized (along dipole axis)

# Total: 1 × 5 × 1 × 1 × 1 = 5 jobs
```

**Spheroid Scattering Study:**
```yaml
sweep_name: "spheroid_aspect_ratio"

parameters:
  wavelength:
    values: [0.65]

  spheroid_equatorial_radius:  # ← Auto-detects spheroid geometry
    values: [0.08]  # a = b

  spheroid_polar_radius:
    values: [0.04, 0.08, 0.12, 0.16, 0.20]  # c (varying)

  propagation_dir:  # ← Auto-detects scattering problem
    values:
      - [0, 0, 1]  # Along symmetry axis
      - [1, 0, 0]  # Perpendicular to axis

  polarization:
    values:
      - [1, 0, 0]

# Total: 1 × 1 × 5 × 2 × 1 = 10 jobs
```

The system automatically:
- Detects problem type (scattering vs antenna) from `propagation_dir` or `amplitude`
- Detects geometry type (ellipsoid, spheroid, dipole, sphere) from parameter names
- Generates all parameter combinations
- Runs simulations in parallel (local or cluster)
- Tracks all parameters, geometry type, and timings
- Aggregates results into CSV with `problem_type` and `geometry_type` fields

## Output & Results

Each simulation produces:

```
results/my_sweep/job_0000/
├── metadata.json     # Full parameter record + results
├── timings.json      # Timing breakdown
└── E_field.vtu       # E-field visualization (if save_solution: true)
```

Post-processing generates:
- `summary.csv` - All jobs in tabular format
- `summary_statistics.json` - Aggregate statistics
- Plots (optional) - Timing analysis, parameter sensitivity

### VTK Visualization

Enable VTK export in your YAML config:

```yaml
output:
  save_solution: true  # Export E-field to VTK
```

**Output:** Each job creates `E_field.vtu` with real and imaginary components of the E-field.

**Open in ParaView:**
1. File → Open → Select `E_field.vtu`
2. Apply
3. View components: `E_real`, `E_imag` (vector fields)
4. Calculate magnitude: Calculator filter with `mag(E_real)` or `mag(E_imag)`

⚠️ **Storage warning:** VTK files are large (~10-100 MB each). Only enable for jobs you need to visualize.

## Performance

Typical timing (coarse mesh, 25k DOFs, λ=0.5):
- Mesh generation: 0.3s
- System assembly: 0.01s
- GMRes solve: 2-3s (200 iterations)
- **Total: ~3 seconds per job**

Scaling:
- **Local**: 1-10 jobs on workstation
- **Cluster**: 1000s of jobs in parallel

### Wavelength and Mesh Resolution Considerations

⚠️ **Important**: Current parameter sweeps use relatively **large wavelengths** (λ = 0.5m to 2.0m) for computational feasibility. This provides good initial results but may miss interesting physics at higher frequencies.

**For more physically interesting results** (smaller wavelengths, higher frequencies):
- **Requirement**: Maintain mesh resolution of **10-15 elements per wavelength** (h_max ≈ λ/15)
- **Trade-off**: Reducing wavelength while maintaining resolution dramatically increases computational cost
  - Smaller λ → finer mesh (h must decrease proportionally)
  - Finer mesh → more elements → more DOFs → longer solve times
  - Example: Halving λ requires ~8× more elements (3D scaling), potentially 10-20× longer solve time

**Example computational scaling:**
```
λ = 1.0m, h = λ/15 ≈ 0.067m  →  ~100k elements  →  ~20s solve
λ = 0.5m, h = λ/15 ≈ 0.033m  →  ~800k elements  →  ~160s solve
λ = 0.1m, h = λ/15 ≈ 0.007m  →  ~100M elements  →  hours/days (impractical)
```

**Recommendations:**
1. **Pilot studies**: Test smaller wavelengths on a few cases to estimate computational cost
2. **Adaptive approach**: Start with coarse mesh (h = λ/10) to verify physics, then refine
3. **Frequency regime selection**: Choose wavelengths based on physics of interest vs. available resources
4. **Consider problem size**: Scatterer/antenna electrical size (ka, L/λ) matters more than absolute wavelength

## Testing

Run all tests with a single command:

```bash
python3 run_all_tests.py
```

**Output:**
```
======================================================================
MaxwellScattererAndAntennaNGSolve Test Suite
======================================================================
Running 5 test suites...

✅ PASS  tests/test_refactored.py
✅ PASS  tests/test_occ_syntax.py
✅ PASS  tests/test_ellipsoid.py
✅ PASS  tests/test_final_syntax.py
✅ PASS  tests/test_geometry_generators.py

✅ All tests passed!
```

**Individual test suites:**
- `tests/test_refactored.py` - Core module imports and structure
- `tests/test_occ_syntax.py` - OpenCascade primitive verification
- `tests/test_ellipsoid.py` - Ellipsoid geometry creation
- `tests/test_final_syntax.py` - Complete geometry workflows
- `tests/test_geometry_generators.py` - All parametric generators (10 tests)

**Verbose output:**
```bash
python3 run_all_tests.py --verbose
```

## Examples

### Basic Scattering Simulation

```bash
./run_scatterer.py --wavelength 0.5 --prop-dir 0 0 1 --polarization 1 0 0
```

### Basic Antenna Simulation

```bash
./run_antenna.py --wavelength 0.5 --polarization 1 0 0 --amplitude 2.0
```

### Scattering Wavelength Sweep (Cluster)

Study scattering as function of wavelength:

```bash
cd cluster
python3 scripts/generate_jobs.py configs/test_local.yaml
python3 scripts/run_local_test.py test_local
```

### Antenna Parametric Sweep (Cluster)

Study antenna radiation patterns:

```bash
cd cluster
python3 scripts/generate_jobs.py configs/antenna_test_local.yaml
sbatch --array=0-1 slurm/array_job.slurm antenna_test_local
```

### Large-Scale Cluster Run

Submit 1000 jobs with rate limiting:

```bash
sbatch --array=0-999%50 slurm/array_job.slurm my_big_sweep
```

(Runs max 50 jobs concurrently)

## Development Workflow

1. **Develop** - Use Jupyter notebook or Python scripts locally
2. **Test** - Run `test_quick.yaml` config locally (~3s)
3. **Validate** - Run small cluster test (10-20 jobs)
4. **Deploy** - Submit full parameter sweep
5. **Analyze** - Post-process results, generate plots

## Dependencies

### Core
- **ngsolve** >= 6.2 (FEM solver)
- **netgen** (mesh generation)
- **numpy** (numerical arrays)
- **scipy** (scientific computing)

### Cluster Workflow
- **pyyaml** (config files)
- **pandas** (result aggregation)
- **matplotlib** (optional, plotting)

### Environment
- Python >= 3.8
- SLURM (for cluster execution)

## Large-Scale Parameter Sweeps

For comprehensive parameter exploration (1000+ cases), use the automated sweep generator:

```bash
# Generate all sweep configurations
python3 cluster/scripts/generate_sweep_configs.py
```

This creates:
- `scattering_1000_sweep.yaml` - 1000 scattering configurations
- `antenna_1000_sweep.yaml` - 1000 antenna configurations
- `pilot_scattering.yaml` - 9 validation cases
- `pilot_antenna.yaml` - 12 validation cases

**Always run pilot studies first:**
```bash
# Local pilot test
python3 cluster/scripts/run_simulation.py pilot_scattering 0

# Cluster pilot (21 jobs total)
sbatch --array=0-8 cluster/slurm/array_job.slurm pilot_scattering
sbatch --array=0-11 cluster/slurm/array_job.slurm pilot_antenna
```

**After validation, run full sweeps:**
```bash
# 1000 scattering cases
sbatch --array=0-999 cluster/slurm/array_job.slurm scattering_1000_sweep

# 1000 antenna cases (with rate limiting)
sbatch --array=0-999%50 cluster/slurm/array_job.slurm antenna_1000_sweep
```

See `cluster/PARAMETER_SWEEPS.md` for complete documentation on parameter selection strategy, validation methodology, and expected physical coverage.

## Documentation

- `README.md` (this file) - Project overview and architecture
- `cluster/README.md` - Detailed cluster documentation
- `cluster/QUICK_START.md` - Quick reference guide
- `cluster/PARAMETER_SWEEPS.md` - Large-scale parameter sweep guide
- `cluster/scripts/README.md` - Script documentation
- `cluster/slurm/README.md` - SLURM job array documentation

## References

- **NGSolve**: https://ngsolve.org
- **Netgen**: https://ngsolve.org/netgen
- **Maxwell equations**: Time-harmonic formulation
- **PML theory**: Perfectly Matched Layers for electromagnetic problems
- **Nedelec elements**: HCurl-conforming finite elements

## License

This project uses NGSolve (LGPL) and Netgen.

## Authors

Developed by Camilo & Maximilian for electromagnetic scattering research.

## Citation

If you use this code in research, please cite:
- NGSolve project
- This repository (if made public)
