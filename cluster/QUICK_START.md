# Quick Start Guide

## First-Time Setup

**Before running any jobs, set up the environment:**

```bash
bash cluster/scripts/setup_environment.sh
```

This creates a virtual environment with all dependencies. See `cluster/README.md` "Environment Setup" section for details.

## Unified Architecture

MaxwellScattererAndAntennaNGSolve supports **both scattering and antenna problems with multiple geometry types** through a unified cluster workflow. The system automatically detects problem type and geometry type based on config parameters:

**Problem Type Detection:**
- **Scattering**: Config has `propagation_dir` parameter
- **Antenna**: Config has `amplitude` parameter

**Geometry Type Detection:**
- **Ellipsoid**: Config has `ellipsoid_semi_axis_a` parameter
- **Spheroid**: Config has `spheroid_equatorial_radius` parameter
- **Dipole**: Config has `dipole_length_factor` parameter
- **Sphere**: Default fallback (backward compatible)

## Local Testing (No SLURM Required)

Test the workflow on your local machine before deploying to the cluster:

**Ellipsoid Scattering Test:**
```bash
cd cluster

# 1. Generate job list from test config
python3 scripts/generate_jobs.py configs/test_quick_ellipsoid.yaml

# 2. Run test locally (1 job, ~2 seconds)
cd ..
python3 cluster/scripts/run_simulation.py \
    --job-file cluster/configs/test_quick_ellipsoid_jobs.json \
    --job-id 0 \
    --output-dir cluster/results/test_quick_ellipsoid \
    --num-threads 4 \
    --local

# 3. Check results
cat cluster/results/test_quick_ellipsoid/job_0000/metadata.json
```

**Spheroid Scattering Test:**
```bash
cd cluster
python3 scripts/generate_jobs.py configs/test_ellipsoid_local.yaml
cd ..
python3 cluster/scripts/run_simulation.py \
    --job-file cluster/configs/test_ellipsoid_local_jobs.json \
    --job-id 0 \
    --output-dir cluster/results/test_ellipsoid_local \
    --num-threads 4 \
    --local
```

**Legacy Scattering Test (Spherical):**
```bash
cd cluster
python3 scripts/generate_jobs.py configs/test_local.yaml
python3 scripts/run_local_test.py test_local
cat results/test_local/summary.csv
```

This runs quick simulations locally to verify everything works for all geometry types.

## Cluster Deployment

Once local testing passes, deploy to the cluster:

**Ellipsoid Scattering Sweep:**
```bash
# 1. Use provided ellipsoid config (or create your own)
# Edit configs/ellipsoid_geometry_sweep.yaml if needed

# 2. Generate full job list
python3 scripts/generate_jobs.py configs/ellipsoid_geometry_sweep.yaml

# 3. Submit to SLURM (48 jobs)
sbatch --array=0-47 slurm/array_job.slurm ellipsoid_geometry_sweep

# 4. Monitor progress
squeue -u $USER
ls results/ellipsoid_geometry_sweep/job_*/metadata.json | wc -l

# 5. Post-process when complete
python3 scripts/post_process.py results/ellipsoid_geometry_sweep
```

**Dipole Antenna Sweep:**
```bash
# Use provided dipole config
python3 scripts/generate_jobs.py configs/dipole_geometry_sweep.yaml
sbatch --array=0-9 slurm/array_job.slurm dipole_geometry_sweep
python3 scripts/post_process.py results/dipole_geometry_sweep
```

**Spheroid Scattering Sweep:**
```bash
# Use provided spheroid config
python3 scripts/generate_jobs.py configs/spheroid_geometry_sweep.yaml
sbatch --array=0-9 slurm/array_job.slurm spheroid_geometry_sweep
python3 scripts/post_process.py results/spheroid_geometry_sweep
```

**Legacy Spherical Sweep:**
```bash
# Backward compatible with old configs
cp configs/template.yaml configs/my_sweep.yaml
python3 scripts/generate_jobs.py configs/my_sweep.yaml
sbatch --array=0-N slurm/array_job.slurm my_sweep
python3 scripts/post_process.py results/my_sweep
```

## Key Files

**Geometry Test Configs:**
- **configs/test_quick_ellipsoid.yaml** - Fast ellipsoid test (1 job, ~2s, large wavelength)
- **configs/test_ellipsoid_local.yaml** - Standard ellipsoid test (1 job, ~20s)
- **configs/test_local.yaml** - Legacy spherical test (2 jobs, coarse mesh)

**Production Sweep Configs:**
- **configs/ellipsoid_geometry_sweep.yaml** - Tri-axial ellipsoid parametric sweep (48 jobs)
- **configs/dipole_geometry_sweep.yaml** - Dipole antenna length study (10 jobs)
- **configs/spheroid_geometry_sweep.yaml** - Spheroid aspect ratio study (10 jobs)
- **configs/template.yaml** - Legacy template for spherical sweeps

**Scripts (unified for all geometries and problem types):**
- **scripts/generate_jobs.py** - Generate job list from config
- **scripts/run_simulation.py** - Run single job (auto-detects geometry & problem type)
- **scripts/post_process.py** - Aggregate results into CSV
- **scripts/run_local_test.py** - Local testing helper
- **slurm/array_job.slurm** - SLURM job array script

## Workflow Summary

```
1. CONFIG   → Write YAML with parameter sweep
2. GENERATE → Create job list JSON
3. TEST     → Run locally (run_local_test.py)
4. DEPLOY   → Submit to SLURM (sbatch)
5. ANALYZE  → Post-process results (post_process.py)
```

## Example: Run 1 Local Job Manually

```bash
cd cluster
python3 scripts/generate_jobs.py configs/test_local.yaml
python3 scripts/run_simulation.py \
    --job-file configs/test_local_jobs.json \
    --job-id 0 \
    --output-dir results/test_local \
    --num-threads 4 \
    --local
```

## VTK Visualization

Enable field visualization in your YAML config:

```yaml
output:
  save_solution: true  # Export E-field to VTK format
```

**What you get:**
- `job_XXXX/E_field.vtu` - VTK file with E-field data
- Contains `E_real` and `E_imag` vector components
- Open directly in ParaView for visualization

**ParaView workflow:**
```bash
# After job completes
paraview cluster/results/sweep_name/job_0000/E_field.vtu

# In ParaView:
# - Apply to load data
# - View E_real and E_imag in dropdown
# - Use Calculator filter: mag(E_real) for magnitude
# - Use Glyph filter for vector arrows
```

⚠️ **Storage:** VTK files are ~10-100 MB each. Only enable for jobs you want to visualize!

## Custom Results Directory

**Default location:** `cluster/results/`

**Custom location (recommended for cluster):**

```bash
# Option 1: Set environment variable before submitting
export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna_results
sbatch --array=0-47 slurm/array_job.slurm ellipsoid_sweep

# Option 2: Set in ~/.bashrc or job submission script
echo "export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna_results" >> ~/.bashrc

# Post-process with custom path
python3 scripts/post_process.py /scratch/$USER/antenna_results/ellipsoid_sweep
```

**Why use custom path?**
- Scratch filesystems have faster I/O
- Better quota management
- Separate from code repository
- Shared project directories

## Tips

- **Start small**: Test with `configs/test_quick_ellipsoid.yaml` (~2s) before running production sweeps
- **Use --validate-only**: Check config without running
- **Monitor resources**: Check logs/ for memory/time usage
- **Local first**: Always test locally before cluster deployment
- **Auto-detection works for both**:
  - **Problem type**: `propagation_dir` → Scattering, `amplitude` → Antenna
  - **Geometry type**: Parameter names automatically select the right geometry generator
    - `ellipsoid_semi_axis_a` → Tri-axial ellipsoid
    - `spheroid_equatorial_radius` → Spheroid (ellipsoid of revolution)
    - `dipole_length_factor` → Cylindrical dipole
    - None of the above → Spherical (backward compatible)

## Large-Scale Parameter Sweeps

For comprehensive parameter exploration (1000+ cases), use the parameter sweep system:

### Generate Sweep Configurations

```bash
# Create all sweep configs (scattering + antenna + pilots)
python3 scripts/generate_sweep_configs.py
```

This creates:
- `scattering_1000_sweep.yaml` - 1000 scattering cases
- `antenna_1000_sweep.yaml` - 1000 antenna cases
- `pilot_scattering.yaml` - 9 validation cases
- `pilot_antenna.yaml` - 12 validation cases

### Run Pilot Study First (Recommended)

Always validate with pilot study before running full sweeps:

```bash
# Local pilot test
python3 scripts/run_simulation.py pilot_scattering 0

# Cluster - all scattering pilots (9 jobs)
sbatch --array=0-8 slurm/array_job.slurm pilot_scattering

# Cluster - all antenna pilots (12 jobs)
sbatch --array=0-11 slurm/array_job.slurm pilot_antenna

# Monitor pilot progress
squeue -u $USER
ls results/pilot_scattering/job_*/metadata.json | wc -l
```

**Pilot study validates**:
- Mesh resolution is adequate (h_max = λ/15)
- Runtime estimates are accurate (20-60s per job)
- Solver converges for all cases
- Physics looks reasonable

### Run Full 1000-Case Sweeps

**Only after pilot validation passes!**

```bash
# Scattering sweep (1000 jobs)
sbatch --array=0-999 slurm/array_job.slurm scattering_1000_sweep

# Antenna sweep (1000 jobs)
sbatch --array=0-999 slurm/array_job.slurm antenna_1000_sweep

# With rate limiting (max 50 concurrent jobs)
sbatch --array=0-999%50 slurm/array_job.slurm scattering_1000_sweep

# Run specific job ranges
sbatch --array=0-99 slurm/array_job.slurm scattering_1000_sweep    # First 100
sbatch --array=500-599 slurm/array_job.slurm antenna_1000_sweep    # Jobs 500-599
```

### Monitor Large Sweeps

```bash
# Check queue status
squeue -u $USER

# Count completed jobs
ls results/scattering_1000_sweep/job_*/metadata.json | wc -l

# Find failed/incomplete jobs
for i in {0..999}; do
  if [ ! -f "results/scattering_1000_sweep/job_$(printf "%04d" $i)/metadata.json" ]; then
    echo "Job $i incomplete"
  fi
done

# Check specific job log
cat logs/scattering_1000_sweep_42.out
```

### What the Sweeps Explore

**Scattering (1000 cases)**: 5 wavelengths × 40 geometries × 5 wave configs
- Wavelengths: 0.5m to 2.0m
- Geometries: Spheres, prolate/oblate spheroids, tri-axial ellipsoids
- Electrical sizes: ka from 0.3 to 3.0 (Rayleigh → Mie regimes)
- Incident directions: 5 angles (normal + oblique)

**Antenna (1000 cases)**: 5 wavelengths × 10 lengths × 4 radii × 5 excitations
- Wavelengths: 0.5m to 2.0m
- Dipole lengths: L/λ from 0.3 to 2.0 (sub-resonant → super-resonant)
- Dipole radii: r/λ from 0.005 to 0.02 (thin wire → thick)
- Orientations: x, y, z axes
- Feed positions: center + offset

See **PARAMETER_SWEEPS.md** for complete documentation on parameter selection strategy, physical coverage, validation methodology, and post-processing.

---

See README.md and PARAMETER_SWEEPS.md for detailed documentation.
