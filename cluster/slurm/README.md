# SLURM Scripts

SLURM batch scripts for HPC cluster execution. Supports all geometry types and problem types through unified workflow.

## array_job.slurm

Main job array script for parametric sweeps.

### Usage

```bash
sbatch --array=0-N slurm/array_job.slurm sweep_name
```

Where:
- `N` = number of jobs - 1 (e.g., 0-47 for 48 jobs)
- `sweep_name` = name of the sweep (matches config filename without .yaml)

### What it does

1. **Environment Setup**
   - Loads required modules (Python, NGSolve)
   - Activates Python virtual environment (if configured)
   - Sets thread count from SLURM allocation

2. **Execution**
   - Runs `scripts/run_simulation.py` for array task ID
   - Automatically detects geometry type from parameters
   - Automatically detects problem type from parameters

3. **Output**
   - Job results: `results/sweep_name/job_XXXX/`
   - Logs: `logs/sweep_name_TASKID.out` (stdout + stderr)

### Example Commands

**Ellipsoid Scattering (48 jobs):**
```bash
sbatch --array=0-47 slurm/array_job.slurm ellipsoid_geometry_sweep
```

**Dipole Antenna (10 jobs):**
```bash
sbatch --array=0-9 slurm/array_job.slurm dipole_geometry_sweep
```

**Spheroid Scattering (10 jobs):**
```bash
sbatch --array=0-9 slurm/array_job.slurm spheroid_geometry_sweep
```

**Rate Limiting (max 50 concurrent):**
```bash
sbatch --array=0-999%50 slurm/array_job.slurm large_sweep
```

### Monitoring

```bash
# Check queue status
squeue -u $USER

# Count completed jobs
ls results/sweep_name/job_*/metadata.json | wc -l

# Check specific job log
cat logs/sweep_name_0.out

# Cancel all jobs in array
scancel <job_id>
```

### Resource Configuration

Default resources (adjust in array_job.slurm):
- **CPUs**: 4 cores per job
- **Memory**: 8GB per job
- **Time**: 1 hour per job
- **Partition**: compute (adjust for your cluster)

Typical timing estimates:
- Ellipsoid (λ=0.65, h_max=λ/15, order=3): ~20s per job
- Dipole (λ=0.65, h_max=λ/15, order=3): ~60s per job
- Spheroid (λ=0.65, h_max=λ/15, order=3): ~15s per job
- Quick test (λ=2.0, coarse mesh, order=2): ~2s per job

### Geometry Support

The script works with all geometry types automatically:

| Geometry | Detection | Example Config |
|----------|-----------|----------------|
| Tri-axial ellipsoid | `ellipsoid_semi_axis_a` | `ellipsoid_geometry_sweep.yaml` |
| Spheroid | `spheroid_equatorial_radius` | `spheroid_geometry_sweep.yaml` |
| Dipole | `dipole_length_factor` | `dipole_geometry_sweep.yaml` |
| Sphere | (default) | `template.yaml` |

No changes to SLURM script needed for different geometries - detection is automatic!

### Custom Results Directory

**Default behavior:** Results saved to `cluster/results/sweep_name/`

**Custom directory (recommended for cluster scratch):**

```bash
# Set environment variable
export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna_results

# Submit jobs (automatically uses custom path)
sbatch --array=0-47 slurm/array_job.slurm ellipsoid_sweep

# Results will be in: /scratch/$USER/antenna_results/ellipsoid_sweep/
```

**Persistent setup:** Add to `~/.bashrc` or job submission wrapper:

```bash
# In ~/.bashrc
export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna_results
```

**Advantages:**
- Faster I/O on cluster scratch filesystems
- Better quota management
- Separate results from code repository
- Easy cleanup after experiments

### VTK Output

Enable VTK field export in YAML config:

```yaml
output:
  save_solution: true
```

**Files created:**
- `job_XXXX/E_field.vtu` - E-field visualization (~10-100 MB)
- `job_XXXX/metadata.json` - Includes `"vtk_output": ["E_field.vtu"]` field

**Typical file sizes:**
- Coarse mesh (h_max = λ/5): ~10 MB per job
- Fine mesh (h_max = λ/15): ~40-80 MB per job

⚠️ **Be selective:** 48 jobs × 50 MB = 2.4 GB. Only enable for visualization jobs.

### Cluster-Specific Setup

Modify `array_job.slurm` for your cluster:

1. **Module loads**: Adjust for your environment
2. **Virtual environment**: Set path to your venv
3. **Partition**: Use your cluster's compute partition name
4. **Resources**: Adjust CPU/memory/time based on your job size
5. **Results path**: Optionally set `ANTENNA_RESULTS_BASE` in the script
