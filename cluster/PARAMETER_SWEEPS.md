# Parameter Sweep Documentation

Comprehensive guide for 1000-case parameter sweeps exploring diverse electromagnetic scattering and antenna configurations.

## Overview

Two major parameter sweeps have been designed to systematically explore:

1. **Scattering Problems**: 1000 configurations across different geometries, electrical sizes, and incident wave directions
2. **Antenna Problems**: 1000 configurations across different dipole dimensions, resonance conditions, and excitations

Both sweeps are designed with:
- Computational feasibility (20-60s per job, manageable mesh sizes)
- Physical diversity (covers key dimensionless parameter regimes)
- Wavelength range: 0.5m to 2.0m (avoids mesh resolution issues)
- Mesh resolution: h_max = λ/15 (10-15 elements per wavelength)

## Quick Start

### 1. Generate Configuration Files

```bash
# Generate all sweep configurations
python3 cluster/scripts/generate_sweep_configs.py
```

This creates:
- `scattering_1000_sweep.yaml` - 1000 scattering cases
- `antenna_1000_sweep.yaml` - 1000 antenna cases
- `pilot_scattering.yaml` - 9 validation cases
- `pilot_antenna.yaml` - 12 validation cases

### 2. Run Pilot Study (Recommended First Step)

```bash
# Local test - single pilot case
python3 cluster/scripts/run_simulation.py pilot_scattering 0

# Cluster - all scattering pilots (9 jobs)
sbatch --array=0-8 cluster/slurm/array_job.slurm pilot_scattering

# Cluster - all antenna pilots (12 jobs)
sbatch --array=0-11 cluster/slurm/array_job.slurm pilot_antenna
```

### 3. Validate Results

Check pilot study results for:
- Mesh quality (elements, DOFs)
- Solver convergence
- Runtime estimates
- Physical reasonableness

### 4. Run Full Sweep (After Validation)

```bash
# Scattering sweep (1000 jobs)
sbatch --array=0-999 cluster/slurm/array_job.slurm scattering_1000_sweep

# Antenna sweep (1000 jobs)
sbatch --array=0-999 cluster/slurm/antenna_1000_sweep

# With rate limiting (max 50 concurrent jobs)
sbatch --array=0-999%50 cluster/slurm/array_job.slurm scattering_1000_sweep
```

---

## Scattering Problems (1000 Cases)

### Strategy

**5 wavelengths × 40 geometries × 5 wave configurations = 1000 cases**

### Wavelength Sweep

```python
wavelengths = [0.5, 0.75, 1.0, 1.5, 2.0]  # meters
```

**Rationale**: Wavelength range chosen to:
- Avoid tiny meshes (λ too small → excessive refinement)
- Keep domain size manageable (R = 2λ)
- Cover range of electrical sizes when combined with geometry sweep

### Geometry Configurations (40 Total)

#### Spheres (5 configurations)

```python
r/λ = [0.1, 0.15, 0.2, 0.25, 0.3]
```

**Electrical size range**: ka ∈ [0.63, 1.88]

Covers:
- Rayleigh scattering (ka < 1): r/λ = 0.1
- Transition region (ka ~ 1): r/λ = 0.15-0.2
- Mie scattering (ka > 1): r/λ = 0.25-0.3

#### Prolate Spheroids - "Cigar" Shapes (10 configurations)

```python
(equatorial_radius/λ, polar_radius/λ) pairs:
(0.1, 0.3), (0.1, 0.4), (0.1, 0.5),
(0.15, 0.3), (0.15, 0.4), (0.15, 0.5),
(0.2, 0.3), (0.2, 0.4), (0.2, 0.5),
(0.15, 0.45)
```

**Aspect ratios**: ~1:3 to 1:5 (elongated along polar axis)

**Physical interest**:
- Orientation-dependent scattering
- Resonances at different ka for different orientations
- Polarization sensitivity

#### Oblate Spheroids - "Pancake" Shapes (10 configurations)

```python
(equatorial_radius/λ, polar_radius/λ) pairs:
(0.3, 0.1), (0.3, 0.15), (0.3, 0.2),
(0.4, 0.1), (0.4, 0.15), (0.4, 0.2),
(0.5, 0.1), (0.5, 0.15), (0.5, 0.2),
(0.35, 0.15)
```

**Aspect ratios**: ~3:1 to 5:1 (flattened)

**Physical interest**:
- Complementary to prolate (tests aspect ratio effects)
- Different resonance structure
- Backscattering enhancement for certain orientations

#### Tri-axial Ellipsoids (15 configurations)

```python
(a/λ, b/λ, c/λ) triplets with aspect ratio families:

Moderate (1:1.5:2):
  (0.1, 0.15, 0.2), (0.15, 0.225, 0.3), (0.2, 0.3, 0.4)

Strong (1:2:3):
  (0.1, 0.2, 0.3), (0.12, 0.24, 0.36), (0.15, 0.3, 0.45)

Mixed (1:1.5:3):
  (0.1, 0.15, 0.3), (0.12, 0.18, 0.36), (0.15, 0.225, 0.45)

Biaxial variations:
  (0.1, 0.1, 0.3), (0.3, 0.3, 0.1), (0.15, 0.2, 0.25), ...
```

**Physical interest**:
- Most general case (no symmetry axes)
- Three principal polarizabilities
- Complex resonance structure

### Wave Configurations (5 per geometry)

**Incident directions**:
1. `(0, 0, 1)` - propagation along z-axis
2. `(1, 0, 0)` - propagation along x-axis
3. `(0, 1, 0)` - propagation along y-axis
4. `(0.7071, 0.7071, 0)` - diagonal in xy-plane (45°)
5. `(0.5774, 0.5774, 0.5774)` - 3D diagonal (equal angles)

**Polarization**: Perpendicular to propagation direction (standard linear polarization)

**Physical coverage**:
- Normal incidence along principal axes (1-3)
- Oblique incidence (4-5)
- Tests angular dependence of scattering

### Dimensionless Parameters Explored

| Parameter | Range | Physical Meaning |
|-----------|-------|------------------|
| ka | 0.3 - 3.0 | Electrical size (2πr/λ) |
| Aspect ratio | 1:1 - 5:1 | Shape anisotropy |
| Incident angle | 0° - 90° | Directional dependence |

### Expected Physical Diversity

- **Rayleigh regime** (ka < 1): Dipole scattering, σ ∝ k⁴a⁶
- **Transition regime** (ka ~ 1-2): Interference effects, complex patterns
- **Mie regime** (ka > 2): Resonances, forward scattering dominance
- **Shape effects**: Orientation-dependent RCS, polarization coupling
- **Resonances**: Size/shape-dependent resonance peaks

---

## Antenna Problems (1000 Cases)

### Strategy

**5 wavelengths × 10 lengths × 4 radii × 5 excitation/orientation = 1000 cases**

### Wavelength Sweep

```python
wavelengths = [0.5, 0.75, 1.0, 1.5, 2.0]  # meters
```

Same rationale as scattering problems.

### Dipole Lengths (10 values)

```python
L/λ = [0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0]
```

**Physical regimes**:

| L/λ | Regime | Impedance | Pattern |
|-----|--------|-----------|---------|
| 0.3 | Short dipole | Capacitive | Nearly isotropic |
| 0.4 | Sub-resonant | Capacitive | Broadening null |
| **0.5** | **Half-wave resonance** | **~73Ω real** | **Classic dipole** |
| 0.6 | Near-resonant | Inductive | Narrowing beam |
| 0.8 | Super-resonant | Inductive | Directive |
| **1.0** | **Full-wave** | **High Z** | **Two lobes** |
| 1.2 | Long dipole | High Z | Multiple lobes |
| 1.5 | 3λ/2 resonance | Moderate Z | Three lobes |
| 1.8-2.0 | Very long | High Z | Multiple lobes |

### Dipole Radii (4 values)

```python
r/λ = [0.005, 0.01, 0.015, 0.02]
```

**Aspect ratios**:
- Minimum: L/r = 15 (L=0.3λ, r=0.02λ) - moderately thick
- Maximum: L/r = 400 (L=2.0λ, r=0.005λ) - very thin wire

**Physical effects**:
- Thin (L/r > 100): Narrow bandwidth, sharp resonance
- Thick (L/r < 50): Wide bandwidth, resistive losses more significant
- Affects input impedance and current distribution

### Excitation and Orientation (5 configurations)

```python
1. z-orientation, I₀ = 1.0, center feed
2. z-orientation, I₀ = 2.0, center feed
3. x-orientation, I₀ = 1.0, center feed
4. y-orientation, I₀ = 1.0, center feed
5. z-orientation, I₀ = 1.0, offset feed (0.25L from center)
```

**Physical diversity**:
- **Orientation** (z/x/y): Different radiation patterns relative to coordinate system
- **Amplitude** (1.0 vs 2.0): Tests linearity, power scaling
- **Feed position**: Center vs offset changes current distribution

### Dimensionless Parameters Explored

| Parameter | Range | Physical Meaning |
|-----------|-------|------------------|
| L/λ | 0.3 - 2.0 | Electrical length |
| r/λ | 0.005 - 0.02 | Electrical radius |
| L/r | 15 - 400 | Aspect ratio (thin wire parameter) |
| ka | 0.05 - 0.25 | Wire electrical radius |

### Expected Physical Diversity

- **Impedance characteristics**: Capacitive → resonant → inductive
- **Radiation patterns**: Isotropic → dipole → multi-lobe
- **Current distribution**: Sinusoidal (resonant) vs exponential (off-resonant)
- **Bandwidth**: Thin (narrow) vs thick (wide)
- **Directivity**: Low (short) → high (resonant) → very high (long arrays)

---

## Pilot Study Configurations

### Purpose

Validate computational approach before running full 1000-case sweeps:

1. **Mesh convergence**: Verify h_max = λ/15 is sufficient
2. **Runtime estimates**: Confirm 20-60s per job predictions
3. **Physics validation**: Check against analytical solutions
4. **Solver stability**: Ensure convergence for all cases

### Scattering Pilots (9 cases)

Representative cases from each electrical size regime:

**Small sphere (ka ~ 0.5)**: 3 incident directions
- Tests Rayleigh limit (σ ∝ k⁴a⁶)
- Should see isotropic scattering

**Medium ellipsoid (ka ~ 1.5)**: 2 incident directions
- Transition regime
- Directional dependence from asymmetry

**Large spheroid (ka ~ 2.5)**: 4 cases (prolate + oblate, 2 directions each)
- Mie regime
- Strong shape/orientation effects

All pilots use:
- λ = 1.0 m (mid-range)
- VTK output enabled (for visualization)
- Direct solver

### Antenna Pilots (12 cases)

Systematic exploration of resonance conditions:

**Short dipole (L=0.3λ)**: 2 radii
- Capacitive input impedance
- Nearly omnidirectional pattern

**Resonant half-wave (L=0.5λ)**: 4 cases
- Most important configuration
- 3 radii (tests aspect ratio effect)
- 1 x-oriented (tests orientation)

**Full-wave (L=1.0λ)**: 2 radii
- Second harmonic resonance
- Two-lobe pattern

**Long dipoles (L=1.5λ, 2.0λ)**: 4 cases
- Multi-lobe patterns
- High impedance

All pilots use:
- λ = 1.0 m
- VTK output enabled
- Direct solver

---

## Computational Estimates

### Per-Job Resources

**Mesh statistics** (typical for λ=1.0m, h_max=λ/15):

| Geometry | Elements | DOFs (order=3) | Memory |
|----------|----------|----------------|--------|
| Sphere | ~80k | ~240k | 4-6 GB |
| Ellipsoid | ~120k | ~360k | 6-8 GB |
| Dipole | ~150k | ~450k | 6-8 GB |

**Runtime estimates**:
- Mesh generation: 2-5s
- Assembly: 5-15s
- Direct solve: 10-40s
- Post-processing: 2-5s
- **Total: 20-60s per job**

### Full Sweep Resources

**1000 jobs in parallel**:
- Wall time: ~1 hour (assumes parallel execution)
- Total CPU hours: 20-60 hours
- Storage (no VTK): ~50-100 MB per job → 50-100 GB total
- Storage (with VTK): ~100-200 MB per job → 100-200 GB total

**Recommendations**:
- Use scratch filesystem for results (faster I/O)
- Disable VTK for production sweeps (enable only for select cases)
- Use rate limiting if cluster policy requires: `--array=0-999%50`

---

## Validation Methodology

### Phase 1: Pilot Study

**Goal**: Validate computational setup

**Steps**:
1. Run all pilot cases (9 scattering + 12 antenna)
2. Check mesh quality:
   - Elements per wavelength: should be 10-15
   - Aspect ratios: should be reasonable (<10:1)
3. Verify solver convergence
4. Compare runtime to estimates
5. Visualize fields in ParaView (VTK output)

**Success criteria**:
- All jobs complete without errors
- Runtimes within 20-60s range
- Fields look physically reasonable
- No mesh pathologies

### Phase 2: Convergence Study

**Goal**: Verify mesh resolution is adequate

**Method**: Select 3-5 representative cases, run with multiple mesh sizes:
- Coarse: h_max = λ/10
- Medium: h_max = λ/12
- Fine: h_max = λ/15
- Very fine: h_max = λ/20

**Metrics**:
- Scattering: RCS at θ=0° (backscattering)
- Antenna: Input impedance, directivity

**Success criteria**:
- < 1% change from h_max=λ/15 to h_max=λ/20
- Monotonic convergence with refinement

### Phase 3: Physics Validation

**Scattering tests**:

1. **Rayleigh limit** (small sphere, ka << 1):
   - σ_total ≈ (8π/3) k⁴a⁶ (theory)
   - Nearly isotropic scattering pattern

2. **Mie theory** (sphere, ka ~ 1-2):
   - Compare to analytical Mie series solution
   - Check resonance peak positions

3. **Reciprocity**:
   - RCS(θ, φ, pol) should satisfy reciprocity relations

**Antenna tests**:

1. **Half-wave dipole** (L = λ/2):
   - Input impedance ≈ 73 + j42.5 Ω (thin wire theory)
   - Directivity ≈ 1.64 (2.15 dBi)
   - Pattern: sin(θ) in E-plane

2. **Full-wave dipole** (L = λ):
   - High input impedance
   - Two-lobe pattern with null at θ=90°

3. **Power conservation**:
   - Radiated power = Input power (no losses in PEC model)

### Phase 4: Dimensionless Parameter Check

**Goal**: Verify solutions depend only on dimensionless parameters

**Test**: Run pairs of cases with:
- Different λ but same ka (scattering)
- Different λ but same L/λ and r/λ (antenna)

**Success criteria**:
- Normalized quantities (RCS/πa², input impedance/Z₀) should be identical
- Confirms proper scaling

---

## Usage Examples

### Local Testing

```bash
# Test single pilot case
python3 cluster/scripts/run_simulation.py pilot_scattering 0

# Test specific job from full sweep
python3 cluster/scripts/run_simulation.py scattering_1000_sweep 42
```

### Cluster Submission

```bash
# Pilot studies
sbatch --array=0-8 cluster/slurm/array_job.slurm pilot_scattering
sbatch --array=0-11 cluster/slurm/array_job.slurm pilot_antenna

# Full sweeps (after validation!)
sbatch --array=0-999 cluster/slurm/array_job.slurm scattering_1000_sweep
sbatch --array=0-999 cluster/slurm/array_job.slurm antenna_1000_sweep

# With rate limiting (max 50 concurrent)
sbatch --array=0-999%50 cluster/slurm/array_job.slurm scattering_1000_sweep

# Select specific job range
sbatch --array=0-99 cluster/slurm/array_job.slurm scattering_1000_sweep  # First 100
sbatch --array=500-599 cluster/slurm/array_job.slurm antenna_1000_sweep  # Jobs 500-599
```

### Monitoring Progress

```bash
# Check queue status
squeue -u $USER

# Count completed jobs
ls cluster/results/scattering_1000_sweep/job_*/metadata.json | wc -l

# Check specific job output
cat logs/scattering_1000_sweep_0.out

# Find failed jobs
for i in {0..999}; do
  if [ ! -f "cluster/results/scattering_1000_sweep/job_$(printf "%04d" $i)/metadata.json" ]; then
    echo "Job $i failed or incomplete"
  fi
done
```

### Post-Processing

```bash
# Collect all results into summary
python3 cluster/scripts/collect_results.py scattering_1000_sweep

# Extract RCS patterns
python3 cluster/scripts/extract_rcs.py scattering_1000_sweep

# Extract antenna impedances
python3 cluster/scripts/extract_impedance.py antenna_1000_sweep

# Visualize parameter space coverage
python3 cluster/scripts/plot_parameter_coverage.py scattering_1000_sweep
```

---

## File Structure

```
cluster/
├── configs/
│   ├── scattering_1000_sweep.yaml    # 1000 scattering cases
│   ├── antenna_1000_sweep.yaml       # 1000 antenna cases
│   ├── pilot_scattering.yaml         # 9 validation cases
│   └── pilot_antenna.yaml            # 12 validation cases
│
├── scripts/
│   ├── generate_sweep_configs.py     # Creates sweep YAML files
│   ├── run_simulation.py             # Single job runner
│   ├── collect_results.py            # Post-processing (TODO)
│   ├── extract_rcs.py                # RCS extraction (TODO)
│   └── extract_impedance.py          # Impedance extraction (TODO)
│
├── slurm/
│   └── array_job.slurm               # SLURM job array script
│
└── results/
    ├── scattering_1000_sweep/
    │   ├── job_0000/
    │   │   ├── metadata.json
    │   │   ├── convergence.json
    │   │   └── E_field.vtu (if enabled)
    │   ├── job_0001/
    │   └── ...
    └── antenna_1000_sweep/
        └── ...
```

---

## Configuration File Format

Each sweep configuration YAML contains:

```yaml
metadata:
  total_jobs: 1000
  problem_type: "scattering"  # or "antenna"
  description: "1000-case parameter sweep"

jobs:
  - job_id: 0
    problem_type: "scattering"
    wavelength: 0.5
    wavenumber: 12.566370614359172
    geometry_type: "sphere"
    sphere_radius: 0.05  # 0.1 * wavelength
    incident_direction: [0, 0, 1]
    incident_polarization: [1, 0, 0]
    domain_radius: 1.0
    pml_width: 0.125
    max_mesh_size: 0.03333333333333333
    fes_order: 3
    solver_method: "direct"
    save_solution: false
    geometry_label: "sphere_r0.10λ"
    wave_label: "z_prop"

  - job_id: 1
    ...
```

The `run_simulation.py` script reads this format and extracts the configuration for the specified job_id.

---

## Troubleshooting

### Job Fails with "Mesh too large"

**Cause**: Wavelength too small or geometry too detailed

**Fix**: Increase mesh size factor or wavelength:
```python
max_mesh_size: wavelength / 10.0  # Instead of /15
```

### Solver Doesn't Converge

**Cause**: Poor mesh quality or too-thin geometry

**Fix**:
- Check mesh aspect ratios
- Increase dipole radius (for antennas)
- Use direct solver instead of iterative

### Runtime Exceeds 1 Hour

**Cause**: DOF count too high

**Fix**:
- Reduce FES order: `fes_order: 2` instead of 3
- Coarsen mesh: `max_mesh_size: wavelength / 10`
- Use smaller domain: `domain_radius: 1.5 * wavelength`

### Storage Quota Exceeded

**Cause**: VTK files are large (~50-200 MB each)

**Fix**:
- Disable VTK for production: `save_solution: false`
- Use scratch filesystem: `export ANTENNA_RESULTS_BASE=/scratch/$USER/antenna`
- Clean up old results

---

## References

**Scattering Theory**:
- Bohren & Huffman, "Absorption and Scattering of Light by Small Particles"
- Mie scattering: Analytical solutions for spheres
- Rayleigh limit: ka << 1 approximation

**Antenna Theory**:
- Balanis, "Antenna Theory: Analysis and Design"
- Half-wave dipole: Fundamental resonance at L=λ/2
- Thin wire approximation: L/r >> 1

**Numerical Methods**:
- Jin, "The Finite Element Method in Electromagnetics"
- PML boundary conditions: Perfectly Matched Layers
- NGSolve documentation: https://docu.ngsolve.org

---

## Next Steps

1. ✅ Generate sweep configurations
2. ⏳ Run pilot studies (9 scattering + 12 antenna)
3. ⏳ Validate convergence and physics
4. ⏳ Implement post-processing scripts
5. ⏳ Deploy full 1000-case sweeps
6. ⏳ Analyze results and generate visualizations

**Do not proceed with full sweeps until pilot validation is complete!**
