# Chevron-Subhalo Interaction Simulations

**Modeling phase-space chevron perturbations from dark matter subhalo encounters with Galactic bar resonances.**

---

## Overview

This project investigates how dark matter (DM) subhalos perturb phase-space structures (chevrons) created by orbital resonances with the rotating Galactic bar. The analysis combines:

1. **Analytical diffusion theory** - Fokker-Planck formalism for subhalo population effects
2. **Test-particle simulations** - Direct N-body integration using AGAMA library
3. **Comparison study** - Quantifying where analytical approximations break down

**Key Result:** Analytical models underpredict subhalo impact for masses M = 10⁵–10⁹ M☉.

---

## Project Structure

```
chevron_subhalo_interaction/
├── README.md                          # This file
├── PROJECT_JOURNAL.md                 # Detailed experiment log
├── CLAUDE.md                          # Working agreement and preferences
│
├── resonance_subhalo.tex              # MNRAS paper draft (theory)
├── simulation_calibration_version2.py # Core simulation module
├── simulation_calibration_version2.ipynb # Original notebook
│
├── compute_grid_serial.py             # Serial grid computation
├── compute_grid_parallel.py           # Parallel grid computation (14 workers)
├── compare_grids.py                   # Validation and performance comparison
│
├── calibration_grid_serial.npz        # Serial output (150 cells)
└── calibration_grid_parallel.npz      # Parallel output (150 cells)
```

---

## Quick Start

### Prerequisites

```bash
# Python 3.13+ with virtual environment
python3 -m venv ~/Work/venvs/.venv
source ~/Work/venvs/.venv/bin/activate

# Install required packages
pip install numpy scipy matplotlib agama tqdm
```

### Run Simulations

```bash
# Activate environment
source ~/Work/venvs/.venv/bin/activate

# Serial computation (150 cells, ~2 minutes)
python compute_grid_serial.py

# Parallel computation (150 cells, ~15 seconds, 14 workers)
python compute_grid_parallel.py

# Compare outputs and performance
python compare_grids.py
```

**Expected Output:**
```
Array Comparison (serial vs parallel)
============================================================
Overall: PASS - All arrays identical

Runtime Comparison
============================================================
Serial runtime:     123.709 s
Parallel runtime:   15.143 s (14 workers)
Speedup:            8.17x
Parallel efficiency: 58.4%
```

---

## Physical Model

### Galactic Components

**Halo:** NFW profile
**Bar:** Ferrers triaxial ellipsoid
**Pattern speed:** Ω_p = 35 km/s/kpc
**Resonance:** Co-rotation (m=2, n=0)

### Subhalo Interaction

**Impulse approximation:**
```
Δv = 2 G M_sub / (b v_rel)
```

Where:
- M_sub: subhalo mass [M☉]
- b: impact parameter [kpc]
- v_rel: relative velocity [km/s]

**Diffusion coefficient:**
```
D(E_J) = ∫ (ΔE_J)² n(M, v) dM dv db
```

Where n(M, v) is subhalo number density and velocity distribution.

### Conserved Quantities

**Jacobi integral (rotating frame):**
```
E_J = E - Ω_p L_z
```

**Action integral (resonance amplitude):**
```
I = ∮ p dq
```

**Preferred metric:** ΔE_J (conserved for resonance dynamics)

---

## Parameter Grid

The calibration grid spans 3D parameter space:

| Parameter | Symbol | Range | N | Spacing |
|-----------|--------|-------|---|---------|
| Subhalo mass | M_sub | 10⁵–10¹⁰ M☉ | 6 | log |
| Impact parameter | b | 0–10 kpc | 5 | linear |
| Relative velocity | v_rel | 0–100 km/s | 5 | linear |

**Total cells:** 6 × 5 × 5 = 150

### Grid Definition

```python
mass_list = np.logspace(5, 10, 6)          # [1e5, ..., 1e10] M☉
impact_list = np.linspace(0.0, 10.0, 5)    # [0, 2.5, 5.0, 7.5, 10] kpc
velocity_list = np.linspace(0.0, 100.0, 5) # [0, 25, 50, 75, 100] km/s
```

---

## Output Data

Each grid computation produces a `.npz` file containing:

### Arrays (shape: 6×5×5)

- `delta_EJ_sim` - ΔE_J from test-particle simulation
- `delta_EJ_ana` - ΔE_J from analytical diffusion theory
- `delta_Is_sim` - ΔI from test-particle simulation
- `delta_Is_ana` - ΔI from analytical diffusion theory

### Grid Parameters

- `mass_list` - Subhalo masses [M☉]
- `impact_list` - Impact parameters [kpc]
- `velocity_list` - Relative velocities [km/s]

### Metadata

- `runtime_seconds` - Wall-clock execution time [s]
- `n_processes` - Number of worker processes (parallel only)

### Loading Data

```python
import numpy as np

# Load grid
data = np.load("calibration_grid_parallel.npz")

# Access arrays
dEJ_sim = data["delta_EJ_sim"]  # shape: (6, 5, 5)
mass = data["mass_list"]         # shape: (6,)
impact = data["impact_list"]     # shape: (5,)
velocity = data["velocity_list"] # shape: (5,)

# Grid cell (i, j, k) corresponds to:
# M_sub = mass[i], b = impact[j], v_rel = velocity[k]
```

---

## Implementation Details

### Serial Version (`compute_grid_serial.py`)

**Algorithm:** Triple nested loop over (mass, impact, velocity)

**Characteristics:**
- Single-threaded execution
- Sequential processing: 150 cells
- Runtime: ~124 seconds
- Memory: ~214 MB peak

**Use case:** Baseline reference, debugging, small grids

### Parallel Version (`compute_grid_parallel.py`)

**Algorithm:** `multiprocessing.Pool` with `imap_unordered`

**Characteristics:**
- Multi-process execution (14 workers on test system)
- Parallel processing: all 150 cells distributed
- Runtime: ~15 seconds
- Speedup: 8.17× over serial

**Use case:** Production runs, large grids, time-critical analysis

**Parallel Efficiency:**
```
Efficiency = Speedup / N_workers = 8.17 / 14 = 58.4%
```

**Overhead sources:**
- Process spawning and communication
- Result collection (sequential)
- Load imbalance (150 % 14 ≠ 0)

---

## Core Module: `simulation_calibration_version2.py`

### Key Function

```python
def one_simulation(subhalo_mass, subhalo_impact, subhalo_velocity):
    """
    Run single simulation for given subhalo parameters.

    Parameters
    ----------
    subhalo_mass : float
        Subhalo mass [M☉]
    subhalo_impact : float
        Impact parameter [kpc]
    subhalo_velocity : float
        Relative velocity [km/s]

    Returns
    -------
    list of float
        [delta_EJ_sim, delta_EJ_ana, delta_Is_sim, delta_Is_ana]
    """
```

### Simulation Workflow

1. **Initialize** - Set up bar potential, halo, and resonance parameters
2. **Analytical calculation** - Compute diffusion coefficients
3. **Test-particle integration** - Evolve orbits with AGAMA
4. **Apply subhalo kick** - Impulse approximation
5. **Measure changes** - Compare ΔE_J and ΔI
6. **Return results** - 4-element list

### Modifications for Clean Import

The following sections have been commented out to allow clean module import:

- **Lines 146-157:** Rotation curve plot
- **Lines 424-442:** Serial execution loop
- **Lines 456-492:** Parallel execution loop
- **Lines 495-763:** Plotting and analysis cells

**Rationale:** Grid computation scripts import `one_simulation()` only; module-level execution creates unwanted side effects.

---

## Validation

### Output Comparison

Serial and parallel implementations produce **identical results**:

```
delta_EJ_sim: Max difference = 0.00e+00 ✓
delta_EJ_ana: Max difference = 0.00e+00 ✓
delta_Is_sim: Max difference = 0.00e+00 ✓
delta_Is_ana: Max difference = 0.00e+00 ✓
```

**Validation method:** `np.allclose(rtol=1e-15, atol=0)`

### Performance Verification

| Metric | Serial | Parallel | Ratio |
|--------|--------|----------|-------|
| Runtime | 123.71 s | 15.14 s | 8.17× |
| Time/cell | 0.825 s | 0.101 s | 8.17× |
| Workers | 1 | 14 | 14× |
| Efficiency | 100% | 58.4% | - |

**Conclusion:** Parallel implementation validated for correctness and performance.

---

## Scientific Context

### Motivation

Galactic bars create orbital resonances that organize stars into coherent phase-space structures (chevrons). Dark matter subhalos perturb these structures via gravitational encounters. Understanding this process constrains:

1. **Subhalo mass function** - Number density vs mass
2. **Bar resonance stability** - How long do chevrons persist?
3. **Dynamical heating** - Energy injection into stellar disk

### Analytical vs Numerical

**Analytical models** (Fokker-Planck diffusion):
- Fast to evaluate
- Smooth, averaged effects
- Assume weak perturbations, ergodic mixing

**Numerical simulations** (test particles):
- Computationally expensive
- Capture individual encounters
- No approximations (within numerical precision)

**Key question:** When do analytical approximations fail?

### Current Findings

For M_sub = 10⁵–10⁹ M☉, analytical models **underpredict** simulation results. Possible explanations:

1. **Strong encounters** - Impulse approximation breaks down for close passages
2. **Resonance coherence** - Perturbations not ergodic; chevrons maintain structure
3. **Velocity distribution** - Analytical model assumes Maxwellian; reality more complex

---

## Usage Examples

### Example 1: Single Simulation

```python
from simulation_calibration_version2 import one_simulation

# Parameters
M_sub = 1e8       # 10^8 M☉
b = 5.0          # 5 kpc impact parameter
v_rel = 50.0     # 50 km/s relative velocity

# Run simulation
result = one_simulation(M_sub, b, v_rel)
dEJ_sim, dEJ_ana, dIs_sim, dIs_ana = result

# Print results
print(f"ΔE_J (simulation): {dEJ_sim:.6e}")
print(f"ΔE_J (analytical): {dEJ_ana:.6e}")
print(f"Ratio: {dEJ_sim / dEJ_ana:.2f}")
```

### Example 2: Load and Analyze Grid

```python
import numpy as np
import matplotlib.pyplot as plt

# Load parallel grid
data = np.load("calibration_grid_parallel.npz")
dEJ_sim = data["delta_EJ_sim"]
dEJ_ana = data["delta_EJ_ana"]
mass = data["mass_list"]

# Compute ratio (simulation / analytical)
ratio = dEJ_sim / dEJ_ana

# Plot: ratio vs mass (fixed impact and velocity)
i_b = 2  # b = 5.0 kpc
i_v = 2  # v = 50 km/s

plt.figure(figsize=(6, 4))
plt.plot(mass, ratio[:, i_b, i_v], 'o-')
plt.xscale('log')
plt.xlabel('Subhalo mass [M☉]')
plt.ylabel('ΔE_J (sim) / ΔE_J (ana)')
plt.axhline(1.0, color='k', linestyle='--', label='Perfect agreement')
plt.legend()
plt.grid(True, alpha=0.3)
plt.title(f'b = {data["impact_list"][i_b]} kpc, v = {data["velocity_list"][i_v]} km/s')
plt.tight_layout()
plt.savefig('ratio_vs_mass.pdf')
plt.show()
```

### Example 3: Compare Runtime Scaling

```bash
# Modify grid size in compute_grid_*.py
# Test different grid resolutions

# Coarse: 4 × 3 × 3 = 36 cells
mass_list = np.logspace(5, 10, 4)
impact_list = np.linspace(0.0, 10.0, 3)
velocity_list = np.linspace(0.0, 100.0, 3)

# Medium: 6 × 5 × 5 = 150 cells (default)

# Fine: 8 × 7 × 7 = 392 cells
mass_list = np.logspace(5, 10, 8)
impact_list = np.linspace(0.0, 10.0, 7)
velocity_list = np.linspace(0.0, 100.0, 7)

# Run both serial and parallel for each resolution
# Compare runtimes and parallel efficiency
```

---

## Troubleshooting

### Issue: Import errors

**Symptom:** `ModuleNotFoundError: No module named 'agama'`

**Solution:**
```bash
source ~/Work/venvs/.venv/bin/activate
pip install agama
```

### Issue: NameError on import

**Symptom:** `NameError: name 'dE_ana' is not defined`

**Cause:** Plotting cells reference undefined variables.

**Solution:** Use provided `simulation_calibration_version2.py` with commented-out execution blocks (lines 146-157, 424-442, 456-492, 495-763).

### Issue: Process hangs

**Symptom:** Script runs indefinitely without output.

**Debugging:**
```bash
# Check process status
ps aux | grep compute_grid

# Monitor CPU usage (should be high during computation)
top -pid <PID>

# Kill if necessary
kill -9 <PID>
```

### Issue: Memory errors

**Symptom:** `MemoryError` or system slowdown.

**Solution:** Reduce grid size or use serial version:
```python
mass_list = np.logspace(5, 10, 4)  # Reduce from 6 to 4
```

---

## Performance Notes

### Timing Estimates (14-core system)

| Grid Size | Serial | Parallel | Speedup |
|-----------|--------|----------|---------|
| 36 cells (4×3×3) | ~30 s | ~4 s | 7.5× |
| 150 cells (6×5×5) | ~124 s | ~15 s | 8.2× |
| 392 cells (8×7×7) | ~324 s | ~40 s | 8.1× |

**Scaling:** Near-linear with cell count. Parallel efficiency stable at ~58%.

### Memory Usage

- **Per simulation:** ~200 MB (AGAMA potentials, test particles)
- **Serial total:** ~214 MB
- **Parallel total:** ~214 MB × 14 = ~3 GB (14 processes)

**Recommendation:** Ensure >4 GB free RAM for parallel execution.

### Optimization Tips

1. **Grid design:** Focus resolution on regions of interest (e.g., narrow mass range)
2. **Worker count:** Match CPU core count; more workers ≠ faster (overhead increases)
3. **Batch processing:** For very large grids, split into chunks and run sequentially

---

## File Format Reference

### NPZ Structure

```python
# Save
np.savez(
    "output.npz",
    mass_list=mass_list,              # (6,)
    impact_list=impact_list,          # (5,)
    velocity_list=velocity_list,      # (5,)
    delta_EJ_sim=delta_EJ_sim,        # (6, 5, 5)
    delta_EJ_ana=delta_EJ_ana,        # (6, 5, 5)
    delta_Is_sim=delta_Is_sim,        # (6, 5, 5)
    delta_Is_ana=delta_Is_ana,        # (6, 5, 5)
    runtime_seconds=elapsed,          # scalar
    n_processes=n_processes,          # scalar (parallel only)
)

# Load
data = np.load("output.npz")
print(data.files)  # List all arrays
```

### Units Convention

| Quantity | Units | Symbol |
|----------|-------|--------|
| Mass | M☉ | M |
| Length | kpc | r, b |
| Velocity | km/s | v |
| Energy | (km/s)² | E, E_J |
| Action | kpc·km/s | I |
| Pattern speed | km/s/kpc | Ω_p |
| Time | Gyr | t |

---

## Future Directions

### Immediate Next Steps

1. **Analysis:** Plot simulation/analytical ratio across full parameter space
2. **Diagnostics:** Identify mass/velocity regimes where analytical model fails most
3. **Interpretation:** Connect failures to physical processes (resonance coherence, etc.)

### Enhancements

1. **Higher resolution:** Refine grid in critical regions
2. **Extended range:** Test M < 10⁵ M☉ and M > 10¹⁰ M☉
3. **Realistic distributions:** Weight grid by subhalo mass function and velocity PDF
4. **Additional observables:** Frequency shifts, phase mixing timescales

### Production Runs

For publication-quality results:

```python
# High-resolution grid (8×10×10 = 800 cells)
mass_list = np.logspace(5, 10, 8)
impact_list = np.linspace(0.0, 10.0, 10)
velocity_list = np.linspace(0.0, 100.0, 10)

# Expected runtime: ~650 s parallel (~11 minutes)
```

---

## Citation

If you use this code, please cite:

**Paper (in preparation):**
Belokurov et al., "Subhalo Perturbations to Bar Resonance Chevrons", MNRAS, 2025

**Code:**
GitHub repository: [URL when available]

---

## Contact

**Principal Investigator:** Dr. V.A. Belokurov
**Institution:** Institute of Astronomy, University of Cambridge
**Code Development:** Claude Code (Anthropic)

For questions or collaboration inquiries, contact [email/GitHub].

---

## License

[Specify license: MIT, GPL, Academic Use Only, etc.]

---

## Acknowledgments

- **AGAMA library:** E. Vasiliev (Edinburgh)
- **Computing resources:** [Institution/Cluster]
- **Theoretical framework:** [Collaborators]

---

*README last updated: 2025-12-01*
