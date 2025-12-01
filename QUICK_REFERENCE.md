# Quick Reference: Chevron-Subhalo Grid Computation

**Last updated:** 2025-12-01

---

## Common Commands

### Run Serial Computation (150 cells, ~2 minutes)
```bash
source .venv/bin/activate
python compute_grid_serial.py
```

### Run Parallel Computation (150 cells, ~15 seconds)
```bash
source .venv/bin/activate
python compute_grid_parallel.py
```

### Validate Results
```bash
source .venv/bin/activate
python compare_grids.py
```

### Expected Output
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

## Grid Parameters

| Parameter | Range | N | Type |
|-----------|-------|---|------|
| Mass | 10⁵–10¹⁰ M☉ | 6 | log |
| Impact | 0–10 kpc | 5 | linear |
| Velocity | 0–100 km/s | 5 | linear |
| **Total** | - | **150** | - |

---

## Output Files

### calibration_grid_serial.npz
```python
mass_list        # (6,)     - Subhalo masses [M☉]
impact_list      # (5,)     - Impact parameters [kpc]
velocity_list    # (5,)     - Relative velocities [km/s]
delta_EJ_sim     # (6,5,5)  - ΔE_J from simulation
delta_EJ_ana     # (6,5,5)  - ΔE_J from analytical model
delta_Is_sim     # (6,5,5)  - ΔI from simulation
delta_Is_ana     # (6,5,5)  - ΔI from analytical model
runtime_seconds  # scalar   - Execution time [s]
```

### calibration_grid_parallel.npz
Same as serial + `n_processes` (number of workers)

---

## Load and Use Data

```python
import numpy as np

# Load
data = np.load("calibration_grid_parallel.npz")

# Extract
dEJ_sim = data["delta_EJ_sim"]  # Simulation results
dEJ_ana = data["delta_EJ_ana"]  # Analytical predictions
mass = data["mass_list"]
impact = data["impact_list"]
velocity = data["velocity_list"]

# Grid cell (i, j, k):
#   M_sub = mass[i]
#   b = impact[j]
#   v_rel = velocity[k]

# Example: Get result for M=10^8 M☉, b=5 kpc, v=50 km/s
i_m = 2  # mass[2] ≈ 1e8
i_b = 2  # impact[2] = 5.0
i_v = 2  # velocity[2] = 50.0
result = dEJ_sim[i_m, i_b, i_v]
```

---

## Modify Grid Size

Edit `compute_grid_*.py` (both files identically):

```python
# Current (6×5×5 = 150 cells, ~2 min serial)
mass_list = np.logspace(5, 10, 6)
impact_list = np.linspace(0.0, 10.0, 5)
velocity_list = np.linspace(0.0, 100.0, 5)

# Coarse (4×3×3 = 36 cells, ~30 sec serial)
mass_list = np.logspace(5, 10, 4)
impact_list = np.linspace(0.0, 10.0, 3)
velocity_list = np.linspace(0.0, 100.0, 3)

# Fine (8×7×7 = 392 cells, ~5 min serial)
mass_list = np.logspace(5, 10, 8)
impact_list = np.linspace(0.0, 10.0, 7)
velocity_list = np.linspace(0.0, 100.0, 7)
```

**Note:** Change in **both** scripts to ensure identical grids.

---

## Performance Estimates (14-core system)

| Grid | Cells | Serial | Parallel | Speedup |
|------|-------|--------|----------|---------|
| 4×3×3 | 36 | ~30 s | ~4 s | 7.5× |
| 6×5×5 | 150 | ~124 s | ~15 s | 8.2× |
| 8×7×7 | 392 | ~324 s | ~40 s | 8.1× |
| 10×10×10 | 1000 | ~825 s | ~101 s | 8.2× |

---

## Troubleshooting

### Problem: Import error
```
ModuleNotFoundError: No module named 'agama'
```
**Solution:**
```bash
source .venv/bin/activate
pip install agama numpy scipy matplotlib tqdm
```

### Problem: NameError
```
NameError: name 'dE_ana' is not defined
```
**Solution:** Use provided `simulation_calibration_version2.py` with commented execution blocks (already done in current version).

### Problem: Script hangs
**Check status:**
```bash
ps aux | grep compute_grid
```
**Kill if needed:**
```bash
killall python
```

### Problem: Memory error
**Solution:** Reduce grid size or use serial version.

---

## File Structure

```
chevron_subhalo_interaction/
├── README.md                    # Full documentation
├── PROJECT_JOURNAL.md           # Detailed experiment log
├── EXPERIMENT_SUMMARY.md        # Results summary
├── QUICK_REFERENCE.md           # This file
│
├── compute_grid_serial.py       # Serial computation
├── compute_grid_parallel.py     # Parallel computation
├── compare_grids.py             # Validation script
│
├── simulation_calibration_version2.py  # Core module
├── calibration_grid_serial.npz         # Serial output
└── calibration_grid_parallel.npz       # Parallel output
```

---

## Key Results from 2025-12-01 Experiment

✓ **Output validation:** PASS (exact match)
✓ **Runtime:** 123.7 s → 15.1 s
✓ **Speedup:** 8.17× with 14 workers
✓ **Efficiency:** 58.4%

**Conclusion:** Parallel implementation validated for production use.

---

## Physics Quick Facts

- **Bar pattern speed:** Ω_p = 35 km/s/kpc
- **Resonance:** Co-rotation (m=2, n=0)
- **Preferred metric:** ΔE_J (Jacobi integral change)
- **Key finding:** Analytical model underpredicts for M = 10⁵–10⁹ M☉

---

## Contact

**PI:** Dr. V.A. Belokurov
**Code:** Claude Code (Anthropic)
**Date:** 2025-12-01

For issues, see full documentation in `README.md` and `PROJECT_JOURNAL.md`.
