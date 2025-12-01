#!/usr/bin/env python3
"""
Compare serial and parallel grid outputs.
"""
import numpy as np

# Load both files
serial = np.load("calibration_grid_serial.npz")
parallel = np.load("calibration_grid_parallel.npz")

# Extract arrays
arrays_to_compare = [
    ("delta_EJ_sim", serial["delta_EJ_sim"], parallel["delta_EJ_sim"]),
    ("delta_EJ_ana", serial["delta_EJ_ana"], parallel["delta_EJ_ana"]),
    ("delta_Is_sim", serial["delta_Is_sim"], parallel["delta_Is_sim"]),
    ("delta_Is_ana", serial["delta_Is_ana"], parallel["delta_Is_ana"]),
]

print("Array Comparison (serial vs parallel)")
print("=" * 60)

all_identical = True
for name, arr_s, arr_p in arrays_to_compare:
    identical = np.allclose(arr_s, arr_p, rtol=1e-15, atol=0)
    max_diff = np.max(np.abs(arr_s - arr_p))

    print(f"\n{name}:")
    print(f"  Shape: {arr_s.shape}")
    print(f"  Identical: {identical}")
    print(f"  Max absolute difference: {max_diff:.2e}")

    if not identical:
        all_identical = False

print("\n" + "=" * 60)
print(f"Overall: {'PASS - All arrays identical' if all_identical else 'FAIL - Arrays differ'}")

# Runtime comparison
t_serial = float(serial["runtime_seconds"])
t_parallel = float(parallel["runtime_seconds"])
speedup = t_serial / t_parallel
n_proc = int(parallel["n_processes"])

print("\nRuntime Comparison")
print("=" * 60)
print(f"Serial runtime:     {t_serial:.3f} s")
print(f"Parallel runtime:   {t_parallel:.3f} s ({n_proc} workers)")
print(f"Speedup:            {speedup:.2f}x")
print(f"Parallel efficiency: {100 * speedup / n_proc:.1f}%")
