#!/usr/bin/env python3
"""
Compute the 3D (mass, impact, velocity) calibration grid in serial.

This script expects that `simulation_calibration_version2.py` is importable
and that it defines:

    def one_simulation(subhalo_mass, subhalo_impact, subhalo_velocity) -> list[float]

which returns a list:
    [delta_EJ_sim, delta_EJ_ana, delta_Is_sim, delta_Is_ana]

The script:
  - defines the active parameter grid
  - loops over all 3D grid points
  - calls `one_simulation` at each point
  - saves the resulting 3D arrays and run time to a .npz file
"""

import time
import numpy as np

from simulation_calibration_version2 import one_simulation


def main():
    # Active parameter grid (as in your original notebook/script)
    mass_list = np.logspace(5, 10, 6)          # 6 values
    impact_list = np.linspace(0.0, 10.0, 5)    # 5 values
    velocity_list = np.linspace(0.0, 100.0, 5) # 5 values

    n_m = len(mass_list)
    n_b = len(impact_list)
    n_v = len(velocity_list)

    # Allocate result arrays: shape (n_m, n_b, n_v)
    delta_EJ_sim = np.zeros((n_m, n_b, n_v), dtype=float)
    delta_EJ_ana = np.zeros((n_m, n_b, n_v), dtype=float)
    delta_Is_sim = np.zeros((n_m, n_b, n_v), dtype=float)
    delta_Is_ana = np.zeros((n_m, n_b, n_v), dtype=float)

    t0 = time.time()

    # Triple nested loop over the 3D grid
    for i_m, m in enumerate(mass_list):
        for i_b, b in enumerate(impact_list):
            for i_v, v in enumerate(velocity_list):
                result = one_simulation(m, b, v)

                if result is None or len(result) != 4:
                    raise RuntimeError(
                        "one_simulation(m, b, v) must return a list of 4 floats: "
                        "[delta_EJ_sim, delta_EJ_ana, delta_Is_sim, delta_Is_ana]. "
                        "Got: {!r}".format(result)
                    )

                dEJ_sim, dEJ_ana, dIs_sim, dIs_ana = result

                delta_EJ_sim[i_m, i_b, i_v] = dEJ_sim
                delta_EJ_ana[i_m, i_b, i_v] = dEJ_ana
                delta_Is_sim[i_m, i_b, i_v] = dIs_sim
                delta_Is_ana[i_m, i_b, i_v] = dIs_ana

    elapsed = time.time() - t0

    # Save everything in a single file
    output_file = "calibration_grid_serial.npz"
    np.savez(
        output_file,
        mass_list=mass_list,
        impact_list=impact_list,
        velocity_list=velocity_list,
        delta_EJ_sim=delta_EJ_sim,
        delta_EJ_ana=delta_EJ_ana,
        delta_Is_sim=delta_Is_sim,
        delta_Is_ana=delta_Is_ana,
        runtime_seconds=elapsed,
    )

    print("Serial grid computation finished.")
    print("  Grid shape: {} x {} x {}".format(n_m, n_b, n_v))
    print("  Total cells: {}".format(n_m * n_b * n_v))
    print("  Runtime: {:.3f} s".format(elapsed))
    print("  Saved to: {}".format(output_file))


if __name__ == "__main__":
    main()
