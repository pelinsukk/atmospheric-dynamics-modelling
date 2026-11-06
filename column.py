"""Idealized dry-column experiments; all inputs are synthetic."""
import argparse
import json
from pathlib import Path

import numpy as np


HEIGHT = 1000.0
SURFACE = 290.0
AMPLITUDE = 8.0
DURATION = 6 * 3600.0


def pressure_change(net_inflow, area, seconds, gravity=9.81):
    """Hydrostatic pressure change (Pa) for constant net mass inflow (kg/s)."""
    return gravity * net_inflow * seconds / area


def potential_temperature(temperature, pressure, reference=100000.0):
    """Dry potential temperature in kelvin; pressure inputs are in pascals."""
    return temperature * (reference / pressure) ** (287.0 / 1004.0)


def simulate(cells, diffusivity, duration=DURATION):
    """Return cell-centred numerical and exact profiles and the flux budget."""
    if cells < 2 or not np.isfinite(diffusivity) or diffusivity <= 0:
        raise ValueError("Use at least two cells and positive finite diffusivity.")
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError("Duration must be positive and finite.")
    dz = HEIGHT / cells
    z = (np.arange(cells) + 0.5) * dz
    wave = np.pi / (2 * HEIGHT)
    # Exact initial cell averages of the sine mode.
    initial = AMPLITUDE * np.sinc(wave * dz / (2 * np.pi)) * np.sin(wave * z)
    anomaly = initial.copy()
    steps = max(1, int(np.ceil(duration * diffusivity / (0.25 * dz**2))))
    dt = duration / steps
    flux = np.zeros(cells + 1)
    boundary_integral = 0.0
    for _ in range(steps):
        flux[0] = -2 * diffusivity * anomaly[0] / dz
        flux[1:-1] = -diffusivity * np.diff(anomaly) / dz
        # Upper boundary is insulating: flux[-1] remains zero.
        anomaly -= dt * np.diff(flux) / dz
        boundary_integral += dt * (flux[0] - flux[-1])
    exact = initial * np.exp(-diffusivity * wave**2 * duration)
    residual = dz * np.sum(anomaly - initial) - boundary_integral
    return z, initial, anomaly, exact, float(residual)


def experiment(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out.mkdir(parents=True, exist_ok=True)
    fig, (profiles, convergence) = plt.subplots(1, 2, figsize=(10, 4.3))
    colors = ["#0072B2", "#D55E00"]
    report = {}
    grids = np.array([20, 40, 80, 160])
    for k, color in zip((2.0, 20.0), colors):
        z, initial, numerical, exact, residual = simulate(80, k)
        if k == 2.0:
            profiles.plot(SURFACE + initial, z, "--", color="0.4", label="Initial")
        profiles.plot(SURFACE + numerical, z, color=color, label=f"K = {k:g} m²/s")
        profiles.plot(SURFACE + exact[::6], z[::6], "o", color=color, ms=3, mfc="none")
        errors = []
        for n in grids:
            _, _, computed, analytical, _ = simulate(n, k)
            errors.append(float(np.sqrt(np.mean((computed - analytical)**2))))
        convergence.loglog(HEIGHT / grids, errors, "o-", color=color, label=f"K = {k:g}")
        report[str(k)] = {
            "grid_cells": grids.tolist(),
            "rms_error_K": errors,
            "last_refinement_order": float(np.log2(errors[-2] / errors[-1])),
            "column_budget_residual_K_m_at_80_cells": residual,
            "fraction_of_initial_column_anomaly_at_6h": float(numerical.sum() / initial.sum()),
        }
    profiles.set(xlabel="Potential temperature [K]", ylabel="Height [m]",
                 title="After 6 h • circles: exact solution")
    convergence.set(xlabel="Vertical grid spacing [m]", ylabel="RMS error [K]",
                    title="Grid refinement against the exact mode")
    for ax in (profiles, convergence):
        ax.grid(alpha=.24)
        ax.legend()
    fig.tight_layout()
    fig.savefig(out / "column.png", dpi=165)
    plt.close(fig)
    report["synthetic_pressure_example_hPa"] = pressure_change(5e6, 1e10, 21600) / 100
    (out / "metrics.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path("outputs"))
    experiment(parser.parse_args().out_dir)
