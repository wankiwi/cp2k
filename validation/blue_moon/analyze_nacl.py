#!/usr/bin/env python3
"""Analyze the full-range PDF benchmark without equating mechanical work to TI."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COULOMB = 14.3996  # eV angstrom, as used in the supplied PDF benchmark.
KB_EV_PER_K = 8.617333262145e-5


def load_log(path):
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != 12 or not np.all(np.isfinite(data)):
        raise ValueError(f"Invalid Blue Moon table: {path}")
    if not np.all(data[:, 2] == 1):
        raise ValueError(f"Expected one distance constraint: {path}")
    return data


def energy(r):
    return COULOMB * (1 / 2.5 - 1 / np.asarray(r))


def radial_free_energy_1k(r):
    # The radial configuration measure is r^2 dr, giving -2 kBT log(r/r0).
    return energy(r) - 2 * KB_EV_PER_K * np.log(np.asarray(r) / 2.5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    ti = []
    for path in args.root.glob("ti_windows/r_*/ti_*-1.blueMoonLog"):
        data = load_log(path)
        tail = data[len(data) // 2 :]
        # Blue Moon estimator: ratio of means, not mean of framewise ratios.
        force = np.mean(tail[:, 8]) / np.mean(tail[:, 6])
        ti.append((np.mean(tail[:, 4]), force))
    ti = np.array(sorted(ti))
    if len(ti) != 31 or not np.allclose(ti[[0, -1], 0], [2.5, 10]):
        raise ValueError("Incomplete 2.5--10 angstrom, 31-window TI benchmark")
    ti_work = np.r_[0, np.cumsum(np.diff(ti[:, 0]) * (ti[:-1, 1] + ti[1:, 1]) / 2)]
    summary = {
        "ti": {
            "windows": len(ti),
            "work_eV": float(ti_work[-1]),
            "analytic_eV": float(energy(10)),
            "error_eV": float(ti_work[-1] - energy(10)),
            "radial_fes_1k_eV": float(radial_free_energy_1k(10)),
            "estimator": "mean(weight*(lambda+GkT))/mean(weight), last 25 of 50 frames",
        }
    }
    style = {
        "font.size": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.dpi": 220,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(2, 2, figsize=(8, 6), layout="constrained")
        force_ax, work_ax, error_ax, temp_ax = axes.flat
        r = np.linspace(2.5, 10, 400)
        force_ax.plot(r, COULOMB / r**2, color="#333333", label="Coulomb force")
        work_ax.plot(r, energy(r), color="#333333", label="Coulomb energy")
        work_ax.plot(
            r, radial_free_energy_1k(r), "--", color="#666666", label="Radial FES (1 K)"
        )
        force_ax.plot(
            ti[:, 0], ti[:, 1], "o", ms=3, color="#009E73", label="TI (weighted mean)"
        )
        work_ax.plot(
            ti[:, 0], ti_work, "o-", ms=3, color="#009E73", label="TI (trapezoid)"
        )
        error_ax.plot(
            ti[:, 0],
            ti_work - energy(ti[:, 0]),
            "o-",
            ms=3,
            color="#009E73",
            label="TI",
        )
        colors = ["#0072B2", "#CC79A7", "#D55E00"]
        styles = ["-", "--", ":"]
        for path, color, linestyle in zip(
            sorted(args.root.glob("sg_*/sg_*-1.blueMoonLog")), colors, styles
        ):
            data = load_log(path)
            name = path.parent.name
            rate = float(name.removeprefix("sg_").replace("m", "-"))
            if abs(data[-1, 4] - 10) > 1e-8:
                raise ValueError(f"Incomplete slow growth range: {path}")
            label = f"SG {rate:g} A/fs"
            force_ax.plot(
                data[:, 4],
                data[:, 5],
                color=color,
                ls=linestyle,
                lw=0.7,
                alpha=0.8,
                label=label,
            )
            work_ax.plot(
                np.r_[2.5, data[:, 4]],
                np.r_[0, data[:, 11]],
                color=color,
                ls=linestyle,
                label=label,
            )
            error_ax.plot(
                np.r_[2.5, data[:, 4]],
                np.r_[0, data[:, 11] - energy(data[:, 4])],
                color=color,
                ls=linestyle,
                label=label,
            )
            ener = np.loadtxt(path.parent / f"{name}-1.ener", comments="#", ndmin=2)
            # Energies are sampled more often than the diagnostic table; map step to scheduled CV.
            temp_ax.plot(
                np.minimum(2.5 + ener[:, 0] * rate, 10),
                ener[:, 3],
                color=color,
                ls=linestyle,
                lw=0.6,
                label=label,
            )
            summary[name] = {
                "steps": int(data[-1, 0]),
                "rate_angstrom_per_fs": rate,
                "work_eV": float(data[-1, 11]),
                "analytic_eV": float(energy(10)),
                "error_eV": float(data[-1, 11] - energy(10)),
                "difference_from_radial_fes_1k_eV": float(
                    data[-1, 11] - radial_free_energy_1k(10)
                ),
                "diagnostic_stride": int(data[1, 0] - data[0, 0]),
                "difference_from_ti_eV": float(data[-1, 11] - ti_work[-1]),
                "max_abs_cv_error_angstrom": float(
                    np.max(np.abs(data[:, 3] - data[:, 4]))
                ),
                "within_0p02_eV_work_tolerance": bool(
                    abs(data[-1, 11] - energy(10)) <= 0.02
                ),
                "mean_temperature_K": float(np.mean(ener[:, 3])),
            }
        for ax in axes.flat:
            ax.set(xlabel="Na--Cl distance (angstrom)", xlim=(2.5, 10))
            ax.grid(alpha=0.15)
        force_ax.set(
            ylabel="Constraint force (eV/angstrom)", title="Instantaneous force and TI"
        )
        work_ax.set(
            ylabel="Energy / mechanical work (eV)",
            title="Common CV range; W(2.5 A) = 0",
        )
        error_ax.axhline(0, color="#777777", lw=0.6)
        error_ax.axhline(0.02, color="#777777", lw=0.6, ls="--")
        error_ax.axhline(-0.02, color="#777777", lw=0.6, ls="--")
        error_ax.set(
            ylabel="Difference from Coulomb energy (eV)",
            title="Integration and finite-rate error",
        )
        temp_ax.axhline(1, color="#777777", lw=0.6)
        temp_ax.set(ylabel="Temperature (K)", title="Nose--Hoover, nominal 1 K")
        force_ax.legend(frameon=False)
        work_ax.legend(frameon=False)
        fig.savefig(args.out / "nacl_full_range.png")
        fig.savefig(args.out / "nacl_full_range.pdf")
        plt.close(fig)
    with (args.out / "ti.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "r_angstrom",
                "weighted_force_eV_per_angstrom",
                "ti_integral_eV",
                "coulomb_eV",
            ]
        )
        writer.writerows(zip(ti[:, 0], ti[:, 1], ti_work, energy(ti[:, 0])))
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    assert len(summary) == 4
    assert all(
        abs(v["error_eV"]) <= 0.02 for v in summary.values()
    ), "Benchmark tolerance exceeded"


if __name__ == "__main__":
    main()
