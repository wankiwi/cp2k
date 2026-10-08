#!/usr/bin/env python3
"""Plot all native samples from three matched-range NaCl growth-rate controls."""

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.axisbelow": True,
    "legend.fontsize": 8.5,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "path.simplify": False,
    "savefig.facecolor": "white",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(fig, name, output):
    files = {}
    for extension in ("png", "pdf", "svg"):
        path = output / f"{name}.{extension}"
        require(not path.exists(), f"Refusing to overwrite {path}")
        fig.savefig(path, dpi=300, facecolor="white", transparent=False)
        if extension == "svg":
            lines = path.read_text().splitlines()
            path.write_text("\n".join(line.rstrip() for line in lines) + "\n")
        files[path.name] = {"sha256": digest(path), "bytes": path.stat().st_size}
    with Image.open(output / f"{name}.png") as raster:
        files[f"{name}.png"].update(
            pixels=list(raster.size), mode=raster.mode, dpi=list(raster.info["dpi"])
        )
    plt.close(fig)
    return files


CASES = (
    ("nacl_growth_0p00025_steps.csv", 0.00025, 10000),
    ("nacl_radial_steps.csv", 0.0005, 5000),
    ("nacl_growth_0p001_steps.csv", 0.001, 2500),
)
RATE_COLORS = ("#0072B2", "#D55E00", "#009E73")
FIXED_DASHES = ((0, (2, 2)), (0, (6, 2)), (0, (6, 2, 1, 2)))
ANALYTIC = "#303030"


def load_data(directory):
    data = []
    for name, growth, steps in CASES:
        raw = np.loadtxt(directory / name, delimiter=",", skiprows=1)
        require(raw.shape == (steps, 8) and np.all(np.isfinite(raw)), f"Invalid {name}")
        require(
            np.array_equal(raw[:, 0], np.arange(1, steps + 1)), f"Frame gaps: {name}"
        )
        require(
            np.allclose(raw[:, 1], 2.5 + growth * raw[:, 0], rtol=0, atol=1e-12),
            f"Wrong target schedule: {name}",
        )
        require(
            np.allclose(raw[:, 4], 14.3996 / raw[:, 1] ** 2, rtol=0, atol=1e-12)
            and np.allclose(
                raw[:, 7],
                np.cumsum(raw[:, 4] * np.diff(np.r_[2.5, raw[:, 1]])),
                rtol=0,
                atol=1e-12,
            ),
            f"Wrong analytic force/work: {name}",
        )
        data.append((name, growth, raw))
    return data


def figure(data):
    fig = plt.figure(figsize=(10.5, 4.8), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1, 0.23))
    force = fig.add_subplot(grid[0, 0])
    work = fig.add_subplot(grid[0, 1])
    legend = fig.add_subplot(grid[1, :])
    legend.set_axis_off()
    masters, fixes = [], []
    reference_handle = None
    for index, (_, growth, raw) in enumerate(data):
        color = RATE_COLORS[index]
        for label, line, force_column, work_column in (
            (f"Master: {growth:g} Å/fs", "-", 2, 5),
            (f"Fixed: {growth:g} Å/fs", FIXED_DASHES[index], 3, 6),
        ):
            for axis, column in ((force, force_column), (work, work_column)):
                handle = axis.plot(
                    raw[:, 1],
                    raw[:, column],
                    label=label,
                    color=color,
                    linestyle=line,
                    linewidth=1.5 if label.startswith("Master") else 1.8,
                )[0]
                if axis is force:
                    (masters if label.startswith("Master") else fixes).append(handle)
        stride = max(1, (len(raw) - 1) // 16)
        # Force is one shared analytic function; discrete work remains rate-specific.
        references = ((force, 4), (work, 7)) if index == 0 else ((work, 7),)
        for axis, column in references:
            handle = axis.plot(
                raw[:, 1],
                raw[:, column],
                color=ANALYTIC,
                linestyle="none",
                marker="o",
                markevery=(0 if axis is force else index * stride // 3, stride),
                markersize=4.5,
                markerfacecolor="white",
                markeredgewidth=1.0,
                label="Analytic Coulomb",
                zorder=4,
            )[0]
            if axis is force:
                reference_handle = handle
    force.set(title="A  NaCl: raw constraint force", ylabel="Constraint force (eV/Å)")
    work.set(title="B  NaCl: accumulated mechanical work", ylabel="Work (eV)")
    for axis in (force, work):
        axis.grid(alpha=0.16)
        axis.set(xlim=(2.5, 5.0), xlabel="Na–Cl target (Å)")
        axis.ticklabel_format(axis="x", style="plain", useOffset=False)
    handles = [
        masters[0],
        fixes[0],
        reference_handle,
        masters[1],
        fixes[1],
        masters[2],
        fixes[2],
    ]
    legend.legend(
        handles=handles,
        loc="center",
        ncol=3,
        handlelength=3.8,
        columnspacing=2.3,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).parent / "data")
    parser.add_argument("--out", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    receipt_path = args.out / "nacl_growth_rate_export_receipt.json"
    require(not receipt_path.exists(), "Refusing to overwrite export receipt")
    data = load_data(args.data_dir)
    with plt.rc_context(STYLE):
        files = export(figure(data), "nacl_growth_rate_comparison", args.out)
    metrics = {}
    for name, growth, raw in data:
        metrics[name] = {
            "growth_A_per_fs": growth,
            "rendered_rows": len(raw),
            "master_force_MAE_eV_per_A": float(np.mean(np.abs(raw[:, 2] - raw[:, 4]))),
            "fixed_force_MAE_eV_per_A": float(np.mean(np.abs(raw[:, 3] - raw[:, 4]))),
            "master_final_work_error_eV": float(raw[-1, 5] - raw[-1, 7]),
            "fixed_final_work_error_eV": float(raw[-1, 6] - raw[-1, 7]),
        }
    receipt = {
        "generated_SGT": datetime.now(ZoneInfo("Asia/Singapore")).isoformat(),
        "script_sha256": digest(Path(__file__)),
        "inputs_sha256": {name: digest(args.data_dir / name) for name, _, _ in data},
        "metrics": metrics,
        "smoothing": False,
        "offset_adjustment": False,
        "excluded_frames": 0,
        "same_axis_scales_per_column": True,
        "presentation": {
            "layout": "Two shared panels with all three growth rates overlaid",
            "master": "Solid line, one color per growth rate",
            "fixed": "Rate-colored short-dashed, long-dashed and dash-dot lines; no markers",
            "analytic": "Dark open-circle markers, no connecting line",
            "analytic_force": "One common Coulomb function sampled on the slow-rate native grid",
            "analytic_work": "Each rate's discrete reference retained; sparse marker positions staggered",
            "analytic_marker_stride": {
                name: max(1, (len(raw) - 1) // 16) for name, _, raw in data
            },
            "all_master_and_fixed_samples_drawn": True,
        },
        "reference_work": "Sum of analytic force_n times target increment_n, including first interval",
        "files": files,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
