#!/usr/bin/env python3
"""Run and analyze Coulomb-only NaCl controls with native CP2K output."""

import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from zoneinfo import ZoneInfo
import numpy as np

CASES = (
    ("growth_0p00025_A_per_fs", 0.00025, 10000, "nacl_growth_0p00025_steps.csv"),
    ("growth_0p0005_A_per_fs", 0.0005, 5000, "nacl_radial_steps.csv"),
    ("growth_0p001_A_per_fs", 0.001, 2500, "nacl_growth_0p001_steps.csv"),
)
# Frozen from src/common/physcon.F at master d59727108bfd58f66da98f7ff772c2be91f08b08.
HARTREE_EV = 27.211383856556296
BOHR_A = 0.5291772085900001


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frames(path):
    lines = path.read_text().splitlines()
    data = {}
    i = 0
    while i < len(lines):
        n = int(lines[i])
        assert n == 2, path
        step = int(re.search(r"\bi\s*=\s*(\d+)", lines[i + 1])[1])
        assert step not in data
        data[step] = np.array(
            [[float(x) for x in line.split()[1:4]] for line in lines[i + 2 : i + 4]]
        )
        i += 4
    return data


def native_force(directory, steps):
    values = []
    for step in range(1, steps + 1):
        file = directory / f"nacl-1_{step}.LagrangeMultLog"
        rows = re.findall(
            r"(?m)^\s*Shake\s+Lagrangian Multipliers:\s+([-+\d.Ee]+)", file.read_text()
        )
        assert len(rows) == 1, file
        values.append(-float(rows[0]) * HARTREE_EV / BOHR_A)
    force = np.array(values)
    assert np.all(np.isfinite(force))
    return force


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--master", type=Path, required=True)
    parser.add_argument("--fixed", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k-data-dir", type=Path)
    parser.add_argument(
        "--disable-aslr",
        action="store_true",
        help="Optional Linux Debug/LSAN workaround.",
    )
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / "data").mkdir()
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    if args.cp2k_data_dir:
        env["CP2K_DATA_DIR"] = str(args.cp2k_data_dir.resolve(strict=True))
    binaries = {
        "master": args.master.resolve(strict=True),
        "fixed": args.fixed.resolve(strict=True),
    }
    metrics = []
    receipts = []
    for name, growth, steps, csv in CASES:
        target = 2.5 + growth * np.arange(1, steps + 1)
        increments = np.diff(np.r_[2.5, target])
        records = {}
        analytic = 14.3996 / target**2
        analytic_work = np.cumsum(analytic * increments)
        for variant, binary in binaries.items():
            source = here / "inputs" / name / variant
            directory = args.out / name / variant
            directory.mkdir(parents=True)
            for filename in ("nacl.inp", "nacl.xyz"):
                shutil.copy2(source / filename, directory / filename)
            command = [str(binary), "-i", "nacl.inp", "-o", "cp2k.out"]
            if args.disable_aslr:
                command = ["setarch", "x86_64", "-R"] + command
            with (directory / "launch.log").open("w") as log:
                subprocess.run(
                    command,
                    cwd=directory,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
            output = (directory / "cp2k.out").read_text()
            assert "PROGRAM ENDED AT" in output and "ABORT" not in output, directory
            position = frames(directory / "nacl-pos-1.xyz")
            velocity = frames(directory / "nacl-vel-1.xyz")
            assert set(position) == set(velocity) == set(range(steps + 1))
            radius = []
            rate = []
            for step in range(1, steps + 1):
                r = position[step][1] - position[step][0]
                distance = np.linalg.norm(r)
                radius.append(distance)
                rate.append(np.dot(velocity[step][1] - velocity[step][0], r / distance))
            position_error = float(np.max(np.abs(np.array(radius) - target)))
            velocity_error = float(np.max(np.abs(np.array(rate) - growth)))
            assert position_error < 1e-7
            if variant == "fixed":
                assert velocity_error < 1e-7
            force = native_force(directory, steps)
            work = np.cumsum(force * increments)
            records[variant] = (force, work)
            receipts.append(
                dict(
                    case=name,
                    variant=variant,
                    steps=steps,
                    max_position_error_A=position_error,
                    max_velocity_error_A_per_fs=velocity_error,
                    input_sha256=digest(directory / "nacl.inp"),
                )
            )
        table = np.column_stack(
            [
                np.arange(1, steps + 1),
                target,
                records["master"][0],
                records["fixed"][0],
                analytic,
                records["master"][1],
                records["fixed"][1],
                analytic_work,
            ]
        )
        np.savetxt(
            args.out / "data" / csv,
            table,
            delimiter=",",
            comments="",
            header="step,target_A,master_lambda_eV_per_A,fixed_lambda_eV_per_A,coulomb_lambda_eV_per_A,master_work_eV,fixed_work_eV,coulomb_discrete_work_eV",
        )
        metrics.append(
            dict(
                growth_A_per_fs=growth,
                steps=steps,
                master_force_MAE_eV_per_A=float(
                    np.mean(np.abs(records["master"][0] - analytic))
                ),
                fixed_force_MAE_eV_per_A=float(
                    np.mean(np.abs(records["fixed"][0] - analytic))
                ),
                master_final_work_error_eV=float(
                    records["master"][1][-1] - analytic_work[-1]
                ),
                fixed_final_work_error_eV=float(
                    records["fixed"][1][-1] - analytic_work[-1]
                ),
            )
        )
        print(name, metrics[-1], flush=True)
    result = dict(
        generated_SGT=datetime.now(ZoneInfo("Asia/Singapore")).isoformat(),
        metrics=metrics,
        cases=receipts,
        executable_sha256={k: digest(v) for k, v in binaries.items()},
        constants=dict(hartree_eV=HARTREE_EV, bohr_A=BOHR_A),
        excluded_frames=0,
        smoothing=False,
        offset_adjustment=False,
    )
    (args.out / "results.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
