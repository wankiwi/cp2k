#!/usr/bin/env python3
"""Exercise the radial moving-constraint bug and fail on the inertial offset."""

import argparse
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    args = parser.parse_args()
    scripts = Path(__file__).resolve().parent
    subprocess.run(
        [
            "python3",
            str(scripts / "reproduce_moving_target.py"),
            "--out",
            str(args.out),
            "--cp2k",
            str(args.cp2k),
        ],
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        check=True,
    )
    results = json.loads((args.out / "moving_target_results.json").read_text())
    failures = []
    for name, data in results.items():
        error = abs(data["excess_lambda_eV_per_angstrom"])
        velocity_error = abs(
            data["final_radial_velocity_angstrom_per_fs"]
            - data["growth_angstrom_per_fs"]
        )
        if error >= 5e-3 or velocity_error >= 1e-7:
            failures.append(
                f"{name}: force error={error:.8g} eV/A, "
                f"velocity error={velocity_error:.8g} A/fs"
            )
    assert not failures, "\n".join(failures)
    print("PASS: moving-constraint force and radial velocity")


if __name__ == "__main__":
    main()
