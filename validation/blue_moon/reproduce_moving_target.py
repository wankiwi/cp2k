#!/usr/bin/env python3
"""Isolate the moving-constraint velocity condition in radial NaCl NVE."""

import argparse
import json
from pathlib import Path
import re

from reproduce import ROOT, generator, run_case


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, default=ROOT / "exe/minimal/cp2k.ssmp")
    args = parser.parse_args()
    summary = {}
    for name, growth, dt, steps in (
        ("fixed_target", None, 1.0, 20),
        ("moving_target_0p0005_angstrom_per_fs_dt_1fs", 0.0005, 1.0, 20),
        ("moving_target_0p001_angstrom_per_fs_dt_1fs", 0.001, 1.0, 10),
        ("moving_target_0p0005_angstrom_per_fs_dt_0p5fs", 0.0005, 0.5, 40),
        ("moving_target_slow_rate_default_tolerance", 0.00001, 1.0, 20),
    ):
        directory = args.out / "nacl_radial_nve_zero_initial_velocity" / name
        directory.mkdir(parents=True, exist_ok=False)
        generator.write_case(
            directory,
            name,
            steps,
            dt,
            0.0,
            2.5,
            growth=growth,
            limit=10.0,
            print_stride=1,
        )
        input_file = directory / f"{name}.inp"
        text = input_file.read_text().replace("ENSEMBLE NVT", "ENSEMBLE NVE")
        if name != "moving_target_slow_rate_default_tolerance":
            text = text.replace(
                "  &CONSTRAINT\n", "  &CONSTRAINT\n    SHAKE_TOLERANCE 1.0E-10\n"
            )
        text = re.sub(
            r"    &THERMOSTAT\n.*?    &END THERMOSTAT\n", "", text, count=1, flags=re.S
        )
        text = text.replace(
            "&END MOTION",
            "  &PRINT\n    &VELOCITIES\n"
            "      UNIT angstrom*fs^-1\n      &EACH\n"
            "        MD 1\n      &END EACH\n"
            "    &END VELOCITIES\n  &END PRINT\n&END MOTION",
        )
        input_file.write_text(text)
        logs = run_case(args.cp2k, directory, input_file)
        last = logs["blueMoonLog"]["rows"][-1]
        velocity_file = directory / f"{name}-vel-1.xyz"
        lines = velocity_file.read_text().splitlines()
        first_velocity = [float(v) for v in lines[-2].split()[1:4]]
        second_velocity = [float(v) for v in lines[-1].split()[1:4]]
        analytic_force = 14.3996 / last[4] ** 2
        summary[name] = {
            "growth_angstrom_per_fs": growth or 0.0,
            "dt_fs": dt,
            "final_target_angstrom": last[4],
            "lambda_eV_per_angstrom": last[5],
            "coulomb_force_eV_per_angstrom": analytic_force,
            "excess_lambda_eV_per_angstrom": last[5] - analytic_force,
            "final_radial_velocity_angstrom_per_fs": second_velocity[0]
            - first_velocity[0],
            "last_blue_moon_row": last,
        }
        print(name, json.dumps(summary[name]))
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "moving_target_results.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
