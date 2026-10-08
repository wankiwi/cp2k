#!/usr/bin/env python3
"""Integration assertions for work, restart, per-CV units and output aliases."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    args = parser.parse_args()
    subprocess.run(
        [
            "python3",
            str(Path(__file__).with_name("reproduce.py")),
            "--out",
            str(args.out),
            "--cp2k",
            str(args.cp2k),
        ],
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
        check=True,
    )
    data = json.loads((args.out / "reproduction_results.json").read_text())
    rows = lambda name, key="blueMoonLog": data[name][key]["rows"]
    every = rows("slow_growth_print_every_step")
    stride = rows("slow_growth_print_every_five_steps")
    alias = rows("slow_growth_print_every_step", "colvarSlowGrowthLog")
    alone = rows("slow_growth_alias_only", "colvarSlowGrowthLog")
    assert every == alias, "print keys must be exact readers of the same cache"
    assert abs(every[-1][11] - stride[-1][11]) < 1e-10, "print stride changes work"
    assert abs(every[-1][11] - alone[-1][11]) < 1e-10, "alias-only work differs"
    assert every[0][10] != 0, "first pulling interval was lost"
    assert abs(sum(row[10] for row in every) - every[-1][11]) < 1e-10
    resumed = rows("slow_growth_restart_next_twenty_steps")
    assert (
        resumed[0][0] == 21 and resumed[0][10] != 0
    ), "first resumed interval was lost"
    assert abs(resumed[-1][11] - every[-1][11]) < 1e-8, "restart resets work history"
    assert (
        data["slow_growth_restart_next_twenty_steps"]["blueMoonLog"]["metadata"][
            "cv_unit"
        ]
        == "angstrom"
    )
    ang = rows("fixed_distance_target_angstrom")[-1]
    meters = rows("fixed_distance_target_meters")[-1]
    assert math.isclose(meters[6], ang[6] * 1e10, rel_tol=1e-10), "meter metric scale"
    assert math.isclose(meters[5], ang[5] * 1e10, rel_tol=1e-6), "meter force scale"
    mixed = rows("fixed_distance_and_angle_input_units")
    assert abs(mixed[-2][4] - 2.5) < 1e-10 and abs(mixed[-1][4] - 90) < 1e-10
    ev = rows("fixed_angle_internal_cv_energy_eV")[-1]
    ha = rows("fixed_angle_internal_cv_energy_hartree")[-1]
    assert math.isclose(
        ev[5] / ha[5], 27.211386245988, rel_tol=1e-6
    ), "internal CV energy conversion"
    assert math.isclose(
        ev[7] / ha[7], 27.211386245988, rel_tol=1e-6
    ), "GkT energy conversion"
    print("PASS: aliases, stride, first interval, restart, mixed CV and SI units")


if __name__ == "__main__":
    main()
