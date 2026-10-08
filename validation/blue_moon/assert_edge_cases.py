#!/usr/bin/env python3
"""Regression checks for clamping, RESPA, REFTRAJ and dimensionless CV units."""

import argparse
from pathlib import Path
import re
import shutil

from assert_restart_and_errors import expect_failure
from reproduce import generator, run_case


def nve_case(root, name, steps=20, growth=0.0005, limit=10.0):
    directory = root / name
    directory.mkdir(parents=True, exist_ok=False)
    generator.write_case(directory, name, steps, 1.0, 0.0, 2.5, growth, limit)
    inp = directory / f"{name}.inp"
    text = inp.read_text().replace("ENSEMBLE NVT", "ENSEMBLE NVE")
    text = re.sub(r"    &THERMOSTAT\n.*?    &END THERMOSTAT\n", "", text, flags=re.S)
    text = text.replace("  &CONSTRAINT\n", "  &CONSTRAINT\n    SHAKE_TOLERANCE 1e-10\n")
    text = text.replace(
        "&END MOTION",
        "  &PRINT\n    &VELOCITIES\n"
        "      UNIT angstrom*fs^-1\n      &EACH\n        MD 1\n"
        "      &END EACH\n    &END VELOCITIES\n  &END PRINT\n&END MOTION",
    )
    inp.write_text(text)
    return directory, inp


def radial_velocity(directory, name, frame=-1):
    lines = (directory / f"{name}-vel-1.xyz").read_text().splitlines()
    frames = [lines[i : i + 4] for i in range(0, len(lines), 4)]
    velocities = [[float(v) for v in line.split()[1:4]] for line in frames[frame][2:]]
    return velocities[1][0] - velocities[0][0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    args = parser.parse_args()

    name = "already_clamped_initial_radial_velocity"
    directory, inp = nve_case(args.out, name, limit=2.5)
    inp.write_text(
        inp.read_text().replace(
            "  &END SUBSYS",
            "    &VELOCITY\n"
            "      -0.0001 0 0\n      0.0001 0 0\n    &END VELOCITY\n  &END SUBSYS",
        )
    )
    rows = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"]
    assert (
        abs(radial_velocity(directory, name, 0)) < 1e-7
    ), "clamped initial projection skipped"
    assert abs(rows[0][5] - 14.3996 / 2.5**2) < 5e-3
    assert all(row[10] == 0 for row in rows)

    for sign in (1, -1):
        name = f"partial_limit_{'forward' if sign == 1 else 'reverse'}"
        directory, inp = nve_case(
            args.out, name, steps=4, growth=sign * 0.0005, limit=2.5 + sign * 0.00075
        )
        rows = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"]
        for frame, rate in enumerate(
            (sign * 0.0005, sign * 0.00025, 0.0, 0.0), start=1
        ):
            assert (
                abs(radial_velocity(directory, name, frame) - rate) < 1e-7
            ), "clamped actual rate mismatch"
        assert (
            abs(rows[-1][11] - rows[1][11]) < 1e-12
        ), "work changed after target stopped"

    name = "respa_five_inner_steps"
    directory, inp = nve_case(args.out, name)
    text = inp.read_text()
    start, end = text.index("&FORCE_EVAL"), text.index("&END FORCE_EVAL") + len(
        "&END FORCE_EVAL"
    )
    force = text[start:end]
    fast = force[: force.index("  &SUBSYS")] + "&END FORCE_EVAL\n"
    fast = fast.replace("CHARGE 1.0", "CHARGE 0.0").replace("CHARGE -1.0", "CHARGE 0.0")
    text = (
        text[:start]
        + "&MULTIPLE_FORCE_EVALS\n  FORCE_EVAL_ORDER 2 1\n&END MULTIPLE_FORCE_EVALS\n"
        + fast
        + text[start:]
    )
    text = text.replace(
        "  &END MD", "    &RESPA\n      FREQUENCY 5\n    &END RESPA\n  &END MD"
    )
    inp.write_text(text)
    rows = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"]
    assert abs(rows[-1][4] - 2.51) < 1e-10, "RESPA advanced target using outer dt"
    assert abs(radial_velocity(directory, name) - 0.0005) < 1e-7
    expected_work = 14.3996 * (1 / 2.5 - 1 / 2.51)
    assert abs(rows[-1][11] - expected_work) < 1e-4, "RESPA inner work not accumulated"

    reference_trajectory = directory / f"{name}-pos-1.xyz"
    name = "reftraj_snapshot_three_no_constraints"
    directory, inp = nve_case(args.out, name, steps=3)
    text = re.sub(
        r"  &CONSTRAINT\n.*?  &END CONSTRAINT\n", "", inp.read_text(), flags=re.S
    )
    text = text.replace(
        "ENSEMBLE NVE",
        "ENSEMBLE REFTRAJ\n    &REFTRAJ\n"
        "      EVAL ENERGY\n      FIRST_SNAPSHOT 3\n      LAST_SNAPSHOT 5\n"
        "      TRAJ_FILE_NAME reftraj.xyz\n    &END REFTRAJ",
    )
    inp.write_text(text)
    shutil.copy2(reference_trajectory, directory / "reftraj.xyz")
    run_case(args.cp2k, directory, inp)

    name = "dimensionless_coordination_incompatible_length"
    directory, inp = nve_case(args.out, name, steps=4, growth=None)
    text = inp.read_text().replace(
        "&DISTANCE\n        ATOMS 1 2\n      &END DISTANCE",
        "&COORDINATION\n        ATOMS_FROM 1\n        ATOMS_TO 2\n"
        "        R0 [angstrom] 2.0\n        NN 6\n        ND 12\n      &END COORDINATION",
    )
    text = text.replace("TARGET [angstrom] 2.50000000", "TARGET 0.2076973784290857")
    text = text.replace("    &BLUE_MOON\n", "    &BLUE_MOON\n      CV_UNIT angstrom\n")
    inp.write_text(text)
    expect_failure(args.cp2k, directory, inp, "CV_UNIT is incompatible")
    print(
        "PASS: already-clamped initialization, RESPA substeps/work, REFTRAJ and dimensionless units"
    )


if __name__ == "__main__":
    main()
