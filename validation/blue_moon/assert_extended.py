#!/usr/bin/env python3
"""Limit, reverse, print-off and joint-metric integration checks."""

import argparse
import math
from pathlib import Path
import re

from reproduce import generator, run_case


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    args = parser.parse_args()
    for name, growth, limit in (
        ("forward_limit_partial_step", 0.0005, 2.50525),
        ("reverse_limit_partial_step", -0.0005, 2.49475),
    ):
        directory = args.out / "nacl_radial_nve" / name
        generator.write_case(directory, name, 16, 1.0, 0.0, 2.5, growth, limit)
        inp = directory / f"{name}.inp"
        text = inp.read_text().replace("ENSEMBLE NVT", "ENSEMBLE NVE")
        text = re.sub(
            r"    &THERMOSTAT\n.*?    &END THERMOSTAT\n", "", text, flags=re.S
        )
        text = text.replace(
            "  &CONSTRAINT\n", "  &CONSTRAINT\n    SHAKE_TOLERANCE 1e-10\n"
        )
        inp.write_text(text)
        rows = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"]
        assert abs(rows[10][4] - limit) < 1e-12
        assert rows[10][10] != 0 and all(row[10] == 0 for row in rows[11:])
        assert abs(rows[-1][5] - 14.3996 / limit**2) < 5e-3

    name = "coupled_distances_equal_masses"
    directory = args.out / "three_ions_joint_metric" / name
    generator.write_case(directory, name, 4, 0.5, 0.0, 2.5)
    inp = directory / f"{name}.inp"
    text = inp.read_text().replace("ENSEMBLE NVT", "ENSEMBLE NVE")
    text = text.replace(
        "  &END SUBSYS",
        "    &KIND Na\n      MASS 1.0\n    &END KIND\n"
        "    &KIND Cl\n      MASS 1.0\n    &END KIND\n"
        "    &COLVAR\n      &DISTANCE\n        ATOMS 2 3\n"
        "      &END DISTANCE\n    &END COLVAR\n  &END SUBSYS",
    )
    text = text.replace(
        "    &BLUE_MOON\n",
        "    &COLLECTIVE\n      COLVAR 2\n"
        "      INTERMOLECULAR\n      TARGET [angstrom] 2.5\n"
        "    &END COLLECTIVE\n    &BLUE_MOON\n",
    )
    inp.write_text(text)
    (directory / "nacl.xyz").write_text(
        "3\nCoupled ion distances\nNa 0 0 0\n"
        "Cl 2.5 0 0\nNa 3.75 2.165063509461097 0\n"
    )
    rows = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"]
    xyz = (directory / f"{name}-pos-1.xyz").read_text().splitlines()
    p = [[float(x) for x in line.split()[1:4]] for line in xyz[-3:]]
    u = [p[1][i] - p[0][i] for i in range(3)]
    v = [p[2][i] - p[1][i] for i in range(3)]
    cosine = sum(a * b for a, b in zip(u, v)) / math.sqrt(
        sum(a * a for a in u) * sum(b * b for b in v)
    )
    expected_weight = 1 / math.sqrt(4 - cosine * cosine)
    assert abs(rows[-1][6] - expected_weight) < 2e-6, "joint cross-atom metric omitted"
    assert (
        rows[-1][6] == rows[-2][6]
    ), "both constraints must use the same joint determinant"
    assert rows[-2][2:3] == [1] and rows[-1][2:3] == [2]
    print(
        "PASS: reverse growth, partial/stopped limit, coupled joint metric and stable IDs"
    )


if __name__ == "__main__":
    main()
