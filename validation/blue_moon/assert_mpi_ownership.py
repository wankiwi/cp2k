#!/usr/bin/env python3
"""Four local molecule constraints coupled to one global CV, on 1/2/4 ranks."""

import argparse
import json
import math
from pathlib import Path

from assert_edge_cases import nve_case
from reproduce import run_case


PSF = """PSF EXT

         1 !NTITLE
NaCl dimer for MPI ownership validation

         2 !NATOM
{atoms}

         1 !NBOND
         1         2

         0 !NTHETA

         0 !NPHI

         0 !NIMPHI

         0 !NDON

         0 !NACC

         0 !NNB
         0         0

         1 !NGRP NST2
         0         0         0
"""


def determinant(a):
    a = [row[:] for row in a]
    result = 1.0
    for i in range(len(a)):
        pivot = a[i][i]
        assert pivot > 0
        result *= pivot
        for j in range(i + 1, len(a)):
            scale = a[j][i] / pivot
            for k in range(i + 1, len(a)):
                a[j][k] -= scale * a[i][k]
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    name = "four_dimers_local_global_coupling"
    directory, inp = nve_case(args.out, name, steps=4)
    text = inp.read_text()
    text = text.replace(
        "    &FORCEFIELD\n",
        "    &FORCEFIELD\n      &BOND\n"
        "        ATOMS Na Cl\n        KIND HARMONIC\n        K [eV*angstrom^-2] 1.0\n"
        "        R0 [angstrom] 2.0\n      &END BOND\n",
    )
    text = text.replace(
        "CONN_FILE_FORMAT OFF",
        "CONN_FILE_FORMAT MOL_SET\n      &MOL_SET\n"
        "        &MOLECULE\n          NMOL 4\n          CONN_FILE_FORMAT PSF\n"
        "          CONN_FILE_NAME dimer.psf\n        &END MOLECULE\n      &END MOL_SET",
    )
    text = text.replace(
        "  &END SUBSYS",
        "    &KIND Na\n      MASS 1.0\n    &END KIND\n"
        "    &KIND Cl\n      MASS 1.0\n    &END KIND\n"
        "    &COLVAR\n      &DISTANCE\n        ATOMS 1 3\n      &END DISTANCE\n"
        "    &END COLVAR\n  &END SUBSYS",
    )
    text = text.replace("      INTERMOLECULAR\n", "      MOLECULE 1\n")
    text = text.replace(
        "    &BLUE_MOON\n",
        "    &COLLECTIVE\n      COLVAR 2\n"
        "      INTERMOLECULAR\n      TARGET [angstrom] 5.0\n"
        "    &END COLLECTIVE\n    &BLUE_MOON\n",
    )
    inp.write_text(text)
    atoms = "\n".join(
        f"{i:10d} {'ION':7s}  {1:8d} {'ION':7s}  {symbol:6s}  "
        f"{symbol:6s}{charge:10.6f}      {1.0:8.3f}           0"
        for i, symbol, charge in [(1, "Na", 1.0), (2, "Cl", -1.0)]
    )
    (directory / "dimer.psf").write_text(PSF.format(atoms=atoms))
    (directory / "nacl.xyz").write_text(
        "8\nFour NaCl dimers\nNa 0 0 0\nCl 2.5 0 0\n"
        "Na 3 4 0\nCl 5.5 4 0\nNa 0 10 0\nCl 2.5 10 0\n"
        "Na 0 15 0\nCl 2.5 15 0\n"
    )
    rows = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"]
    assert len(rows) == 20
    assert [r[2] for r in rows[-5:]] == [1, 2, 3, 4, 5]
    assert all(
        abs(r[5]) > 1e-3 for r in rows[-5:-3]
    ), "remote-owned multipliers missing"
    assert all(
        abs(r[4] - 2.502) < 1e-10 for r in rows[-5:-1]
    ), "target advanced per molecule"
    assert all(
        abs(r[3] - r[4]) < 1e-9 for r in rows
    ), "global correction invalidated local CV"
    lines = (directory / f"{name}-pos-1.xyz").read_text().splitlines()[-8:]
    positions = [[float(v) for v in line.split()[1:4]] for line in lines]
    gradients = []
    for a, b in [(0, 1), (2, 3), (4, 5), (6, 7), (0, 2)]:
        vector = [positions[b][i] - positions[a][i] for i in range(3)]
        norm = math.sqrt(sum(v * v for v in vector))
        g = [[0.0] * 3 for _ in positions]
        g[a] = [-v / norm for v in vector]
        g[b] = [v / norm for v in vector]
        gradients.append(g)
    metric = [
        [
            sum(gi[a][d] * gj[a][d] for a in range(8) for d in range(3))
            for gj in gradients
        ]
        for gi in gradients
    ]
    expected_weight = 1 / math.sqrt(determinant(metric))
    assert all(
        abs(r[6] - expected_weight) < 2e-7 for r in rows[-5:]
    ), "local/global metric cross terms missing"
    if args.reference:
        reference = json.loads(args.reference.read_text())
        assert len(reference) == len(rows)
        assert (
            max(abs(a - b) for ra, rb in zip(rows, reference) for a, b in zip(ra, rb))
            < 1e-9
        )
    result = args.out / "ownership_rows.json"
    result.write_text(json.dumps(rows, indent=2) + "\n")
    print(
        "PASS: four molecule owners, one target advance, local/global cross metric and MPI equivalence"
    )


if __name__ == "__main__":
    main()
