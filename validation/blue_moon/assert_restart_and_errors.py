#!/usr/bin/env python3
"""Print-off, EXT_RESTART, split restart and explicit failure-path assertions."""

import argparse
import os
from pathlib import Path
import re
import shutil
import subprocess

from reproduce import make_case, run_case


def restart_work(path):
    payload = re.search(
        r"(?m)^\s*STATE\s+([\s\S]*?)^\s*&END CONSTRAINT", path.read_text()
    ).group(1)
    values = [float(x) for x in payload.replace("\\", " ").split()]
    assert len(values) == 9
    return values[1]


def expect_failure(executable, directory, inp, message):
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    launcher = (
        ["mpirun", "-n", env["BLUE_MOON_MPI_RANKS"]]
        if "BLUE_MOON_MPI_RANKS" in env
        else []
    )
    result = subprocess.run(
        launcher + [str(executable), "-i", inp.name, "-o", "cp2k.out"],
        cwd=directory,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=120,
    )
    output = result.stdout + (directory / "cp2k.out").read_text()
    assert result.returncode != 0 and message in output, output[-2000:]
    (directory / "expected_failure.log").write_text(output)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    args = parser.parse_args()
    directory, inp = make_case(args.out, "work_continuous_reference")
    reference = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"][-1]
    reference_work = restart_work(directory / "work_continuous_reference-1.restart")
    directory, inp = make_case(args.out, "work_print_keys_disabled")
    inp.write_text(
        inp.read_text()
        .replace("&BLUE_MOON\n", "&BLUE_MOON OFF\n")
        .replace("&COLVAR_SLOW_GROWTH\n", "&COLVAR_SLOW_GROWTH OFF\n")
    )
    assert not run_case(args.cp2k, directory, inp)
    assert (
        abs(
            restart_work(directory / "work_print_keys_disabled-1.restart")
            - reference_work
        )
        < 1e-12
    )

    # Ordinary constraints must not acquire joint-metric checks or restart state.
    directory, inp = make_case(
        args.out, "absent_diagnostics_duplicate_fixed_cv", steps=4, growth=None
    )
    text = inp.read_text().replace("PRINT_LEVEL LOW", "PRINT_LEVEL HIGH")
    text = text.replace(
        "    &BLUE_MOON\n",
        "    &COLLECTIVE\n      COLVAR 1\n"
        "      INTERMOLECULAR\n      TARGET [angstrom] 2.5\n"
        "    &END COLLECTIVE\n    &BLUE_MOON\n",
    )
    for key in ("BLUE_MOON", "COLVAR_SLOW_GROWTH"):
        text, count = re.subn(rf"    &{key}\n.*?    &END {key}\n", "", text, flags=re.S)
        assert count == 1
    inp.write_text(text)
    assert not run_case(args.cp2k, directory, inp)
    assert (
        "BLUE_MOON_RESTART"
        not in (
            directory / "absent_diagnostics_duplicate_fixed_cv-1.restart"
        ).read_text()
    )

    first, first_inp = make_case(args.out, "ext_restart_first_twenty", steps=20)
    run_case(args.cp2k, first, first_inp)
    restart = first / "ext_restart_first_twenty-1.restart"
    directory, inp = make_case(args.out, "ext_restart_next_twenty", steps=20)
    inp.write_text(
        "&EXT_RESTART\n  RESTART_FILE_NAME "
        + str(restart.resolve())
        + "\n  RESTART_DEFAULT T\n&END EXT_RESTART\n"
        + inp.read_text()
    )
    final = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"][-1]
    assert final[0] == 40 and abs(final[11] - reference[11]) < 1e-8

    directory, inp = make_case(args.out, "ext_restart_reset_counters", steps=20)
    inp.write_text(
        "&EXT_RESTART\n  RESTART_FILE_NAME "
        + str(restart.resolve())
        + "\n  RESTART_DEFAULT T\n  RESTART_COUNTERS F\n&END EXT_RESTART\n"
        + inp.read_text()
    )
    final = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"][-1]
    assert final[0] == 20 and abs(final[11] - reference[11]) < 1e-8

    first, first_inp = make_case(args.out, "split_restart_first_twenty", steps=20)
    first_inp.write_text(
        first_inp.read_text().replace(
            "&END MOTION",
            "  &PRINT\n    &RESTART\n"
            "      SPLIT_RESTART_FILE T\n    &END RESTART\n  &END PRINT\n&END MOTION",
        )
    )
    run_case(args.cp2k, first, first_inp)
    directory, inp = make_case(args.out, "split_restart_next_twenty", steps=20)
    for path in first.glob("*restart*"):
        if path.is_file():
            shutil.copy2(path, directory / path.name)
    restart = first / "split_restart_first_twenty-1.restart"
    text = re.sub(
        r"(?m)^(\s*PROJECT_NAME\s+).+$", r'\1"split_continuation"', restart.read_text()
    )
    binary = first / "split_restart_first_twenty-1.restart.bin"
    inp.write_text(
        "&EXT_RESTART\n  RESTART_FILE_NAME "
        + str(restart.resolve())
        + "\n  BINARY_RESTART_FILE_NAME "
        + str(binary.resolve())
        + "\n  RESTART_DEFAULT T\n&END EXT_RESTART\n"
        + text
    )
    final = run_case(args.cp2k, directory, inp)["blueMoonLog"]["rows"][-1]
    assert final[0] == 40 and abs(final[11] - reference[11]) < 1e-8

    directory, inp = make_case(args.out, "mismatched_restart_id", steps=20)
    inp.write_text(
        re.sub(r"(?m)^(\s*ID\s+)1 0 1$", r"\g<1>99 0 1", restart.read_text())
    )
    for path in first.glob("*restart*"):
        if path.is_file():
            shutil.copy2(path, directory / path.name)
    expect_failure(args.cp2k, directory, inp, "Blue Moon restart constraint ID changed")

    directory, inp = make_case(
        args.out, "dependent_duplicate_constraints", steps=4, growth=None
    )
    inp.write_text(
        inp.read_text().replace(
            "    &BLUE_MOON\n",
            "    &COLLECTIVE\n      COLVAR 1\n"
            "      INTERMOLECULAR\n      TARGET [angstrom] 2.5\n"
            "    &END COLLECTIVE\n    &BLUE_MOON\n",
        )
    )
    expect_failure(args.cp2k, directory, inp, "Dependent")

    directory, inp = make_case(
        args.out, "incompatible_cv_print_unit", steps=4, growth=None
    )
    inp.write_text(
        inp.read_text().replace(
            "    &BLUE_MOON\n", "    &BLUE_MOON\n      CV_UNIT deg\n"
        )
    )
    expect_failure(args.cp2k, directory, inp, "CV_UNIT is incompatible")
    print(
        "PASS: print-off work, EXT/split restarts, restart mismatch and invalid metric/units"
    )


if __name__ == "__main__":
    main()
