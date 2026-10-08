#!/usr/bin/env python3
"""Run short, isolated reproductions against the currently built CP2K."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
GENERATOR = Path(__file__).with_name("prepare_inputs.py")
spec = importlib.util.spec_from_file_location("nacl_input_generator", GENERATOR)
generator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generator)


def parse_log(path):
    metadata = {}
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            key, sep, value = line[1:].partition(" = ")
            if sep:
                metadata[key.strip()] = value.strip()
            continue
        fields = line.split()
        if fields:
            rows.append([float(value) for value in fields])
    return {"metadata": metadata, "rows": rows}


def run_case(executable, directory, input_file):
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    with (directory / "console.log").open("w") as console:
        result = subprocess.run(
            (
                ["mpirun", "-n", os.environ["BLUE_MOON_MPI_RANKS"]]
                if "BLUE_MOON_MPI_RANKS" in os.environ
                else []
            )
            + [str(executable), "-i", input_file.name, "-o", "cp2k.out"],
            cwd=directory,
            env=env,
            stdout=console,
            stderr=subprocess.STDOUT,
            timeout=120,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"CP2K exit {result.returncode}: {directory / 'cp2k.out'}")
    logs = {}
    for extension in ("blueMoonLog", "colvarSlowGrowthLog"):
        candidates = sorted(directory.glob(f"*-1.{extension}"))
        if candidates:
            logs[extension] = parse_log(candidates[0])
    return logs


def make_case(outdir, name, steps=40, stride=1, growth=0.005):
    directory = outdir / "nacl_point_charges_1k" / name
    directory.mkdir(parents=True, exist_ok=False)
    generator.write_case(
        directory,
        name,
        steps,
        1.0,
        1.0,
        2.5,
        growth=growth,
        limit=10.0,
        print_stride=stride,
    )
    return directory, directory / f"{name}.inp"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, default=ROOT / "exe/minimal/cp2k.ssmp")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    results = {}

    for name, stride in (
        ("slow_growth_print_every_step", 1),
        ("slow_growth_print_every_five_steps", 5),
    ):
        directory, input_file = make_case(args.out, name, stride=stride)
        results[name] = run_case(args.cp2k, directory, input_file)

    name = "slow_growth_alias_only"
    directory, input_file = make_case(args.out, name)
    text = re.sub(
        r"    &BLUE_MOON\n.*?    &END BLUE_MOON\n",
        "",
        input_file.read_text(),
        count=1,
        flags=re.S,
    )
    input_file.write_text(text)
    results[name] = run_case(args.cp2k, directory, input_file)

    name = "slow_growth_restart_first_twenty_steps"
    directory, input_file = make_case(args.out, name, steps=20)
    results[name] = run_case(args.cp2k, directory, input_file)
    restart = directory / f"{name}-1.restart"
    restart_text = restart.read_text()
    continuation = (
        args.out / "nacl_point_charges_1k/slow_growth_restart_next_twenty_steps"
    )
    continuation.mkdir(parents=True, exist_ok=False)
    restart_text = re.sub(
        r"(?m)^(\s*PROJECT_NAME\s+).+$", r'\1"continuation"', restart_text
    )
    continuation_input = continuation / "continuation.inp"
    continuation_input.write_text(restart_text)
    (continuation / "nacl.xyz").write_text((directory / "nacl.xyz").read_text())
    results["slow_growth_restart_next_twenty_steps"] = run_case(
        args.cp2k, continuation, continuation_input
    )

    for name, target in (
        ("fixed_distance_target_angstrom", "[angstrom] 2.50000000"),
        ("fixed_distance_target_atomic_units", "4.72431531156543"),
        ("fixed_distance_target_meters", "[m] 2.50000000E-10"),
    ):
        directory, input_file = make_case(args.out, name, steps=4, growth=None)
        input_file.write_text(
            input_file.read_text().replace("[angstrom] 2.50000000", target)
        )
        results[name] = run_case(args.cp2k, directory, input_file)

    name = "fixed_distance_and_angle_input_units"
    directory, input_file = make_case(args.out, name, steps=4, growth=None)
    text = input_file.read_text().replace(
        "  &END SUBSYS",
        "    &COLVAR\n      &ANGLE\n        ATOMS 1 2 3\n"
        "      &END ANGLE\n    &END COLVAR\n  &END SUBSYS",
    )
    text = text.replace(
        "    &BLUE_MOON\n",
        "    &COLLECTIVE\n      COLVAR 2\n"
        "      INTERMOLECULAR\n      TARGET [deg] 90.0\n"
        "    &END COLLECTIVE\n    &BLUE_MOON\n",
    )
    input_file.write_text(text)
    (directory / "nacl.xyz").write_text(
        "3\nNa-Cl-Na right angle\nNa 0.0 0.0 0.0\nCl 2.5 0.0 0.0\nNa 2.5 2.5 0.0\n"
    )
    results[name] = run_case(args.cp2k, directory, input_file)

    for energy_unit in ("eV", "hartree"):
        name = f"fixed_angle_internal_cv_energy_{energy_unit}"
        directory, input_file = make_case(args.out, name, steps=4, growth=None)
        text = input_file.read_text().replace(
            "&DISTANCE\n        ATOMS 1 2\n      &END DISTANCE",
            "&ANGLE\n        ATOMS 1 2 3\n      &END ANGLE",
        )
        text = text.replace("TARGET [angstrom] 2.50000000", "TARGET 1.5707963267948966")
        text = text.replace("ENERGY_UNIT eV", f"ENERGY_UNIT {energy_unit}")
        input_file.write_text(text)
        (directory / "nacl.xyz").write_text(
            "3\nNa-Cl-Na right angle\nNa 0.0 0.0 0.0\nCl 2.5 0.0 0.0\nNa 2.5 2.5 0.0\n"
        )
        results[name] = run_case(args.cp2k, directory, input_file)

    report = args.out / "reproduction_results.json"
    report.write_text(json.dumps(results, indent=2) + "\n")
    for name, logs in results.items():
        print(name)
        for extension, data in logs.items():
            rows = data["rows"]
            print(
                extension,
                "unit",
                data["metadata"].get("cv_unit"),
                "rows",
                len(rows),
                "first",
                rows[0],
                "last",
                rows[-1],
            )
    print(report)


if __name__ == "__main__":
    main()
