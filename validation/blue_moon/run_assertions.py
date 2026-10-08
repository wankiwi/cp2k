#!/usr/bin/env python3
"""Run the short scientific assertions with an explicit binary and provenance."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    binary = args.cp2k.resolve(strict=True)
    scripts = [
        "assert_moving_constraint",
        "assert_diagnostics",
        "assert_extended",
        "assert_restart_and_errors",
        "assert_edge_cases",
        "assert_mpi_ownership",
    ]
    env = {
        **os.environ,
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    if args.mpi_ranks is not None:
        if args.mpi_ranks < 1:
            parser.error("--mpi-ranks must be positive")
        env["BLUE_MOON_MPI_RANKS"] = str(args.mpi_ranks)
    personality = Path("/proc/self/personality")
    results = {
        "binary": str(binary),
        "sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "mpi_ranks": args.mpi_ranks,
        "assertions": {},
        "process_personality": personality.read_text().strip()
        if personality.exists()
        else None,
        "lsan_options": env.get("LSAN_OPTIONS", ""),
        "private_library_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in binary.parent.glob("libcp2k.so*")
            if p.is_file()
        },
    }
    for script in scripts:
        command = [
            sys.executable,
            str(Path(__file__).with_name(f"{script}.py")),
            "--out",
            str((args.out / script).resolve()),
            "--cp2k",
            str(binary),
        ]
        with (args.out / f"{script}.log").open("w") as log:
            result = subprocess.run(
                command, env=env, stdout=log, stderr=subprocess.STDOUT
            )
        results["assertions"][script] = result.returncode
        (args.out / "results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(f"{script}: {'PASS' if result.returncode == 0 else 'FAIL'}", flush=True)
        if result.returncode:
            raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
