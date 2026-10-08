#!/usr/bin/env python3
"""Execute the PDF NaCl sweep and TI windows, preserving all individual logs."""

import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
from zoneinfo import ZoneInfo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cp2k", type=Path, required=True)
    args = parser.parse_args()
    args.cp2k = args.cp2k.resolve(strict=True)
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}
    cases = sorted(args.out.glob("sg_*/*.inp"), reverse=True) + sorted(
        args.out.glob("ti_windows/*/*.inp")
    )
    assert len(cases) == 34, f"Expected three rates and 31 TI windows, got {len(cases)}"
    binary_hash = hashlib.sha256(args.cp2k.read_bytes()).hexdigest()
    provenance = {
        "started_SGT": datetime.now(ZoneInfo("Asia/Singapore")).isoformat(),
        "binary": str(args.cp2k),
        "binary_sha256": binary_hash,
        "inputs_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in cases
        },
    }
    (args.out / "execution_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    for inp in cases:
        assert (
            hashlib.sha256(args.cp2k.read_bytes()).hexdigest() == binary_hash
        ), "Binary changed during benchmark"
        directory = inp.parent
        print(f"START {directory}", flush=True)
        start = time.monotonic()
        with (directory / "console.log").open("w") as console:
            result = subprocess.run(
                [str(args.cp2k), "-i", inp.name, "-o", "cp2k.out"],
                cwd=directory,
                env=env,
                stdout=console,
                stderr=subprocess.STDOUT,
            )
        assert result.returncode == 0, f"CP2K failed: {directory}"
        output = (directory / "cp2k.out").read_text()
        assert "PROGRAM ENDED AT" in output, f"Incomplete calculation: {directory}"
        print(f"DONE {directory}, elapsed={time.monotonic() - start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
