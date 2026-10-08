#!/usr/bin/env python3
import argparse
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent

TEMPLATE = """&GLOBAL
  PROJECT {project}
  RUN_TYPE MD
  PRINT_LEVEL LOW
&END GLOBAL

&FORCE_EVAL
  METHOD FIST
  &MM
    &FORCEFIELD
      &CHARGE
        ATOM Na
        CHARGE 1.0
      &END CHARGE
      &CHARGE
        ATOM Cl
        CHARGE -1.0
      &END CHARGE
      &SPLINE
        RCUT_NB [angstrom] 20.0
      &END SPLINE
      &NONBONDED
        &LENNARD-JONES
          ATOMS Na Na
          EPSILON [eV] 0.0
          SIGMA [angstrom] 1.0
          RCUT [angstrom] 20.0
        &END LENNARD-JONES
        &LENNARD-JONES
          ATOMS Na Cl
          EPSILON [eV] 0.0
          SIGMA [angstrom] 1.0
          RCUT [angstrom] 20.0
        &END LENNARD-JONES
        &LENNARD-JONES
          ATOMS Cl Cl
          EPSILON [eV] 0.0
          SIGMA [angstrom] 1.0
          RCUT [angstrom] 20.0
        &END LENNARD-JONES
      &END NONBONDED
    &END FORCEFIELD
    &POISSON
      &EWALD
        EWALD_TYPE NONE
      &END EWALD
    &END POISSON
  &END MM
  &SUBSYS
    &CELL
      ABC 40.0 40.0 40.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME nacl.xyz
      CONN_FILE_FORMAT OFF
    &END TOPOLOGY
    &COLVAR
      &DISTANCE
        ATOMS 1 2
      &END DISTANCE
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL

&MOTION
  &MD
    ENSEMBLE NVT
    STEPS {steps}
    TIMESTEP {dt}
    TEMPERATURE {temp}
    &THERMOSTAT
      REGION GLOBAL
      TYPE NOSE
      &NOSE
        LENGTH 3
        MTS 2
        TIMECON 50.0
        YOSHIDA 3
      &END NOSE
    &END THERMOSTAT
  &END MD
  &CONSTRAINT
    &COLLECTIVE
      COLVAR 1
      INTERMOLECULAR
      TARGET [angstrom] {target:.8f}
{growth_block}    &END COLLECTIVE
    &BLUE_MOON
      &EACH
        MD {print_stride}
      &END EACH
      TIME_UNIT fs
      ENERGY_UNIT eV
    &END BLUE_MOON
    &COLVAR_SLOW_GROWTH
      &EACH
        MD {print_stride}
      &END EACH
      TIME_UNIT fs
      ENERGY_UNIT eV
    &END COLVAR_SLOW_GROWTH
  &END CONSTRAINT
&END MOTION
"""


def write_case(
    path, project, steps, dt, temp, target, growth=None, limit=None, print_stride=1
):
    growth_block = ""
    if growth is not None:
        growth_block = (
            f"      TARGET_GROWTH [angstrom*fs^-1] {growth:.10f}\n"
            f"      TARGET_LIMIT [angstrom] {limit:.8f}\n"
        )
    text = TEMPLATE.format(
        project=project,
        steps=steps,
        dt=dt,
        temp=temp,
        target=target,
        growth_block=growth_block,
        print_stride=print_stride,
    )
    path.mkdir(parents=True, exist_ok=True)
    (path / f"{project}.inp").write_text(text)
    (path / "nacl.xyz").write_text((ROOT / "nacl.xyz").read_text())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "generated"))
    ap.add_argument("--full-range", action="store_true")
    ap.add_argument("--slow-growth-print-stride", type=int, default=1)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    for growth in (5.0e-4, 1.0e-4, 1.0e-5):
        name = f"sg_{growth:.0e}".replace("-", "m")
        steps = 1500
        if args.full_range:
            steps = math.ceil((10.0 - 2.5) / (growth * 1.0))
        write_case(
            out / name,
            name,
            steps=steps,
            dt=1.0,
            temp=1.0,
            target=2.5,
            growth=growth,
            limit=10.0,
            print_stride=args.slow_growth_print_stride,
        )
    for i, target in enumerate([2.5 + 0.25 * j for j in range(31)]):
        write_case(
            out / "ti_windows" / f"r_{target:.2f}",
            f"ti_{i:03d}",
            steps=50,
            dt=1.0,
            temp=1.0,
            target=target,
        )


if __name__ == "__main__":
    main()
