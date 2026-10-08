# Moving COLVAR: Coulomb-only NaCl controls

These artifacts validate the moving-constraint fix against official master
`d59727108bfd58f66da98f7ff772c2be91f08b08` and the bug-only implementation
`bc3e6119ccde643b08f5098429c95a3fd8ab8df0`. They are kept on a personal evidence branch, outside the
official PR diff.

The model uses `METHOD FIST`, charges +1 and -1, zero Lennard-Jones `EPSILON` for every pair, and
`CONN_FILE_FORMAT OFF`. Thus only Coulomb interactions are enabled. `EWALD_TYPE NONE` selects the
direct-space electrostatic interaction in a 40 Å periodic XYZ cell; the pair cutoff is 20 Å.

All cases use NVE, zero initial velocity, a 1 fs timestep and `SHAKE_TOLERANCE 1e-10`. The distance
target grows from 2.5 to 5 Å:

| Growth (Å/fs) | Steps | Master force MAE (eV/Å) | Fixed force MAE (eV/Å) | Master final work error (eV) | Fixed final work error (eV) |
| ------------- | ----: | ----------------------: | ---------------------: | ---------------------------: | --------------------------: |
| 0.00025       | 10000 |                0.722889 |            0.000176328 |                    +1.807222 |                +0.000440819 |
| 0.0005        |  5000 |                1.445774 |            0.000349123 |                    +3.614435 |                +0.000872808 |
| 0.001         |  2500 |                2.891544 |            0.000694715 |                    +7.228861 |                +0.001736787 |

Every master/fixed pair uses identical input bytes. The input `TARGET_LIMIT` is 10 Å; each
simulation ends at 5 Å because of its prescribed step count. The 5000-step example is in
`inputs/growth_0p0005_A_per_fs/fixed/nacl.inp`, with `nacl.xyz` in the same directory.

The force is the negative native SHAKE multiplier, converted from atomic units using the frozen CP2K
constants in `PROVENANCE.json`. Its analytic reference is `14.3996 / r²` eV/Å. Mechanical work is
the discrete sum `Σ F_n (r_target,n − r_target,n−1)`, including the first interval. The analytic
work uses the same native target grid and sum at each rate. It is not an equilibrium free-energy
estimate. All 35000 master/fixed MD frames are retained; no force offsets, smoothing or frame
exclusions are applied.

The two-panel figure overlays all three rates. Master uses solid lines, Fixed uses short dashes,
long dashes and dash-dot lines, and the analytic references use dark open circles. All master/fixed
samples are drawn; analytic markers are displayed at sparse positions for readability. PNG, PDF and
SVG exports are in `figures/`.

## Reproduce

Build unmodified master and the bug-only implementation as separate CP2K executables. Python needs
NumPy, and figure export also needs Matplotlib and Pillow. Use a new output directory:

```sh
python reproduce_nacl.py --master /path/to/master/cp2k.ssmp \
  --fixed /path/to/fixed/cp2k.ssmp --out /path/to/nacl-controls
python figures/gen_nacl_growth_rate_figure.py \
  --data-dir /path/to/nacl-controls/data --out /path/to/nacl-figures
```

The parser validates completed calculations, all native position/velocity frames, the target
schedule and the fixed radial velocity. It reads one native SHAKE record per NVE step and writes the
CSVs consumed by the plotting script. The supplied parser was checked against all six original
trajectories (`analysis_validation.json`). `--cp2k-data-dir` can set CP2K's data directory;
`--disable-aslr` applies the Linux Debug/LSAN startup workaround if needed. No new CP2K print key is
required.

The short official regression uses existing `LAGRANGE_MULTIPLIERS` output and checks fixed, forward,
reverse and slow/default-tolerance cases. The fixed case passes unmodified master; all three moving
cases fail its native force matcher. Additional native controls exercise initialization,
partial/final target limits, RESPA inner intervals, ROLL, shared molecule kinds, overlapping
local/global constraints and text restart. Numerical records and build provenance are in
`PROVENANCE.json`; MPI numerical checks do not imply MPI leak-detector coverage.
