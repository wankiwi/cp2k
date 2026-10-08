# NaCl moving-constraint growth-rate controls

These inputs were executed against pristine master
`8a1319df473a130350026b9c8f4392d2a6d93b6e` and the fixed solver.
Both use the same two-ion FIST Coulomb model, zero initial velocities, NVE,
a 1 fs timestep and SHAKE tolerance 1e-10. The Na-Cl target advances from
2.5 to 5.0 angstrom.

| Growth (angstrom/fs) | Steps | Duration (ps) | Input directory |
| ---: | ---: | ---: | --- |
| 0.00025 | 10000 | 10 | growth_0p00025_A_per_fs |
| 0.0005 | 5000 | 5 | growth_0p0005_A_per_fs |
| 0.001 | 2500 | 2.5 | growth_0p001_A_per_fs |

All three exact executed input pairs are supplied. The middle-rate result is
reused from the completed 5000-step control. Only PROJECT, STEPS and
TARGET_GROWTH vary between rates for a given solver.

Each directory contains the common geometry and the exact executed input
for each solver. NaCl_master.inp uses existing native constraint printing;
NaCl_fixed.inp additionally prints Blue Moon and slow-growth diagnostics.
Their physical settings are identical.

For each variant, copy its input and nacl.xyz into a separate fresh working
directory, then run the corresponding CP2K executable there:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /path/to/cp2k.ssmp \
  -i NaCl_master.inp -o cp2k.out
```

Use NaCl_fixed.inp with the fixed executable. Supply the executable's own
shared-library directory through LD_LIBRARY_PATH if required by the build.
Do not run directly in these frozen input directories.

Acceptance checks verify normal termination, contiguous position, velocity,
energy and SHAKE-multiplier records, and the target constraint residual.
The fixed radial velocity must match TARGET_GROWTH within 1e-7 angstrom/fs;
master reproduces the near-zero radial velocity. Native and diagnostic
forces/work are independently compared.

Force is the negative native SHAKE multiplier converted to eV/angstrom.
Work is the discrete sum of force times target increment, including the
first interval from 2.5 angstrom. The analytic force is 14.3996/r^2
eV/angstrom; reference work uses each rate's own target increments.
All frames enter the comparison, without smoothing or fitted offsets.
The numerical tables and regeneration script are in ../figures.
