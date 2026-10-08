# Moving COLVAR constraints and Blue Moon diagnostics

The short checks exercise the imposed velocity, moving target limits, shared molecule-kind targets,
joint multi-CV mass metrics, units, work accumulation, restarts, RESPA and MPI ownership. Run each
invocation in a fresh output directory.

Build CP2K with its CMake build system and pass the resulting executable explicitly:

```bash
python3 validation/blue_moon/run_assertions.py \
  --out /path/to/validation/checked_debug \
  --cp2k /path/to/build/bin/cp2k.sdbg
python3 validation/blue_moon/run_assertions.py \
  --out /path/to/validation/mpi_2_ranks \
  --cp2k /path/to/build/bin/cp2k.pdbg --mpi-ranks 2
```

Repeat the MPI invocation with 1 and 4 ranks when checking ownership invariance. The driver records
the executable SHA-256 and each of the six suite exit codes. The native regression inputs are in
`tests/Fist/regtest-blue-moon/`.

A fixed distance constraint imposes `grad(CV).v = 0`. A moving distance constraint imposes the
actual target increment divided by the integration interval, including a partial final interval at
TARGET_LIMIT. The NaCl point-charge controls compare the emitted radial velocity with that rate and
the negative SHAKE multiplier with `14.3996/r**2` eV/angstrom. The force acceptance tolerance is
0.005 eV/angstrom.

BLUE_MOON and COLVAR_SLOW_GROWTH serialize the same per-step state. Their twelve columns are
`step time constraint value target lambda inv_sqrt_z gkt blue_moon_force dxi dW W`. The constraint
index maps to a stable input/molecule/local identity in the header. CV units inherit each TARGET
unless CV_UNIT overrides them. Headers specify per-constraint units for heterogeneous sets.

The joint metric includes all active COLVAR constraints and shared-atom cross terms. `inv_sqrt_z` is
`det(Z)**(-1/2)` and `blue_moon_force` is the weighted numerator `inv_sqrt_z*(lambda+gkt)`. The
equilibrium TI estimator is the ratio of the mean numerator to the mean weight. W is mechanical
pulling work, summed over all integration intervals, including the first, and preserved in restart
state. Printing less often or using both aliases does not change that state. Explicit OFF keys
retain bookkeeping; absent diagnostic keys and restart state leave the extension inactive.

Moving targets and explicitly requested extension bookkeeping enable joint local/global convergence
checks. This includes explicit OFF print keys and Blue Moon restart state. Ordinary fixed-target
calculations retain the legacy iteration. Enabling the extension can change mixed local/global
fixed-constraint trajectories when global corrections invalidate previously converged local CVs.

The three-rate NaCl reproduction inputs are in `nacl_growth_rate_controls/`. The corresponding force
and mechanical-work comparison, numerical tables and regeneration script are in `figures/`.
