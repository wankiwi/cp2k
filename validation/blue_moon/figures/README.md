# NaCl moving-COLVAR constraint comparison

The frozen numerical tables compare master and the fixed solver against the
analytic Coulomb reference. Plotting does not rerun CP2K.

`nacl_growth_rate_comparison.png` overlays all three growth rates in two shared
panels: raw constraint force and accumulated mechanical work. Every measured
frame is drawn, without smoothing, fitted offsets or exclusions. All cases use
an integration interval of 1 fs and targets from 2.5 to 5.0 angstrom:

| Growth (angstrom/fs) | Steps | Duration (ps) |
| ---: | ---: | ---: |
| 0.00025 | 10000 | 10 |
| 0.0005 | 5000 | 5 |
| 0.001 | 2500 | 2.5 |

Color identifies the growth rate. Master uses solid lines; Fixed uses
short-dashed, long-dashed and dash-dot lines, respectively. The analytic
reference uses dark open circles without a connecting line. Reference markers
sample native target values for readability; the measured curves retain every
frame. The force reference is `14.3996/r**2` eV/angstrom. Reference work is each
rate's own discrete sum of force times target increment, including the first
interval.

The CSV files in `data/` contain the force and work samples, three-rate metrics
and short NaCl controls. `PROVENANCE.json` records source hashes and
transformations. Exact executed inputs for all three rates are under
`../nacl_growth_rate_controls/`.

Regenerate PNG, PDF and SVG exports into a fresh output directory:

```bash
python3 validation/blue_moon/figures/gen_nacl_growth_rate_figure.py \
  --out /path/to/fresh/figures
```

`nacl_growth_rate_export_receipt.json` records script, source-table and export
hashes, numerical metrics, dimensions and presentation settings. PNG uses an
opaque white background at 300 dpi; PDF and SVG are vector exports.
