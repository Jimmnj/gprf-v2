# Galaxy dynamics from Newtonian gravity and a non-baryonic mass set by the baryons

Source, analysis scripts and pre-registrations for

> J. R. Marsen (2026), *Galaxy dynamics from Newtonian gravity and a non-baryonic mass set by the baryons*,
> v2.0.0, Zenodo, [doi:10.5281/zenodo.23748074](https://doi.org/10.5281/zenodo.23748074)
> (all versions: [10.5281/zenodo.23748073](https://doi.org/10.5281/zenodo.23748073)).

## Layout

| Path | Contents |
|---|---|
| `paper/` | The posted PDF and its LaTeX source; `Marsen_2026_v2.0.0_annotated.tex` is the same text with the provenance tags (`[RESULT_...]`, `[script.py]`) that link each number to the record below; `figs_v2/` the figures |
| `notebooks/` | Scripts (`*.py`), their outputs (`*.csv`), the result records (`RESULT_*.md`), the pre-registrations (`PREREG_*.md`, `RESULT_DWARF_PREREG_*.md`) and one derivation note |
| `data/` | Small published tables used in the environment and lensing tests (see `data/README.md`) |
| `tools/fetch_sparc.py` | Downloads SPARC (not redistributed) into `sparc_galaxy_data/` |

## Reproducing

Python 3.10+, with `numpy`, `scipy`, `pandas`, `matplotlib`.

```bash
python3 tools/fetch_sparc.py          # SPARC into sparc_galaxy_data/
cd notebooks
python3 pass3_relation3_20261022.py   # e.g. the acceleration relation with relation (3)
```

Scripts are run from `notebooks/` and read `../sparc_galaxy_data/` and `../data/`. Several import shared code from
`gate2c_robustness_20260925.py`, `n9_dwarfs_20261003.py` or `pass2_sparc_full_20261022.py`.
The pre-registered runs take an argument, e.g. `python3 n16_closure_20261028.py S1` (parts S1, S3, S2, S2k) and
`python3 sform_family_20261028.py univ` (or `rel3`). Figure scripts write to `notebooks/figs_v2/`; the paper's copies are in `paper/figs_v2/`.
Two runs were re-executed from this layout before release and reproduce the archived reports exactly
(`n16_closure_20261028.py S1`, `pass3_relation3_20261022.py`).

## Where each section's numbers come from

| Section | Records and scripts |
|---|---|
| 2.1 Statement | `RESULT_N10`, `RESULT_SFORM_FAMILY_20261028` |
| 3.1 A cored isothermal response with BTFR normalisation reproduces the acceleration relation; the transition shape follows from the self-consistent solve | `RESULT_PASS3_RELATION3_20261022`, `RESULT_PASS2_SPARC_FULL_20261022`, `RESULT_N12_EMERGENT_NU`, `RESULT_SFORM_FAMILY_20261028`, `RESULT_N12` |
| 3.2 Individual rotation curves, and the comparison with the acceleration relation | `pass2_table2_full_20261022.py`, `RESULT_MW_DECLINE_20261021` |
| 3.3 The core surface density | `RESULT_TEST1`, `fig4_mu0_mond_20261027.py`, `RESULT_FIG4_MU0_MOND_20261027` |
| 4.1 Surface-density dependence at fixed mass | `pass4_oman_drift_20261022.py`, `RESULT_PASS4_OMAN_DRIFT_20261022`, `RESULT_PASS5_SIGMA_FREE_OVERSHOOT_20261022`, `RESULT_PASS6_LOCAL_SIGMA_20261022`, `RESULT_PASS2_SPARC_FULL_20261022` |
| 4.2 Crater II | `RESULT_RHOC_SIGMA_B_OOS`, `RESULT_DWARF_PREREG`, `RESULT_RHOC`, `RESULT_CRATER2_TIDES_20261021` |
| 4.4 Environment | `RESULT_N10` |
| 6.1 Wide binaries | `pass9_wb_rhoc_sweep_20261025.py`, `RESULT_SIMROUND_C1-C4`, `RESULT_WB_R1PC`, `RESULT_FRM3_WB_FROZEN_LAW` |
| 6.5 Tests in reach | `pass11_lensing_forecast_20261025.py` |
| 7. Open list | `RESULT_CLUSTER_RELEASE_20261022`, `lagrangian_checks_20261021.py`, `DERIVATION_ACTION_20261021` |
| A.1 The Tully-Fisher relation and curve features | `RESULT_OWED_VFLAT_RENZO`, `RESULT_TEST1`, `RESULT_MW_DECLINE_20261021`, `mw_decline_20261021.py` |
| A.2 Dwarf spheroidals and globular clusters | `RESULT_N9` |
| B.1 What sets rho_c | `RESULT_N16_CLOSURE_20261028`, `RESULT_RHOC_SIGMA_B_OOS` |
| B.2 Dwarf dispersions, pre-registered | `RESULT_N16_CLOSURE_20261028`, `RESULT_CRATER2_TIDES_20261021` |
| Appendix C. Lensing: compatibility, truncation, satellites and the mass budget | `RESULT_N10`, `pass10_lensing_fixedVf_20261025.py`, `RESULT_SIMROUND_C1-C4`, `RESULT_LENSING_FIT_SATELLITES`, `RESULT_N11_PASS1_20261026`, `PREREG_N11_20261026` |
| Appendix D. The environment test in full | `RESULT_N8B`, `RESULT_N8D` |
| E.1 Cloud-9 | `RESULT_BODY_RULE_20261022` |

The `RESULT_*.md` files are the working records written when each run was made. They name the script,
the inputs and the output files, and they are kept as written. A `[UNRATIFIED]` mark means the record was
written before the author ruled on it; the paper states only what was ruled. Some records refer to internal
files not included here (a ruling log and session notes); the numbers they report come from the scripts and
outputs that are included.

## Pre-registrations

`PREREG_N16_20261028.md` (closure of the central density through MOND's core surface density; FAIL),
`PREREG_SFORM_20261028.md` (shape of the saturation; 'immaterial' FAILS at 0.04 dex),
`PREREG_N11_20261026.md` (lensing field test), and the dwarf protocol in `RESULT_DWARF_PREREG_20261020.md` /
`dwarf_prereg_20261020.py`. Each was fixed before its run; the outcomes are in the matching `RESULT_*.md`.

## AI usage

See the paper's AI-usage statement. Several records name the assistant that ran them.

## Licence

CC0 1.0 (see `LICENSE`). Third-party data keep their own terms; cite the original papers.
