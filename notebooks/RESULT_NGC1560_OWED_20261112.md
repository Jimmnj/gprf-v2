# RESULT_NGC1560_OWED_20261112 -- owed items (a)-(c) of SESSION_HANDOFF_20261111 sec 3.4 [RATIFIED JM 2026-11-12 01:25 (OW-1..OW-4); post hoc, not pre-registered, not scored]
Written 2026-11-12 01:08 EST (clock read). Queue item 2 of SESSION_HANDOFF_20261112. JM go 01:04.
Script: ngc1560_owed_20261112.py (expectations in its header, written 01:06 before the run); report ngc1560_owed_20261112_report.txt; run 01:06:23-01:06:25 EST in the device VM.
Inputs: ngc1560_corr_*_20261111.csv (ERRATUM_NGC1560_COLUMNS_20261111) and, for validation, ngc1560_curves_20261108.csv. POST-PROCESSING ONLY: no solver, nothing touches the law (rule 3a not triggered).
Recipes = ngc1560_feature_diag_20261111.py H4 (transmission) and H5 (shelf, step), unchanged.

## 1. (a) Transmission and shelf, CORRECTED curves
Transmission = OLS slope of detrended log V on detrended log V_bar, 2-8 kpc (45 points). Obs errors: quoted-error perturbation sd / jackknife.
| D, M/L_I, M_* (at 1.43) | obs | eq. (1) universal | eq. (1) fitted | acceleration relation |
|---|---|---|---|---|
| 3.45, 1.43, 7.0e8 (primary) | 0.354 (0.139 / 0.036) | 0.210 | 0.255 | 0.560 |
| 3.45, 1.43, 5.7e8 / 7.9e8 | same | 0.216 / 0.206 | 0.255 | 0.560 |
| 3.45, 1.00, 7.0e8 | 0.299 (0.118 / 0.031) | 0.186 | 0.211 | 0.554 |
| 2.94, 0.98, 7.0e8 | 0.296 (0.117 / 0.031) | 0.190 | 0.187 | 0.554 |
| (old, swapped M_* 1.3e9) | 0.354 | 0.188 | 0.256 | 0.560 |
Slopes in km/s per kpc of each run's radius scale; segments 4.0-6.0 (shelf) and 6.0-6.9 (step) kpc at 3.45 Mpc. Pull = (model - obs) / weighted obs error.
| D, M/L_I, M_* | shelf obs | acc. rel. | eq. (1) univ. | eq. (1) fitted | step obs | acc. rel. | eq. (1) univ. | eq. (1) fitted |
|---|---|---|---|---|---|---|---|---|
| 3.45, 1.43, 7.0e8 | 1.89 +- 0.48 | 2.34 (+0.9) | 3.60 (+3.6) | 4.31 (+5.1) | 5.96 +- 0.84 | 9.32 (+4.0) | 4.51 (-1.7) | 5.78 (-0.2) |
| 3.45, 1.43, 5.7e8 / 7.9e8 | same | 2.34 | 3.41 (+3.2) / 3.72 (+3.9) | 4.13 / 4.42 | same | 9.32 | 4.40 / 4.58 | 5.60 / 5.89 |
| 3.45, 1.00, 7.0e8 | 1.89 +- 0.48 | 2.49 (+1.3) | 3.62 (+3.6) | 4.40 (+5.3) | 5.96 +- 0.84 | 10.26 (+5.1) | 4.53 (-1.7) | 5.64 (-0.4) |
| 2.94, 0.98, 7.0e8 | 2.21 +- 0.56 | 2.71 (+0.9) | 4.25 (+3.7) | 4.07 (+3.3) | 7.00 +- 0.98 | 11.16 (+4.2) | 5.30 (-1.7) | 5.08 (-2.0) |
| (old, swapped) | 1.89 +- 0.48 | 2.34 (+0.9) | 4.36 (+5.2) | 4.90 (+6.3) | 5.96 +- 0.84 | 9.32 (+4.0) | 4.96 (-1.2) | 6.44 (+0.6) |

## 2. (b) DM-20: gas/star handling in feature_axisym_20261109.py -- the SWAP IS PRESENT
- Line 14: `r,vobs,lo,hi,vg,vd=d[:,0..5]` reads the mislabelled csv as labelled. So "vg" (reconstructed as GAS: h 0.2 kpc, soft mass constraint 1.10e9) is the STELLAR curve, and "vd" (reconstructed as the stellar DISC, scaled to M_* = 1.3e9) is the GAS curve. The fidelity gate's "gas peak 3.39 -> 3.51 kpc" checked the stellar curve's peak. M_b = 2.40e9 (sigma about 6% high), as in the spherical pre-registered run.
- WHAT SURVIVES [T3]: verdict (iii) GEOMETRY IRRELEVANT (FX-1). A1 (disc) and R1 (sphere) used the same in-plane V_bar (symmetric in the swap at M/L 1.43), the same M_b and the same sigma; the two scale heights are near equal (code: 0.20 gas, 0.196 R_d^0.633 = 0.23 kpc stars), so swapping which curve sits at which height barely moves the 3-D potential. The ratio A(A1)/A(R1) = 1.02 is like for like. The ~1/3 damping is still the law's.
- WHAT DOES NOT: (i) the A3 bracket rows (M_* 1.0e9 / 1.6e9) scale the GAS curve, not the stars: INVALID as labelled. (ii) FX-4 side results (disc fit 10^7.15; rms 0.061; chi2 1033 at 7.45) were computed at the swapped masses; they are superseded in the same way as the spherical 7.10 -> 7.17. Not recomputed. (iii) The gas-mass soft constraint was applied to the stellar curve (F2 flag: unconstrained mass beyond the last point).
- The tex A.3 still quotes "10^{7.15} in disc geometry" (l.702). The handoff 3.4 rewrite drops it. Keep it dropped unless recomputed.

## 3. (c) "0.013" vs 0.014
- Source of "0.013" (tex l.106 Summary limits row; l.702 A.3): RESULT_FEATURE_AXISYM_20261109 sec 5 wording, from A0 0.0128 / A1 0.0127 -- disc-geometry runs at the SWAPPED masses. The spherical pre-registered value was 0.0124 ("0.012" in RESULT_NGC1560).
- Corrected spherical, eq. (1) universal, peak-to-trough: 0.0137 (primary); 0.0135-0.0141 across the M_* bracket; 0.0144 at M/L 1.00; 0.0149 at (2.94, 0.98). Fitted: 0.0146-0.0172. Acceleration relation 0.0356 (M/L 1.43), 0.042 (1.00, 0.98). Observed 0.0344 in every row.
- Disc geometry at the corrected masses: NOT computed. Carried only as the ratio 1.02 (sec 2).
- So: "0.014" replaces "0.013" (two figures, rule (cm)). "In disc geometry alike" must be worded as the ratio, or the disc run redone (ballot OW-3).

## 4. Against the expectations (script header)
X1 RIGHT (H4 and H5 reproduce exactly from the old csv; the corrected 1.3e9 row reproduces them too). X2 RIGHT (obs 0.354 and acceleration relation 0.560 unchanged at M/L 1.43; shelf 3.60 / 4.31 / 2.34 = plot report). X3 RIGHT (0.210; 0.255). X4 RIGHT on the obs fall (0.30), the acceleration relation (0.55) and eq. (1) (0.19); jackknife RIGHT (0.031); perturbation sd 0.117-0.118 just BELOW the stated 0.12-0.17. X5 RIGHT (0.0137) and the source of "0.013" confirmed.

## 5. Reading [T3]
- The correction weakens the shelf result against eq. (1) but does not remove it: 3.2-3.9 sigma (universal density) at every normalisation tried, was 5.2. The acceleration relation draws the shelf (within 1.3 sigma) everywhere. On the quoted errors, which are conservative (side-to-side differences, about 4 x the point scatter, H2).
- The acceleration relation overshoots the rise after 6 kpc by 4.0-5.1 sigma at every normalisation; eq. (1) undershoots it by about 1.7.
- Transmission: the data sit between the two at every normalisation. At M/L 1.43 the data are nearer eq. (1) (1.0 sigma quoted) than the acceleration relation (1.5). At M/L 1.00 and at Gentile's MOND-fit values the data FALL to 0.30 and the gap widens: eq. (1) 0.9 sigma, acceleration relation 2.2 sigma (jackknife: 3.4-3.6 and 8.2-8.3).
- So the two post-hoc shape statistics lean opposite ways. The shelf favours the acceleration relation; the overall transmission and the step favour eq. (1). Neither is pre-registered. The honest sentence reports both and calls neither decisive. Same principle as the 17:24 ruling: what survives every normalisation is the shape, and the shape is split.
- Handoff 3.4 draft A.3 numbers CONFIRMED: rms 0.054 vs 0.073; chi2 722 vs 879; 515 vs 385; 92 vs 41; fitted 10^7.17, 120; amplitudes 0.014-0.015 vs 0.036-0.042, observed 0.034; shelf 1.9 +- 0.5 against 2.3 and 3.6. ADD (DM-11): transmission 0.35 (+-0.14 quoted, +-0.04 jackknife) against 0.21 and 0.56; at M/L 1.00, 0.30 against 0.19 and 0.55; the step overshoot.

## 6. Also found (tex, for the tex pass)
- l.702 A.3: "McGaugh (2026c) names it". 2026c is part 2; part 1 (2026b) is where NGC 1560 is named (CHECK_NP17 sec 8). Fixed by the 3.4 rewrite (DM-12, ruled). Whether part 2 also names it: not checked.
- l.702 A.3 still carries the pre-erratum numbers (0.067 vs 0.073; 1293; 10^{7.10}; 10^{7.15}; 0.013). Replaced by the 3.4 rewrite.

## 7. Ballot -- RULED JM 2026-11-12 01:25 "rec's accepted": OW-1 KEEP, OW-2 KEEP, OW-3 = (i), OW-4 KEEP
OW-1 Adopt secs 1 and 3 as the post-hoc NGC 1560 shape numbers for DM-11 and the A.3 rewrite, reported both ways (shelf against eq. (1); transmission and step against the acceleration relation). [KEEP]
OW-2 DM-20 outcome: FX-1 verdict (iii) stands as a like-for-like ratio; A3 bracket rows INVALID as labelled; FX-4 side results superseded (computed at swapped masses); pointer written into RESULT_FEATURE_AXISYM_20261109. [KEEP]
OW-3 Disc geometry at corrected masses: (i) word it as the ratio ("disc geometry changes the amplitude by 2 per cent"), no run; or (ii) rerun A1 only with the columns assigned correctly (solver run: rule 3a read first, expectations first, time estimate first). [Rec (i) for v2.0.2; (ii) optional]
OW-4 Tex l.106: "0.013" -> "0.014" and "spherical or disc" -> per OW-3. [KEEP; tex pass]
# END
