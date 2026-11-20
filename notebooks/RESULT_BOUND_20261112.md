# RESULT_BOUND_20261112 -- BOUND-1: is the rho_c search bound arbitrary? [UNRATIFIED]
> **JM 2026-11-12 23:43:** "0289 is very close to 10.5 if it is arbitrary"
Script bound_check_20261112.py (verbatim chi_fn of errprop_globalfit_20261112; expectations E1-E4 and the hypothesis written in its header before each run). Runs 23:44-23:49. Outputs bound_check_20261112_part.csv, bound_check_first_20261112.csv.
NOTE: several stamps in the script header were first written by guess (23:49-23:58) and corrected after reading the clock (23:47).

## 1. The bounds
- Fixed arm: bounds=(4.5,10.5) first appear in pass2_table2_full_20261022.py l.21. No rationale recorded anywhere. ARBITRARY. Inherited by axisym_r1, clean_conc, errprop, feature_transmission, h2_refit, sigmafree.
- Marginalised arm A (marg_ns_20261030.py): flat prior on log rho_c over [6,10]. Also arbitrary, and tighter.

## 2. Widening the fixed-arm search to 13 (all 147)
- 3 galaxies find a lower chi^2 above 10.5: NGC 5907 10.50 -> 10.84 (6.13 -> 5.61); NGC 2998 8.05 -> 11.35 (18.03 -> 16.26); UGC 11914 8.63 -> 12.21 (25.68 -> 16.93).
- NGC 0289 stays at 10.37 (interior minimum; chi^2 rises past it). E1 held. E3 (<= 3 movers) held at its limit, but two mid-range fits jump by 3.3-3.6 dex: NOT anticipated.

## 3. chi^2 oscillates in log rho_c (E4 wrong in kind: no flat floor)
- Above ~10, chi^2(log rho_c) is a damped oscillation. Minimum-to-minimum spacing: NGC 0289 10.35 -> 14.40 (4.05), NGC 5907 10.85 -> 14.90 (4.05), NGC 2998 11.35 -> 15.50 (4.15), UGC 11914 12.20 -> 16.35 (4.15). Max-to-max 4.05-4.10. Later dips shallower.
- Identification [T3; math T2]: the isothermal sphere's homology oscillation. Lane-Emden asymptote psi ~ ln(xi^2/2) + A xi^(-1/2) cos((sqrt7/2) ln xi + d): period 4 pi/sqrt7 = 4.750 in ln r = 2.063 dex in r. With r0 = sigma/sqrt(4 pi G rho_c), PREDICTED period 4.13 dex in log rho_c. Measured 4.05-4.15. The first (pre-asymptotic) gap is shorter (3.35 dex).
- Meaning: on the second branch the core is far smaller than the data (NGC 0289 at 10.35: r0 ~0.12 kpc, King radius ~0.36 kpc; innermost fitted point 7.6 kpc). chi^2 there responds to the oscillation phase, not to a resolved central density. At the first minimum (6.92) r0 ~6 kpc.

## 4. What it does to the headline numbers (all 147, fixed arm)
| Variant | median | robust spread | std | max | fits >= 10 |
| V1 grid [4.5,10.5] global (table of record after EP-1) | 7.415 | 0.497 | 0.628 | 10.50 | 2 |
| V2 grid [4.5,13] global | 7.415 | 0.497 | 0.805 | 12.20 | 4 |
| V3 first (lowest-rho_c) minimum | 7.376 | 0.495 | 0.528 | 9.57 | 0 |
- V3 changes only NGC 0289 (10.37 -> 6.92; chi^2 13.3 -> 21.7) and NGC 5907 (10.50 -> 7.44; 6.1 -> 14.4). NGC 2998 and UGC 11914 already sit on their first minimum in V1.
- Median and robust spread are stable (7.38-7.42; 0.50). std is not (0.53-0.81): do not quote std.

## 5. Consequences
- EP-2's "one fit (NGC 5907) sits at the upper bound of the search" is an artefact of an arbitrary bound. HOLD EP-2.
- EP-4's "std 0.68 -> corrected 0.63": the std depends on the variant. The note should quote the robust spread only.
- The marginalised arm A's [6,10] prior silently imposes the first-minimum choice on all four galaxies.
- ERRPROP perturbed fits used the [4.5,10.5] global grid. tau must be redone under whichever variant is chosen (12 affected galaxies x 5 fits, ~2 min).
- The bimodality (9/147, now 12 with >1 minimum or a lower minimum above 10.5) is explained by the homology oscillation, not by an optimizer quirk alone.

## 6. Ballot (UNRULED)
BD-1 Record both bounds as arbitrary (sec 1). [KEEP]
BD-2 Treatment of the second branch: (a) V2 widen everywhere + rerun arm A NS with a wider prior; (b) V3 define rho_c as the first-minimum (resolved-core) value, consistent with arm A's prior, the four second-branch galaxies flagged with their chi^2 in a footnote; (c) keep [4.5,10.5] -- not defensible. [Claude recommends (b), with the V2 numbers shown alongside]
BD-3 After BD-2: redo ERRPROP tau for the 12 affected galaxies; then redraft EP-2. [YES]
BD-4 Homology-oscillation identification (sec 3) as paper material: framing ledger (cp). [JM call]

- 2026-11-13 00:27: precision fix to sec 3: first gap (first minimum -> second) is 3.30-3.57 dex from the refined first-minimum fits (0289 3.43, 5907 3.41, 2998 3.30, UGC 11914 3.57), not 3.35. Gap 2->3 is 4.05-4.15 vs predicted 4.125 (residuals -0.08 to +0.02). Core scale for NGC 0289: r0 6.14 kpc (first min), 3.34 (universal), 0.12 (second min). Spreadsheet: https://claude.ai/artifact/LUSNunAJAPuycG6nCckKK3

BD-5 (proposed 2026-11-13 00:56) [UNRATIFIED]: adopt the name Joseph-Lundgren oscillation (JLO) for the BOUND-1 effect; framing ledger (cq). Claude recommends YES.
- 2026-11-13 19:42: BD-2 RULED (b) by JM ("for 1. b"): V3 first-minimum definition adopted; V2 numbers shown alongside; four second-branch galaxies footnoted with chi^2. Record: BALLOT RULING BD-2.
