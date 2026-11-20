# RESULT_ERRPROP_20261112 -- ERRPROP-1: catalogue errors through the fixed-M/L fit of rho_c, and an optimizer failure found on the way [UNRATIFIED]
Written 2026-11-12 22:12 EST. Prereg PREREG_ERRPROP_20261112.md (secs 1-6 written 21:51 before any perturbed fit; sec 7 addenda 21:55 and 22:01, each before its run). Rule 3a (new text, stamped 21:49) followed: canon lines read 21:50; solver = pass2 solve_gal verbatim; validation 147/147 exact (max |diff| 0.000).
Scripts / outputs: errprop_20261112.py (+part.csv), errprop_analysis_20261112.py (+report), errprop_grid_20261112.py (+report), errprop_globalfit_20261112.py (+part.csv), errprop_globalcmp_20261112.py (+report), errprop2_20261112.py (+part.csv), errprop2_analysis_20261112.py (+report), errprop_table2_fix_20261112.py (+report). All in the device VM, 21:51-22:11.

## 1. OPTIMIZER FAILURE in the fixed-arm table (pass2_table2_part.csv)
The bounded Brent fit (minimize_scalar, [4.5, 10.5]) is not a global fitter here: chi^2(log rho_c) is bimodal for 9 of 147 galaxies (second minimum near 10.3-10.4; cf. sec 2.2 of the note, several solutions beyond a central contrast of ~14). A 0.05-dex grid + refine finds that pass2 MISSED the global minimum in 3 galaxies:
| galaxy | pass2 | global | chi2 pass2 -> global |
|---|---|---|---|
| NGC 4013 | 10.49 (bound) | 7.27 | 39.3 -> 9.7 |
| UGC 1281 | 10.49 (bound) | 7.42 | 0.5 -> 0.4 (flat: unconstrained) |
| NGC 0289 | 6.92 | 10.37 | 21.7 -> 13.3 |
Fixed arm, all 147, global fits: median 7.42 (unchanged), std 0.63 (was 0.68), robust half-width 0.50 (was 0.52). ONE fit at the bound (NGC 5907, i 88 deg, global minimum at 10.50), not three. Table 2 "law, rho_c free" rms 0.0813 -> 0.0806 (0.081 at the quoted precision; RAR function, M/L free, 0.079): unchanged as quoted.

## 2. ERRPROP (global fitter; one-at-a-time +-1 sigma; 143 unflagged galaxies)
Median per-galaxy error on log rho_c: 0.42 dex (16-84: 0.18-0.63) -- distance 0.32, M/L 0.10, inclination 0.08. Robust spread of fitted log rho_c 0.48 (std 0.57). Maximum-likelihood intrinsic scatter tau (errors in quadrature): 0.48 dex. Direct-distance galaxies (SPARC f_D != 1, N 65): median error 0.27, robust spread 0.45, tau 0.53. Hubble-flow (N 78): error 0.51, spread 0.49, tau 0.37. Error size does not track deviation from the median (Spearman -0.08, p 0.33). Flagged at a bound in some perturbation: 4 (NGC 5371, NGC 5907, NGC 6674, UGC 11557).

## 3. Against the expectations (prereg sec 5 / 7)
X1 RIGHT (validation exact). X2 median error 0.42: just ABOVE the stated 0.2-0.4; distance-dominated RIGHT. X3 tau 0.48: just ABOVE 0.3-0.45. X4 direct-distance error 0.27 (stated 0.1-0.2: WRONG, higher); tau within 0.1 of Hubble-flow: WRONG (0.53 vs 0.37 -- the better-measured galaxies scatter MORE). X5 (revised 1-4 at a bound): RIGHT (4 flagged; 1 at the bound at catalogue values). Grid G1 (NGC 4013 artefact) RIGHT; G2 (NGC 5907 real) RIGHT; G3 (UGC 1281 flat) RIGHT. "2-8 missed minima" RIGHT (3). Not foreseen: NGC 0289 missed in the opposite direction.

## 4. Reading [T3] (prereg sec 6 rule: tau >= 0.3 -> "a real spread of about tau dex remains")
- A real spread of about 0.5 dex in rho_c remains after the catalogue errors are propagated. The catalogue errors are as large as the spread (median 0.4 dex, through the D^-5 lever), but they do not explain it: galaxies with direct distances, whose errors are smaller, scatter as much or more, and the size of a galaxy's error does not predict its deviation. B.1's "it is not measurement error" stands, now with numbers.
- Caveats: linear one-at-a-time propagation (no nuisance covariance); Gaussian likelihood; the Hubble-flow tau (0.37) may be low if SPARC's Hubble-flow distance errors are generous; fixed arm only (the marginalised arm samples, so it is not affected by the optimizer issue).

## 5. RETRACTIONS (by this result)
- NP-20's premise "three fits sit at the fitter's upper bound" (Claude, 2026-11-12 01:34; ruled 18:37): WRONG for NGC 4013 and UGC 1281 -- optimizer failures. Only NGC 5907 is at the bound.
- The edge-on reasoning built on the three (18:30, 18:39-18:41: "5.9 robust sigma", "lower limits", the edge-on sentence): void for NGC 4013 and UGC 1281. The marginalised-arm evidence comparison of the three (BALLOT 18:41) is UNAFFECTED (nested sampling).
- "std 0.68": a fixed-arm number from a table with three wrong fits; corrected 0.63.

## 6. Ballot (UNRULED)
EP-1 Adopt the global-fit table as the fixed arm of record (errprop_globalfit_20261112_part.csv); pass2_table2_part.csv kept, with a pointer. [KEEP]
EP-2 NP-20 sentence, replacing the 18:37 ruling's premise: "The fitted densities have median 10^{7.42} and a robust spread of 0.5 dex; one fit (NGC 5907) sits at the upper bound of the search. Propagating the catalogue errors on distance, inclination and mass-to-light ratio gives a median error of 0.4 dex per galaxy, mostly from distance; the spread that remains is about 0.5 dex, and galaxies with direct distances spread as much as the rest." [KEEP; JM wording]
EP-3 FIT-1 audit, owed before the tex pass: every number in the note that comes from a bounded scalar fit of rho_c (fixed-arm correlations in B.1 and sec 7 incl. "0.65 of 0.68", RESULT_CLEAN_CONC FIX0, rho_c rules, NGC 1560 R2 fits 7.17 / 7.23 / 7.49, NGC 6946 7.29, the Milky Way 7.88) re-checked with the grid + refine fitter. [YES]
EP-4 Retractions sec 5 [KEEP].

## 7. EP-1 RULED (JM 2026-11-12 22:26 "accepted"): fixed-arm table corrected IN PLACE 22:26 (pass2_globalfix_20261112.py; 3 rows; backup kept). Post-fix: median 7.42, std 0.63, robust 0.50, Table 2 law rms 0.0806.

## 8. NGC 0289 high branch: post hoc diagnosis (22:37) [UNRATIFIED]
Script ngc0289_diag_20261112.py + report (verbatim solver via errprop_globalfit chi_fn). Asked by JM; diagnosis, not a test.
- Expectation written first: high branch bought by inner points. WRONG: the six points at r < 7 kpc are masked (ev/v >= 0.1).
- The high branch is bought by the curve's shape at r > 7.5 kpc: peak ~185 at 10-14 kpc, dip to ~156 at 20-26, rise to 182 at 40-44, fall to ~158 at 62-71. The universal-rho_c model runs flat at 183-197. Raising rho_c to 10.37 packs the nb mass inward and the model falls from 200 (10.6 kpc) to 169 (71 kpc). Chi2 13.3 vs 21.7 (6.92) vs 37.0 (7.45).
- FRAGILE: any 1-sigma nudge in D, i or Y_disk moves the fit to 6.4-7.8 (errprop2). The propagated error is several dex, so the tau MLE already gives it almost no weight. No exclusion needed.
- The driving features are the ones the observers flag: Walsh et al. 1997 model the 10-30 kpc dip with a warp and spiral density waves; small bar with oval distortion (Pence & Blackman 1984).
- Outer level: V_f 163 vs V_BTFR 184 (0.88) at D 20.8 Mpc. Gas-rich (M_gas ~3.7e10 of M_b 7.26e10), HI to 71 kpc. Every model in the marginalised arm pulls D shorter (U zD -1.24, MOND -1.7) and all tie once nuisances are free (lnZ -83.8 to -86.9).
- Distance is NOT a fix by assertion: independent distances do not point shorter (Walsh 21.7, CF3 ~23.3). The 13.7 Mpc TF distance would be circular here.
- Light profile: SBeff 1588 vs SBdisk 70 and V_disk peaking at 2.5-5 kpc (not 2.2 R_d ~ 15 kpc): a bright barred inner disc plus a faint outer disc forced into one R_d. Does not drive the fit (inner points masked).
- Reading [T3]: a kinematics-plus-distance case, not evidence of a high rho_c. Candidate referee answer for the high tail (with NGC 5907).
- Checked: NGC 0289 sigma_tot = 5.3 dex (sD 0.71, sI 3.78, sU 3.63), so its weight in the tau MLE is negligible. Side issue: sD = |lrcDp - lrcDm|/2 understates the swing when both nudges leave a branch (here 0.71 vs a 4.0 dex jump). Affects at most the 9 bimodal galaxies; tau unlikely to move, not rerun.
- CORRECTION 23:36 [UNRATIFIED]: (1) Walsh et al. 1997 ADOPT 21.0 Mpc (Table 1; redshift with Local Group + Virgocentric corrections), not 21.7 (21.7 is their value before the Virgocentric step). Redshift distances in their Table 1 span 19.4-22.4. (2) Tully-Fisher: Bottinelli et al. 1986 give 16.7 Mpc; Walsh's own corrected TF gives 13.7 (not adopted). Both TF values point shorter; both are circular for a BTFR test, and Walsh et al. note large LSB discs are overluminous for their rotation. D at which V_f = V_BTFR is 16.2 Mpc. (3) Marginalised arm is NOT a strict tie: ln Z spans 3.1 (NFW -83.8 best; NFW leads U by 2.3, MOND by 3.1). 'All models tie' in the 22:37 answer and in sec 8 was loose. Spreadsheet: https://claude.ai/artifact/4yqEjWpfaM2cX2PSsUM4wA

## 9. BD-3: ERRPROP under the V3 first-minimum definition (2026-11-13 19:50-20:10) [UNRATIFIED]
Trigger: BD-2 RULED (b), JM 2026-11-13 19:42. Scripts: errprop3_first_20261113.py (fitter: scan 4.5-13.0 at 0.1, first interior local minimum, refine +-0.1; chi_fn verbatim, rule 3a), errprop3_analysis_20261113.py (+ _report.txt), pass2_v3_preview_20261113.py (+ _report.txt; read-only on the table).
- Scope: the 12 multi-minimum galaxies (bound_check_first_20261112.csv) PLUS 5 found by the X4 screen (any errprop2 perturbed fit > 1 dex from baseline): DDO 064, NGC 3521, NGC 5371, NGC 6674, NGC 6946. A branch jump is >= ~3.3 dex, so the 1-dex screen catches all. Three of the five had second-branch PERTURBED fits: DDO 064 D- 9.58 -> 8.31; NGC 5371 D- 10.50 -> 8.40; NGC 6674 D+ 10.50 -> 7.63. NGC 3521, NGC 6946 unchanged.
- Expectations: X1 HOLDS (baseline vs lrc_first, max 0.001 dex). X2 HOLDS (NGC 0289 sigma_tot 5.29 -> 0.80; NGC 5907 3.57 -> 0.59). X3 HOLDS (tau 0.479 -> 0.475). X4 WRONG (5 galaxies outside the 12 needed checking; 3 changed).
- tau (MLE, errors in quadrature): all 0.475 (N 147, was 0.479 on N 143); direct distance 0.523 (N 66); Hubble flow 0.365 (N 81). Median sigma_tot 0.42 (D 0.32, I 0.08, M/L 0.10). Spearman sigma_tot vs |deviation| -0.07 (p 0.41). No galaxy is flagged now (was 4: NGC 5371, NGC 5907, NGC 6674, UGC 11557).
- UGC 11557: its inclination-perturbed fit has NO low-branch minimum (the low minimum disappears), so V3 returns 10.41; sigma_tot 2.10 -> 3.98. It enters tau with negligible weight. A genuine first minimum, not a branch choice.
- Fixed-arm table under V3 (preview, table NOT changed): median 7.38 (was 7.42), std 0.53 (was 0.63), robust 0.50 (unchanged), none at the bound (was 1), law rms 0.0810 (was 0.0806). NGC 0289 10.37 -> 6.92 (chi^2 13.3 -> 21.7; galaxy rms 0.063 -> 0.086); NGC 5907 10.49 -> 7.44 (6.2 -> 14.4; 0.049 -> 0.080).
- Reading [T3]: BD-2 (b) moves nothing that carries the argument. The intrinsic spread ~0.5 dex stands; the std falls to the robust value once the second branch is excluded (the 0.63-0.81 instability was the branch, not the data). The cost is a worse fit for two galaxies, to be shown in the footnote.
- Toolchain: first run wrote 12/13-column rows (lrcUb only for bulge galaxies); fixed, bad file in _to_delete/. argv overwrite by the exec'd chain recurred; captured before exec.
- 2026-11-13 22:06: BD-6 RULED YES (JM 22:04) and APPLIED: table of record now V3 (pass2_v3fix_20261113.py; md5 82dce260; backup kept). Fixed arm: median 7.38, std 0.53, robust 0.50, law rms 0.0810.
- 2026-11-13 22:43: EP-2 RULED (JM 'yes'): 20:10 redraft adopted as worded. See BALLOT.
- 2026-11-13 22:59: EP-4 RULED KEEP (JM "yes") with sec 5 UPDATED after BD-2/BD-6 (sec 5 text above kept as the 11/12 record):
  (1) NP-20 premise "three fits sit at the fitter's upper bound" WITHDRAWN: NGC 4013, UGC 1281 optimizer failures; NGC 5907 a JLO minimum; after BD-2/BD-6 no fit sits at a bound.
  (2) Edge-on reasoning ("5.9 robust sigma", "lower limits", the edge-on sentence ruled 18:37) WITHDRAWN; EP-2 sentence replaces it. Marginalised-arm evidence comparison (BALLOT 18:41) unaffected (nested sampling).
  (3) "std 0.68" WITHDRAWN: three wrong fits + two JLO fits. Corrected std 0.53 (interim 0.63 after EP-1); robust 0.50. The note quotes the robust spread (EP-2).
