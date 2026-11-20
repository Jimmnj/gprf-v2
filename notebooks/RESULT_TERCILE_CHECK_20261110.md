# RESULT_TERCILE_CHECK_20261110 -- Sigma_h terciles of the fitted rho_c; sample audit of the B.1 slopes [RATIFIED JM 2026-11-10 17:43 "accepted": secs 1-2 and TC-1..TC-4; APPLIED 17:55; consistency check on existing fits, not a pre-registered test]
Written 2026-11-10 17:50 EST. Script tercile_check_20261110.py -> tercile_check_20261110_report.txt (numpy only; ~5 s). Input rhoc_corr_20261105_table.csv (fixed arm lrc = pass2_table2_part, N 147; marginalised arm lrc_marg = marg_fit_params arm A, MAP; M_b, Sigma_h, f_gas at CATALOGUE values in both). Asserted Rd(NGC2841) = 3.64, Rd(DDO154) = 0.37.
> **JM 2026-11-10 17:30:** "run the tercile check"
Question (ledger (ck), BALLOT 17:15 OWED): can B.1's surface-density slope and its "class medians within 0.15 dex of 7.45" both hold in the fixed arm?

## 1. Numbers (147 galaxies; bootstrap 95%)
| | FIXED arm | MARG arm (MAP) |
|---|---|---|
| median / sd of log rho_c | 7.416 / 0.677 | 7.412 / 0.596 |
| OLS (M_b, Sigma_h): Sigma_h slope | +0.464 [+0.22, +0.71]; resid sd 0.648 | +0.317 [+0.05, +0.58]; resid 0.572 |
| OLS (M_b, Sigma_h, f_gas): M_b / Sigma_h / f_gas | -0.381 / +0.374 [+0.11, +0.61] / -0.808; resid 0.634 | -0.300 / +0.275 [+0.03, +0.51] / -0.377; resid 0.568 |
| variance removed by Sigma_h beyond M_b | 7.3% | 4.5% |
| raw Sigma_h thirds (median log Sigma_h 1.40 / 1.95 / 2.43): median log rho_c | 7.29 / 7.51 / 7.59; high-low +0.30 [+0.03, +0.57] | 7.33 / 7.53 / 7.41; +0.08 [-0.17, +0.36] |
| thirds at fixed M_b | 7.11 / 7.31 / 7.65; +0.54 [+0.22, +0.76] | 7.22 / 7.35 / 7.64; +0.42 [+0.01, +0.63] |
| thirds at fixed M_b and f_gas | 7.17 / 7.36 / 7.50; +0.33 [+0.05, +0.54] | 7.17 / 7.42 / 7.51; +0.34 [+0.03, +0.61] |
Sigma_h 16-84 range 1.07 dex; at fixed M_b 0.78 dex.
Validation: the fixed-arm three-variable fit reproduces the pass7 full-sample refit (10.934, -0.381, +0.374, -0.808; pass9_wb_rhoc_sweep_20261025.py l.14) to three decimals.

## 2. Findings [T3]
(a) ANSWER: the two B.1 statements coexist only because surface density is not among the listed classes. In the fixed arm the Sigma_h thirds at fixed mass run 7.11 / 7.31 / 7.65 (0.54 dex, interval excludes zero). "Class medians within 0.15 dex of 7.45" is true for bulge, gas, r_peak/R_d and mass classes (CLEAN_CONC) and FALSE for surface density. B.1 must name the exception.
(b) SAMPLE MISMATCH IN THE NOTE (pre-existing, found here): B.1's "+0.37 [+0.14, +0.62] with mass as the only other variable" and rule (3)'s +0.268 are the N=120 TEST1 fits of 2026-10-18/19 (RESULT_RHOC_SIGMA_B_OOS_20261019 step 1; RESULT_DWARF_PREREG_20261020 l.2). On the 147-galaxy pass2 fits used everywhere else in v2: +0.46 (mass only) and +0.37 (mass and gas fraction). The note gives no N beside these slopes.
(c) "UNCHANGED" IS WRONG LIKE-FOR-LIKE: RESULT_MARG caveat 3 and B.1 say the Sigma_h slope survives marginalisation "at the same value (+0.255 vs +0.268)". That compares the 147-galaxy marginalised fit with the 120-galaxy fixed fit. On the same 147 galaxies and the same controls: +0.37 fixed -> +0.26 to +0.28 marginalised (+0.255 with fitted M_b, Sigma_h, RESULT_MARG M3; +0.275 with catalogue values, here). The slope falls by about a quarter and stays positive (lower bound +0.03).
(d) CLAUDE ERROR 2026-11-10 17:15-17:20: the "like-for-like" correction to AN-1/AN-2 (+0.27 fixed, +0.26 marginalised) inherited mismatch (c); it is now IN THE TEX (Summary rho_c row; sec 7 is unaffected, it quotes no slope). Correct pair on 147: +0.37 fixed, +0.26 marginalised.
(e) What stands: "no baryon property predicts rho_c" -- Sigma_h removes 7% (fixed) / 5% (marginalised) of the variance beyond mass; spread 0.68 -> 0.65 and 0.60 -> 0.57. In the marginalised arm the RAW Sigma_h thirds are flat (+0.08); the trend appears only at fixed mass.
(f) Rule (3) itself is NOT to be changed: it is the published, pre-registered rule (N=120 fit) used in Appendix B, sec 3.1, N16, SIGMAD, CDR. Only its description needs the N.

## 3. Text (APPLIED 2026-11-10 17:55; BALLOT same stamp)
TC-1 Summary rho_c row: "+0.27 fixed, +0.26 marginalised" -> "+0.37 fixed, +0.26 marginalised, 147 galaxies".
TC-2 B.1 opening: say rule (3) and the +0.37 [+0.14, +0.62] mass-only slope are the 120-galaxy fits with a flat part; add "on the full 147: +0.46 with mass only, +0.37 with gas fraction included". Replace "the surface-density slope is unchanged (+0.26, [+0.03, +0.50])" with "the surface-density slope falls from +0.37 to +0.26 ([+0.03, +0.50]) on the same 147 galaxies".
TC-3 B.1 class-median sentence, add: "Surface density is the exception: at fixed mass the thirds in $\Sigma_h$ have medians 7.11, 7.31 and 7.65."
TC-4 RESULT_MARG caveat 3: annotate "(+0.255 vs +0.268)" as a 147-vs-120 comparison (record fix only).
## 4. Files: tercile_check_20261110.py, tercile_check_20261110_report.txt, this file. UNCOMMITTED.
> **JM 2026-11-10 17:43:** "accepted"
# END
