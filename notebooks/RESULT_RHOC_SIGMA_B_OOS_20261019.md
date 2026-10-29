# RESULT_RHOC_SIGMA_B_OOS_20261019 [UNRATIFIED] -- does the SPARC rho_c(Sigma_b) trend predict Crater II and NGC1052-DF2/DF4? Pre-registered.
Script rhoc_sigma_b_oos_20261019.py; report rhoc_sigma_b_oos_report_20261019.txt. Frozen law (TEST1 tanh2 rho_c, N=120); dwarfs via the N9 solver
(Plummer stars, BTFR sigma, host ambient; uncapped exponent -- cap irrelevant, x_b < 0.5 for all four). Dwarf inputs [mem, VERIFY] (script header).
Declared PASS criterion: predicted log rho_c inside the required band for BOTH Crater II and DF2; And XIX must stay a pass.

## Step 1 (SPARC): log rho_c = 8.35 - 0.15 log M_b + 0.37 log Sigma_h   (c: [+0.14, +0.62]; b: [-0.32, +0.01]); resid 0.65 dex; median 7.45.
   Sigma_h = M_b / 2 pi R_h^2, R_h = 1.678 R_d. SPARC log Sigma_h spans 0.6..2.9 (median 2.0).

## Steps 2-3 (required vs predicted log rho_c; prediction error = SPARC resid (+) extrapolation)
| object    | log Sigma_h | REQUIRED        | PREDICTED     | at median 7.45 -> sigma | verdict |
| And XIX   |  -1.58      | 6.85 [6.25, >8] | 6.81 +/- 0.74 | 8.3 (obs 7.8)           | PASS (control) |
| Crater II |  -1.44      | 6.13 [6.00,6.27]| 6.98 +/- 0.74 | 5.2 (obs 2.7)           | HALF: trend supplies 0.47 of the 1.32 dex; predicted sigma ~4.6, 1.1 sigma(pred) above band |
| DF2       |   0.82      | 5.77 [<5, 6.35] | 7.38 +/- 0.66 | 20.8 (obs 8.5)          | FAIL: 1.6 sigma(pred); Sigma_h is INSIDE the SPARC range -> not an extrapolation; DF2 needs ~1.6 dex less rho_c than SPARC galaxies of its M_b and Sigma_h |
| DF4       |   0.97      | <5.96 (1 sigma) | 7.45 +/- 0.66 | 17.4 (obs 4.2)          | FAIL: observed sigma below the stars-only floor (7.8); only the +4.4 error admits any nb |

## Verdict: FAIL by the declared criterion. The surface-density trend is real in SPARC (c = +0.37, 95% excludes 0) and moves Crater II a third of
the way and And XIX exactly onto its band, but it does not reach Crater II and does nothing for DF2/DF4. Consequence (framing (as)): 'mostly' = FOUR
concessions; DF2/DF4 named beside Crater II; MOND's EFE account (Famaey, McGaugh & Milgrom 2018) cited as the competing explanation.
## Sharper statement for the open list: DF2 sits at a SPARC-typical (M_b, Sigma_h) and needs rho_c ~1.6 dex below SPARC galaxies there. Under the
frozen law that is not a trend; it is an outlier. The law has no lever (no EFE, no stripping, R7 makes stars always source). Candidate levers,
all outside the frozen law: (i) distance (13 Mpc would lower M_b 3x and r_h; recompute if adopted), (ii) stellar non-equilibrium (tides mild, r_J ~5 kpc),
(iii) a rho_c dependence on something DF2 lacks (gas? rotation?) -- SPARC low-Sigma_h galaxies are gas-rich rotating discs; DF2/DF4 are gas-free spheroids.
(iii) is testable: split SPARC rho_c residual by gas fraction. NOT RUN.
## Caveats: DF2/DF4 parameters from memory; sigma_los within r_h, isotropic; host ambient from projected 80 kpc; rho_c definitions (TEST1 central
density vs N9 rho_amb e^{phi0/s2}) agree to leading order.

## ADDENDUM 19:25 -- gas-fraction split (rhoc_fgas_split_20261019.py). Joint fit f_gas coefficient -0.93 [-1.71, -0.14] dex/unit: gas-rich galaxies have LOWER
rho_c at fixed (M_b, Sigma_h); binned residuals non-monotonic. Wrong direction for gas-free DF2/DF4 (+0.37 dex from median). Lever (iii) DEAD. DF2/DF4 remain outliers.

## ADDENDUM 2026-10-20 20:05 -- VERIFIED DF2/DF4 inputs; rerun (rhoc_sigma_b_oos_df_rerun_20261020.txt). DF2: stellar sigma 8.4 +/- 2.1 (Danieli+2019, KCWI), r_h 2.7 kpc,
M* ~2e8 (van Dokkum+2018; GC sigma < 10.5 at 90%). DF4: STELLAR sigma 8.0 +2.3/-1.9 (Shen+2023, KCWI; M_dyn(<r_h) 8 +6/-4 e7; stars-only expectation 7 +/- 1) --
the earlier 4.2 +4.4/-2.2 was the 7-GC value (van Dokkum+2019). CORRECTION: 'DF4 below the stars-only floor' WITHDRAWN; DF4 is a stars-only object like DF2.
Rerun: DF2 required log rho_c 5.82 [<5, 6.22] vs predicted 7.31 -> still a ~1.5 dex outlier; DF4 required 5.44 [<5, 6.52] vs predicted 7.45 -> ~1 dex outlier at 1 sigma.
Verdict unchanged (FAIL; outliers), wording changed: both are 'consistent with stars alone; the law predicts 17-22 km/s at the sample density'.
MOND (verified): isolated ~20; with NGC 1052's field 13.4 +4.8/-3.7, 9.7-18.2 over x2 in M/L (Famaey, McGaugh & Milgrom 2018); Haghi+2019 agree if D ~20 Mpc and
separation < 150 kpc. Observed 8.4 +/- 2.1 is ~1.2 sigma below the EFE central value -- 'accommodated' remains the fair word for MOND.
