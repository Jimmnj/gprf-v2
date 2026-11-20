# PREREG_BIGSPARC_20261108 -- predictions of the saturated isothermal relation for BIG-SPARC, fixed before the mass models exist
Written 2026-11-08 17:05 EST. [UNRATIFIED] until JM stamps. Trigger: JM 16:26 "are we ready or able to make predictions about big sparc - can we target it?"; 16:55 "go". No BIG-SPARC data have been touched. All numbers below come from the 147-galaxy SPARC analysis of v2.0.1 (Zenodo 10.5281/zenodo.23815559) and the ratified RESULT files cited.
> **JM 2026-11-08 16:26:** "are we ready or able to make predictions about big sparc - can we target it?"

## 0. Target
BIG-SPARC (Haubner, Lelli, Di Teodoro, Duey, McGaugh, Schombert, Hess et al.; IAU S392 proceedings, online 2026-07-14): ~4000 galaxies, HI datacubes from public archives (APERTIF, ASKAP, ATCA, GMRT, MeerKAT, VLA, WSRT), WISE near-IR photometry, homogeneous rotation curves, surface brightness profiles and mass models. First chunk: Haubner et al. 2026, A&A 714, A250 (568 Apertif galaxies; rotation curves, HI surface densities, kinematic inclinations, CF4 distances; NO stellar mass models; 3-23 rings, mode 4; 30" beam at median 60 Mpc ~ 8.7 kpc per beam). This prereg scores against the BIG-SPARC MASS-MODEL release (V_gas, V_disk, V_bul at a stated M/L), not against the Apertif-568 alone. Nothing is fitted to the Apertif-568 before that release.

## 1. The model, as frozen (v2.0.1, canon R1)
rho_nb = rho_c exp[S(x_b) - S(x_b0) - dphi_nb/sigma^2] - rho_amb, x_b = phi_b/sigma^2, S saturating at x_b = 2 (g_b = a0); sigma^2 = sqrt(G a0 M_b)/2 (BTFR input, a0 = 1.2e-10); rho_c the one fitted constant, universal value 10^7.45 Msun kpc^-3 (median of the SPARC fit). Nothing in this file changes the model. Solver: the v2.0.1 archive (10.5281/zenodo.23995058), rule 3a (SNAP_20260922_001 sec 0c) before any new code.

## 2. Nuisance protocol (fixed)
M/L: BIG-SPARC's own stated values (WISE W1), applied identically to the relation and to the comparison function. Distances and inclinations: catalogue values with catalogue errors; the primary scoring is at FIXED catalogue nuisances (as v2.0.1 sec 3.1); a marginalised arm (D, i, M/L free under the catalogue priors, as PREREG_MARG_20261030) is supplementary. Comparison function: the empirical RAR function of McGaugh et al. (2016) with g_dagger = 1.2e-10, same nuisances. Points: g_bar > 0, V_obs error < 10% of V_obs (SPARC cut), galaxies with >= 5 points for the RAR items.

## 2a. Compute plan (fixed; JM 18:02 "accepted")
Fixed-nuisance arm on the FULL mass-model sample (~5 s/galaxy, ~6 CPU-h). Marginalised arm (dynesty, as PREREG_NS_20261030) on a 500-galaxy random subsample drawn with seed 20261108 from the galaxies meeting the P1 point cut, ONE sigma_int = 0.03, 16-32 workers on the 5950X (~110 CPU-h, one overnight). Solver = v2.0.1 archive unchanged; a cached-phi_b variant may be used ONLY after reproducing the canon test (pb_solve_canon vflat 0 -> Q1 clean +0.158 [+0.065, +0.275]) and the 10/30 marginalised medians on 10 galaxies to 0.02 dex. No GPU (FP64, sequential sampler; see SNAP_20261108_002 sec 2).

## 3. Predictions and criteria (numbers fixed now)
### Resolution-robust (score on the full mass-model sample)
P1 MEAN ACCELERATION RELATION at the universal density. Binned mean log g_obs vs log g_bar, 0.2-dex bins, g_bar in [1e-12, 1e-9] m s^-2: eq. (1) at rho_c = 10^7.45, no per-galaxy freedom, lies within 0.10 dex of the data mean in every bin with >= 30 points. Fallback (BTFR): the sample's BTFR at the flat part has slope 4 +- 0.3 and normalisation G a0 M_b within 0.1 dex at their M/L (eq. 2 is the BTFR by construction; a slope outside [3.7, 4.3] with their M/L is a FAIL of the input, to be stated as such). PASS: all bins within 0.10. FAIL: any bin with >= 30 points off by > 0.15. Between: report.
P2 ENVIRONMENT (the sec 6.2 test). Per galaxy: mean residual log(g_obs/g_RAR) over points with g_bar < a0/3, against the Newtonian environmental field log e_N (Chae et al. 2021 method, maximum-clustering estimate, recomputed for the BIG-SPARC footprint from a public redshift catalogue; the recipe and catalogue to be named BEFORE the run), at fixed log M_b and log Sigma_h (OLS with controls; bootstrap CI). Relation predicts slope 0. External-field effect predicts -0.08 dex per dex (v2.0.1 sec 6.2, SPARC-like sample; to be re-derived for the BIG-SPARC mass range before the run and the number frozen then). KILL of the relation's environmental prediction: slope <= -0.08 with 95% CI excluding 0 and excluding -0.02 (i.e. -0.08 +- 0.03 or tighter). PASS: CI includes 0 and excludes -0.08. Between: inconclusive. Forecast precision: SPARC gives +-0.16 on 60 galaxies; at 1000 galaxies with field values ~ +-0.04, at 3000 ~ +-0.023 (sqrt-N scaling, same per-galaxy scatter) -> the test is decisive at 2-3 sigma on the full sample.
P3 THE CENTRAL-DENSITY DISTRIBUTION. Per-galaxy fitted log rho_c (fixed nuisances; chi^2 minimum, flat prior [6,10]), on galaxies whose posterior/profile width is < 0.5 dex: median 7.45 +- 0.15; 16-84 width 1.0 +- 0.25 dex (SPARC marginalised: 1.01; fixed-nuisance: 0.68 -> state which arm is scored: FIXED, so the prediction is width 0.7 +- 0.2 for the fixed arm, 1.0 +- 0.25 for the marginalised arm); Spearman with log SB_disk |rho| < 0.2; with f_gas after the log M_b trend |rho| < 0.2; with log M_b between -0.35 and 0 (SPARC -0.27). FAIL (B.1 was a small-sample artefact): width < 0.4 dex on >= 500 constrained galaxies. FAIL (a rule exists): any of SB, f_gas-after-mass with |rho| >= 0.4, CI excluding 0.2. [RESULT_RHOC_CONSTRAINED_20261108 secs 1, 4; RESULT_RHOC_CORR_20261105]
### Resolution-dependent (score ONLY on the subsample meeting the floor; the floor is part of the prediction)
Floor F: >= 7 rings, innermost ring < 1.5 R_d (WISE disc scale length), beam < 2 kpc at the galaxy's distance. Galaxies not meeting F are UNDECIDED for P4-P6, not failures.
P4 SHAPE OFFSETS of rho_c at fixed mass: galaxies with a WISE bulge component have fitted log rho_c lower than pure discs by 0.16 +- 0.12 dex (Mann-Whitney p < 0.05 on >= 30 bulge galaxies); galaxies without a flat part (no V_f by the Lelli 2016 flag or its BIG-SPARC equivalent) lower by 0.3 +- 0.15. Sign and magnitude both predicted; opposite sign at > 2 sigma = FAIL. [RESULT_RHOC_CONSTRAINED secs 4, 6 -- UNRATIFIED at writing; this item is scored only if those are ratified]
P5 INTRINSIC SCATTER at low acceleration: eq. (1) at the universal density predicts 0.07-0.11 dex of intrinsic scatter about the mean relation for g_bar < a0/3 (profile variety); MOND predicts 0. Scored only where the published error model allows the intrinsic term to be separated to +- 0.03 dex (their stated per-point errors + distance + inclination + 0.1 dex M/L propagated as in sec 3.1). PASS: 0.04-0.14. FAIL: < 0.03 with CI excluding 0.07. Otherwise undecided.
P6 FEATURE CORRELATION (Renzo's rule): detrended curve features (log V minus a per-galaxy quadratic in log r, as v2.0.1 sec 7) of eq. (1)'s curves correlate with the observed at r = 0.33 +- 0.08 over the pooled points, against 0.24 +- 0.08 for the RAR function and 0.27 for the baryons alone. PASS: relation > RAR function by >= 0.05 with bootstrap CI excluding 0. FAIL: relation < RAR function by >= 0.05, CI excluding 0.

## 4. Not predicted / excluded from scoring
Individual-curve chi^2 at the BTFR sigma (v2.0.1 says the relation does not describe individual curves to the data's precision; median chi^2 25.5 vs 6.0 free); the Milky Way class (compact massive discs, not a SPARC-like sample); clusters; anything beyond the last measured point.

## 5. Failure branches, stated in advance
P1 FAIL -> the universal-density claim of sec 3.1 is withdrawn for the larger sample; per-galaxy arm reported alone. P2 KILL -> sec 6.2 and the source rule of sec 2.3 fall (as v2.0.1 says); the relation keeps no environmental prediction. P3 width FAIL -> B.1 rewritten: rho_c tight, the SPARC spread was sample-limited (a GAIN for the model, to be stated plainly). P3 rule FAIL -> that property becomes a candidate rule, canon ruling needed. P4 opposite sign -> the shape reading of RESULT_RHOC_CONSTRAINED sec 6 is dead. P5 FAIL -> the profile-variety mechanism of sec 3.1 (ii) is wrong at the stated level. P6 FAIL -> Renzo's-rule print-through is not a strength; sec 7 sentence withdrawn.

## 6. Expectation (Claude, for the record): P1 PASS at 0.1 dex with the low-acceleration offset in slowly rising galaxies (sec 3.1 (iii)) visible; P2 undecided unless >= 1500 galaxies carry field values, then PASS or KILL at 2-3 sigma -- this is the item that matters; P3 width PASS (spread real), no rule; P4 PASS on sign, magnitude uncertain; P5 undecided (error model); P6 undecided at BIG-SPARC resolution.

## 7. Ballot
BS-1 stamp this prereg (keep / amend / kill). BS-2 the sec 7 "Big SPARC sentence" of v2.0.2 points to this file by name and date. BS-3 archive this file in the v2.0.2 software release.
## 8. STAMP
> **JM 2026-11-08 18:02:** "accepted"  -> BS-1 KEEP as written: P2's -0.08 re-derived for the BIG-SPARC mass range and frozen BEFORE the run (recipe + catalogue named then); P3 FIXED-nuisance arm primary; P4 scored only if RESULT_RHOC_CONSTRAINED secs 4/6 are ratified; sec 2a compute plan. BS-2, BS-3 follow in v2.0.2.
# END

## ADDENDUM 2026-11-14 18:40 -- BS-AMEND-1 (i)-(ii), per RULING NP4-8 (JM 2026-11-14 18:38 "accept all"; BALLOT same stamp)
Appended after the stamped text, which is unchanged above. No BIG-SPARC data have been touched.
- P3 width statistic: the robust spread, half the 16-84 range of the per-galaxy fitted log rho_c (EP-4), on galaxies with at least three fitted points
  (EP3-4), with no fit at a scan edge.
- Fixed arm: fits by BD-2 (b) (first interior minimum of chi^2 in log rho_c, scan 4.5-13.0 at 0.1, refined; V3), with the NP4-6 floor clause. Under
  this rule no SPARC fit sits at a bound.
- SPARC values in this statistic: fixed arm 0.47 dex (141 galaxies, V3 table of record md5 82dce260; full 16-84 range 0.93; median 7.42);
  marginalised arm 0.50 dex (113 of the 141 with posterior width < 0.5 dex). [RESULT_EP3_20261114 sec 7; RESULT_NP4_20261114 sec 1]
- The band in P3 ("0.7 +- 0.2" fixed, "1.0 +- 0.25" marginalised) was stated in other statistics (std 0.68; full width 1.01). Its restatement in
  this statistic is OPEN: BS-AMEND-2, for JM, before any BIG-SPARC data. Claude's proposal: fixed 0.47 +- 0.12, marginalised 0.50 +- 0.12 (the
  original tolerances halved with the widths).
- Not covered here: BS-AMEND-1 (iv), edge-on galaxies (Claude 2026-11-12: do not exclude; report as a class).

## ADDENDUM 2026-11-14 21:39 -- BS-AMEND-2 and BS-AMEND-1 (iv), per RULING 2026-11-14 21:33 (JM "opens accepted then batch 8"; BALLOT same stamp)
Appended after the stamped text and the 18:40 addendum, both unchanged above. No BIG-SPARC data have been touched.
- BS-AMEND-2: the P3 band in the statistic of the 18:40 addendum (robust spread) is fixed arm 0.47 +- 0.12 dex, marginalised arm 0.50 +- 0.12 dex
  (Claude's 18:40 proposal, accepted as written). It replaces P3's "0.7 +- 0.2" and "1.0 +- 0.25" for scoring.
- BS-AMEND-1 (iv): edge-on galaxies (i >= 80 deg) are reported separately as a class (count, median, robust spread). No galaxy is cut for its
  inclination.
