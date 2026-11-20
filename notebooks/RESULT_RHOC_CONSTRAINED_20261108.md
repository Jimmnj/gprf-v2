# RESULT_RHOC_CONSTRAINED_20261108 -- is the per-galaxy rho_c spread real, or prior wandering? [UNRATIFIED; check, not a pre-registered test]
**RATIFICATION 2026-11-09 01:06 (JM "accepted"; BALLOT same stamp): secs 1, 4, 5 (record), 6 bulge offset, 8 KEEP; sec 6 shape reading PARK (lead); sec 7 concentration PARK pending the clean-data check (clean_conc_20261109.py).**
2026-11-08 01:40 EST. Inputs: marg_ns_20261030.jsonl arm A (per-galaxy rho_c), sig_int 0.03, posterior 16/50/84 percentiles of log rho_c; marg_fit_params_20261030 for M_b; SPARC table for V_f, SB_d, R_d, MHI. No new fits. Script: inline (device shell), reproduced in this file's sec 3.
> **JM 2026-11-08 01:27:** "sounds like maybe we're not so far after off after all - where does the .68 'failure' come in?"
Hypothesis tested (Claude): the 0.68-dex spread of fitted rho_c is largely unconstrained fits wandering over the flat prior U[6,10]; keeping only galaxies that pin rho_c would shrink it. RESULT: FALSE.

## 1. Numbers
- Posterior width (p84-p16) of log rho_c: p10 0.15, p25 0.19, p50 0.28, p75 0.45, p90 0.96 dex (flat prior would give ~2.7). 114 of 147 galaxies constrain it to < 0.5 dex; 77 to < 0.3.
- Spread of the medians: all 147: 16-84 = 1.10 dex, std 0.56. Width < 0.5 (N=114): 16-84 = 1.01, std 0.50, range 6.09-8.49. Width < 0.3 (N=77): 16-84 = 0.89, std 0.50. The constrained subset is NOT narrower.
- Deviation from the universal 7.45 in posterior half-widths: 77 galaxies > 2, 52 galaxies > 3 (of 147).
- Correlations on the 114 (Spearman, slope dlog rho_c / dx): log M_b -0.27 (p 4e-3), slope -0.14; log V_f -0.30 (p 2e-3), -0.69; log SB_d -0.02 (n.s.); log R_d -0.27, -0.40; f_gas +0.19 (p 0.04), +0.30; Hubble T +0.15 (n.s.). Residual after the M_b trend: std 0.48 vs 0.50 raw; residual vs SB_d +0.14 (n.s.), vs f_gas -0.01.
- Low end (6.1-6.9): NGC 5371, 6674, 801, 2998, 6015, IC 2574, NGC 6195, UGC 2885, 6614, NGC 2841, 4013 -- massive spirals (and IC 2574). High end (8.1-8.5): UGC 7125, 731, 2259, NGC 24, UGC 5764, 8550, IC 4202, UGC 191, 8490, 4325, NGC 2915, UGC 5721, NGC 3877 -- small/intermediate systems. A weak mass trend, giants low.
- NOTE: the paper's "0.68 dex" is the spread at fixed catalogue nuisances (sec 3.1/B.1); the marginalised-nuisance spread here is larger (1.0 dex 16-84). Different statistic; both say the same thing.

## 2. Reading [T3]
The per-galaxy rho_c is a measurement the curves make, to +-0.1-0.25 dex in most galaxies, and it varies by a decade with at most a weak dependence on mass and none on surface brightness or gas fraction. This confirms B.1 and RHOC_CORR on the subset where the parameter is actually determined, and removes one escape a referee might otherwise offer or attack: the spread is not prior volume. Implication for the text (v2.0.2, B.1): add one sentence -- "114 of 147 galaxies constrain log rho_c to better than 0.5 dex (median posterior width 0.28); the spread among them is 1.0 dex (16-84), so the variation is required by the curves, not by the prior." Implication for the figure of 11/08: NGC 2841 (6.85, width 0.18) is a clean case of a giant demanding a low central density; DDO 154 (7.43, width 0.13) a clean case at the universal value.
What it does NOT say: why. The rule is still open. (The giants-low trend, if it is physical, runs against a naive "compact = dense" and toward a core that does not scale with the baryon potential depth -- cf. B.1's remark that rho_c r_c sits near MOND's surface density; worth a scan of rho_c r_c vs M_b on the constrained set next.)

## 3. Reproduce
python3 - <<'PY'
import json,numpy as np
A=[json.loads(l) for l in open('marg_ns_20261030.jsonl')]; A=[r for r in A if r['arm']=='A' and r['sig_int']==0.03]
med=np.array([r['post']['lrc'][1] for r in A]); w=np.array([r['post']['lrc'][2]-r['post']['lrc'][0] for r in A])
for cut in [0.3,0.5,9]: k=w<cut; print(cut,k.sum(),np.percentile(med[k],84)-np.percentile(med[k],16),med[k].std())
PY

## 4. ADDENDUM 2026-11-08 15:20 -- nuisance and shape diagnostic (device-shell script, same inputs + SPARC_Lelli2016c.mrt Vflat/SBdisk) [UNRATIFIED]
> **JM 2026-11-08 15:02:** "let's do the owed first - what's with low g galaxiex ? not enough rho_c"
> **JM 2026-11-08 15:08:** "gas fraction seems to throw things off too?"
- Flat part (V_f>0), constrained N=100: median log rho_c 7.57, 16-84 7.06-8.00, f_gas 0.33. No flat part, N=14 (27 total): 7.24, 6.92-8.00, f_gas 0.63. Paper sec 3.1 (iii) said 10^7.19 at fixed nuisances; 7.24 marginalised -> (iii) STANDS under both nuisance treatments. Apparent sign conflict with sec 1 of this file RESOLVED: the M_b trend is a bad summary of a U-shape (both ends low, see below).
- f_gas: Spearman +0.19 (p 0.04) raw; -0.01 (p 0.91) after the M_b trend. f_gas is a mass proxy, non-monotone (most gas-rich subset sits LOW). Should not enter a rule; consistent with RESULT_MARG D2 (error-driven) and rule (3) OOS failure.
- SB_disk: -0.02 (p 0.83). Nothing.
- Distance: rho_c vs zD -0.30 (p 0.001); after M_b detrend -0.27 (p 0.004); slope -0.12 dex per sigma_D. Real trade, small: +-2 sigma_D buys +-0.25 dex of the 1.0 dex spread. M/L (lud) -0.12 (p 0.19), inclination -0.14 (p 0.15): n.s.
- Low rho_c (<7.0, N=19): zD +0.06, f_gas 0.18 (giants). High (>8.0, N=19): zD -0.18, f_gas 0.47 (small flat-curve discs). Plus the no-flat dwarfs low at 7.24.
- READING [T3]: the driver is CURVE SHAPE. Slowly rising gas-rich dwarfs want less central nb mass; compact bulge giants (NGC 2841, 6.85) want the response not to track the bulge potential; ordinary flat-curve discs sit high. U-shape in M_b. Same statement as fig_threeway_20261108 in numbers. Supports the rho_c*r_c (core surface density) scan as the next lead: a constant surface density with varying r_c would give exactly a U in rho_c.
- Fix to sec 1: the "weak mass trend, giants low" line is superseded by the U-shape reading; the -0.14 dex/dex slope is reported but not a trend to quote.

## 5. SCAN 2026-11-08 16:20 -- core surface density Sig_c = rho_c * r_c, r_c = sqrt(9 sigma^2 / 4 pi G rho_c), sigma^2 = sqrt(G a0 M_b)/2 (R1), a0 = 1.2e-10, M_b = MAP from marg_fit_params [UNRATIFIED; exploratory scan, not pre-registered]
> **JM 2026-11-08 16:14:** "go"
- 114 constrained: log Sig_c 16-84 0.62, std 0.29 (log rho_c: 1.01, 0.50). Sig_c vs log M_b: Spearman +0.50 (p<1e-3), slope +0.18 dex/dex. r_c vs M_b: +0.80, slope +0.32. Bins (M_b 7-8.5 / 8.5-9.5 / 9.5-10.5 / 10.5-12): Sig_c 7.79 / 8.14 / 8.34 / 8.43; r_c 2.6 / 3.1 / 5.0 / 12.0 kpc. Median Sig_c 2.1e8 Msun/kpc^2 = 212 Msun/pc^2 (MOND Sigma_M 138; Donato+2009 ~140).
- DEAD END, WITH REASON: with r_c defined by the law, log Sig_c = 0.5 log rho_c + 0.25 log M_b + const. The spread halving and the +0.18 slope (= 0.5*(-0.14) + 0.25) are arithmetic. No new information; NOT a universal surface density (rises as M_b^0.18).
- U-shape check: quadratic in log M_b lowers rho_c residual std 0.481 -> 0.446 (coef -0.185). Present, weak. Bin medians rho_c 7.39 / 7.61 / 7.69 / 7.33.
- The informative version needs a MEASURED core radius: rho_c,fit x r_c(Burkert, Li+2020) on the 109 matched galaxies of sec 3.3. Queued as a scan; not run.

## 6. PATTERN CHECK 2026-11-08 16:30 -- is the spread a shape effect? [UNRATIFIED; exploratory]
> **JM 2026-11-08 16:23:** "but is it pointing to a pattern"
- Bulge flag (jsonl 'bul'), constrained 114: bulge N=29 median log rho_c 7.25 (16-84 6.88-7.83); no bulge N=85 7.61 (7.11-8.07); Mann-Whitney p 0.03. Residual after the M_b trend: bulge -0.10, no bulge +0.06 (-0.16 dex at fixed mass). M_b>10.3: bulge 7.22 (N 28) vs no bulge 7.39 (N 19).
- [RETRACTED 23:30 -- COLUMN ERROR: Reff, Rd, RHI were read from the wrong MRT fields (e_L, SBeff, MHI). Correct values in sec 7. Bulge flag, M_b, SB_disk (sec 4), V_flat unaffected.] ~~Reff/Rd -0.11, RHI/Rd -0.04, Hubble T +0.15: nothing.~~
- PATTERN [T3]: rho_c is low where phi_b departs from an exponential disc at either end -- centrally peaked (bulge, 7.25; NGC 2841 6.85) and too shallow (slowly rising gas-rich dwarfs, 7.24); flat-curve discs sit at the universal value (7.6). Candidate reading: the fitted rho_c absorbs a SHAPE error of the saturating response S(x_b) (tuned, with BTFR sigma, to disc potentials), not only a density. Consistent with fig_threeway (NGC 2841 inner hump high, outer low).
- CANON NOTE: a fix would touch S or x_b = a change to the frozen law (FRM-4-FREEZE). Not proposed. Lead only.
- PRE-REGISTRABLE TEST (not written, needs JM stamp as PREREG): bulge galaxies refit with bulge and disc entering x_b separately (saturation on the disc phi_b; bulge handled as stated in the prereg). Criterion fixed before the run: bulge-sample median rho_c within 0.1 dex of the disc median AND rms not worse -> spread was shape; else density. One evening on the 5950X. Rule 3a applies.

## 7. CONCENTRATION SCAN 2026-11-08 23:30 (columns verified: NGC 2841 Rd 3.64, DDO 154 Rd 0.37) [UNRATIFIED; exploratory]
> **JM 2026-11-08 23:21:** "so does the baryon distribution of a disk - deviation from exponential track rho_c ?"
- 112 constrained galaxies with rotmod + Rd. Spearman vs log rho_c | vs residual after the M_b trend:
  log(r_peak,bar/Rd) +0.34 (p<1e-3) | +0.16 (p 0.10);  Reff/Rd +0.19 (p 0.05) | +0.14 (p 0.13);  rms dev of V_bar from a Freeman disc -0.05 | +0.02;  RHI/Rd +0.07 | +0.11;  Hubble T +0.17 | -0.04.
- By r_peak,bar/Rd: <1.5 (N 30, 70% bulge) rho_c 7.43, resid -0.01; 1.5-3 (N 31) 7.44, -0.08; >3 (N 51, 2% bulge, 18% no-flat) 7.72, +0.14.
- READING [T3]: departure from exponential per se does NOT track rho_c; baryon CONCENTRATION does -- concentrated phi_b (bulge, compact disc) -> lower fitted rho_c; extended baryons (LSB/gas discs) -> higher. One monotone axis; mass was partly its proxy (giants are concentrated); ~12% of variance; residual at fixed mass marginal (p 0.10). Supersedes the "U-shape in M_b" description of sec 4 as the primary summary. Consistent with the shape-absorption reading (sec 6). Not a rule; the right variable for the bulge-split prereg and a candidate P4 line for PREREG_BIGSPARC ("rho_c rises with r_peak,bar/Rd, Spearman +0.3 +- 0.1").
- PRECONDITION (JM 2026-11-09 00:00 "it was the ratio"): the 9/23 r/Rd RAR-residual signal (RESULT_RAR_RRD_DEPENDENCE) died in the clean data (RESULT_RAR_RRD_ROBUSTNESS 9/25: Q=1, innermost ~2 points dropped -> slope ~0). r_peak,bar is often the innermost point, so sec 7 inherits that systematic. Before +0.34 is quoted anywhere: redo with Q=1 only, innermost two points dropped from both the index and the rho_c fit, and a smoother index (half-baryonic-mass radius / Rd). Survives -> result; collapses -> the 9/23 ghost.
- CLEAN-DATA CHECK RUN 2026-11-09 01:20: RESULT_CLEAN_CONC_20261109 -- the +0.34 is a MARGINALISED-ARM artefact (nuisance D/i moves co-vary with concentration and mass); fixed-M/L arm gives 0.00 with or without the inner points, Q=1 or all. Inner points: 0.00 dex effect (not the 9/23 mechanism). Recommendation KILL as physics, keep as record.
- Caveat: r_peak is often the innermost rotmod point for bulge galaxies (coarse); a smoother concentration index (e.g. r at half the baryonic enclosed mass / Rd) should replace it before any prereg.
## 8. PRIOR ART IN THE REPO (JM 2026-11-08 23:59 "didn't we notice something like this weeks ago")
- RESULT_G2_SATURATING_20261004 sec 9 (JM 10/4 18:53 "do bulges make things worse for gprf"), pre-FRM-4 law: all 31 bulge galaxies in the high-mass half (G4's mass trend = bulge split); bulge halo piled inward, outer starved; cap 7.5 (bulge) vs 8-8.5 (bulgeless); suspects (a) response over-concentrates in a compact round potential, (b) Upsilon_bul/decomposition; "bulges are easy for MOND, hard for GPRF" (g_N > a0 there). Today on the frozen law with Upsilon_bul marginalised: offset survives (-0.16 at fixed mass), continuous in concentration -> suspect (b) partly excluded, (a) = the sec 6 shape-absorption reading. Same diagnosis, better footing.
- RESULT_RHOC_CORR_20261105 (d): "bulge fraction ... null" on the FIXED-nuisance rho_c vs bulge FRACTION. Today: marginalised rho_c vs bulge FLAG and vs r_peak/Rd finds it. Different arm and statistic; both true; B.1 must not say "no morphology dependence" without this qualifier.
- Lesson: the lead sat in an exploratory section for five weeks. Carry sec 7 into the bulge-split prereg explicitly.
# END [UNRATIFIED]
