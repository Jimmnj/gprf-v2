# RESULT_CLEAN_CONC_20261109 -- clean-data check of the concentration-vs-rho_c scan (RESULT_RHOC_CONSTRAINED_20261108 sec 7) [RATIFIED 2026-11-09 12:55 -- JM "accepted": CC-1..CC-7 as recommended; BALLOT_20260926 12:55 line]
Written 2026-11-09 01:20 EST. Script: clean_conc_20261109.py (per-galaxy table clean_conc_part_20261109.csv). Trigger:
> **JM 2026-11-09 01:06:** "1 - done ; 2 accepted ; 3 go"
Precondition being tested (sec 7, JM 00:00 "it was the ratio"): does the +0.34 Spearman (log r_peak,bar/Rd vs log rho_c, 112 constrained galaxies) survive Q=1, dropping the innermost two points from both the index and the rho_c fit, and a smoother index?

## 1. Setup
- Sample: the 113 galaxies with marginalised posterior width < 0.5 dex (marg_ns_20261030.jsonl, arm A, sig_int 0.03), Rd > 0, >= 5 rotmod points. Q=1: 78; Q=2: 35.
- Two rho_c estimates per galaxy: (i) MARG = median of the marginalised posterior (Upsilon_disc, Upsilon_bul, D, i free under SPARC priors; the sec 7 quantity); (ii) FIX = fixed-M/L 1D refit (pass2 method: Upsilon 0.5/0.7, catalogue D and i, 10% velocity weights, Lelli-style mask), run twice: all points (fix0) and innermost two points masked out (fix2).
- Two indices: r_peak,bar/Rd (sec 7) and the half-baryonic-mass radius / Rd (M(<r) = V_bar^2 r/G), each from all points and from r[2:].
- Columns asserted: Rd(NGC2841) = 3.64, Rd(DDO154) = 0.37.

## 2. Numbers (Spearman of log index vs log rho_c | partial after log M_b)
| arm | sample | index | rho_c | rho (p) | at fixed M_b (p) |
|---|---|---|---|---|---|
| A0 reproduce | 113 | r_peak/Rd, all pts | MARG | +0.33 (<0.001) | +0.21 (0.023) |
| A0 | 113 | r_half/Rd, all pts | MARG | +0.30 (0.001) | +0.18 (0.053) |
| control | 112 | r_peak/Rd, all pts | FIX0 | **-0.00 (0.98)** | -0.06 (0.54) |
| A1 Q=1 | 78 | r_peak/Rd, all pts | MARG | +0.45 (<0.001) | +0.18 (0.115) |
| A1 Q=1 | 78 | r_half/Rd, all pts | MARG | +0.32 (0.005) | +0.02 (0.83) |
| A2 Q=1 | 78 | r_peak/Rd, r[2:] | FIX2 | -0.01 (0.91) | -0.16 (0.15) |
| A3 Q=1 | 78 | r_half/Rd, r[2:] | FIX2 | +0.01 (0.92) | -0.06 (0.63) |
| A3 Q=1 | 78 | r_half/Rd, r[2:] | MARG | +0.31 (0.005) | +0.03 (0.78) |
| all Q | 112 | r_peak/Rd, r[2:] | FIX2 | -0.01 (0.88) | -0.11 (0.24) |
- Dropping the innermost two points moves the fixed-M/L rho_c by 0.00 dex (median |fix2 - fix0|): the inner points do not drive the fit. This is NOT the 9/23 mechanism.
- Bins by r_peak/Rd (r[2:]), Q=1, FIX2 median log rho_c: <1.5 (N 20) 7.59; 1.5-3 (N 21) 7.44; >3 (N 37) 7.49. Flat. (Sec 7 MARG bins were 7.43 / 7.44 / 7.72.)

## 3. Where the signal lives: the nuisance parameters
- FIX0 vs MARG rho_c: Spearman only 0.38; delta = MARG - FIX0 has median -0.01 but 16-84 range -0.41..+0.45 dex. The per-galaxy rho_c is arm-dependent at the +-0.4 dex level, even in the "constrained" subset. The SPREAD is not: FIX0 std 0.59 (16-84 1.07), MARG std 0.50 (16-84 1.02) -- sec 1 (R-1, ratified) stands in both arms.
- delta vs log(D_fit/D_0): rho -0.67 (p<1e-3); vs inclination shift: -0.37 (p<1e-3); vs log Upsilon_disc: -0.22 (0.02); vs log r_peak/Rd: +0.26 (0.005); vs log M_b: -0.05 (n.s.).
- The index itself: log r_peak/Rd vs log M_b -0.67 (giants are concentrated), vs log(D/D_0) -0.22 (0.018), vs inclination shift -0.19 (0.045).
- Reading: in the marginalised fit, a galaxy pulled to a smaller distance or a larger inclination gets a higher rho_c; those pulls correlate (weakly) with concentration and (strongly, through M_b) with mass. Hold the nuisances at catalogue values and the concentration correlation is zero to two decimals.

## 4. Verdict [T3]
- Sec 7's "+0.34" does NOT survive. Not because of inner points (the 9/23 ghost: ruled out, 0.00 dex effect) but because it is a property of the marginalised arm only: nuisance moves (D, i) co-vary with concentration and mass. In the fixed-M/L arm, with or without the inner points, Q=1 or all, both indices: rho = 0.00 +- 0.1.
- In the marginalised arm itself, the fixed-mass partial dies in Q=1 (+0.02/+0.03 for the half-mass index); only the raw correlation survives, and it is carried by mass.
- Recommendation: R-5 (sec 7 concentration) KILL as a physical rho_c correlation. KEEP sec 7 as the record of a dead end, with this file as the reason. The bulge-split prereg (R-4b, parked lead) should NOT use concentration as its variable.
- Not a rule, not quotable, not a P4 line for PREREG_BIGSPARC. BIG-SPARC P4 keeps the ratified secs 4/6 content only.

## 5. Flags a referee would raise (new, from this check)
- F1 ARM DEPENDENCE: per-galaxy rho_c differs between fixed and marginalised arms by +-0.4 dex (16-84). B.1 must say which arm each number comes from; the v2.0.2 B.1 sentences (handoff sec 4) should quote the arm. The 1.0 dex spread is arm-independent (safe).
- F2 BULGE OFFSET (R-4a, ratified KEEP 01:06): bulge galaxies (N 29) MARG median 7.25 vs no-bulge 7.62 (-0.37 raw); FIX0 7.48 vs 7.43 (+0.05). The -0.16-at-fixed-mass result is a marginalised-arm result (as its label says); in the fixed arm it is null (RHOC_CORR primary, consistent). Same mechanism as sec 7? Not tested here (Upsilon_bul free in MARG). Recommend: keep R-4a as ratified but add "marginalised arm only; null at fixed M/L" wherever it is quoted; the bulge-split prereg must state both arms.
- F3 MASS TREND (R-2 (iii), "stands marginalised"): MARG vs log M_b -0.25 (p 0.007); FIX0 -0.14 (p 0.14, n.s.) on this sample. Weaker than quoted in the fixed arm. Check against RHOC_CORR primary before any B.1 mass-bin sentence is applied.
- F4 the marginalised posterior "width < 0.5 dex" selects galaxies where D and i are well constrained too; it does not make the MARG median closer to the FIX value (0.38 rank correlation). Posterior width is not accuracy against the fixed-nuisance answer.

## 5b. ADDENDUM 2026-11-09 01:27 -- the sec 6 shape reading (R-4b) split by arm
> **JM 2026-11-09 01:19:** "you said it looked like rho_c tracked "shape" of rar?"  / 01:24 "ok" (run it)
Sec 6's pattern was: rho_c low where phi_b departs from an exponential disc at either end (bulges 7.25; slowly rising gas-rich dwarfs 7.24), universal for flat discs (7.6). Medians of log rho_c, constrained 113, MARG | FIX (fixed M/L, catalogue D and i):
  ALL 7.55 | 7.43;  bulge (N 29) 7.25 | 7.48;  no-bulge flat discs lMb>=9.5 (N 43) 7.57 | 7.42;  dwarfs lMb<9.5 (N 38) 7.73 | 7.46;  gas-rich dwarfs (N 32) 7.66 | 7.40;
  no V_f slow risers (N 13) 7.25 | 7.15 (Mann-Whitney vs flat discs p 0.21 FIX, 0.16 MARG);  no V_f & lMb<9.5 (N 9) 7.24 | 7.45;  gas-rich any mass (N 45) 7.69 | 7.33;  gas-poor no-bulge (N 14) 7.39 | 7.48.
  Spearman vs log M_b: MARG -0.25 (p 0.007) | FIX -0.14 (p 0.14);  vs f_gas: MARG +0.20 (p 0.03) | FIX -0.11 (p 0.26).
- In the fixed arm every subclass median sits within 0.1 dex of 7.45 except the 13 slow risers (7.15, not significant, and 7.45 for the low-mass ones among them). Both ends of the sec 6 pattern are marginalised-arm features. The "shape of phi_b" reading has no support once D, i and Upsilon are held at catalogue values.
- What stands in both arms: the within-class spread (~1 dex). What moves between arms: the class medians, by up to 0.4 dex, in the direction the nuisances pull.
- Recommendation: R-4b (shape-absorption lead) KILL. Queue item 5 (bulge-split prereg) loses its motivation; drop it, or re-motivate it from the axisym run (item 4) if disc geometry changes the bulge medians. R-4a (bulge offset -0.16, ratified) stands as labelled: a marginalised-arm statement; null at fixed M/L.
- Framing consequence (for the ledger, not yet written): "the spread in rho_c is real; its correlations with galaxy type are not robust to the treatment of distance, inclination and M/L" is the safe sentence. The universal value 7.45 is the median of every subclass in the fixed arm.

## 5c. ADDENDUM 2026-11-09 01:35 -- what the arm split hints at: the lever arms [UNRATIFIED]
> **JM 2026-11-09 01:30:** "interesting pattern - does it hint at anything"
- FIX arm, Spearman vs log rho_c (N 100-112): BTFR residual (log V_f - 0.25 log(G a0 M_b)) **+0.72** (p<1e-3); r_last/Rd +0.16 (0.08); npts +0.06; inc -0.03; log D -0.19 (0.05); error floor -0.07. Linear slope: **6.2 dex of rho_c per dex of V_f offset**; BTFR-residual std 0.062 dex -> 0.38 dex of rho_c scatter; FIX std 0.59 -> 0.45 after removing the BTFR trend.
- MARG arm: BTFR residual +0.06 (n.s.) -- the free D, i, Upsilon absorb the V_f/M_b offset and rho_c no longer carries it; instead rho_c carries the nuisance pulls (sec 3).
- delta(MARG - FIX) vs log(D_fit/D_0): slope -5.2 -> effectively **rho_c ~ D^-5** at fixed angular data (sigma^2 ~ sqrt(M_b) ~ D rises with D; the fit lowers rho_c to push the nb rise outward). Pull std 0.066 dex -> 0.34 dex.
- HINT [T3]: with sigma pinned to the BTFR, rho_c is the one knob that reconciles the BTFR-implied asymptote with the observed velocities, and it does so through steep levers: ~ (V_f/V_BTFR)^6 and ~ D^-5. A universal rho_c is then close to the statement "zero BTFR scatter"; the per-galaxy spread is BTFR scatter + distance/inclination error, each magnified five- to six-fold. Same content as N14 (free sigma "just reads V_f"): rho_c at fixed sigma does the same job.
- CONSEQUENCE either way: (i) B.1 must state the lever arms -- a 10% distance error is 0.2 dex in rho_c; (ii) the universal-vs-fitted comparison is dominated by how D and i are handled; (iii) OWED CHECK (cheap, pre-register first): propagate the CATALOGUE errors (eD, e_inc, 0.1 dex Upsilon) through the fixed-M/L fit (pass2 15b style, ~20 realisations x 112 galaxies, ~10-15 min chunked, or WSL) and compare the predicted rho_c spread with the observed 1.0 dex. If errors alone give >~0.7 dex, R-1 "spread required by the curves" weakens to "spread consistent with measurement error" -- which would FAVOUR universality, not hurt it. If they give <~0.4 dex, R-1 stands and the residual 0.45 dex is physical. Caveat: the marginalised arm already lets D, i move under SPARC priors and still finds 1.0 dex spread with 0.28 dex widths (sec 1); the check must explain why, or the hint fails.
- Referee flag F5: a parameter with D^-5 sensitivity is weakly testable as a universal constant from rotation curves alone; the note's 0.68/1.0 dex spread statements need the lever-arm sentence beside them.

## 5d. ADDENDUM 2026-11-09 01:40 -- is the distance sensitivity a systematic? [UNRATIFIED]
> **JM 2026-11-09 01:32:** "i'm surprised shifting with distance is a factof i would guess it shouldn't matter unless it hints at a possible systematic -"
- Distance enters every rotation-curve fit: at fixed angular data V_obs is unchanged while V_bar^2 ~ D (M ~ D^2, r ~ D). Standard; SPARC MOND fits adjust D for the same reason (Li+2018). What is specific to the law is the steepness: sigma^2 ~ sqrt(M_b) ~ D pins the asymptote, so rho_c -- which has leverage only in the core -- must undo both the baryon rise (~D) and the asymptote rise (~D) at the sampled radii. Inner-core scaling alone gives rho_c ~ D^-2; the saturated form steepens it to the measured ~D^-5. Structural, not a systematic by itself.
- Systematic test = does the fitted rho_c know the distance? FIX arm: vs log D raw -0.19 (p 0.047); at fixed log M_b -0.15 (p 0.11, n.s.); at fixed BTFR residual -0.25 (p 0.012). log D vs log M_b = +0.68 (Malmquist: distant = luminous), so the distance and mass trends are not separable on this sample. Distance METHOD: Hubble-flow (fD=1, N 63) FIX median 7.42 vs direct (N 49) 7.43, Mann-Whitney p 0.96 -- clean (RHOC_CORR primary agreed: -0.13, p 0.25; the MARG arm shows -0.27, p 0.02 -- again the nuisance arm).
- Verdict [T3]: no distance systematic detectable beyond the (weak, giants-low) mass trend in the fixed arm; the steep lever is a property of the law. It makes rho_c a poor universal-constant test from rotation curves (F5) and makes the error-propagation check (5c (iii)) the right next step. If a distance systematic exists it would show in BIG-SPARC (PREREG_BIGSPARC) as a rho_c trend with D at fixed M_b -- add as a diagnostic there, not a prediction.

## 5e. ADDENDUM 2026-11-09 01:50 -- does the rho_c spread shrink with data quality? (the "Hossenfelder test") [UNRATIFIED]
> **JM 2026-11-09 01:45:** "i think sabine said something about as the data gets better , mond looks weaker - mcgaugh claims the opposite"
Spread of log rho_c (16-84 / std) by data quality, constrained 113, FIX | MARG:
  Q=1 (78) 0.94/0.52 | 0.86/0.50;  Q=2 (35) 1.11/0.73 | 0.97/0.51;  eD/D < 0.12 (41) 1.13/0.55 | 1.07/0.50;  eD/D >= 0.25 (35) 1.15/0.60 | 0.92/0.46;
  direct D (50) 0.95/0.63 | 0.89/0.44;  Hubble-flow D (63) 1.12/0.56 | 1.05/0.54;  Q=1 & eD/D<0.12 & e_inc<=5 (19) 0.86/0.55 | 0.97/0.54;  npts>=15 (71) 0.97/0.62 | 0.88/0.48.
- The spread does NOT shrink with better distances (best-distance third 1.13 vs worst third 1.15), barely with Q, and the cleanest 19 galaxies still span 0.86 dex. If the spread were distance error through the D^-5 lever, the best-distance galaxies would be visibly tighter. They are not.
- CONSEQUENCE: this pre-empts the 5c (iii) hypothesis on real data -- measurement error cannot be the bulk of the 1.0 dex. R-1 (spread real, ratified) is REINFORCED. The error-propagation check is now a quantification (how much of the 1.0 dex is error; expected ~0.3-0.4) rather than a test that could flip R-1. Demote 5c (iii) from "owed" to "optional quantification"; keep F5 (lever-arm sentence) as is.
- Context: Hossenfelder's claim (video "New data spell trouble...", 2026) is that for rotation curves and wide binaries alike "the higher the uncertainty of the data, the better MOND seems to work"; McGaugh's reply (Triton Station) is that MOND's residual scatter in credible galaxy data is fully accounted for by known errors; Cookson+2026 (quality framework, WBs) find apparent MOND signals diminish as rigour improves. For the law, the analogous statement is the opposite of an artefact: a real per-galaxy quantity that better data does not erase. The same statement is bad news for universality: better distances will not rescue a universal rho_c.
> **JM 2026-11-09 01:49:** "is systematics consistent that a median value works better (spread wise) - but your test disagrees"
- Agreement, not disagreement: under pure systematics a universal (median) rho_c would fit each galaxy about as well as its own noisy value (shrinkage), and the best-measured galaxies would sit nearest the median. Observed: A beats U by ~600 in log Z (CDM-arms table 10/31) and 52/147 galaxies sit > 3 posterior half-widths from 7.45; and (5e) the best-distance galaxies are no tighter. Both tests say the per-galaxy rho_c is real; the levers inflate the scatter by a few tenths of a dex, they do not create it.

## 6. Reproduce
python3 clean_conc_20261109.py   (~70 s; resumable; prints the table) ; sec 3 diagnostic: join clean_conc_part_20261109.csv with marg_fit_params_20261030.csv (arm A) on gal; delta = lrc_marg - lrc_fix0; dD = log10(D/D0); di = inc - i0; Spearman as listed.

## 7. Files: clean_conc_20261109.py, clean_conc_part_20261109.csv (113 rows), GPRF_CleanConc_ArmSplit_20261109.xlsx (class medians, arm table, per-galaxy; JM 01:27), this file. UNCOMMITTED (after 077d698).
# END -- secs 2-5 [UNRATIFIED]
