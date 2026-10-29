# RESULT TEST 1 -- FRM-4-FREEZE law vs SPARC (2026-10-18 15:45, Claude, [UNRATIFIED])
Script: test1_capped_baryon_20261018.py (md5 at write 420338480f42...). Output: test1_capped_baryon_20261018.csv, test1_report_20261018.txt. Run in desktop VM, 2 workers, ~2 min.
Law: BALLOT PROPOSED FRM-4-FREEZE 15:30. sigma fixed by BTFR (s2 = sqrt(G a0 M_b)/2), rho_c the only free parameter, cap on the baryon exponent only.
Pass criteria declared in the script header BEFORE the run.

## Validation
lin law (no cap) reproduces N14 lrhoc_btfr / chi_btfr: d log rho_c median 0.000, max 0.009; chi ratio 1.000. Machinery is the N14 machinery.

## Numbers (N = 120; Q=1 N = 83 gives the same picture)
x_b0 = phi_b0/s2 at BTFR sigma: median 2.50 (16-84: 1.05-7.83); 59% of galaxies have x_b0 > 2 -> the cap acts on most of the sample.
| law | median chi | std log rho_c | resid vs pred (dex) | r(fit, pred) | self-term / s2 |
| lin (uncapped) | 28.0 | 1.13 | 1.54 | 0.56 | 11.6 |
| tanh2 (MAIN) | 13.8 | 0.67 | 0.67 | 0.13 | 12.2 |
| tanh1 | 15.6 | 0.73 | 0.72 | 0.07 | 12.4 |
| tanh4 | 17.9 | 0.50 | 0.54 | 0.42 | 11.9 |
| rat2 | 15.4 | 0.76 | 0.75 | 0.13 | 12.5 |
Reference: N14 free sigma, no cap: median chi 5.6. Diagnostic tanh2 with sigma free: chi 6.1, fitted sigma = BTFR to +0.01 dex (median).

## Verdict against the declared criteria
P1 PASS: capping halves chi at BTFR sigma (28.0 -> 13.8). Saturation is SUPPORTED by the rotation curves, not only required by GCs/stars.
P2 FAIL: residual 0.67 dex, bar 0.30. Better than the uncapped baseline at the same sigma (1.54) only because the cap COMPRESSES rho_c
  (std 1.13 -> 0.67); the capped baryon factor no longer PREDICTS rho_c (r = 0.13). rho_c is a universal constant to 0.67 dex.

## What survives (candidate note content, all UNRATIFIED)
1. Self-term = 12.2 sigma^2, stable to +/-0.4 across all four saturation shapes and uncapped (11.6). A universal number, still unexplained.
2. BTFR sigma is the right CENTRE (free-sigma fit lands on it, +0.01 dex), but per-galaxy scatter about it costs chi 13.8 vs 6.1.
3. Under the frozen law: rho_c ~ rho_amb e^{~2} e^{12.2} with 0.67 dex scatter. Not predictive per galaxy.

## Caveats
- N9b's 0.45 dex was at FREE sigma; not directly comparable to 0.67 at BTFR sigma. Like-for-like (free sigma, capped): 0.84.
- tanh4 has the lowest scatter (0.54) but worse chi (17.9). NOT adopted: shape choice after the run would be tuning. Reported as sensitivity only.
- Velocity-space cost of 0.67 dex in rho_c not computed.

## Addendum 16:30 -- P2 anatomy (JM 16:22 "tell me more about p2")
- rho_c is TIGHT per galaxy: +/-0.2 dex raises chi 5-40%; +/-0.67 dex doubles to quadruples it (DDO064, NGC3198, NGC2403, UGC02953). The 0.67 dex scatter is real, not fit noise.
- tanh2: ln(rho_c/rho_amb) = 13.78 + 0.01 x_b0 (scatter 0.67 dex). rho_c is INDEPENDENT of the baryon potential under the cap.
- lin: self-term = 13.63 - 0.66 x_b0 (r = -0.80). The fitted self-term CANCELS two thirds of the baryon term; data want total exponent slope 0.34 on x_b0,
  not 1 (lin) and not ~0 (cap). P2 cannot pass for any law in which rho_c = rho_amb exp(baryon term + universal self-term) with slope 1.
- Residual (self-term) vs galaxy variables: log Mb r -0.27, log Rd -0.29, f_gas +0.15, log SBd +0.06. No strong second variable.
- REFEREE ITEM: core surface density rho_c r_King (r_King = 3 sigma / sqrt(4 pi G rho_c)): tanh2 median 175 Msun/pc^2 (16-84: 98-374), 0.39 dex, slope +0.24 on log Mb;
  lin 279, 0.68 dex, +0.53. Compare Donato+2009 constant halo surface density log(rho_0 r_0) ~ 2.15 +/- 0.2 [mem, VERIFY; Burkert r_0 != King radius].
  Under the cap GPRF lands near the observed value with ~2x the scatter and a residual mass trend. Constant surface density is a MOND prediction (Milgrom 2009 [mem]).

## Addendum 16:55 -- Burkert-equivalent core surface density (JM 16:29 "verify Donato+2009 and use a Burkert-equivalent radius")
Donato+2009 VERIFIED (arXiv:0904.4054; MNRAS 397, 1169): log(r0 rho0) = 2.15 +/- 0.2 Msun/pc^2 (140 +80/-30), Burkert rho0 r0^3/((r+r0)(r^2+r0^2)),
-8 > M_B > -22, from RC decompositions, URC, weak lensing, dSph Jeans. +/-0.2 described as variation across data types.
Scripts: test1b_burkert_surface_20261018.py (md5 8008ed19; model halos), test1c_burkert_data_20261018.py (CONTROL: Burkert fit to SPARC data, same 120,
Yd 0.5, same cuts/weights). CSVs alongside. One units slip (+3) caught and fixed before any number was recorded.
| quantity (N=109, r0 < rmax) | median log mu0 | std | slope vs log Mb |
| SPARC data, Burkert fit | 2.29 | 0.34 | +0.21 (95%: +0.15..+0.27) |
| GPRF capped (tanh2), Burkert fit to model halo | 2.27 | 0.34 | +0.20 |
| GPRF uncapped (lin) | 2.31 | 0.81 | +0.36 |
| universal rho_c (10^7.45) + BTFR sigma, NO per-galaxy parameter | data - pred +0.04 | resid 0.29 | +0.25 (predicted) |
READING:
1. The 16:30 "model has a mass trend that should not be there" is WITHDRAWN: SPARC itself, by Burkert fit, has slope +0.21 and scatter 0.34 in this sample.
   Donato's flat relation is excluded in OUR sample/method (slope > 0 in 2000/2000 bootstraps). Depends on Yd 0.5 and method; literature on a mu0-luminosity
   trend exists [mem, VERIFY before citing].
2. Per-galaxy model-vs-data mu0: scatter 0.12 dex, r 0.93. PARTLY BUILT IN (rho_c fitted to the same curves; same caveat as N15). Not a prediction.
3. The non-trivial piece: a single universal rho_c plus the BTFR sigma predicts mu0 = rho_c r_core ~ sigma ~ Mb^(1/4): slope +0.25 vs observed +0.21
   (inside the 95% band, ~1.3 sigma high); normalisation offset +0.04 dex (one global constant); residual 0.29 vs data spread 0.34 (r 0.52).
   Structural origin: isothermal core r_c = 3 sigma/sqrt(4 pi G rho_c) -> rho_c r_c = 3 sigma sqrt(rho_c/4 pi G). Inputs: BTFR sigma (a0) + one constant.
4. P2 restated in the surface-density observable: rho_c's 0.67 dex enters mu0 as rho_c^(1/2) -> ~0.33 dex, matching the 0.29 residual. Same failure, smaller face.
  XLSX (2026-10-18 16:45): GPRF_TEST1_Burkert_Surface_20261018.xlsx (a99429fa) -- Summary / Per galaxy / Constants / Reading; stats are live formulas. Correction: chat 16:55 'uncapped' row was N=120 (2.31/0.81/+0.36); matched N=109 = 2.35/0.71/+0.37.

## Addendum 17:00 -- is a mass-independent rho_c real? (JM 16:50 "surprised ... constant core density for different mass galaxies")
DATA (Burkert fit to SPARC, N=109, log Mb 8.0-11.4): log rho0 slope vs log Mb -0.06 (95% -0.19..+0.06), median 7.68, scatter 0.52 dex; log r0 slope +0.27 (+0.21..+0.33).
MODEL capped log rho_c: slope -0.02 (-0.17..+0.13), median 7.46, scatter 0.68. MODEL uncapped: slope +0.56 (+0.38..+0.76) -- REJECTED by the data slope.
Binned medians (data rho0 / model rho_c): 7.61/7.54, 7.81/7.63, 7.69/7.42, 7.43/7.39, 7.61/7.48 for log Mb bins 8-9, 9-9.5, 9.5-10, 10-10.5, 10.5-11.5.
WHY in the law: ln(rho_c/rho_amb) = S(x_b0) + Phi_self/s2. Both terms are dimensionless and (capped) mass-independent; rho_amb is universal -> rho_c carries no mass.
  Mass enters only via sigma, which sets the core SIZE r_c = 3 sigma/sqrt(4 pi G rho_c) ~ Mb^(1/4) (data r0: +0.27). Uncapped, x_b0 rises 1 -> 8 with mass -> rho_c rises.
  Universal Phi_self/s2 (~12) = every galaxy's nb truncated at the same r_t/r_c (~500, N9b) -- the unexplained piece.
CAVEATS: "constant" = no trend, with 0.5-0.7 dex scatter. Yd fixed 0.5 (a mass-dependent M/L would tilt data rho0). r0 < rmax selection. No UFD/dSph extension.
  REFEREE ITEM: Kormendy & Freeman (2004/2016) report halo core density DECREASING with luminosity [mem, VERIFY]; Donato-constant mu0 also implies rho0 ~ 1/r0, falling.
  Our SPARC/Yd 0.5 sample shows neither. Must reconcile before the note claims it.

## Addendum 17:25 -- Yd robustness (JM 17:07 "yes"; variants declared in script header before run)
Script test1d_yd_robust_20261018.py; per-variant CSVs test1d_yd_{0.3,0.5,0.7,mdep}_20261018.csv; report test1d_report_20261018.txt. Validation: Yd 0.5 reproduces TEST 1/1c exactly (max 0.000).
mdep = Yd 0.3 + 0.4 (log L - 8)/3, clipped 0.3-0.7. Slopes vs log M_b, [95% bootstrap]:
| Yd | N(r0<rmax) | model chi | MODEL log rho_c | DATA log rho0 | DATA log r0 | DATA log mu0 | MODEL log mu0 |
| 0.3 | 114 | 20.4 | +0.16 [+0.03,+0.31] | +0.09 [-0.04,+0.20] | +0.22 | +0.31 [+0.25,+0.38] | +0.31 |
| 0.5 | 109 | 13.8 | -0.02 [-0.16,+0.13] | -0.06 [-0.18,+0.06] | +0.27 | +0.21 [+0.14,+0.27] | +0.20 |
| 0.7 | 97 | 13.7 | -0.26 [-0.38,-0.15] | -0.17 [-0.32,-0.04] | +0.30 | +0.13 [+0.03,+0.21] | +0.11 |
| mdep | 100 | 16.1 | -0.30 [-0.40,-0.18] | -0.19 [-0.31,-0.05] | +0.31 | +0.12 [+0.04,+0.19] | +0.09 |
Medians: model log rho_c 7.72 / 7.45 / 7.21 / 7.28; data log rho0 7.85 / 7.68 / 7.55 / 7.61.
VERDICT: "rho_c flat in mass" is NOT robust -- it is a Yd 0.5 result. The sign of the trend is set by the M/L assumption; at Yd 0.7 and mass-dependent
  M/L, rho0 FALLS with mass (Kormendy & Freeman direction [mem]) and the flat slope is excluded. RETRACT the 17:00 / 17:05 "no trend with mass" wording.
ROBUST: (1) model tracks data in every variant (slopes agree within errors; built in by the fit, so consistency only); (2) r0 ~ M_b^(0.22-0.31);
  (3) |slope of rho0| <= 0.3 per dex: rho0 changes by <~1 dex over 3.4 dex of M_b; (4) SPARC median 10^7.2-7.9 vs dSph 10^7.95-8.0 (dSph values are
  dispersion-based, no Yd) -- the near-common density across ~8 dex survives at the 0.1-0.8 dex level.
CONSEQUENCE for (ao): the universal-rho_c prediction mu0 ~ M_b^0.25 matches the data only at Yd 0.3-0.5 (0.31, 0.21); at Yd 0.7 / mdep the data give
  0.12-0.13 and rho0 falls, so universal rho_c is excluded there. (ao) is M/L-CONDITIONAL.
  XLSX (2026-10-18 17:30): GPRF_TEST1d_Yd_Robustness_20261018.xlsx (76144783) -- Summary (live SLOPE/MEDIAN over 4 variant tabs) / Yd 0.3 / 0.5 / 0.7 / mass-dep / Verdict.
