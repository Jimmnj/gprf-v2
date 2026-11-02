# PREREG_CDR_20261030 -- Central density relation test  [STAMPED: JM 2026-10-30 19:38 "go"]

Written 2026-10-30 19:10 before any run. Trigger: McGaugh, Lelli, Li & Schombert 2019 (arXiv:1909.02011) list five regularities; the note scores four. The fifth, the central density relation (CDR; Lelli, McGaugh, Schombert & Pawlowski 2016, ApJL 827, L19, arXiv:1607.02145), is a statement about the centre, where rho_c lives.

## Question
Q1. Does the relation reproduce the CDR -- the run of central dynamical surface density Sigma_dyn(0) with central stellar surface density Sigma_*(0), in particular its low-Sigma branch -- (a) with per-galaxy rho_c, (b) with universal rho_c = 10^7.45, (c) with rho_c from relation (3)?
Q2. Is the CDR's low-Sigma slope evidence that rho_c varies between galaxies (i.e. does universal rho_c fail it)?

## Recipe (identical for observed and model curves)
- Sigma_dyn(0) = (1+q0)/(2 pi G) * [ sum_j V_j^2/R_j^2 dR_j + V_N^2/R_N ], Lelli+2016 eq. 6 (Toomre 1963), rings dR_j from midpoints with R_{-1/2}=0, flat-V tail beyond the last point; q0 = 0.2 constant (Lelli+2016 state constant q0 changes little; it cancels in model-vs-observed). Units M_sun/pc^2. All SPARC radii used (no error cut), since the inner points dominate.
- Sigma_*(0) = 0.5 SBdisk + 0.7 SBbul at the first rotation-curve radius (rotmod columns 7,8; Upsilon as Lelli+2016). DEVIATION recorded: Lelli+2016 average the [3.6] profile within the first radius; we have the rotmod sampling only. The same x-axis is used for observed and model, so the model-vs-observed tests are unaffected; the external comparison with beta = 0.62+-0.04 carries this caveat.
- Sample: harness galaxy list (Q<=2, i>=30, vbar2>0) further cut to i<=85 as Lelli+2016.
- Fit: eq. 10 of Lelli+2016, Sigma_dyn = Sigma_0 (1+Sigma_*/Sigma_crit)^(alpha-beta) (Sigma_*/Sigma_crit)^beta with alpha=1 fixed, orthogonal distance regression (scipy.odr) in log-log with equal 0.14 dex errors on both axes (their average formal errors), starting values Sigma_0=1468, Sigma_crit=1571, beta=0.62.
- Arms (primary, nominal nuisance parameters Upsilon_d=0.5, Upsilon_b=0.7, D=D_SPARC, i=i_SPARC): OBS (observed V); A0 (rho_c fitted per galaxy, 1-D bounded minimisation of the harness chi^2, 6<=log rho_c<=10); U0 (log rho_c=7.45); R0 (relation 3); B0 (RAR fitting function, McGaugh+2016, as the empirical comparator). Model V from marg_fit_20261030.Gal.model, i.e. the pass2/pass3 solver unchanged.
- Arms (secondary, for the record): A/U/R/B at the marginalised MAP parameters of marg_fit/marg_supp_params_20261030.csv, with the observed curve rescaled identically (radii x D'/D, V x sin i0/sin i).

## Pre-set decision rules (primary arms)
C1 per-galaxy: Delta = log Sigma_dyn,arm - log Sigma_dyn,OBS. PASS if rms(Delta) <= 0.10 dex and |median(Delta)| <= 0.05 dex.
C2 slope: beta_arm vs beta_OBS from the same recipe. PASS if |beta_arm - beta_OBS| <= 2 sqrt(s_arm^2 + s_OBS^2) (ODR standard errors). Reported alongside: Lelli+2016 beta = 0.62+-0.04; MOND (Milgrom 2016, arXiv:1607.05103) low-Sigma slope 1/2.
C3 scatter: rms residual of arm points about their own eq.-10 fit <= 1.5 x that of OBS.
D1 (Q1): an arm "reproduces the CDR" if it passes C1, C2 and C3.
D2 (Q2): "CDR requires varying rho_c" if U0 fails C2 with beta_U0 < beta_OBS - 2 sqrt(s^2+s^2) while A0 passes C2. If U0 passes C2, the CDR does not by itself discriminate and the note must say so.

## Stated expectations (before running)
A0 passes C1, C2, C3 (it fits the inner curve; little credit). U0 fails C2 flat (central nb column ~ sqrt(rho_c) sigma ~ M_b^1/4 only). R0 between, likely fails C2. B0 passes C2 against OBS and against 0.62+-0.04 within the Sigma_*(0) caveat.

## Outputs
cdr_test_20261030.py; cdr_points_20261030.csv (per galaxy, per arm: Sigma_*, Sigma_dyn, Delta); cdr_fits_20261030.csv (per arm: Sigma_0, Sigma_crit, beta, errors, rms); RESULT_CDR_20261030.md.
No note edits until stamped and the result is written.
