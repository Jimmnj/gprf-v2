# RESULT_PASS5_SIGMA_FREE_OVERSHOOT_20261022 [UNRATIFIED] -- 2026-10-22 21:55 [time label corrected 22:03]
Handoff 10/22 sec 4 item 2. Diagnostics suggested in SIMULATED reviews (4c: no credit to persons); tests and numbers ours.
Scripts: pass5_sigma_free_20261022.py (+ _fits.csv, _report.txt); pass5b_overshoot_20261022.py (+ _report.txt). Same 147 Lelli-style galaxies, law and solver as pass 2.
PROCESS NOTE: first pass5 run DISCARDED -- scipy Nelder-Mead default initial step from 0 is 0.00025, so 61/147 galaxies never left log q = 0 (median exactly 0.000 in every quartile gave it away). Fixed with an explicit simplex (0.15 dex); after the fix 5/147 sit near 0. Bugged files in _to_delete.
## A. sigma free per galaxy (sigma^2 = q sigma^2_BTFR; rho_c also free; catalogue errors, 3% floor)
- log q median -0.006 (16-84: -0.19..+0.19) -> the BTFR sigma is right on average; per-galaxy scatter ~ +-0.09 dex in sigma (~ +-24%). 5 at the bound |log q| = 1.
- Freeing sigma: median Delta chi2 8.9 for 1 dof; 89/147 improve by > 4.
- Sigma_h dependence WEAK: Spearman +0.14 (p 0.1); quartile medians log q -0.062 / -0.056 / +0.034 / +0.003 (low -> high Sigma_h) = sigma ratio 0.93 / 0.94 / 1.04 / 1.00. Univariate slope sigma^2 ∝ Sigma_h^0.04.
- Multivariate (bootstrap): log q = +0.059(+-0.053) log Sh -0.488(+-0.148) f_gas -0.152(+-0.049) log M_b. The driver is GAS FRACTION (3.3 sigma), not surface density (n.s.); plus a residual M_b slope (3 sigma) at fixed f_gas. Collinear regressors (f_gas, Sigma_h, M_b).
- Size check: LSB sigma 6-7% low <-> outer speed ∝ sqrt(2) sigma 6-7% low = the 3-8% LSB residual. Consistent.
## B. isothermal overshoot (residual at R_last vs R_last/r0, r0 = 3 sigma/sqrt(4 pi G rho_c))
- R_last/r0 median 2.3 (univ), range ~1.4-3.9 (16-84). Correlates with Sigma_h (+0.35 univ, +0.58 rel3).
- Universal rho_c: quartile medians V_model/V_obs 0.980 / 1.032 / 1.036 / 0.989 at R/r0 1.3 / 2.0 / 2.8 / 4.2 -> a HUMP. Quadratic fit: x +0.243(+-0.070), x^2 -0.274(+-0.088) (3.1 sigma), peak R/r0 = 2.8 (F11 9/26 said isothermal v^2 overshoots flat by 27% at ~3 r0 -- same feature). Sigma_h term survives: -0.044(+-0.011) (4 sigma).
- Relation (3): x^2 -0.408(+-0.146) (2.8 sigma), peak 2.3; Sigma_h -0.026(+-0.013) (2 sigma).
## Reading
1. The overshoot is in the data at the predicted place (~3 r0). It is a real shape error of the isothermal form where last points sample 2-3 r0, and it is independent of the Sigma trend.
2. The LSB residual is not a Sigma-modulated sigma in the simple sense: once f_gas is in, Sigma_h carries nothing. The data prefer sigma lowered in gas-rich galaxies (and a slightly different M_b slope). Candidates: (i) the gas weighting in the BTFR sigma (1.33 M_HI; He/H2 treatment), (ii) M/L of gas-poor discs (fixed 0.5), (iii) a genuine local-sigma effect (FKP-type: sigma from M_b(<r) is lower where gas extends beyond the last point -- gas-rich galaxies are exactly those). (iii) is testable: compute sigma from M_b(<R_last) instead of total M_b and refit -- OWED, cheap.
3. Handoff wording 'Sigma-modulated sigma named as candidate' (sec 4.1) should become 'a gas-fraction-dependent or locally set sigma'. JM to rule.
VERIFY: none external. Caveat: catalogue errors with 3% floor (not 10% weights) -- chi2 values not comparable to Table 2's.
