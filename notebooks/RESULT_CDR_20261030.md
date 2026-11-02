# RESULT_CDR_20261030 -- Central density relation test  [UNRATIFIED]
Prereg: PREREG_CDR_20261030.md (stamped 19:38). Script cdr_test_20261030.py; outputs cdr_points_20261030.csv, cdr_fits_20261030.csv, cdr_report_20261030.txt. Run 19:39, 52 s, N=129 (i<=85).

## Recipe check
OBS: beta = 0.578+-0.048, Sigma_crit 1263, residual rms 0.215 dex. Lelli+2016: 0.62+-0.04, ~1571, 0.21. The recipe reproduces the published CDR despite the Sigma_*(0) sampling deviation.

## Scores (primary arms, nominal parameters)
| arm | beta | C1 rms Delta | C2 | C3 | D1 |
|---|---|---|---|---|---|
| A0 per-galaxy rho_c | 0.626+-0.035 | 0.167 FAIL | PASS | PASS | FAIL (C1) |
| U0 universal 10^7.45 | 0.276+-0.071 | 0.214 FAIL | FAIL | PASS | FAIL |
| R0 relation (3) | 0.373+-0.048 | 0.215 FAIL | FAIL | PASS | FAIL |
| B0 RAR function | 0.496+-0.025 | 0.213 FAIL | PASS | PASS | FAIL (C1) |
Secondary (MAP parameters, observed curve rescaled): A 0.626, U 0.542, R 0.460, B 0.606; all pass C2, all fail C1 (rms 0.15-0.17).

## Decisions
D1: no arm passes; C1 fails for every arm including the empirical RAR function (0.213 dex). The 0.10 dex threshold was set without a calibration: the Toomre sum is dominated by V^2/R^2 at the innermost point, where the velocity errors are largest, so no smooth model can match the observed Sigma_dyn(0) galaxy by galaxy to 0.10 dex. C1 is recorded as failed and as miscalibrated; the per-arm rms values stand as the measurement.
D2: YES. Universal rho_c gives beta = 0.28+-0.07 against 0.58+-0.05 observed, 4.3 sigma low, while per-galaxy rho_c gives 0.63+-0.04. The CDR's low-Sigma branch requires rho_c to vary with the baryons.
Not pre-registered, reported as observation: relation (3) also fails C2 (0.37+-0.05, 3.0 sigma low). Its Sigma_h exponent is too weak to carry the CDR. The RAR function gives 0.50+-0.03, Milgrom's (2016) MOND value, against 0.58 observed. At the MAP parameters universal rho_c recovers beta = 0.54: the distance and inclination shifts that buy the RAR fit (RESULT_MARG) also buy the CDR.

## Scratch implication (not a result)
At low Sigma the model's central column is ~ (pi/2)(1+q0) rho_c r_c ~ sqrt(rho_c) sigma. Matching beta ~ 0.58 needs sqrt(rho_c) M_b^1/4 ~ Sigma_*^0.58, i.e. rho_c ~ Sigma_*^1.2 M_b^-0.5 before accounting for the sample's Sigma-M_b covariance. Relation (3) has Sigma_h^0.27 M_b^-0.37. The CDR is therefore a constraint on relation (3), or the source of a better one: fitting rho_c to the CDR directly is a candidate pre-registered test.

## For the note (pending ratification)
- Summary table: CDR row -- per-galaxy rho_c: slope reproduced (0.63 vs 0.58); universal rho_c: fails (0.28); relation (3): fails (0.37).
- Sec 3.1 / B.1: the CDR is published evidence that rho_c is not universal; it joins the fitted-rho_c pattern and (if it holds) the NS evidence.
- Sec 7: relation (3) must be revised to carry the CDR slope; the CDR as a route to the Sigma-dependence.
- Sec 1: cite McGaugh, Lelli, Li & Schombert 2019 as the list of regularities; add Lelli+2016 CDR and Milgrom 2016.
