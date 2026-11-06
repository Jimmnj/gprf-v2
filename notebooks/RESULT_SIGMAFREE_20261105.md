# RESULT_SIGMAFREE_20261105 -- sigma free under a BTFR-scatter prior (PREREG_SIGMAFREE_20261105) [RATIFIED JM 2026-11-05 21:31 "accepted" (text applied 21:30)]
Run 2026-11-05 17:25 EST: sigmafree_20261105.py -> sigmafree_part.csv (147 galaxies; fixed M/L, D, i; prior width 0.06 dex on log sigma; also unpenalised).

## 1. Numbers
                      rms log rho_c  half-width  median   chi^2 median / sum
sigma fixed (eq. 2)        0.679        0.517      7.42      6.0 / 2352
sigma prior 0.06 dex       0.678        0.525      7.50      2.3 / 1268
sigma free                 0.673        0.517      7.52      1.9 / 1203
d log sigma (prior): mean -0.004, rms 0.056; slope vs log M_b +0.004 [-0.005, +0.013] dex/dex -> implied BTFR log V slope 0.254.
d log sigma (free):  mean +0.010, rms 0.140; slope vs log M_b -0.012 [-0.046, +0.020] -> 0.238.
Spearman(log rho_c, BTFR residual): fixed 0.65 -> prior 0.54 -> free 0.47. Spearman(d log sigma_free, BTFR residual) = +0.59; vs log M_b = +0.03.

## 2. Criteria
K1 FAIL: the rho_c spread does NOT fall (0.68 -> 0.68 -> 0.67); expectation 0.58-0.64 wrong.
K2 PASS: no trend of the fitted sigma with M_b (slope 0.004 +- 0.005; implied BTFR slope 0.254 vs 0.25) and no mean offset (-0.004 dex). The BTFR slope and the a0 normalisation used in eq. (2) are NOT the culprit.
K3 FAIL: the correlation with the BTFR offset weakens (0.65 -> 0.47) but stays well above 0.3.
K4: a free sigma halves the per-galaxy chi^2 (6.0 -> 2.3 median), as N14 found; the fitted sigma reads V_f (Spearman 0.59 with the BTFR residual).

## 3. Reading [UNRATIFIED] -- RETRACTION of the 17:02 mechanism
(a) The 17:02 chat claim that rho_c "absorbs the BTFR offset because sigma is rigid" is WRONG as a mechanism: letting sigma move does not shrink the rho_c spread and only partly decouples rho_c from the offset. sigma and rho_c move together, not in exchange.
(b) Consistent reading: at fixed D and M/L the galaxies above the BTFR line need BOTH a higher sigma and a higher rho_c, i.e. more nb mass throughout. Marginalisation removes the correlation (0.21) because a higher M/L or D moves M_b, sigma, the baryon curve and the needed rho_c all at once, in the direction that erases the offset. So the rho_c-BTFR correlation at fixed nuisances is the signature of nuisance error (M/L, D) acting on the whole galaxy, not of sigma rigidity. The deflation argument of 17:02 point 1 survives only in this form: part of the fixed-M/L rho_c spread is nuisance error; the marginalised 0.60 dex is the real number. Unchanged conclusion, different mechanism.
(c) BTFR as culprit (JM 17:09): NO. K2 shows the relation's velocity scale is right on average and in slope. What remains is the BTFR's own scatter (0.06 dex rms in log V here, mostly nuisance-driven), which the relation at fixed nuisances cannot absorb into sigma and so puts partly into rho_c -- but only partly (K1), which is the surprise of this run.
(d) The 0.6-0.68 dex rho_c spread survives sigma freedom at fixed nuisances and survives nuisance freedom at fixed sigma. It is robust. B.1's "rho_c carries information none of the baryon properties supplies" stands, stronger than before.
(e) Text (HELD): replace the CORR-2 B.1 paragraph's mechanism sentence with: "At fixed distance and mass-to-light ratio the fitted rho_c tracks the galaxy's offset from the BTFR (Spearman 0.65); the correlation falls to 0.47 if sigma is also freed and to 0.21 when distance and mass-to-light ratio are free, while the spread of rho_c stays at 0.6-0.7 dex in every case. The offset is nuisance-driven, and rho_c's spread is not an artefact of it."
## 4. Ballot: SF-1 scoring; SF-2 retraction (a) accepted into the record; SF-3 text (e) replacing CORR-2's mechanism sentence.
# END
