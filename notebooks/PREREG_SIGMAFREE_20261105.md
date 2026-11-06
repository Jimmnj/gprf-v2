# PREREG_SIGMAFREE_20261105 -- separate the two jobs of rho_c: free sigma under a BTFR-scatter prior
Written 2026-11-05 17:15 EST BEFORE code. [UNRATIFIED]. Trigger: JM 17:09 "sure - run it - could btfr itself be the culprit?"
## 1. Question. At fixed D, i, M/L the fitted rho_c tracks the BTFR offset (RESULT_RHOC_CORR). If sigma is allowed to move within the BTFR scatter, does the rho_c spread fall to the marginalised 0.60 dex with NO distance/M-L freedom? And does the fitted sigma show a systematic trend with M_b (which would indict the BTFR slope or normalisation used in eq. 2)?
## 2. Method. Per galaxy (147, Lelli-style points, 10% weights, fixed M/L 0.5/0.7, catalogue D and i): minimise chi^2 + (log sigma - log sigma_BTFR)^2/w^2 over (log rho_c, log sigma), w = 0.06 dex (the observed rms BTFR offset in log V over 120 galaxies; = log sigma since sigma = V/sqrt 2). Solver: the frozen law (pass2 solve_gal) with sigma^2 passed explicitly. Nelder-Mead from (7.45, log sigma_BTFR). Also record the unpenalised (w -> infinity) optimum from the same start (N14 check).
## 3. Criteria (fixed now)
K1 spread: rms of log rho_c (sigma prior 0.06) vs 0.68 (sigma fixed) and 0.60 (marginalised, RESULT_MARG). Expectation: 0.58-0.64.
K2 BTFR culprit: OLS slope of (log sigma_fit - log sigma_BTFR) on log M_b. |slope| < 0.03 dex/dex with CI incl. 0 -> the BTFR slope 1/4 and a0 normalisation are NOT the culprit; |slope| >= 0.03 with CI excl. 0 -> a slope/normalisation mismatch feeds rho_c (report the implied BTFR slope 1/4 + slope). Also the mean offset: |mean| < 0.02 dex -> normalisation fine.
K3 decoupling: Spearman(new log rho_c, BTFR residual) < 0.3 (was 0.65).
K4 fit: median chi^2 improvement vs rho_c-only; report (Occam not scored here; this is a diagnostic).
Decision: K1+K3 pass -> B.1 sentence: "with sigma free within the BTFR scatter and nuisances fixed, the spread of rho_c is X dex and no longer tracks the BTFR offset; the remainder is the quantity the relation does not explain." K2 fail -> report the implied slope; no model change; a BTFR-slope item enters sec 7.
# END
