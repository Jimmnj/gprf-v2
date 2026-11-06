# PREREG_RHOC_ENV_20261105 -- does the fitted central density track environment?
Written 2026-11-05 16:35 EST BEFORE any code. [UNRATIFIED] until JM stamps. Trigger: JM 16:30 "can you run environment test" (proposed 11/04 23:38).

## 1. Question
rho_c is the one per-galaxy quantity no baryonic property predicts (App B.1). An environmental boost of the galaxy's own response (external potential, capped) would appear as a POSITIVE correlation of fitted log rho_c with environmental density. Eq. (1) as frozen (invariant under phi -> phi + const; environment as density only) predicts NO correlation. MOND's EFE acts on the curves, not on a halo density; translated into our fit it would appear as fitted rho_c FALLING with environment (suppression read as a lower core density), i.e. a negative sign.

## 2. Inputs (fixed)
- Fitted log rho_c at fixed M/L: pass2_table2_part.csv column lrc (147 galaxies; the B.1 values). Supplementary: marginalised per-galaxy rho_c from marg_ns_20261030_pergal.csv if a per-galaxy rho_c column exists (S1).
- Environment: Chae et al. 2021 Table 3, data/chae2021/chae2021_table3_eNenv.csv: log e_N,env (maxclust) primary; noclust as a check. Field values exist for ~60-100 SPARC galaxies.
- Controls: log M_b, log Sigma_h, f_gas (B.1 covariates), from pass2 outputs.

## 3. Statistics
Q1 Spearman rho(log rho_c, log e_N,env) with bootstrap 95% CI (1000 resamples).
Q2 OLS slope of log rho_c on log e_N,env alone and with log M_b, log Sigma_h, f_gas as controls; bootstrap CI.
Q3 Tercile contrast: median log rho_c in the high- vs low-field tercile, with bootstrap CI.

## 4. Criteria (fixed now)
K1 DETECTION of an enhancing effect: Spearman >= +0.30 with 95% CI excluding 0, AND the controlled slope positive with CI excluding 0.
K2 NULL: |Spearman| < 0.20, or CI including 0 for both Q1 and Q2.
K3 MOND-sign: Spearman <= -0.30 with CI excluding 0 (and controlled slope negative, CI excluding 0).
Between 0.20 and 0.30 in |rho|: INCONCLUSIVE; report, no text change.
Decision: K2 -> environment joins the B.1 list of properties that do not set rho_c; one sentence in B.1; the invariance of eq. (1) stands as a prediction consistent with data. K1 -> a capped external-potential response becomes a CANDIDATE rule for rho_c (sec 2.4 / sec 7), not adopted; canon ruling needed before any text beyond "candidate". K3 -> report as a MOND-sign environmental dependence of the fitted core density; App D gains the number; no model change without JM.
Expectation (Claude): K2. Caveat: e_N,env is a Newtonian external field estimate from a flux-limited catalogue; its errors (0.25-0.3 dex) and the fixed-M/L rho_c (anticorrelated with fitted M/L, -0.37) both dilute any correlation; a null here is a weak bound (|slope| <~ 0.3 dex per dex at 95%, to be computed).

## 5. Not done here
Marginalised rho_c vs environment with full covariance; region (underdense/overdense) assignments (not available to us).
# END
