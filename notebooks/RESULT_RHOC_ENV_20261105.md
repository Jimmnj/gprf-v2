# RESULT_RHOC_ENV_20261105 -- fitted rho_c vs environment (PREREG_RHOC_ENV_20261105) [RATIFIED JM 2026-11-05 16:48 "accepted": ENV-1, ENV-2]
Run 2026-11-05 16:45 EST on the VM: rhoc_env_20261105.py -> rhoc_env_20261105.csv. 90 SPARC galaxies with a fixed-M/L fitted rho_c (pass2_table2_part) and a Chae+2021 Table 3 environmental field. Field range 1 dex (log e_N,env -2.72..-1.74); rho_c spread 1.0 dex (16-84).

## 1. Numbers (maxclust field; noclust in brackets)
Q1 Spearman(log rho_c, log e_N,env) = -0.146, p = 0.17, 95% CI [-0.37, +0.08]  (-0.119 [-0.36, +0.11])
Q2 slope alone -0.22 dex/dex [-0.80, +0.41]; with log M_b, log Sigma_h, f_gas controlled -0.45 [-1.02, +0.13]  (-0.22 [-0.86, +0.37]; -0.48 [-1.00, +0.11])
Q3 tercile contrast (high - low field) median log rho_c -0.27 [-0.58, +0.21]  (-0.19 [-0.50, +0.32])

## 2. Criteria
K1 (enhancing detection, Spearman >= +0.30): NO.
K2 (null, |Spearman| < 0.20 or CIs include 0): YES -- all CIs include 0; |rho| = 0.15.
K3 (MOND-sign, Spearman <= -0.30): NO. Point estimates are negative throughout (MOND's sign) but not significant; the controlled slope -0.45 misses zero by 0.13 at 95%.
Verdict: K2 NULL, as expected. Bound on an ENHANCING environmental effect on the core density: slope < +0.4 dex per dex (alone) / < +0.13 (controlled) at 95% over a 1-dex field range. A capped external-potential boost of the galaxy's own response at the e^2 level (0.87 dex) over this field range is excluded; a boost of <~ 0.1-0.4 dex per dex is not.

## 3. Reading [UNRATIFIED]
(a) Environment joins M_b, Sigma_h, f_gas, the central baryonic potential and the central stellar surface density on the list of properties that do not set rho_c (App B.1). The invariance of eq. (1) under a constant external potential is consistent with the data; it is not demonstrated beyond the bound above.
(b) The negative point estimates, if real, would be MOND's sign: denser fields, lower fitted core density (suppression read as a lower rho_c). Same sign as App D's residual slope (-0.12 +- 0.16). Not significant in either; noted, not claimed.
(c) Power is limited by the 1-dex field range, the 0.25-0.3 dex field errors, and the fixed-M/L rho_c (anticorrelated with fitted M/L). Marginalised rho_c would tighten it; not run.
(d) Text (HELD): App B.1, after the sentence listing what rho_c does not correlate with: "nor with the environmental field of Chae et al. (2021) (Spearman -0.15 +- 0.11 over 90 galaxies; an enhancing dependence steeper than +0.4 dex per dex is excluded) [RESULT_RHOC_ENV_20261105]."

## 4. Ballot: ENV-1 criteria scored (K2 null). ENV-2 B.1 sentence (d). 
# END
