# RESULT_RHOC_CORR_20261105 -- what the fitted rho_c correlates with (PREREG_RHOC_CORR_20261105) [RATIFIED JM 2026-11-05 21:31 "accepted" (text applied 21:30)]
Run 2026-11-05 17:15 EST: rhoc_corr_20261105.py -> rhoc_corr_20261105_{primary,marg}.csv, _table.csv. 13 variables, Spearman + bootstrap CI, Holm over 13, partial Spearman after regressing out log M_b, log Sigma_h, f_gas.

## 1. PRIMARY (fixed-M/L rho_c, pass2_table2_part; N = 147, 120 with V_f > 0)
Holm-significant: (1) BTFR residual log V_f - 0.25 log(G a0 M_b): Spearman +0.65 [0.53, 0.74], partial +0.66 (p < 1e-4). (6) log SB_disk: +0.27 [0.11, 0.43] raw, partial +0.08 (p = 0.36) -> this is the known Sigma_h correlation (B.1), nothing new.
Not significant after Holm: log D -0.22 (partial 0.00); r_last/R_d +0.15; R_HI/R_d +0.14; inner gas fraction -0.12; distance method (Hubble flow vs direct) median diff -0.13 dex (p = 0.25); N points +0.09 (partial +0.21, p = 0.011, not Holm-sig); error floor -0.08; T -0.08; Q -0.08 (partial -0.18, p = 0.027); B/T +0.05; inclination +0.04.
Size: OLS log rho_c = 7.40 + 5.4 x (BTFR residual); rms spread of log rho_c 0.68 -> 0.59 dex after removing the BTFR term; 16-84 half-width 0.48 -> 0.35 dex. The BTFR residual itself has rms 0.061 dex (120 galaxies).

## 2. S1 (marginalised MAP rho_c, marg_fit_params arm A; D, i, M/L free)
Nothing Holm-significant. BTFR residual +0.21 [0.02, 0.39] (p = 0.021; computed with CATALOGUE V_f and M_b, not the fitted D/M/L -- caveat). log D -0.19; distance method -0.27 dex (p = 0.024); all others |rho| < 0.14.

## 3. Reading [UNRATIFIED]
(a) Expectation (1) CONFIRMED: at fixed nuisances the fitted core density tracks the galaxy's offset from the BTFR. Mechanism in the fit: sigma is fixed by M_b (eq. 2), so a flat speed above the BTFR value can only be met by raising rho_c (and below it by lowering). Spearman 0.65 means ~40% of the rank variance of rho_c is this channel; in dex it removes 0.09 of the 0.68 rms, more of the core (half-width 0.48 -> 0.35): the tails of the rho_c distribution are the BTFR outliers.
(b) With D and M/L free the channel largely closes (0.65 -> 0.21): the offset is absorbed by the nuisances, as RESULT_MARG found for the spread (0.68 -> 0.60). The remaining 0.6 dex is NOT BTFR scatter and remains unexplained. So the fixed-M/L B.1 correlations are partly this effect; the marginalised ones are the cleaner statement, and the note already quotes both.
(c) Expectation (7) (short curves -> high rho_c): weak (+0.15, not significant). No measurement systematic is Holm-significant. log D and distance method are the only systematics with |rho| ~ 0.2; both vanish on controlling M_b (dwarfs are nearby).
(d) Physics variables other than the BTFR offset: bulge fraction, Hubble type, HI extent, inner gas fraction -- all null. rho_c does not know about morphology or gas distribution beyond what M_b, Sigma_h, f_gas carry.
(e) Text (HELD), App B.1 after the environment sentence: "At fixed mass-to-light ratio and distance the fitted rho_c does track the galaxy's offset from the BTFR (Spearman 0.65 over the 120 galaxies with a flat part): sigma is fixed by M_b, so a flat speed above the BTFR value is met by a higher core density. With distance and mass-to-light ratio free the correlation falls to 0.21; the offset is absorbed by the nuisances, and the 0.60 dex spread that remains is not BTFR scatter. Bulge fraction, Hubble type, HI extent, inner gas fraction, inclination, quality, distance method and the radial extent of the data show no correlation after correction for multiple tests [RESULT_RHOC_CORR_20261105]."
## 4. Ballot: CORR-1 scoring above. CORR-2 B.1 text (e). Caveat for S1: recompute the BTFR residual at the fitted D and M/L before quoting 0.21 as more than indicative.
# END
