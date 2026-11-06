# PREREG_RHOC_CORR_20261105 -- what else does the fitted rho_c correlate with? (batch, exploratory, corrected)
Written 2026-11-05 17:05 EST BEFORE any code. [UNRATIFIED]. Trigger: JM 16:48 "anything else it could correlate with?" / 16:58 "you go".

## 1. Inputs
Primary rho_c: pass2_table2_part.csv lrc (fixed M/L, 147). Secondary (S1): marg_fit_params_20261030.csv arm A lrc (MAP with D, i, M/L free). SPARC Table 1 (T, D, f_D, inc, L, Reff, SBeff, Rd, SBdisk, MHI, RHI, Vflat, Q) and rotmod files.
## 2. Variables (13), Spearman with log rho_c, bootstrap 95% CI, Holm correction over the 13 at alpha = 0.05; then partial Spearman after regressing out log M_b, log Sigma_h, f_gas from both.
Physics: (1) BTFR residual log Vf - 0.25 log(G a0 M_b) [Vf > 0 only]; (2) bulge fraction proxy 0.7 vb^2/(0.5 vd^2 + 0.7 vb^2) at r_last; (3) Hubble type T; (4) R_HI/R_d; (5) inner gas fraction vg^2/(vg^2 + 0.5 vd^2 + 0.7 vb^2) at the point nearest R_d; (6) log SB_disk (central disc surface brightness).
Systematics: (7) r_last/R_d; (8) N points; (9) inclination; (10) Q; (11) log D; (12) distance method f_D (Hubble flow = 1 vs direct 2-5; Mann-Whitney); (13) median e_V/V.
## 3. Expectations (Claude, fixed now): (1) positive, significant after Holm (rho_c absorbs sigma error); (7) possibly negative (short curves -> high rho_c); all others null. Any Holm-significant systematic (7-13) is reported as a measurement effect and triggers a check that the B.1 correlations survive controlling it.
## 4. Decision: report every number; Holm-significant physics variables become B.1 sentences (as correlations, no mechanism); nothing changes the model. If (1) is significant, the B.1 reading "rho_c carries information none of M_b, Sigma_h, f_gas supplies" gains the clause "part of which is the galaxy's offset from the BTFR".
# END
