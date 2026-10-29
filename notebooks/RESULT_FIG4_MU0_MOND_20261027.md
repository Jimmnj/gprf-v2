# RESULT_FIG4_MU0_MOND_20261027 -- Fig 4 redo: MOND per-galaxy core surface density  [UNRATIFIED]
Scripts: fig4_mu0_mond_20261027.py (figure + CSV fig4_mu0_mond_20261027.csv), fig4_mu0_regress_20261027.py.
Ruling: BALLOT 'RULING 2026-10-27 22:36 -- Fig 4 / sec 4.1 split'.

## 1. MOND formula (VERIFIED 2026-10-27 from Milgrom 2009, MNRAS 398, 1023, arXiv:0909.5184, read in full via ar5iv)
- eq (6)/(8): Sigma_M = a0/2piG = 138 (a0/1.2e-8 cgs) Msun/pc^2, log 2.14.
- eq (7): HSB value Sigma_c* = gamma Sigma_M, gamma = 2 lambda/pi, lambda = int_0^inf L dx, L = x mu'/mu; gamma 0.7-2 over the mu forms Milgrom lists (2/pi limiting form; 1 standard form).
- eq (10): X0 mu(X0) = (2/3) Sigma_b/Sigma_M, with Sigma_b = rho R (mean density x extent) -- Milgrom's definition.
- eq (11): Sigma_c* ~ Sigma_M {(6/pi) X0 [1-mu(X0)] + int_0^X0 L dx}.
- eq (12) LSB limit: ~2.4 (Sigma_b/Sigma_M)^1/2 Sigma_M. Milgrom calls the estimate rough.
- Sigma_c* is defined to compare with Donato's rho_0 r_0 (Burkert) -- same quantity as Fig 4's y axis.

## 2. Our choices (state in caption/text)
- mu = inverse of the RAR nu (McGaugh+16), the nu used elsewhere in the note. Gives lambda = 3.29, gamma = 2.09 (top of Milgrom's 0.7-2 range): MOND's HSB value with this mu is log 2.46, NOT log 2.14. The dashed line on the old Fig 4 (a0/2piG labelled 'MOND value') was therefore neither MOND's per-galaxy prediction nor, for this mu, its HSB value.
- Sigma_b: R = R_h = 1.678 R_d (SPARC), rho = (M_b/2)/(4pi R_h^3/3) -> Sigma_b = (3/4) Sigma_h, Sigma_h = M_b/2piR_h^2 (the note's Sigma_h). Sensitivity Sigma_b x0.5 / x2 shown.
- Sample: the 120-galaxy TEST1 set (the old Fig 4 set); Burkert data = Li+2020 fits with r0 < r_max (109 matched).

## 3. Numbers
Mass slope (log mu0 vs log M_b), medians, sd:
| set | slope | median | sd |
| this work (fitted rho_c, SPARC M/L) | +0.240 | 2.24 | 0.39 |
| MOND per galaxy | +0.145 | 2.31 | 0.20 |
| MOND, Sigma_b x0.5 / x2 | +0.155 / +0.130 | 2.19 / 2.42 | |
| SPARC Burkert (data) | +0.209 | 2.29 | 0.34 |
Per galaxy vs data (N=109): this work median offset -0.05, rms 0.26; MOND -0.06, rms 0.30.
Two-variable fit log mu0 = a + b log M_b + c log Sigma_h (bootstrap 2000; corr(logM_b, logSigma_h) = 0.62):
| data | b +0.14 +- 0.04 | c +0.20 +- 0.07 |
| this work | b +0.16 +- 0.05 | c +0.22 +- 0.07 |
| MOND | b 0.00 | c +0.40 +- 0.01 (by construction: mu0 is a function of Sigma_b only) |

## 4. Reading (Claude, for JM)
(a) Mass trend: JM's 21:20 eye is right. Data +0.21; this work +0.24; MOND per galaxy +0.15 (0.13-0.16 over x4 in Sigma_b).
(b) Decomposed: the data carry BOTH a mass term (+0.14, 3.5 sigma) and a surface-density term (+0.20, 3 sigma). MOND has the Sigma term at twice the data's coefficient and no mass term. This work matches both coefficients.
(c) CAVEAT 1 -- shared data: this work's mu0 comes from rho_c fitted per galaxy to the same rotation curves the Burkert fits use; agreement is partly the same data seen twice, not a prediction. MOND has no per-galaxy freedom (given M/L and mu). This is 'comparable description, different parameter count' (MU0-FRAMING 21:31), not a win.
(d) CAVEAT 2 -- M/L: this work's mass slope falls to ~0.12 at disc M/L 0.7 (RESULT_TEST1 addendum 17:25). MOND's at M/L 0.7 NOT recomputed.
(e) CAVEAT 3 -- Milgrom's eq (11) is an order-of-magnitude estimate (his word: rough); the mu choice moves MOND's normalisation (gamma), not its slopes; Sigma_b definition moves both a little (table).
(f) Consequence for sec 4.1: in mu0, the fitted rho_c DOES carry a Sigma dependence (+0.22). The old sentence 'the dependence on surface density at fixed mass ... is present in the data and not in the relation' is false for mu0 with fitted rho_c. It is true for the tail AMPLITUDE (sigma set by M_b alone; outer-speed residual). The 3.3/4.1 split puts each claim where it holds.
(g) Not claimed: that the mu0 comparison discriminates. Both describe it; one with a0 alone, one with a density per galaxy.
