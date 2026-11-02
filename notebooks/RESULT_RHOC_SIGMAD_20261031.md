# RESULT_RHOC_SIGMAD_20261031 -- rho_c from Milgrom's central-surface-density bound: Stage 1  [UNRATIFIED]

Written 2026-10-31 19:55. Prereg PREREG_RHOC_SIGMAD_20261031.md (stamped "go" 19:38). Script rhoc_sigmaD_20261031.py -> rhoc_sigmaD_20261031.csv/.jsonl, rhoc_sigmaD_run.log. VM, nominal nuisances. Sigma_M = 137 Msun/pc^2 (log 2.14). 145/147 galaxies have a root (2 none in [4,11]).

## Metrics (gamma = 1 primary; target = arm A posterior-median log rho_c at sig_int 0.03)
P1 A - rule: median +0.39 dex, rms 0.65, sd 0.55. (gamma = 2/pi: +0.68 / 0.84; no-gas Sigma_b0: +0.44 / 0.67.) The RULE IS LOW by 0.4 dex, not high.
P2 Spearman(A, rule) = 0.24; slope 0.57. For comparison A's own sd about its mean is 0.57 dex: the rule's rms (0.65) is WORSE than a constant; after removing the median offset, 0.55 = no better than a constant.
P3 A vs log M_b: Spearman -0.21, slope -0.16. Rule vs log M_b: slope -0.33 (closed form says -0.5; the low-Sigma branch flattens it). A regressed on (log M_b, log Sigma_b0): coefficients (-0.245, +0.127), rms 0.54, Spearman 0.18. The fitted rho_c carries almost no M_b or Sigma_b0 information.
P4 per-point RAR rms (dex in g), nominal nuisances, pooled: rule 0.183, universal 0.183, A 0.145 (medians 0.135 / 0.111 / 0.090). Rule beats universal in 65/145 galaxies. Tie.
P5 CDR slope: not computed (Stage 2 / RESULT_CDR recipe).

## Decisions
Q1 rule tracks the fit: NO (Spearman 0.24 < 0.6; rms 0.65 > 0.5).
Q2 no free amplitude: NO (offset -0.39 at gamma 1; -0.68 at 2/pi).
Q3 rule fits curves at least as well as universal: YES, marginally (tie).
-> Stage 2 (NS arm M) NOT triggered by the prereg rule.

## Against stated expectations
Direction right: RIGHT (rule falls with M_b). Spearman 0.4-0.7: WRONG (0.24). Offset rule HIGH +0.2..+0.5: WRONG SIGN (rule low by 0.4). Q3 marginal: RIGHT. Q1 borderline: WRONG (clear NO).

## Reading (UNRATIFIED)
1. Milgrom's bound, converted through the solution's own central column density, does not predict the per-galaxy rho_c the curves want. The failure is not amplitude (a global offset would fix 0.4 dex); it is scatter: the rule and the fitted values are nearly uncorrelated.
2. The deeper finding is P3: the fitted rho_c correlates with neither M_b nor Sigma_b0 (Spearman -0.21 / 0.18; a two-variable fit explains none of the 0.57 dex spread). Combined with RESULT_RHOC_RULES (no information in M_b, Sigma_h, f_gas, x_b0, Sigma_*(0)) and with P4 (rule rms = universal rms), the per-galaxy rho_c is absorbing galaxy-specific curve features and nuisance errors, not a population law. That is the "error-model escape" reading of RESULT_NS D3, now supported from the other side.
3. Consequence for the CDM comparison: A ties DC14 because both have a per-galaxy freedom that absorbs the same idiosyncrasies. The population-level model is the universal one; it ties the RAR function (+110) and loses to per-galaxy haloes exactly as MOND does. The honest framing may be: universal rho_c is the model; per-galaxy rho_c is a diagnostic of what the curves want beyond the law, comparable to halo-fitting freedom, not a second model.
4. Tension to resolve: the CDR requires rho_c to vary with Sigma at the population level (slope 0.58), yet the fitted rho_c shows almost no Sigma_b0 dependence (+0.13 per dex). Either the CDR slope is carried by the background normalisation (RESULT_RHOC_RULES: BG rule gives 0.53 with RAR rms unchanged) -- a rule that is not a per-galaxy fit -- or the CDR test and the per-galaxy fit measure different things (central Sigma_dyn within R_d vs the whole curve). P5 on the Milgrom rule would say which; owed.
5. Caveat: P1-P3 compare rule values at NOMINAL nuisances with A medians fitted with nuisances free. Nuisance freedom (D, i, Upsilon) shifts the fitted rho_c; part of A's scatter is nuisance-driven by construction. A fairer Stage-1b: compare the rule with A's rho_c at A's own posterior-median nuisances (available in the jsonl). Cheap; proposed.

## Proposed next
(a) Stage 1b (nuisance-matched comparison), 10 min. (b) P5: CDR slope under the Milgrom rule via the RESULT_CDR recipe. (c) Note: B.1 gains one sentence -- Milgrom's central-surface-density bound tested as a rho_c rule; tie on curves, no correlation with the per-galaxy values; cite Milgrom 2009, Donato+2009, Gentile+2009, Kaplinghat+2020. (d) Framing decision for JM: universal rho_c as THE model, per-galaxy as diagnostic (reading 3).

## Addendum 21:15 -- Stage 1b (nuisance-matched) and P5 (CDR slope)  [UNRATIFIED]
Script rhoc_sigmaD_1b_20261031.py -> rhoc_sigmaD_1b_20261031.csv, rhoc_sigmaD_cdr_20261031.csv, rhoc_sigmaD_1b_run.log.
Stage 1b: rule evaluated at A's posterior-median nuisances (Upsilon_d, Upsilon_b, D, i; Sigma_b0 scaled by Upsilon). N = 146. A - rule: median +0.37, rms 0.62, sd 0.53; Spearman 0.34 (from 0.24), slope 0.69. Still Q1 NO (0.34 < 0.6; 0.62 > 0.5), Q2 NO. The nuisances move the rule itself by only 0.01 +- 0.10 dex.
Degeneracy: the residual (A - rule) correlates with A's own nuisances: Spearman -0.37 with log Upsilon_d, -0.29 with the distance shift z_D. Higher fitted M/L or distance -> lower fitted rho_c. Part of the per-galaxy rho_c spread is the Upsilon/D trade-off, which supports reading 2 (rho_c absorbs nuisance freedom).
P5 CDR (RESULT_CDR recipe, nominal nuisances, i <= 85, N = 127): OBS beta 0.568 +- 0.050; Milgrom-rule arm M0 beta 0.246 +- 0.059, Delta vs OBS rms 0.217, C2 FALSE. The rule FAILS the CDR slope exactly as universal (U0 0.276) and relation (3) (R0 0.373) did; the low-Sigma branch's +0.17 dex/dex is far too weak. For reference A0 0.626 (pass), B0 (RAR fn) 0.496 (pass), BG-normalisation rule (RESULT_RHOC_RULES) 0.53 (pass).
Answer to the owed question: the CDR and the curve fits are different constraints. The BG rule passes the CDR without improving the curves; the Milgrom rule ties the curves and fails the CDR; only the per-galaxy fit does both, and it does so partly by trading rho_c against Upsilon and D. No rule tested so far does both.
Scorecard: Stage 1b expectation (correlation rises, maybe to pass): rose 0.24 -> 0.34, did not pass. P5 expectation 0.4-0.6: WRONG (0.25).
