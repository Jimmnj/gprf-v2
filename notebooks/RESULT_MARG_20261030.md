# RESULT_MARG_20261030 [UNRATIFIED]

Run 2026-10-30 12:40-13:30 under PREREG_MARG_20261030.md (STAMPED JM 12:20). Scripts marg_fit_20261030.py (fit), marg_metrics_20261030.py (metrics). Data marg_fit_params/points_20261030.csv (arms A, B), marg_supp_params/points_20261030.csv (supplementary arms U, R). Report marg_fit_report_20261030.txt.
Harness check: at nominal nuisances the new solver reproduces pass2 vmod to 3e-9. External check: arm B (RAR function, Ud/Ub/D/i free) gives per-point rms 0.056, matching Li et al. (2018)'s 0.057 for the same exercise.
MAP fits (Nelder-Mead, two starts), not posterior sampling. Nuisance priors: log Ud N(log 0.5, 0.1), log Ub N(log 0.7, 0.1), D N(D, e_D), i N(i, e_i).

## Verdicts by the frozen rules
| | result |
|---|---|
| D1 scatter tension | **RESOLVED** by the rule (M1 A 0.045 <= B 0.056 + 0.01; max window intrinsic estimate A 0.000 <= 0.04) -- but see caveat 1 |
| D2 rho_c patterns | **ERROR-DRIVEN by the rule** (f_gas interval includes zero); Sigma_h pattern SURVIVES (+0.255 [+0.03, +0.50]) |
| D3 nuisance guard | ok (pooled rms z: A 1.27, B 1.77) |

## Numbers
M1 per-point rms log(g_obs/g_mod), 2696 points: A (relation, rho_c free) 0.045 | B (RAR function) 0.056. Fixed-nuisance baselines: universal 0.155, relation (3) 0.140, RAR fn 0.133.
Window (-1.75..-0.25) residual sd: A 0.035-0.046 | B 0.049-0.057.
M2 intrinsic estimate: 0 in every bin for BOTH arms: the injected measurement-only scatter (0.10-0.13 per bin, from D, i and velocity errors at the prior widths) exceeds the post-marginalisation residual, so the subtraction floors at zero. Consistent with Li et al.'s finding of intrinsic scatter consistent with zero, but the method cannot separate 0 from 0.03.
M3 fitted log rho_c (arm A): sd 0.598 (fixed nuisances: 0.68); after removing M_b and Sigma_h trends 0.570 (0.65). OLS: -0.236 log Mb [-0.50, +0.02] +0.255 log Sigma_h [+0.03, +0.50] +0.013 f_gas [-0.83, +0.84].
M4 nuisance shifts (rms in prior sigmas): A -- Ud 1.82, D 0.98, i 0.76 (|z|>2: 27/7/6 of 147). B -- Ud 2.50, D 1.40, i 1.10 (50/16/10).
Chi^2 sums: A 8865 (4-5 par/gal) | B 11050 | U 12232 | R 13853 (3-4 par/gal). Median chi^2 per point: A 0.82, B 1.67, U 1.42, R 1.36.

## Supplementary arms (not pre-registered; added after the A/B run; no verdict)
U = relation at UNIVERSAL rho_c 10^7.45 with the same nuisances: rms **0.057** (fixed-nuisance 0.155). R = rho_c from relation (3) at the fitted M_b, Sigma_h, f_gas: rms **0.057** (fixed 0.140).
So with D, i, Upsilon free within their published errors, the relation at one global density reaches the same per-point rms as the RAR function with the same freedom (0.057 vs 0.056); penalised objective: U better than B in 84 of 147 galaxies, median difference -0.7.
BUT the universal arm pays in distances: |z_D| rms 2.01 (B 1.40, A 0.98); 14 galaxies beyond 3 sigma in D (B 7, A 3), including galaxies with accurate distances: IC 2574 -10.2 sigma (e_D/D 5%), NGC 5055 -6.1 (3%), NGC 0801 -6.6, IC 4202 +7.4, UGC 05764 +7.2. The RAR arm has its own outliers (NGC 5907 -7.7, NGC 0801 -6.2). Arm A (rho_c free) needs almost no distance shifts (max 5.8, F571-8, which every arm moves).

## Caveats
1. D1 'RESOLVED' is partly an artefact of the design. (a) M2 floors at zero for both arms, so its <= 0.04 condition is vacuous. (b) Arm A has one more free parameter per galaxy than arm B; A <= B + 0.01 was always likely. The informative comparisons are U vs B (same parameter count) and the nuisance-shift distributions, neither of which was in the decision rule. By the rule: RESOLVED. By the evidence: the relation's excess scatter at universal rho_c is of a kind that nuisance errors within SPARC's quoted uncertainties can absorb (rms 0.057 vs 0.056), at the cost of distance shifts 2 sigma rms and a tail to 10 sigma that the RAR function does not need to the same degree.
2. MAP, not marginalisation proper. Each galaxy has ~18 points and 3-5 parameters; MAP residuals understate the scatter for every arm. Li et al. (2018) used MCMC; a posterior-sampling version (the 5950X job) is the proper follow-up and is what a referee will ask for.
3. The f_gas slope of rho_c disappears once nuisances are free (it was weakly constrained before: [-1.71, -0.14]). The Sigma_h slope survives at the same value (+0.255 vs +0.268). The mass slope weakens (-0.24, interval touching zero).
4. rho_c's spread is NOT measurement error: 0.68 -> 0.60 dex. The co-variation with the baryons (the fine-tuning statement in sec 8) stands.

## What this changes in the note (proposed, JM to rule)
- sec 3.1: 'A universal rho_c is therefore in tension with the data; a marginalised fit ... would decide' -> report: at fixed nuisances 0.155 vs 0.133; with Ud, D, i free within SPARC's uncertainties (MAP), 0.057 vs 0.056, i.e. the same rms as the empirical function, but with distance shifts of 2 sigma rms against 1.4 for the empirical function and a tail to 10 sigma; a posterior analysis is owed. The 'tension' sentence becomes 'the tension at fixed nuisances is of the size SPARC's quoted distance and inclination errors can absorb; whether the required shifts are credible is the open question.'
- B.1: f_gas term reported as not robust to nuisance marginalisation; Sigma_h term robust.
- Abstract/sec 8 shortfall sentence: 'the scatter about the RAR exceeds the measured intrinsic scatter' -> qualify: 'at fixed distances and mass-to-light ratios'.
