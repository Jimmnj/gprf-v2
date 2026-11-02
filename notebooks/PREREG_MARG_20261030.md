# PREREG_MARG_20261030 [STAMPED JM 2026-10-30 12:20 "accepted - i wouldn't know what to alter"; FROZEN]

Drafted 2026-10-30 11:45 (JM 11:30 "apply and yes"). The marginalised scatter test named in sec 3.1 ("a marginalised fit of the model, which we have not done, would decide") and by all three cold reviews.

## 1. Question
How much of the relation's excess scatter about the RAR survives when each galaxy's nuisance parameters are fitted within their published uncertainties, treated identically for the relation and for the empirical RAR function? And does the fitted rho_c keep its 0.65 dex spread and its Sigma_h / f_gas patterns once those errors are absorbed?

## 2. Sample and harness (frozen)
147 Lelli-style SPARC galaxies, 2696 points (Q<=2, inc>=30, ev/v<0.1), harness pass2_sparc_full_20261022.py (solver unchanged). Baryon models as in the note: disc, bulge, gas (1.33 x HI).

## 3. Free parameters per galaxy and priors (Li et al. 2018 style)
- log Upsilon_disc ~ N(log 0.5, 0.1 dex); log Upsilon_bul ~ N(log 0.7, 0.1 dex) (bulge galaxies only); 3.6 micron.
- Distance D ~ N(D_SPARC, e_D) (SPARC columns), applied to radii, velocities' baryon models and M_b consistently.
- Inclination i ~ N(i_SPARC, e_i), clipped to [20, 90] deg; v_obs rescaled by sin(i_SPARC)/sin(i).
- Arm A (the relation): log rho_c free, flat on [6, 10]; sigma from eq. (2) at the FITTED M_b (Upsilon and D applied). No other freedom.
- Arm B (empirical RAR function, g_dagger = 1.2e-10 fixed): the same nuisance parameters, no rho_c.
- Objective: chi^2 with the SPARC velocity errors plus the Gaussian prior penalties (MAP, as Li et al.). Optimiser: Nelder-Mead from the prior means and the fixed-nuisance best rho_c; two restarts; the better kept.

## 4. Metrics
M1 per-point rms of log(g_obs/g_model) after marginalisation, each arm; and the binned (0.25 dex) scatter in the window log g_bar/a0 = -1.75..-0.25.
M2 intrinsic scatter estimate: the per-bin residual sd minus, in quadrature, the propagated measurement scatter of the fitted model (the error-injection method of sec 3.1 applied to the marginalised curves); each arm.
M3 fitted rho_c: sd of log rho_c across galaxies (before 0.65 dex after removing M_b, Sigma_h); and the partial slopes on log Sigma_h and f_gas with their bootstrap 95% intervals (before: +0.27 and -0.93).
M4 nuisance shifts: median and spread of d log Upsilon_disc, dD/e_D, di/e_i in each arm (a sanity check that the relation is not buying its fit with implausible nuisance values; Li et al. found these consistent with the priors for the RAR).

## 5. Decision rules (fixed now)
D1 The scatter tension is RESOLVED if M1(A) <= M1(B) + 0.01 dex AND M2(A) in the window is <= 0.04 dex. CONFIRMED if M1(A) > M1(B) + 0.03 or M2(A) > 0.06 in two or more window bins. Otherwise PARTLY.
D2 The rho_c patterns SURVIVE if the Sigma_h and f_gas partial slopes keep their sign with 95% intervals excluding zero after marginalisation; they are ERROR-DRIVEN if either interval includes zero.
D3 The relation's nuisance shifts must not exceed the RAR arm's by more than 50% in rms (M4); if they do, D1 is reported as PARTLY regardless, because the fit is being bought with nuisance parameters.
D4 Nothing is adopted on the result; the note's sec 3.1 and B.1 sentences are updated to report it, whichever way it falls.

## 6. Pre-stated expectation
Marginalisation lowers both arms' rms (Li et al. 2018 reach ~0.057 observed for the RAR). The relation's excess shrinks because rho_c absorbs part of the D and Upsilon errors (Spearman 0.83 between low-g offsets and rho_c misfit). Expectation: PARTLY -- excess reduced, not removed; rho_c spread falls to ~0.4-0.5 dex; the Sigma_h slope survives, the f_gas slope becomes uncertain. The older G7 result (a different law family, Upsilon only free) did not rescue that family; it does not bear on this one.

## 7. Procedure
Pilot: 10 galaxies (fixed list: first 10 alphabetically with Q=1), both arms, to time the fit and check convergence; no metrics read. Then the full run in the background (~2-3 h on 2 cores at ~0.07 s per solve and ~400 evaluations per galaxy per arm). Outputs: marg_fit_20261030.py, marg_fit_params_20261030.csv (per galaxy, both arms), marg_fit_points_20261030.csv, marg_fit_report_20261030.txt, RESULT_MARG_20261030.md [UNRATIFIED].
