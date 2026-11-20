# RESULT_AXISYM_R1_20261109 -- the frozen law solved in disc geometry (SPARC, fixed-M/L arm) [RATIFIED 2026-11-09 13:09 -- JM "accepted": AX-1..AX-7 as recommended; BALLOT 13:09 line]
Written 2026-11-09 13:50 EST. Script: axisym_r1_20261109.py (modes C0/C1/full; outputs axisym_r1_<mode>_{gals,points}_part.csv). Queue item: SESSION_HANDOFF_20261109 addendum 11:40 item 3 ("AXISYM R1 code"). Rule 3a: SNAP_20260922_001 sec 0c read first; the file implements canon R1 (2026-10-18) + the 2026-10-21 clarification, not the 9/22 form.
ARM: every number here is the FIXED-M/L arm (Yd 0.5, Yb 0.7, SPARC distances and inclinations as given; rho_amb 35.5). Marginalised-arm numbers differ per galaxy by +-0.4 dex (RESULT_CLEAN_CONC sec 5c).

## 1. What was solved
rho_nb(R,z) = rho_c exp[S(x_b) - S(x_b0) - dphi_nb/sigma^2] - rho_amb, S = 2 tanh(x/2), x_b = phi_b(R,z)/sigma^2, sigma^2 = sqrt(G a0 M_b)/2.
- C0: spherical phi_b from SPARC V_bar (the pass2 construction). Same geometry as every result since TEST 1.
- C1: flattened phi_b from ring sums (stellar disc Sigma = Yd*SBdisk, h_z = 0.196 Rd^0.633; gas exponential disc M = 1.33 MHI with Rg fitted to V_gas, h_g = 0.2 kpc; spherical bulge from V_bul); nb response monopole only.
- full: flattened phi_b; nb response with multipoles l = 2..12 (Green's function, damped fixed point; converged in <= 10 iterations for every galaxy, criterion 1e-3 sigma^2).
ASSUMPTION (as in gate2b): the in-plane model velocity is SPARC V_bar^2 + the nb part; the ring-sum baryon model only sources the response. Ring-sum check: disc mass ratio median 1.01 (IQR 0.94-1.06); in-plane v^2(disc+gas) ring-sum/SPARC median 0.95 (IQR 0.90-0.98).
Fit per galaxy: log rho_c by bounded Brent on [4.5, 10.5], pass2 mask (dV/V < 0.1, V_bar^2 > 0), weights (0.1 V_obs)^-2. Also evaluated at universal 10^7.45.

## 2. Validation (C0 vs the same-law spherical solver)
- VALIDATION TARGET SUBSTITUTED. The handoff named pb_solve_canon's Q1-clean r/Rd slope +0.158 [+0.065, +0.275] and the 10/30 marginalised medians. Both are inapplicable: the first is the 10/04 law (no saturation), the second the nuisance arm. Replacement: pass2 solve_gal and its fixed-arm refits lrc_fix0 (clean_conc_part_20261109.csv). Ballot item AX-1.
- C0 reproduces pass2 solve_gal at 7.45 to <= 2e-4 in v (NGC 2403, DDO 154, NGC 2841); fitted log rho_c within 0.007 dex of lrc_fix0 for all 112 overlapping Q=1 galaxies; medians 7.428 both. PASS.

## 3. Results
(a) Universal rho_c survives disc geometry. Fitted log rho_c, 141 galaxies with fits in all three modes:
    C0 median 7.425 (Q=1, n=91: 7.451), std 0.66; C1 7.440 (7.483), std 0.60; full 7.444 (7.486), std 0.65.
    Shift full - C0: median +0.033 dex, 16-84 +0.003..+0.058, max 0.097. Weakly mass-dependent (Spearman +0.22 with log M_b; by bin: <8.5 +0.008, 8.5-9.5 +0.028, 9.5-10.5 +0.042, >10.5 +0.034). The sign: a flattened baryon potential of the same mass is higher in the plane, so the in-plane baryon term is slightly stronger; the fit compensates with slightly more rho_c (reading, not derived). The shift is below the per-galaxy spread by a factor 20.
(b) The rho_c spread is NOT a geometry artefact: std 0.66 (C0) -> 0.65 (full). R-1 (spread real) stands in disc geometry.
(c) RAR at universal 10^7.45, 2696 Lelli-style points, 147 galaxies: rms data-model C0 0.1551, C1 0.1548, full 0.1560; median residual +0.006 -> +0.011. RAR function 0.1328 (same points). Per-galaxy fit: 0.081 in all modes. Disc geometry changes nothing at the 0.001 level; the fixed-arm gap to the RAR function (0.155 vs 0.133) is the same gap as in the spherical harness and is the marginalised arm's business (RESULT_MARG_20261030: with nuisances free the universal law matches the RAR function, 0.057 vs 0.056).
(d) Shape of the response: the nb isodensity axis ratio at 2 Rd is 0.83 (16-84: 0.74-0.90), 0.89 at 3.2 Rd, 0.94 at 5 Rd, while the baryon isopotential has q = 0.56 at 2 Rd and 0.75 at 5 Rd. The response is much rounder than its source: saturation of the baryon term plus the self-term round it off. The spherical harness is therefore a good approximation for rotation curves (effects <= 0.1 dex in rho_c, 0.03 systematic).
(e) Multipole content: the l>=2 terms change the fit by < 0.01 dex relative to C1 (median full - C1 +0.004). The monopole of the flattened potential carries nearly all of the geometric effect.

## 4. Anomaly
NGC 5907: chi^2(log rho_c) has two minima; C0 and full land at the upper bound (10.49, chi^2 5-6), C1 at 7.55 (chi^2 13.5). Brent artefact, not geometry. Three galaxies at a bound in each mode (C0: 3, C1: 2, full: 3). Excluded from nothing above except that bound values enter the std.

## 5. Referee flags
F1 The baryon model behind the flattened potential is a parametric reconstruction (constant h_z/Rd, exponential gas disc, spherical bulge); the 5% in-plane v^2 deficit (vbratio 0.95) is its systematic. A first-order fix (scaling the ring sums to match SPARC V_bar) was NOT applied.
F2 The in-plane model uses SPARC's V_bar for the baryon force; only the response is flattened-geometry. Consistent with how the RAR itself is built, but a referee may ask for the fully self-consistent version.
F3 Fixed arm only. The 10/30 marginalised arm was not rerun in disc geometry (cost: ~900 CPU-h per arm; not planned).
F4 The +0.03 dex systematic shift is below any quoted precision on 7.45 (the Q=1 median moves 7.45 -> 7.49); if a v2.0.2 sentence quotes 7.45, say "7.45 (spherical) / 7.49 (disc)" or round to 7.5.
F5 Convergence check DONE 14:05: NR 400 -> 800 and criterion 1e-3 -> 1e-4 sigma^2 change the fitted log rho_c by < 0.001 dex and q(2Rd) by <= 0.001 (NGC 2403, DDO 154, NGC 2841). Grid and criterion are adequate.

## 6. What this answers (referee's first question, handoff item 4(i))
"You solved a spherical law for disc galaxies." Answer: solved in disc geometry with l <= 12, the universal rho_c moves by +0.03 dex, the spread by -0.01 dex, and the RAR rms by +0.001; the response is rounder (q ~0.83 at 2 Rd) than the disc that sources it. The spherical harness was not the problem anywhere it has been used.

## 7. Next
- PREREG_FEATURE_AXISYM on NGC 1560 (handoff item 4(ii)): expectation A(disc) ~ 0.02 vs MOND 0.036; criterion A(disc) >= 0.8 A(RAR fn) => spherical harness was the problem for the feature test. Needs the NGC 1560 baryon model in ring-sum form (its data are not in SPARC: use RESULT_NGC1560 extraction + an exponential disc). JM stamp before running.
- v2.0.2: one sentence in sec 3.1 or App B (geometry check) + F4 wording.
# END -- secs 3-6 [UNRATIFIED]

## 8. RATIFICATION BALLOT
AX-1 validation target substitution (sec 2) [KEEP -- the handoff's targets were the wrong law/arm].
AX-2 sec 3(a)-(b) universal rho_c and its spread survive disc geometry (+0.03 dex) [KEEP].
AX-3 sec 3(c) RAR unchanged at the 0.001 level in disc geometry [KEEP].
AX-4 sec 3(d) response rounder than source (q 0.83 at 2 Rd); spherical harness adequate for rotation curves [KEEP].
AX-5 sec 6 referee answer wording [KEEP as framing; ledger entry (cf) if kept].
AX-6 F5 convergence/grid check -- done, passes [KEEP as record].
AX-7 PREREG_FEATURE_AXISYM (NGC 1560) [WRITE next; stamp before run].
