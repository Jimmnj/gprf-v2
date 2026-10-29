# RESULT -- N9b: IS rho_c PREDICTABLE FROM THE BARYONS ALONE?
# gprf\notebooks\RESULT_N9B_RHOC_VS_PHIB_20261017.md
# 2026-10-17 22:37-22:50 (Fable session). Status: ACCEPTED by JM 2026-10-17 23:22 (ballot; note wording: "works for rho_c for dwarfs but .45 dex unexplained"). T2 numbers, T3 reading.
# Trigger: N9 open item N9b; N15 caveat (d); JM 22:37 quoted the bar back: "<= 0.3 dex scatter -> MA5 gate redundant; 0.7-1.1 -> fits only".
# Script: n9b_rhoc_vs_phib_20261017.py -> .csv. Inputs: N14 free fits (lrhoc_free, s2_free; BTFR variant too), SPARC rotmod, Yd 0.5 Yb 0.7.
# Phi_b,0 = int_0^inf vbar^2/r dr (solid-body inside r1, Keplerian beyond rmax; outer share median 0.23). Spherical-equivalent proxy.

## 0. VERDICT (one line)
Partial. rho_c = rho_amb * exp(Phi_self/sigma^2) * exp(Phi_b,0/sigma^2) with Phi_self ~ 12.4 sigma^2 nearly universal and the baryon
factor tracking the fitted rho_c at r = 0.78, slope 0.82, residual 0.45 dex (raw spread 0.73). Not the 0.3-dex bar; not the 0.7-1.1 "knob".

## 1. NUMBERS (N = 120)
 log rho_c (free) vs log Phi_b,0 : slope +0.42, r 0.45, scatter 0.65 dex | vs log SBd: 0.71 | vs log Mb: 0.69 | vs log sigma^2: 0.71
 log rho_c (BTFR) vs log Phi_b,0 : slope +0.78, r 0.53, scatter 0.95 dex (raw 1.13)
 Direct canon definition, pred = log rho_amb + Phi_b,0/(sigma^2 ln10):
   free sigma: pred - fit median -5.38 (16-84 -5.81..-4.80), std 0.47; fit vs pred slope 0.82, r 0.78, scatter 0.45 dex
   BTFR sigma: pred - fit median -5.10 (-5.62..-4.15), std 1.39; slope 0.39, r 0.59, scatter 0.91
 Q=1 only (N=83): slope 0.79, r 0.75, scatter 0.46 (same).
 Phi_b,0 range log 2.72..5.50 km^2/s^2 -> Phi_b,0/sigma^2 ~ 1-3; fits need phi_0/sigma^2 ~ 14.5.

## 2. READING
(a) The 5.4-dex constant: the central depth in the Boltzmann factor is the nb halo's SELF-potential, ~12.4 sigma^2, same for every galaxy.
    Isothermal: depth core->truncation = 2 sigma^2 ln(r_t/r_c) -> r_t/r_c ~ 500; r_c few kpc -> r_t ~ Mpc = the ambient truncation (N10).
    So the normalisation of rho_c is the law's own cosmic-baseline truncation, not a galaxy property. This is N14's "definition, not output"
    made quantitative: rho_c/rho_amb = exp(Phi_self/sigma^2) x exp(Phi_b,0/sigma^2).
(b) The per-galaxy part: exp(Phi_b,0/sigma^2) explains ~60% of the variance of log rho_c (0.73 -> 0.45 dex). Baryon well matters,
    at the right slope (0.82 vs 1 expected), but 0.45 dex is left over.
(c) Against the bar (JM 22:37): 0.45 > 0.3. MA5's inner suppression is MOSTLY a consequence of the phi-law (N15 + N9b), not fully.
    The MA5 gate is not redundant yet; it is reduced to "whatever sets the 0.45-dex residual".
(d) Candidates for the residual (not run): Yd fixed at 0.5 (a 0.3-dex Yd spread alone would do it); spherical Phi_b,0 proxy for discs
    (flat-disc central potential is deeper by O(1) factor, galaxy-dependent); bulge/no-bulge; distance errors. A Yd-free rerun is the next
    cheap step; if the residual drops below 0.3 at free Yd (within SPARC's 0.3-0.8 range), the bar is met.
(e) Honest sentence for the note: "The central nb density is the cosmic baseline raised by the halo's own potential (universal) and by the
    baryon well (per galaxy); the latter accounts for most but not all of the observed spread (0.45 dex residual at fixed Yd)."

## 3. STATUS LINES
  N9b-1  rho_c normalisation = law's own truncation (Phi_self ~ 12.4 sigma^2, r_t/r_c ~ 500 ~ Mpc scale) (T2/T3)
  N9b-2  baryon-well factor tracks fitted rho_c: slope 0.82, r 0.78, residual 0.45 dex (from 0.73) (T2)
  N9b-3  bar not met (0.45 > 0.3); MA5 mostly consequence, not fully; Yd-free rerun queued (T3)

## 4. ADDENDUM 22:55 -- Yd-FREE RERUN (N14b: n14b_yd_free_20261017.py/.csv; n9b_ydfree_20261017.csv; JM 22:44 'go')
 Refit (log s2, log rho_c, Yd in 0.3-0.8), Yb 0.7, from the N14 free point; 120 galaxies, Pool(2), <1 min (scipy now in the VM).
 Yd: median 0.38; 51 at 0.3, 27 at 0.8, 42 interior (disc-halo degeneracy). chi improves 21% median. Raw spread log rho_c 0.69.
 fit vs canon pred (Yd free): slope 0.77, r 0.73, residual 0.47 dex (Q=1: 0.47; interior-Yd only: 0.44). Offset -5.47 (16-84 -5.97..-4.92).
 => Yd is NOT the residual. Remaining suspects: spherical Phi_b,0 proxy (thin exponential disc centre: 2 pi G Sigma_0 Rd, galaxy-dependent
    O(1) factor vs the spherical integral) -- next and last cheap check; bulge treatment; distance. Or the residual is real: a second
    galaxy variable the law does not contain. Bar (<= 0.3) still not met; MA5 ~60% consequence of the phi-law.
  N9b-4  Yd-free: residual 0.47 (unchanged); Yd degenerate with core (78/120 at bounds); not the cause (T2)
  N9b-5  exact thin-disc Phi_b,0 queued; if residual stays ~0.45, the law needs a second variable (T3)

## 5. ADDENDUM 23:15 -- EXACT THIN-DISC Phi_b,0 (N9c: n9c_thindisc_phib_20261017.py/.csv; JM 23:08)
 Phi_disc(0) = 2 pi G Yd int SBdisk dR from the rotmod SB profiles (inner constant, outer exp tail: share ~0). Bulge/gas on spherical proxy.
 Exact/spherical disc ratio median 1.06 (16-84 0.90..1.14); disc share of Phi_b,0 median 0.66.
 fit vs canon pred: slope 0.81, r 0.78, residual 0.45 dex -- IDENTICAL to sec 1. Offset -5.38, std 0.47.
 => Proxy cleared. Residual 0.45 dex is intrinsic at this level: survives Yd free, exact disc geometry, Q=1 selection.
 VERDICT (firm): rho_c = (universal self-potential factor) x (baryon-well factor, ~60% of the spread) x (0.45 dex unexplained).
 FRM-2 should be reframed from "gate vs phi-law" to "what second galaxy variable sets the residual". Candidates (not cheap, not run):
 central baryon acceleration g_b,0 (= MOND's answer; would restore FRM-2's motivation if it correlates); SPARC distance factor f_D; HI/stellar
 concentration. A correlation scan of the residual against these is the next step.
  N9b-6  exact thin-disc potential: residual 0.45 unchanged; proxy cleared (T2)
  N9b-7  residual intrinsic; second variable needed; scan residual vs g_b,0 / f_D / concentration queued (T3)
