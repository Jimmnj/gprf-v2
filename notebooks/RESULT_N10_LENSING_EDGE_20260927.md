# RESULT -- N10: LENSING FLATNESS vs GPRF R_edge, THE C4 CLOSURE, AND THE MASS BUDGET
# gprf\notebooks\RESULT_N10_LENSING_EDGE_20260927.md
# 2026-09-27 (Fable session). Status: [UNRATIFIED]. Script: n10_edge_bounds_20260927.py.
# Run under Claude's recommended ballot rulings (E5 kill; C4 restate; rest keep). JM has NOT ruled.

## 0. VERDICT (one line)
"Flat to 1 Mpc" (Mistele+2024) bounds the isolated ambient density to rho_amb <~ 150 Msun/kpc^3
(~1.1 rho_crit), 7x tighter than the 0c(f) bound; the C4 closure value passes; GPRF predicts a
mass-ordered decline onset at 2-4 Mpc, just beyond current data; and at these ambient densities
R_edge spheres OVERLAP, so the isolated-halo picture (and its depletion trough) does not extend to
the population -- a continuous-field calculation is required, not a kill.

## 1. DATA USED
Mistele, McGaugh, Lelli, Schombert, Li 2024 (arXiv:2406.09685), KiDS-1000 weak lensing,
deprojected circular velocities around ISOLATED galaxies (no neighbor with >10% stellar mass
within 4 Mpc/h70). Four baryonic mass bins, Table 2 (50-300 kpc): log Mb = 10.10, 10.66, 10.96,
11.29; Vflat = 137.3+/-11.5, 182.1+/-10.8, 204.7+/-8.2, 260.2+/-6.2 km/s. "No clear hint of a
decline out to 1 Mpc"; systematics dominate beyond. The "downturn in the last 3-4 points"
discussed on tritonstation is the low-a_b tail of Brouwer+2021's RAR (Mistele+2023 Fig 7: LTG
points at ~0.7-0.8 x MOND below a_b ~ 1e-13), flagged unreliable by the isolation criterion.
=> N10 as a FIT to a downturn is premature. What the data give is a hard lower bound R_edge > 1 Mpc.

## 2. METHOD
Exterior solved as the self-gravitating isothermal-with-floor (canon 0c law, baryons negligible
at Mpc): psi'' + 2psi'/r = 4 pi G [rho_c exp(-psi/s2) - rho_amb], s2 = Vflat^2/2, core r0 = 5 kpc
(sets rho_c only; v_flat is independent of rho_c). R_edge = radius where excess M(<r) peaks.
Decline onset = radius where v^2 falls to 0.9 (and 0.7) of the flat value.
Constants: H0 = 70, rho_crit = 136.0 Msun/kpc^3, Omega_m = 0.315 -> rho_m = 42.8.

## 3. RESULTS
(a) Analytic bound R_edge = sqrt(s2/(2 pi G rho_amb)) > 1 Mpc:
      logMb 10.10: rho_amb < 349 (2.56 rho_crit) | 10.66: < 614 | 10.96: < 775 | 11.29: < 1253
    Numerical 10%-decline bound (R(0.9 vf^2) > 1 Mpc), lowest bin: rho_amb <~ 150.
    => rho_amb(isolated) <~ 150 Msun/kpc^3 ~ 1.1 rho_crit ~ 3.5 rho_m. 0c(f)'s 1e3 is SUPERSEDED.
(b) C4 closure: rho_bg = a0^2/(G c^2) = 35.5 Msun/kpc^3 = 0.26 rho_crit = 0.83 rho_m.
    PASSES (a) in every bin (one-sided). NOTE THE NUMBER: the closure density is the cosmic MEAN
    MATTER DENSITY to ~20%, i.e. a0 = c sqrt(G rho_m) (H0=70, Omega_m=0.315 gives 1.3e-10).
    0c(f) quoted "0.2 rho_crit"; the cleaner statement is rho_bg = rho_m. If the medium's mean
    density IS the cosmic matter density, that is the natural reading, not a coincidence
    between unrelated numbers. Still imposed, not derived. [T3]
(c) Decline-onset prediction (10% drop) at rho_amb = closure:  2.2 / 2.9 / 3.3 / 4.1 Mpc for the
    four bins; R_edge = 3.7 / 4.9 / 5.5 / 7.0 Mpc. At rho_amb = 100: 1.3 / 1.8 / 2.0 / 2.5 Mpc.
    At rho_amb = 300: 0.7 / 1.0 / 1.2 / 1.5 (excluded for bins 1-2 by (a)).
    => PREDICTION: no decline inside 1 Mpc for isolated galaxies; a ~10% decline between ~1.5
       and ~4 Mpc, ONSET RADIUS INCREASING WITH MASS as Mb^(1/4) (sigma ~ Mb^(1/4), R_edge ~ sigma).
       MOND: no decline ever (isolated). LCDM: 1-halo -> 2-halo transition, not mass-ordered
       this way. Euclid/Rubin target. Mistele's 3 Mpc bins sit exactly here but with large
       systematics.
(d) M_nb(<R_edge)/M_b at closure: 873 / 555 / 394 / 376. (The data themselves imply ~300 at 1 Mpc
    for bin 1; this is not GPRF-specific, every theory owes it.)
(e) MASS BUDGET (0922 sec 5 item 5, first numbers). Sum over the GSMF (Baldry+2012 double
    Schechter, Mb = 1.3 M*) of M_nb(<R_edge) with BTFR sigma:
      rho_amb =  35.5: 159 Msun/kpc^3 = 3.7 x rho_m = 4.5 x rho_amb
      rho_amb = 100  :  95            = 2.2 x rho_m = 0.9 x rho_amb
      rho_amb = 300  :  55            = 1.3 x rho_m = 0.2 x rho_amb
      rho_amb = 1000 :  30            = 0.7 x rho_m
    At the closure value the halos would need 4.5x more medium than exists, and the floor
    (excess >= -rho_amb) means troughs cannot compensate. BUT: R_edge = 3-7 Mpc exceeds the mean
    galaxy separation (~2-5 Mpc), so R_edge spheres overlap and the sum double-counts the same
    medium. The isolated-halo picture does not extend to R_edge for the population. This is NOT
    a kill; it says the Mpc-scale medium is a continuous field (the regime LCDM handles with the
    2-halo term) and the "depletion trough" is never realized around real galaxies (0922 sec 1
    already said neighbors fill it). Required calculation: GPRF galaxy-matter correlation
    xi_gm(r) from the medium's response to the galaxy field, compared to KiDS at 1-10 Mpc.
    The GSMF integral is also dominated by dwarfs, many of which are satellites truncated by
    their host's rho_env; the sum is an overestimate. Crude; direction is what matters.

## 4. WHAT THIS DOES TO THE BALLOT / FLAGS
- C4 (restated): closure survives N10 one-sided; the closure density = rho_m is a cleaner
  statement than 0.2 rho_crit. Two-sided test still needs rho_env(isolated) independently.
- C5: 0c(c) demo numbers were at rho_amb = 3e4; now excluded by ~200x, not 30x.
- NEW F17: 0c(f) "rho_amb(isolated) <~ 1e3" superseded by <~ 150 (this file).
- NEW F18: mass budget forces a field treatment at Mpc scale; single-halo R_edge and the trough
  are population-inconsistent at rho_amb <~ 150. The v2.0.0 "cutoff/edge" language must be
  scoped to "isolated galaxies inside ~R_edge/2" and the trough dropped from predictions.
- NEW N11: GPRF xi_gm(r) at 1-10 Mpc (field calculation), the real Mpc-scale test.
- Prediction for the paper (c): mass-ordered decline onset at 1.5-4 Mpc. Zero free parameters
  once rho_amb is fixed; one if rho_env is fit.

## 5. RATIFICATION
  N10a  rho_amb(isolated) <~ 150 Msun/kpc^3 from flatness to 1 Mpc (T2)
  N10b  closure rho_bg = a0^2/(G c^2) = 0.83 rho_m passes; restate closure as a0 = c sqrt(G rho_m) (T3)
  N10c  mass-ordered decline onset 1.5-4 Mpc as a v2.0.0 prediction (T3)
  N10d  mass budget: isolated-halo picture population-inconsistent beyond ~R_edge/2; field
        treatment required; trough dropped (T3)
  N10e  N11 xi_gm calculation added to the queue (T3)
# END -- [UNRATIFIED]
