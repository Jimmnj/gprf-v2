# RESULT -- N9: SATELLITE dSph DISPERSIONS FROM THE CANON LAW (And XIX, Crater II)
# gprf\notebooks\RESULT_N9_DWARFS_20261003.md
# 2026-10-03 (Fable session). Status: [UNRATIFIED]. Script: n9_dwarfs_20261003.py (device sandbox, ~60 s).

## 0. VERDICT (one line)
Ambient truncation does almost nothing to satellite dispersions (3-15% at And XIX, 4% at Crater II): R_edge is
2-4 r_half but the halo mass inside r_half is barely touched. The lever is rho_c, which the law does not set.
With the SPARC-median core (log rho_c ~ 7.8) GPRF gives And XIX sigma_los ~ 8.3 km/s vs observed 7.8 +1.7/-1.5
(Collins+2020): a PASS where isolated MOND (6.0) is 1.2 sigma low and MOND-with-EFE (~2) fails. Crater II (obs
2.7 +- 0.3) needs log rho_c ~ 6.0, 60x below the SPARC median: GPRF can FIT it, not predict it. The forecast
"GPRF ~ isolated MOND for dwarfs" is replaced by "GPRF dwarf sigma is set by rho_c; truncation is negligible."

## 1. SETUP
Law: rho_nb = rho_amb[exp(phi_loc/s2) - 1], floor at -rho_amb. Baryons Plummer (a = r_half/1.305). s2 = BTFR
input. rho_amb = host v_c^2/(4 pi G D^2) + rho_bg (M31 v_c 230 at D = 115 kpc -> 7.4e4; MW 220 at 117 -> 6.6e4).
Tracers follow the light; isotropic Jeans in the total potential; sigma_los luminosity-weighted inside r_half.
Scan log rho_c = 5.5..7.5 (SPARC range 7.3-9.0, median 7.8, N14). Inputs: And XIX M_b 1.6e6 (M_V ~ -10.1,
M/L* 1.6; CHECK), r_half 3.1 kpc; Crater II M_b 2.6e5 (L_V 1.6e5), r_half 1.07 kpc, D 117 kpc (Caldwell+2017).
MOND comparators are MY estimates: isolated sigma^4 = 4 G a0 M/81 (Milgrom 1984); EFE quasi-Newtonian with
G_eff = G a0/g_ext, g_ext = v_c^2/D, Wolf estimator. Not McGaugh-Milgrom 2013 values; verify before citing.

## 2. NUMBERS (sigma_los km/s; iso = rho_amb 35.5, trunc = host ambient)
And XIX (law sigma 8.9, R_edge 6.3 kpc = 2.0 r_half; MOND iso 5.95, MOND EFE 2.1, stars-only 0.74):
  log rho_c  5.5: 4.03 / 3.46 | 6.0: 5.76 / 5.41 | 6.5: 7.32 / 7.06 | 7.0: 8.29 / 8.04 | 7.5: 8.57 / 8.32
  M_nb(<r_half)/M_b at 7.0: 90 / 86  -> M/L_half ~ 280 (observed 278 +146/-198). Consistent.
  Low-mass/far variant (M_b 0.8e6, D 130): 7.0 -> 7.10 / 6.87. Still inside the observed range.
Crater II (law sigma 5.7, R_edge 4.3 kpc = 4.0 r_half; MOND iso 3.8, MOND EFE 1.5, stars-only 0.51):
  log rho_c  5.5: 1.61 / 1.44 | 6.0: 2.50 / 2.41 | 6.5: 3.59 / 3.53 | 7.0: 4.60 / 4.55 | 7.5: 5.24 / 5.20
  Observed 2.7 needs log rho_c ~ 6.1. SPARC-median core gives ~5.2: 2x high.

## 3. READING
(a) The LIT sec 8b forecast ("And XIX truncation SUBSTANTIAL, R_edge/r_half 1.5-2") was right about the ratio and
    wrong about the consequence: the Jeans sigma inside r_half sees the enclosed mass, and the floor removes mass
    only outside ~R_edge/2. 5-15% is the whole truncation effect. The ambient mechanism does NOT produce an
    EFE-like lowering anywhere in the Local Group. Said in advance (LIT 8/8a); now a number.
(b) The 2020 And XIX dispersion (7.8) is ABOVE isolated MOND and far above MOND-EFE. McGaugh's 2018 reading
    ("MOND via EFE") used the 2013 value 4.7. With the new value And XIX is a case where the EFE hurts MOND and
    GPRF with a normal core lands on the data. That is the first satellite where GPRF is ahead. One object; the
    error bar is wide (+1.7/-1.5); M_b is uncertain by ~2x (M/L_half 278 and the tidal elongation both argue the
    light may not trace the mass). Do not oversell.
(c) Crater II stays MOND's. GPRF reaches 2.7 only with a core 60x below the SPARC median. Nothing in the law
    prefers that value. MOND-EFE gives ~1.5-2 (McGaugh's prediction was ~2.1 before the 2.7 measurement; the
    measured value is between MOND-EFE and MOND-isolated). GPRF's honest statement: Crater II is a fit, not a
    prediction, and its low core is unexplained. Possible physics: tidal stripping of the medium by the MW
    (real-mass tides, LCDM's story) or a real rho_c-environment dependence (F23 direction: cores track the baryon
    bump; Crater II's bump is the shallowest in the sample).
(d) The dSph problem in GPRF is therefore the SAME problem as N14/F23: what sets rho_c. If rho_c tracks the
    central baryon potential (F23 slope), then Crater II (lowest Sigma_*) should have the lowest rho_c, which is
    the right direction. Test: compute rho_c(Phi_b,0) from the SPARC N14 fits and extrapolate to Crater II's
    central potential; see whether 10^6 comes out. [N9b]

## 4. FLAGS
F24  dSph sigma is rho_c-limited, not truncation-limited; the "ambient truncation = GPRF's EFE" line in LIT
     sec 4c/8 must be scoped to R_edge < r_half/2 objects, of which there are none in the Local Group. (T2)
F25  And XIX M_b, D(M31) and the MOND comparators need literature values (Martin+2016 photometry; McGaugh &
     Milgrom 2013b prediction; Collins+2020 Sec. 5). (T3)
N9b  rho_c vs central baryon potential from the N14 fits, extrapolated to Crater II. (queued)

## 5. RATIFICATION
  N9a  truncation negligible for Local Group satellites (3-15%) (T2)
  N9b  And XIX: GPRF with SPARC-median rho_c ~ 8.3 vs 7.8 +1.7/-1.5: pass; isolated MOND 6.0; MOND-EFE ~2 (T3)
  N9c  Crater II: GPRF needs log rho_c ~ 6.1 (60x below median): fit not prediction; MOND-EFE closer (T2)
  N9d  dwarf problem = rho_c problem = N14/F23 (ruling)
# END -- [UNRATIFIED]
