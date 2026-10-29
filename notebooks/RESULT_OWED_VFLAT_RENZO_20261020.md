# RESULT_OWED_VFLAT_RENZO_20261020 [UNRATIFIED] -- two owed checks (cold review A9; debate rows 3-4). Script owed_vflat_renzo_20261020.py; report .txt. From the N12 per-point table.
(1) BTFR scatter in model V_flat (outer-third mean) at fixed M_b, N=120: data 0.059 dex (incl. errors) | GPRF fitted rho_c 0.044 | GPRF universal 0.029 | MOND RAR fn 0.026.
    PASS (bar data + 0.02). Slopes 0.26-0.27 all. The 0.65 dex rho_c spread adds ~0.03 dex to V_flat. Outline 3.2 may quote 0.044 vs 0.059.
(2) Renzo per-feature (features = log V minus per-galaxy quadratic in log r; 2789 points, galaxies >= 8 pts): Pearson r model-vs-obs: GPRF fitted 0.329 | universal 0.311 |
    MOND RAR fn 0.241 | baryons-only 0.272. Slope(obs on model) 0.40 / 0.38 / 0.30 / 0.24. nb-only feature rms 0.016 vs MOND phantom part 0.010 (obs features 0.022).
    READING: the potential response carries baryon features at least as well as MOND's acceleration response (r 0.33 vs 0.24). Absolute correlations are low
    (features ~0.02 dex, errors ~0.04): label stays 'potential explanation', now with a number. Owed before any 'addressed': a per-galaxy feature test on the
    canonical Renzo cases (NGC 1560, NGC 6946, IC 2574) with the dip/bump amplitude compared.

## ADDENDUM 18:55 -- canonical cases (renzo_cases_20261020.py / report). NGC 1560 NOT in SPARC; NGC 7793 outside the TEST1 sample.
IC 2574 (Sancisi 2004): dip at 4.3 kpc obs -30 mdex; MOND -30; GPRF fitted -12 (40%, FAILS the declared 50% bar); universal 0; r(fit,obs) -0.24. NGC 2403: r(fit) -0.39 vs MOND +0.23.
NGC 6946 / NGC 3198: inner structure, all models r ~0.8-0.97 (non-discriminating). DDO 154: outer dip obs -31, GPRF -13, MOND -16; r 0.93 vs 0.83. NGC 2841: r 0.83 vs 0.71.
VERDICT: 'addressed' NOT earned; 'potential explanation' stays with IC 2574 named. Physics: response to the potential (density integrated twice) smooths sharp
features; MOND responds to acceleration (once). The pooled r = 0.33 was carried by smooth-inner-structure galaxies. Method caveat: quadratic-in-log-r detrend
is crude; a smooth-model-residual test is the referee-grade version; it will not flip IC 2574. NOTE: do not claim Renzo's rule; say 'features in the baryon
potential are carried into the response; the sharpest density features are partly smoothed (IC 2574: 40% of the observed dip)'.
