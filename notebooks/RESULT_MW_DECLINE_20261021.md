# RESULT_MW_DECLINE_20261021 [UNRATIFIED]
Pre-registered 2026-10-21 16:20 (JM 'run'); script mw_decline_20261021.py; report mw_decline_20261021_report.txt; csv mw_decline_20261021.csv.
Trigger: Triton Station comment thread 2026-10-18 (McGaugh): "The RC of the MW is declining and it should be declining in MOND because that is what the mass distribution says it should do."

## Inputs (verified from arXiv:2303.12838 PDF, fetched 2026-10-21)
Ou, Eilers, Necib & Frebel 2024 (MNRAS 528, 693): Table 1 circular velocities, 37 bins R=6.27-27.31 kpc (234.1 at 7.86 -> 173.0 at 27.31; they note systematics >15% beyond 22 kpc from neglected asymmetric drift). Table 2 baryons (de Salas+2019; Misiriotis+2006): Hernquist bulge 1.55e10 Msun a=0.70 kpc; double-exponential stellar disc 3.65e10, Rd 2.35, hz 0.14; HI 8.2e9 Rd 18.24; H2 1.3e9 Rd 2.57. Total M_b 6.15e10 (gas He factor as in their table; not altered).

## Method
Baryon v^2(R): Hernquist + Freeman thin-disc formula for each disc component. nb response: TEST1 spherical solver, tanh2 cap, BTFR sigma = 125.1 km/s (v_flat 176.9), rho_amb 35.5, universal rho_c 10^7.45 and fitted. MOND comparison: spherical RAR nu on the same baryons. Error model: table errors in quadrature with a 5 km/s floor. Declared pass: law(universal) rms <= 2x MOND rms AND both slopes (10-25 kpc) negative and within 50% of the data slope.

## Result
| model | rms 8-27 kpc (km/s) | slope 10-25 (km/s/kpc) | v(8.2) | v(27.3) | chi2 (N=37) |
|---|---|---|---|---|---|
| Newton (baryons only) | 86.7 | -4.11 | 173.0 | 95.9 | 8441 |
| MOND spherical | 17.2 | -1.67 | 218.0 | 188.0 | 354 |
| law, rho_c = 10^7.45 | 19.7 | -1.88 | 215.0 | 181.6 | 478 |
| law, rho_c fitted 10^7.88 | 8.6 | -3.03 | 237.6 | 184.6 | 70 |
| DATA | - | -2.50 | 232.8 | 173.0 | - |
PASS on both declared criteria. Law (universal) and MOND within 3 km/s of each other at every radius; both decline at ~70-75% of the observed slope; both ~15 km/s low at the Sun under this baryon model (known: MOND fits of the MW use a heavier disc, McGaugh 2019 -- not applied here to keep the comparison at fixed baryons). Fitted rho_c is +0.43 dex above the SPARC median, inside the 0.65 dex spread.

## Beyond the data
Law curves pass below the BTFR speed at 40-100 kpc (universal: 169 at 100 kpc; fitted: 160) before converging -- the isothermal sphere's oscillatory approach to sqrt(2) sigma. Prediction: MW v_c beyond 30 kpc runs 5-10% below 177 km/s. MOND (RAR nu) converges monotonically from above (182 at 100 kpc). A near-Keplerian decline beyond 20 kpc (Jiao et al. 2023) would be a problem for MOND and the law alike.

## Caveats
Spherical nb response (as in the SPARC pipeline); thin-disc baryons (hz ignored, <2% at R>6 kpc); Ou's Jeans curve not the only MW curve (Jiao+2023 steeper; outer tracers differ); no bar/spiral features. Not a MOND-vs-law discriminator: the two agree to within the baryon-model uncertainty.

## For the note
Section 4 (or App. A) line: the law predicts the declining Gaia DR3 Milky Way curve from the Ou et al. baryons at the universal rho_c to 20 km/s rms, the same as MOND on the same baryons; the Sun sits at 1.4 MOND radii, so the decline is the transition, not a Keplerian fall. Label: addressed by construction (shared with MOND).

## Addendum 16:45 -- is the MW a one-off? SPARC massive galaxies vs R/r_MOND (JM 16:34/16:38)
From n12_emergent_nu_20261020.csv (universal rho_c 10^7.45), galaxies with M_b > 3e10 (51; 1545 points), r_M = sqrt(G M_b/a0). Median log(g_obs/g_model):
0-1 r_M: +0.016 (N 758, 51 gals); 1-2: -0.015 (337, 51); 2-3: -0.008 (225, 45); 3-4: +0.010 (124, 32); 4-6: +0.012 (85, 23); 6-10: -0.052 (16, 4: NGC 289, UGC 9133, NGC 3198, UGC 6614).
Reading: to 6 r_M the law tracks massive discs to ~1% in g with no trend; this is the band where the MW law curve is already near v_flat (182->172 km/s at 3-6 r_M). The predicted dip (minimum ~10 r_M, 5-10% in v = 0.04-0.09 dex in g) is untested: 16 points beyond 6 r_M sit 0.05 dex below the law (dip direction) but too few to call. NGC 2841's outer excess (Fig. 3) is individual, not a class property. Jiao+2023's ~150 km/s at 25-27 kpc is below the law's minimum (169) -- not excluded given their errors, not 'compatible'. Mistele+2024 flat lensing at 1 Mpc is the same curve's tail: declining at tens of kpc + flat at Mpc is a shape the law predicts and RAR-nu MOND does not (monotone from above).
