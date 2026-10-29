# RESULT_FRM3_WB_FROZEN_LAW_20261020 [UNRATIFIED] -- wide binaries under the FROZEN LAW (R1 + NEST-1 C). Script frm3_wb_frozen_law_20261020.py; report .txt.
Trigger: JM 14:39 'CHECK OWED from 10/18' (recompute R2 under the total-g gate). ANSWER TO THE CHECK: R2 (wb_conserve_r1pc_20261015.py) has GE = 1.2: the
Galactic field was already in. R2 = QUMOND phantom density with the EFE, negatives floor-clipped (conserving). It is NOT the frozen law and is withdrawn as a GPRF number.

## Calculation (hand-level, exact enough: enclosed-mass argument)
Binary M = 2 Msun. Sub-system sigma^2 = sqrt(G a0 M)/2 = 0.089 (km/s)^2 (sigma 0.30 km/s). Baryon exponent x_b = G M/(s sigma^2) = 2.8 / 1.4 / 0.7 at 7 / 14 / 28 kAU;
cap radius (x_b = 2) 0.048 pc. Own-population density <= e^2 rho_c. MW population: rho ~0.01 Msun/pc^3, boost exp(phi_b/sigma_host^2) = 1.00000x.
gamma - 1 ~ M_nb(<s)/M:
 s      | own, rho_c 10^7.45 | 10^8.0  | 10^8.5  | MW through-pass | MOND w/ EFE (RESULT_WB_R1PC) | R2 (withdrawn)
 7 kAU  | 1.7e-5             | 6.0e-5  | 1.9e-4  | 8e-7            | 0.28                         | 0.13
 14 kAU | 6.1e-5             | 2.2e-4  | 6.8e-4  | 7e-6            | 0.35                         | 0.21
 28 kAU | 1.9e-4             | 6.7e-4  | 2.1e-3  | 5e-5            | 0.36                         | 0.26
gamma = 1.1 at 28 kAU would need 19 Msun/pc^3 = 10^10.3 Msun/kpc^3 inside the orbit: 2 dex above every fitted rho_c and above the GC bound (TEST5).

## Reading
1. PREDICTION: gamma = 1.000 (+ < 0.003) for disc and halo binaries alike, under host OR sub-system sigma. The cap (TEST1-required, GC-confirmed) bounds the
   density; a sub-pc orbit encloses too little volume for a bounded density to matter. MOND's phantom density is unbounded near a point mass; hence 0.3.
2. CORRECTIONS: (a) ballot 10/18 00:10 'gamma ~ 1.1 from 0.2 Msun in 1 pc' used 1 pc^3; the binary encloses ~1e-3 pc^3 -> WRONG, retracted. (b) FRM-4 PROS/CONS
   'MOND escape lost under sub-system sigma' -> REVERSED by the cap. (c) FRM-4 SCALE 'halo WBs show the signal' -> dropped. (d) C13 R2 numbers 1.12/1.20/1.24 ->
   superseded; R2 was a QUMOND construction. (e) FRM-3 BOUND (tau ~ L/sigma) -> moot: no bump to relax into.
3. KILL RULE: confirmed gamma >= 1.1 (Chae 2023/2024 standing after DR4) KILLS the frozen law. Banik+2024 Newton = PASS. [RETRACTED 15:35, WB-LIT-VERIFY: 'Pflamm-Altenburg 1.05-1.11 mild tension' -- P-A 2025 is theory, no such numbers.]
   Newton results: Cookson+2026 (arXiv:2602.24035), Banik+2024, Pittordis+2025 (2504.07569); detection: Chae 2025 (2502.09373), Hernandez. Opposite-direction discriminator from MOND (1.36).
4. With EFE-1 this gives the note two kill tests with opposite MOND expectations: WB (1 vs 1.36) and environment (0 vs dip).
## RECOMMENDATION: rule FRM-3 = NEWTON (gamma 1.000 +<0.003); un-suspend the WB sentence (framing (bb)); C13 restated: prediction 1.00, kill at >= 1.1, DR4.
