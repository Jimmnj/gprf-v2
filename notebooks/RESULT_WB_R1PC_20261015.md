# RESULT_WB_R1PC_20261015 -- conserving floor-clip (R2) wide-binary gamma at R = 1 pc  [UNRATIFIED]
2026-10-15 11:00-11:35 EDT, Claude (Fable 5.1), cloud sandbox (2 vCPU). Trigger: JM 10:56 "run wb r = 1pc" (queue item 2, HANDOFF_20261015).
Script: notebooks\wb_conserve_r1pc_20261015.py = wb_conserve_20261011.py with RS = (1,2,4,8,16,24,32) and a timing line; physics unchanged.
Builds on RESULT_WB_CLIP_GAMMA_20261011.md sec 5c (R <= 16, "still falling"). Purpose: pin the R2 headline number.

## 1. Setup (as 10/11 sec 5c)
QUMOND, RAR nu, g_N,e = 1.2 a0, equal-mass binary, s in r_M (7,030 AU per r_M at 1 Msun). Floor rho_host = 0.01 Msun/pc^3.
R2 = deficit inside sphere R about the COM removed from the positive phantom in the same sphere, uniformly in proportion (f_R = D_R/P_R).
Outer box half-width W = 34 r_M (holds the R = 32 sphere), N = 408 -> h = 0.167 r_M (same h as the 10/11 s = 4 run).
Fine box: s = 4 -> b = 1.6 about each star, Nf = 128; s <= 2 -> single box bb = s/2 + 1.2, Nf = 192. 3-pt GL in cos(theta): 27.5 / 60 / 83.5 deg, weights 5/8/5 (/18).
Runtime 63-66 s per (s, theta); one W = 50, N = 600 run (248 s) for R = 40, 48 at s = 4, theta = 60.

## 2. Validation against 10/11
- s = 4 theta = 60: W = 16 vs W = 34 give R16 = 1.2558 vs 1.2559 (box-independent at matched h).
- Iso MOND 1.281 / 1.354 / 1.361 at s = 1/2/4: identical to 10/11.
- Iso R4 (s = 1) 1.190 = 10/11; iso R8 (s = 2) 1.248 = 10/11; iso R16 (s = 4) 1.270 = 10/11.
- CORRECTION to 10/11 sec 5c: the s = 1 R8/R16 (1.167/1.166) and s = 2 R16 (1.240) entries were BOX-TRUNCATED (W = 6 and 8 < R;
  the sphere extended past the grid). Correct full-sphere values: s = 1 R8 1.161, R16 1.141; s = 2 R16 1.217. The 10/11 s = 4 column
  (W = 16 = R) is unaffected.

## 3. RESULT TABLE (iso gamma_r; per-angle values in sec 6)
 s/r_M  kAU(1Msun) | MOND  | R=16  R=24  R=32  | boost cut at R32 | aniso 0/90 at R32 (MOND)
 1.0     7.0       | 1.281 | 1.141 1.134 1.133 |  47%             | 0.973 (0.960)
 2.0    14.1       | 1.354 | 1.217 1.211 1.210 |  41%             | 1.001 (1.049)
 4.0    28.1       | 1.361 | 1.270 1.257 1.255 |  29%             | 1.006 (1.084)
 (boost cut = 1 - (gamma_R2 - 1)/(gamma_MOND - 1). f_R = 0.57-0.59 at R = 32 for all s.)
 R-convergence: R16 -> R24 -> R32 increments 0.013 -> 0.002 (s = 4); 0.006 -> 0.001 (s = 2); 0.007 -> 0.001 (s = 1). CONVERGED at R = 32.

## 4. Beyond R = 32: the rule stops being physical
 s = 4, theta = 60, W = 50: R16 1.2559, R24 1.2431, R32 1.2415, R40 1.2476, R48 1.2568; f = 0.56, 0.58, 0.58, 0.57, 0.55.
 Mdef (total deficit in box) 0.551 (W = 34) vs 0.554 (W = 50): the deficit integral is SATURATED by ~32 r_M (~1.1 pc at 1 Msun),
 consistent with |rho_ph| = rho_host at ~1.4 pc (10/11). Past that, P_R keeps growing (far positive lobes) while D_R does not, so the
 uniform-proportional rule removes less from the near cocoon and gamma drifts back up. That drift is an artefact of the rule, not transport.
 => The physical R for R2 is the deficit's own extent, R ~ 32 r_M. The minimum gamma is the R2 number. No free length remains in R2.

## 5. READING [Claude, UNRATIFIED]
 (a) R2 HEADLINE: gamma = 1.13 / 1.21 / 1.26 at 7 / 14 / 28 kAU (1 Msun). MOND: 1.28 / 1.35 / 1.36. R0: 1.00.
 (b) Shape: R2 keeps RISING with separation where MOND saturates (1.21 -> 1.26 vs 1.35 -> 1.36). Same discriminator as 10/11, now below MOND.
 (c) Anisotropy erased at 14-28 kAU (1.00-1.01 vs MOND 1.05-1.08). At 7 kAU R2 keeps MOND's mild perpendicular excess (0.97 vs 0.96).
 (d) vs data (refs VERIFY): Chae 2023 1.43 +/- 0.06 at g_N <~ 1e-10 (s >~ 1.1 r_M): R2 1.21-1.26 is 3-4 sigma LOW. Banik 2024 alpha_grav
     = -0.02 +0.07/-0.05 on a Newton-QUMOND line: 1.21-1.26 sits at alpha ~ 0.6-0.7 on that line, excluded by their posterior as stated.
     The 10/11 sec 5c(c) verdict stands: R2 is in the middle and fails both camps as they currently stand. R0 (1.00) is on Banik's number.
 (e) The number is now PINNED (sec 3-4). C13 can be ruled on it. Recommendation unchanged from 00:04: R2 as prediction, R0 as stated
     fallback -- OR, given (d), R0 as prediction and R2 as the "if a 20% residual survives DR4" statement. JM's call.
 (f) Mass dependence: all numbers scale with r_M; at 1.5 Msun the three rows are 8.6 / 17 / 34 kAU.

## 6. Per-angle values (theta = 27.5 / 60 / 83.5 deg)
 s = 4: MOND 1.4306/1.3436/1.3201; R16 1.2858/1.2559/1.2777; R24 1.2715/1.2431/1.2640; R32 1.2695/1.2415/1.2623; Mdef 0.524/0.551/0.571
 s = 2: MOND 1.3905/1.3494/1.3261; R16 1.2223/1.2113/1.2204; R24 1.2110/1.2009/1.2099; R32 1.2095/1.1996/1.2088; Mdef 0.531/0.556/0.578
 s = 1: MOND 1.2519/1.2853/1.3035; R16 1.1235/1.1420/1.1554; R24 1.1172/1.1350/1.1479; R32 1.1163/1.1341/1.1471; Mdef 0.531/0.550/0.566
 (R1, R2 columns at s >= 2 are ill-posed, R < s/2 + 1; ignore.)

## 7. CAVEATS (carried from 10/11 sec 6, unchanged)
 1. QUMOND one-shot; AQUAL feedback of the removed/filled mass not solved (expect O(10-30%), queue item 3).
 2. Uniform-proportional removal is one rule; removal weighted toward the collar would lower gamma further (10/11 5c(d)).
 3. Equal masses, RAR nu only, isotropic 3D weighting; not Chae's projected statistic. Mock-Gaia pass owed before a data comparison is quoted.
 4. Single-star self-acceleration (10/11 sec 4) still unchecked against data.

## 8. NEXT
 1. Rule C13 (R0 vs R2) with sec 3 in hand. 2. AQUAL single-star check (queue 3). 3. Note draft uses sec 5(a) as the R2 line.
