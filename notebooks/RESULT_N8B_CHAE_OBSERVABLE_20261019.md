# RESULT_N8B_CHAE_OBSERVABLE_20261019 [UNRATIFIED] -- Chae's observable vs Chae+2021 e_N,env on SPARC. Pre-registered (header of script).
Script: n8b_chae_observable_20261019.py; per-galaxy table n8b_chae_observable_20261019.csv; raw report n8b_report_20261019.txt. Device VM, scipy 1.15.3.
Sample: SPARC Q<=2, inc>30, >=6 points, in Chae+2021 Table 3 -> 88 galaxies; 82 with >=3 points at g_bar < a0/3.
Declared before run: GPRF (frozen law + NEST-1 C) -> coefficient on log e_N = 0; MOND EFE -> negative on both observables.

## Results (coefficient on log e_N,env max-clustering; 95% bootstrap; controls log M_b + log Sigma_b)
(A) outer RAR residual (mean lgo - lgmond at g_bar < a0/3):   -0.120 [-0.279, +0.044] dex/dex   (e_N alone -0.146 [-0.322, +0.018]; no-clust e_N -0.148 [-0.312, +0.006]; Q=1 only -0.076 [-0.299, +0.123])
(B) outer RC log-slope (outer third of points):                -0.029 [-0.239, +0.178]           (Q=1 only -0.044 [-0.278, +0.182])
Circular check (Chae's RC-fitted e-tilde, NOT a test): (A) -0.49 [-0.86, -0.08]; (B) +0.13 [-0.32, +0.58].

## Reading
- (B) NULL, centred on zero: the outer slope does not track e_N.
- (A) NOT a clean null: point estimate has MOND's sign in every variant, CI touches zero (~1.5-1.9 sigma). Over the 5-95% e_N range (0.58 dex) the implied residual change is ~0.07 dex.
- POWER: CI half-width ~0.16 dex/dex. The MOND-EFE expectation at these e_N (0.002-0.018) and g_bar ~ a0/3 is itself of order -0.1 dex/dex [estimate, not computed here]. So THIS SAMPLE CANNOT SEPARATE 0 (GPRF) FROM MOND'S EXPECTED SLOPE at 95%. N8b is inconclusive, not a null.
- Chae+2020's detection came from per-galaxy EFE fits (e-tilde); our circular row reproduces that e-tilde correlates with the residual, as it must. The environmental (e_N) test is the non-circular one, and it is underpowered on SPARC.

## Consequences for the note (framing (ag') third item)
- MAY SAY: "On SPARC the low-acceleration RAR residual shows no significant dependence on the environmental field (slope -0.12 +/- 0.16 dex/dex); the point estimate has MOND's sign and the test does not yet separate the two predictions."
- MAY NOT SAY: "our own test found no signal" (as if a null were established) or "GPRF passes".
- EFE-1 kill rule unchanged: a confirmed e_N correlation kills the environment sector. Current status: open, underpowered.
- Side result (MISTELE-CHAE check): Chae Table 3 e_N (max clust) spans 0.0019-0.018, median 0.005, 5th pct 0.0022. McGaugh's 'few % a0 almost everywhere' (0.03) exceeds Chae's MAXIMUM. Mistele's per-bin bounds (0.005-0.018) are consistent with Chae's table. The squeeze is on McGaugh's Gallardo argument, not on Chae: either e_N ~0.03 is typical and Chae's e_N scale is low by ~6x, or the Gallardo EFE argument needs a different field.
## Owed: a MOND-EFE forward model of observable (A) at Chae's e_N values (QUMOND or Chae's fitting function) to turn the power statement into a number.

## ADDENDUM 2026-10-20 16:50 -- MOND-EFE FORWARD MODEL (n8c_mond_efe_forward_20261020.py / report / csv; cold review A6; the 'owed' item). Chae+2021 Eq. (qumond1d)
verbatim from RCLSSr1.tex (simple IF nu0), each galaxy at its own e_N (max clust), same 82 galaxies and low-acceleration points as (A). EFE dip median -0.087 dex
(range -0.25..-0.03) at e_N median 0.005. PREDICTED slope of (A) vs log e_N with M_b, Sigma_b controls: MOND -0.080 [-0.127, -0.034]; GPRF 0; DATA -0.120 [-0.283, +0.044].
IF-mismatch leakage 0.000. READING: the point estimate sits AT MOND's value; GPRF's zero is inside the 95% CI near its upper edge; separation 0.08 vs half-width 0.16
= ~0.5 sigma discrimination; rough likelihood odds ~3:1 for the EFE sign, not significant. NOTE WORDING (replaces the MAY-SAY line): 'the sample does not separate
the two predictions (-0.08 for the external-field effect, 0 for the law); the point estimate, -0.12 +/- 0.16, lies at the external-field value.' DO NOT write
'consistent with zero' alone. KILL NUMBER: a slope of -0.08 +/- 0.03 or tighter on a larger sample kills the law's environment sector and favours MOND.
