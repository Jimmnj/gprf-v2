# RESULT_N8D_CHAE_METHOD_20261021 [UNRATIFIED] -- Chae+2021's statistic applied to the observed curves and to the frozen law's own curves. Pre-registered (script header).
Script n8d_chae_method_20261021.py; csv; report. Model: Chae+2021 eq. (nueN) with nu_eN(y) = 1/2 + sqrt(D^2/4 + D/y) - etilde C/y; nuisance Y_disk (0.5 +/- 0.1 dex), Y_bulge (0.7 +/- 0.1),
D (SPARC e_D), inc (SPARC e_inc); SPARC velocity errors; Nelder-Mead from four etilde starts. Selection x_{0,3} < -10.6 from Chae's Table 2 (102 of 119 galaxies).
| curves | median etilde | frac > 0 | N | sign-test p |
| Chae+2021 tabulated fits | +0.025 | 0.61 | 102 | 0.037 |
| observed, our refit | +0.016 | 0.55 | 102 | 0.37 |
| LAW curves, rho_c fitted | +0.022 | 0.60 | 102 | 0.059 |
| LAW curves, rho_c 10^7.45 | +0.012 | 0.57 | 102 | 0.20 |
Our refit vs Chae's table: Spearman +0.877, median difference -0.010 (machinery validated). corr(etilde, log e_N,env max-clust): ours -0.257 (N=60); Chae's table -0.097.
## Reading
1. The law's curves, with NO external field, return positive etilde under Chae's model: half his median at universal rho_c, all of it at fitted rho_c (partly circular -- fitted curves
   track the data). The cored-isothermal outer shape mimics an external field in his parametrisation.
2. The median-etilde detection therefore does not discriminate suppression from this law. His 4.4 sigma (full sample 143) vs our subsample p = 0.04 (102): same direction, weaker here.
3. The discriminating statistic is etilde vs the environmental field: on this subsample weak and NOT of the EFE sign (-0.26, N=60); consistent with Chae's own sec 4.2 (r = 0.057 at
   current errors). His regional contrast (underdense ring ~0 vs overdense ~13 sigma) is the remaining test -- untested here (needs his region assignments).
4. Consequence for the note: EFE-1 'second kill' is NOT currently live against the law via the median statistic; the environment kill must be stated on the e_N dependence
   (slope -0.08 +/- 0.03 or tighter) or the regional contrast. Framing (au') amended in the draft 4.4.
Caveats: our subsample (119 fitted, 102 selected) vs his 162/143; Nelder-Mead vs MCMC; gas ratio fixed at 1.33; priors approximate his.
## ADDENDUM 03:00 -- tercile contrast by environmental field (n8d_tercile_report_20261021.txt), Chae selection with e_N, N=60 (20/20 in the outer terciles):
etilde median low-eN / high-eN / contrast / MWU p: Chae table +0.088/+0.025/-0.064/0.34 | observed refit +0.051/-0.001/-0.052/0.03 | law 7.45 +0.020/-0.032/-0.053/0.02 |
law fitted +0.076/-0.038/-0.114/0.00. obs - law(uni): low +0.025, high 0.000, contrast -0.025, p 0.76. READING: on this subsample etilde is HIGHER in low-field galaxies (opposite
to the EFE), and the law's curves show the same contrast -> the field dependence here is a property of which galaxies sit where, reproduced without an environmental term; no
residual field dependence beyond the law. Chae's regional (cosmic-web) contrast uses different galaxies and region assignments not on file -- still untested.
