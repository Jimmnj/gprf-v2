# RESULT_LENSING_FIT_SATELLITES_20261021 [UNRATIFIED] -- (A) likelihood fit of the frozen law's truncated tail to Mistele+2024 Table 1; (B) satellite-body correlated term under R7.
Scripts: mistele_fit_20261021.py (+ corrected part B appended to report), satellite_body_term_20261021.py; reports .txt. Trigger: Mistele-persona review (1a, 1b, 1c).

## (A) Fit. Model beyond 100 kpc: rho_excess = sigma^2/(2 pi G r^2) - rho_amb to r_t = sigma/sqrt(2 pi G rho_amb), zero beyond; v_c from the enclosed excess; V_f free per bin
(BTFR normalisation marginalised); errors stat (+) syst in quadrature; 41 points, 'All' sample, 4 mass bins. Fit in v_c space (their deprojected quantity), NOT Delta Sigma -- stated caveat.
| bin | N | V_f(fit) | rho_amb best / 68% / 95% upper (Msun/kpc^3) | chi2 at 35.5 / 150 / 1000 / 3000 (dof) |
| 1 | 11 | 140 | 224 / 501 / 891 | 5.0 / 4.7 / 8.9 / 18.0 (9) |
| 2 | 11 | 190 | 1122 / 1995 / 3548 | 10.7 / 8.8 / 3.7 / 6.4 (9) |
| 3 | 10 | 212 | 3 / 282 / 631 | 7.6 / 7.8 / 15.8 / 42.8 (8) |
| 4 | 9 | 258 | 28 / 447 / 891 | 3.7 / 3.8 / 8.2 / 32.3 (7) |
JOINT: best 316; 68% upper 447; 95% upper 631 Msun/kpc^3; chi2_min 24.3 (33 dof). Closure value 35.5: Delta chi2 = 2.7. REPLACES the hand bound '<~ 150' (N10 10%-decline criterion).
Delta Sigma(truncated)/Delta Sigma(SIS), bin 2: rho_amb 35.5 (r_t 4.2 Mpc): 1.00/1.00/1.00/0.99/0.95/0.71 at R = 0.3/0.5/1/2/3/5 Mpc | 150 (2.0): .../0.99/0.86/0.57/0.34 | 630 (1.0): .../0.84/0.42/0.28/0.17.
## (B) Satellite-body term. Each satellite (m < 0.1 M_host) carries its own halo with sigma_s^2 = sqrt(G a0 m)/2, truncated where its tail meets the host's: r_t,s = sigma_s r/sigma_h,
M_nb,s = 2 sigma_s^3 r/(G sigma_h). Population: dN/dlog m ~ m^-0.3 on [1e7, 0.1 M_host], uniform in a 1 Mpc sphere, N(>1e8) = 10 (plausible) or 30 (upper) [counts ROUGH, mem].
Satellite nb / host nb inside 1 Mpc: N=10: 0.46 / 0.28 / 0.25 / 0.17 for hosts V_f 137/182/205/260; N=30: 1.36 / 0.86 / 0.76 / 0.50. Inside 0.5 Mpc: 0.02-0.06 (N=10), 0.06-0.17 (N=30).
## Reading
1. The compatibility bound is rho_amb < 630 (95%), ~4x looser than the hand value; the closure value is comfortably inside. The truncation is a 1-5 Mpc effect in Delta Sigma.
2. Under R7 satellites contribute 20-45% of the host's enclosed nb at 1 Mpc (N=10) and are negligible inside 0.5 Mpc. The term is positive and mass-ordered opposite to the host's
   (larger fraction for low-mass hosts) -- it masks a truncation and loosens the rho_amb bound further; 'flat to 1 Mpc' is partly satellite bodies in this law.
3. Note wording (sec 5): 'A fit of the truncated profile to the Mistele et al. circular velocities (statistical and systematic errors in quadrature) gives rho_amb < 630 Msun kpc^-3
   at 95% (joint over four mass bins, chi2 24 for 33 degrees of freedom), with the closure value well inside; the satellites of an isolated galaxy, each a body under the source
   rule, add ~20-45% to its enclosed non-baryonic mass at 1 Mpc and little inside 0.5 Mpc, so the stacked profile beyond ~0.5 Mpc is not the host's halo alone.'
Caveats: v_c-space fit, not Delta Sigma; satellite counts and radial distribution are rough; satellite halos' own truncation by each other ignored; host tail taken as SIS beyond 100 kpc.

## ADDENDUM 03:10 -- (C) Delta Sigma-space fit (mistele_dsigma_fit_20261021.py; report). Forward model host tail + satellite term -> Sigma -> Delta Sigma at their radii -> deprojected
with the same beyond-last-bin 1/R assumption their systematic column encodes (Sbar from dSbar/dR = -2 DS/R; Abel inversion with the R = sqrt(r^2+t^2) substitution) -> v_c vs their
values with STATISTICAL errors. Pipeline validated on an exact SIS: recovered 250-261 for 260 (2-4% low at the innermost and outermost radii). Grid 0.2 dex in rho_amb, V_f step 8.
| satellite term | joint best rho_amb | 68% upper | 95% upper | chi2_min (41 pts) |
| none | ~100 | ~160 | ~160-250 | 44.2 |
| N(>1e8)=10 level (f_sat 0.46/0.28/0.25/0.17) | ~400 | ~400 | ~630 | 60.2 |
READING: the Delta Sigma-space bound without satellites (~160-250) is close to the original hand value (150) and tighter than the v_c-space fit (630): the deprojection with a 1/R tail
assumption penalises an inner truncation more than a direct v_c comparison. Adding the satellite term at the R7 estimate WORSENS the fit (Delta chi2 +16) unless the truncation moves
inward: the stacks do not want 20-45% extra mass at 1 Mpc on an isothermal tail; the satellite term and r_t trade off, and the data prefer less satellite mass than the rough count.
Closure value 35.5 inside the bound in both cases. (D) Tercile contrast (n8d_tercile_report): see RESULT_N8D addendum.

## ERRATUM 2026-10-22 16:50 -- TABLE PARSER BUG (found by the simulated Mistele review)
All four lensing scripts (mistele_fit, mistele_dsigma_fit, fig7_truncation, euclid_forecast) parsed Mistele+2024 Table 1 with whitespace split(); the rows at 1836.51 and 2496.47 kpc have blank bin-1 (and bin-2) columns, so their values were assigned to the wrong mass bins. Fixed with a byte-position parser (notebooks/mistele_table.py); pre-fix scripts in _to_delete/*_preparsefix.py.
| fit | before | after |
|---|---|---|
| v_c space, R >= 100 kpc | best 316, 95% < 631, chi2 24.3 | best 355, 95% < 631, chi2 21.5 |
| v_c space, 100 kpc <= R <= 1 Mpc (review request) | -- | best 251, 95% < 708, chi2 19.5 |
| Delta Sigma, no satellites | best 100, 95% < 158, chi2 44.2 | best 398, 95% < 631, chi2 30.2 |
| Delta Sigma, satellites (s=1) | best 398, 95% < 631, chi2 60.2 | best 631, 95% < 1000, chi2 30.5 |
WITHDRAWN: 'Delta Sigma tightens the bound to <~160-250 without satellites'; 'satellite term disfavoured (Delta chi2 +16)' (now +0.3, no preference); 'fit value rho_amb ~150'. REPLACEMENT: fitted rho_amb ~350-400 (v_c and Delta Sigma agree), 95% < 631 (v_c) / < 631-1000 (Delta Sigma), and the data within 1 Mpc bound rho_amb only from above. Preference for truncation over an untruncated tail: v_c chi2 at 35.5 = 28.0 vs 21.5 min (Delta chi2 ~6.5), driven by points beyond 1 Mpc where isolation and deprojection are least secure -- NOT a detection.
Euclid forecast rerun (euclid_forecast_20261022_fixed.txt) at fiducial 355: KiDS-level 355-vs-flat Delta chi2 10.4 (V_f fixed); f = 0.25: 166; closure 35.5-vs-flat still needs 5 Mpc (9.4 at f = 0.25). Fig 7 regenerated with 35.5 / 355 / 630.
Downstream to update: sec 6.3, App C, Fig 7 caption, truncation-acceleration line (g_bar(r_t) at the fit value), cluster estimate case 'rho_amb 150' (rerun at 355), abstract r_t range.
