# RESULT_SIMROUND_C1-C4_20261025 [UNRATIFIED] -- 2026-10-25 00:48 [time label corrected] (run unattended; JM asleep, PC on)
> **JM 2026-10-25 00:37:** "yes - and i'll hit the sack but leave my pc running"

Four 'easy' analyses from the 10/24 simulated review round (reviews simulated; tests and numbers ours; 4c). No note text changed.

## C1 Cross-validated relation (3) -- pass7_cv_relation3_20261025.py, pass7_cv_results.csv
Refit on the full 147-galaxy per-galaxy rho_c (pass2): log rho_c = 10.934 - 0.381 log M_b + 0.374 log Sigma_h - 0.808 f_gas (Sigma_h in Msun/pc^2); in-sample rms 0.139 (published constants: 0.140).
K-fold by galaxy, 20 repeats each: 5-fold 0.1424 +- 0.0015, 10-fold 0.1424 +- 0.0010. Out-of-sample cost 0.003 dex. (universal 0.155; RAR fn 0.133.)
READING: relation (3) is not an in-sample artefact; the 'fitted on the same data' objection costs 0.003 dex.

## C2 Galaxy-level bootstrap -- pass8_bootstrap_20261025.py (+ _report.txt), 2000 resamples of whole galaxies
rms differences (A - B), 95% interval: law rho_c free - RAR M/L free +0.002 [-0.010,+0.013] (P(A<B) 0.35: 'equal within the noise' CONFIRMED);
universal - RAR fixed +0.022 [+0.010,+0.035] (worse, significant); relation (3) - RAR fixed +0.007 [-0.001,+0.015] (P 0.04, marginal); relation (3) - universal -0.015 [-0.023,-0.008].
Galaxy-level trends, universal law (N 147): mean residual vs log Sigma_h +0.25 (95% [+0.09,+0.40]) SIGNIFICANT; vs f_gas -0.16 [-0.31,+0.01] marginal; vs log M_b n.s.;
within-galaxy radial slope median -0.055 [-0.092,-0.015] (real, modest; replaces the per-point Spearman -0.17 as the defensible statistic).

## C3 Wide binaries vs the binary's own rho_c -- pass9_wb_rhoc_sweep_20261025.py (+ _report.txt); model of frm3_wb_frozen_law_20261020.py
gamma - 1 is linear in rho_c. log rho_c(binary) needed for gamma-1 = 0.10 / 0.05 / 0.02 / 0.01: at 28 kAU 10.17 / 9.87 / 9.47 / 9.17; 14 kAU 10.67/10.37/9.97/9.67; 7 kAU 11.22/...; 3 kAU 12.32/...
Relation (3) extrapolated to a 2 Msun binary (Sigma_h 30-130 Msun/pc^2, f_gas 0): log rho_c 11.35-11.61 -> gamma(28 kAU) = 2.5-3.7 -- EXCLUDED by present Newtonian data.
GC bound already in the note: log rho_c(GC) <~ 7.2 -> satisfies the gamma < 1.1 condition by ~3 dex.
READING (pre-post item 11): 'wide binaries Newtonian' is CONDITIONAL on rho_c(binary) <~ 10^10.2 (gamma < 1.1 at 28 kAU); compact stellar systems (GCs) are known to sit far below; relation (3) is a disc calibration that must NOT be extended to stellar systems. Conversely a 1%-level Newtonian result bounds rho_c(binary) < 10^9.2.
DRAFT WORDING (Table 1 row + sec 6.1): "Newtonian, provided the binary's own central density is below ~10^10 Msun kpc^-3 (gamma - 1 < 0.1 at 28 kAU) -- the globular-cluster bound (<~10^7.2) clears this by three decades; relation (3), calibrated on discs, does not extend to stellar systems. Conversely, Newtonian data at 1% precision would bound it below 10^9.2."

## C4 Lensing with the amplitude fixed by the law -- pass10_lensing_fixedVf_20261025.py (+ _report.txt)
M_b,eff from Mistele et al. 2024b Table 2 ('All', R<300 kpc; gas included): 10.10/10.66/10.96/11.29. Law V_f = (G a0 M_b)^(1/4) = 119/164/195/236 km/s vs lensing V_flat 137/182/205/260 -> ratios 0.87/0.90/0.95/0.91 (bin 4 -3.9 sigma).
Joint chi2 (41 points, no free amplitude): nominal M_b 73.8, rho_amb 95% upper 158; M_b +0.1 dex 33.4, upper 355; M_b -0.1 dex 144.
READING: with the amplitude fixed, the law under-predicts the lensing flat speeds by 5-13% at the stacked masses; a +0.1 dex mass shift (within stellar-mass/gas systematics) fits well. The amplitude test is currently systematics-limited, and V^4 = G a0 M_b is MOND's asymptote too, so the offset is SHARED with MOND (not a discriminator). The ambient-density bound depends on the amplitude choice: <160-355 instead of <630.
DRAFT WORDING (App C): "With the amplitude fixed by the law rather than fitted (V_f^4 = G a0 M_b, also MOND's asymptote), the lensing flat speeds lie 5-13% above the prediction at the stacked baryonic masses (chi2 = 74 for 41 points); a 0.1 dex increase in M_b, within the stellar-mass and gas systematics, gives chi2 = 33. The ambient density is then bounded by 160-355 (95%)."
