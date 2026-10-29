# RESULT_PASS4_OMAN_DRIFT_20261022 [UNRATIFIED] -- 2026-10-22 21:15
Script pass4_oman_drift_20261022.py; report pass4_oman_drift_20261022_report.txt. Handoff 10/22 sec 4 item 1. The check was suggested in a SIMULATED review (4c: no credit to persons); the test and numbers are ours.
QUESTION: is the LSB outer-speed over-prediction (V_model > V_obs at R_last, Spearman -0.30 vs Sigma_h) uncorrected pressure support in the HI?
METHOD: V_c^2 = V_obs^2 + n sigma^2, n = 2R/R_HI (exponential outer HI disk, h_g = R_HI/2, Sigma(R_HI) = 1 Msun/pc^2; fallback h_g = 1.5 R_d for 4 galaxies without R_HI). Brackets sigma 6/8/12 km/s, n x0.5/1/1.5. Same 147 galaxies / 2696 points, no refit (universal rho_c 7.45 and relation (3) as in pass 3). Last-point ratios reproduce pass 2 to 4e-16.
## Results (ratio = V_model/V_c at R_last; Sigma_h quartiles low -> high)
| case, drift | median | Spearman(Sh) | low-Sh quartile | high-Sh quartile | low-g bins law-data | RAR fn-data |
| univ, none | 1.017 | -0.30 (2e-4) | 1.085 | 0.981 | +0.014 +0.084 +0.040 | -0.008 +0.031 +0.028 |
| univ, s8 n1 (fiducial) | 1.016 | -0.27 (8e-4) | 1.055 | 0.979 | -0.007 +0.073 +0.031 | -0.029 +0.020 +0.020 |
| univ, s12 n1.5 (extreme) | 0.991 | -0.17 (0.04) | 1.024 | 0.974 | -0.050 +0.048 +0.012 | -0.072 -0.005 +0.000 |
| rel3, none | 1.045 | -0.29 (4e-4) | 1.112 | 1.008 | +0.045 +0.076 +0.039 | (same as above) |
| rel3, s8 n1 | 1.028 | -0.26 (1e-3) | 1.105 | 1.005 | +0.024 +0.064 +0.030 | |
| rel3, s12 n1.5 | 1.011 | -0.17 (0.04) | 1.051 | 0.997 | -0.019 +0.039 +0.011 | |
(low-g bins: log g_bar/a0 in [-2,-1.75], [-1.75,-1.5], [-1.5,-1.25]; mean log g_model - log g_obs.)
Fiducial drift at R_last: median dV/V 0.6% (all), 1.1% (low-Sh quartile), 0.2% (high).
OLS on the UNCORRECTED log ratio, bootstrap (2000 galaxy resamples):
 univ: (10/V)^2 n +0.146 +- 0.115; log Vmax +0.018 +- 0.046; log Sh -0.039 +- 0.016 (2.4 sigma).
 rel3: (10/V)^2 n +0.188 +- 0.118; log Vmax -0.045 +- 0.038; log Sh -0.013 +- 0.014 (n.s.).
 Regressors collinear: drift proxy ~ log Sh -0.60; log Vmax ~ log Sh +0.67.
 Pure-drift expectation for the drift coefficient at sigma = 8: 0.217 x 0.64 = 0.139. Measured 0.146 +- 0.115: consistent with drift at the expected size, and also with zero.
## Reading
1. Drift is present at the expected size but SMALL: ~1% at R_last in LSBs at fiducial values. At universal rho_c it removes about a third of the low-Sh quartile excess (8.5% -> 5.5%); only the extreme bracket (sigma 12, steep gas) removes most (-> 2.4%). Relation (3): little change (11.2% -> 10.5% fiducial).
2. At universal rho_c the Sigma term survives drift + Vmax control at 2.4 sigma. Under relation (3) Sigma and V cannot be separated (collinear; Vmax is the strongest single correlate, -0.34). The residual follows 'low V / low Sigma / gas-rich' jointly.
3. Drift shifts the data, not the models: law-minus-RAR-function offsets in the low-g bins are unchanged by construction (+0.022/+0.053/+0.012 universal). Drift cannot help the law relative to MOND in the binned RAR; it only moves both toward/below the data.
4. Combined with the 18:00 quality addendum (clean subsample ~4%): the LSB residual is ~3-5% after drift and quality cuts, MOND's sign, 1-2.4 sigma depending on sample and rho_c rule. Not an artefact of drift; not robustly significant either. The draft's '4-8%' stands as a range; candidate wording '3-8%, of which about a third at universal rho_c is pressure support in the fiducial case' -- JM to rule.
## Caveats / VERIFY
- Whether SPARC curves already include asymmetric-drift corrections for some dwarfs (source papers, e.g. LITTLE THINGS / Oh et al. 2015) -- if so, this double-corrects those galaxies and overstates the drift. VERIFY in Lelli et al. 2016 sec 2 before quoting.
- Outer HI dispersion 6-12 km/s range from memory (typical ~8 in dwarfs; Iorio et al. 2017) -- VERIFY.
- n prescription crude (exponential through Sigma(R_HI) = 1 Msun/pc^2); real outer HI profiles may fall faster. Drift-corrected curves (Iorio 2017, Oh 2015) not used: no network in the VM. Owed: direct comparison for the SPARC overlap.
- Medians by quartile are fragile (rel3 vs univ quartile response differs by which galaxy sits at the median).

## ADDENDUM 2026-10-24 22:15 -- VERIFY items resolved from local PDFs (fetch block 3)
- HI DISPERSION: Iorio et al. 2017 (arXiv:1611.03865), Table 2 col. 9, median best-fit sigma_v per galaxy for 17 LITTLE THINGS dwarfs: 4.8-12.6 km/s, plus 21.0 for the starburst NGC 1569; median of the 17 = 8.0 km/s. SPARC (Lelli+2016, sec 1): "The H I velocity dispersion is typically of ~10 km s^-1, so pressure support becomes significant only in very low-mass galaxies with rotation velocities of ~20 km s^-1". -> pass 4 bracket 6-12, fiducial 8 SUPPORTED; cite Iorio+2017 Table 2 and Lelli+2016 sec 1.
- SPARC DRIFT TREATMENT: the SPARC paper describes NO asymmetric-drift correction of its own (no occurrence of 'asymmetric drift' in the text; only the sec 1 remark above). SPARC compiles literature curves, so individual source curves may carry their authors' corrections (e.g. LITTLE THINGS-type analyses). Double-correction risk is therefore limited to galaxies whose source curves were drift-corrected -- per-galaxy check against SPARC's reference column OWED if the sec 4.1 sentence is to be exact.
