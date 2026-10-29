# RESULT_PASS6_LOCAL_SIGMA_20261022 [UNRATIFIED] -- 2026-10-22 22:01 [time label corrected 22:03]
Script pass6_local_sigma_20261022.py; report pass6_local_sigma_20261022_report.txt; pass6_local_sigma_chi2.csv. Owed test (iii) from RESULT_PASS5. Motivation for a local velocity scale: Famaey, Khoury & Penco 2018 (v^4 ~ a0 G M_b(<r)).
TEST: sigma^2 = sqrt(G a0 M_b(<R_last))/2 instead of total M_b. M_b(<R_last) = stars (SPARC SBdisk/SBbul integrated, M/L 0.5/0.7; recovers 0.96 of 0.5 L) + gas M_g[1-(1+x)e^-x], x = R_last/h_g, h_g = R_HI/2 (pass 4; fallback 1.5 R_d). Median M_b(<R_last)/M_b = 0.79 (16-84 0.63-0.93).
PROCESS NOTE: first version used M = R V_bar^2/G -- over-states a disc's enclosed mass (median 1.17 M_b), so 'local' sigma came out LARGER; DISCARDED (not diagnostic).
## Results (one global rho_c each, chosen by per-point rms of log g)
| sigma rule | global log rho_c | rms(log g) | last-pt median | Spearman vs log Sh | vs f_gas | Sh quartiles (low->high) | f_gas quartiles (low->high) |
| BTFR (total M_b) | 7.5 | 0.153 | 1.025 | -0.31 (1e-4) | +0.25 (3e-3) | 1.114 1.047 0.983 0.994 | 0.982 1.002 1.053 1.078 |
| local (M_b(<R_last)) | 7.6 | 0.147 | 0.990 | -0.21 (1e-2) | +0.09 (0.3) | 1.042 1.011 0.969 0.985 | 0.985 0.990 0.999 1.007 |
(reference: RAR function 0.133; relation (3) 0.140)
Per-galaxy rho_c fits (catalogue errors, 3% floor): median chi2 BTFR-sigma 25.2, local 32.6, sigma free 7.5; local better in 74/147 (a wash).
## Reading
1. The local proxy REMOVES the gas-fraction trend (+0.25 -> +0.09, n.s.) and cuts the low-Sigma excess from 11% to 4%; the Sigma trend weakens (-0.31 -> -0.21) but stays at p 0.01. Global rms improves 0.006 dex. This supports reading (iii) of pass 5: the LSB/gas-rich residual is largely 'sigma set by mass not yet enclosed'.
2. Per galaxy it is a wash (better in 74/147; median chi2 worse) -- one sigma per galaxy taken at R_last is not a radius-dependent law; shapes inside are not improved.
3. NOT A LAW AS WRITTEN: sigma from M_b(<R_last) depends on where the data end (survey depth). A referee kills that immediately. The physical version is a sigma(r) varying with radius (FKP-type), which changes the law (sec 7 sentence, 4c). This result is a motivation for that research item (handoff item 5), not a fix for this note.
4. Sensitive to the gas extrapolation (h_g = R_HI/2): the gas beyond R_last is the crux. Owed if pursued: bracket h_g, and a true sigma(r) solve.
## For the note (JM to rule)
- sec 4.1: residual tracks gas fraction; a velocity scale set by the enclosed rather than the total baryon mass removes the gas-fraction trend and most of the LSB excess in a proxy test, but a radius-dependent sigma is a different law (sec 7). Or: leave out of this note, keep as research.
