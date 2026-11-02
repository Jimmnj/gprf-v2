# RESULT_CDM_20261031 -- NFW halo arms under the PREREG_NS evidence machinery  [UNRATIFIED]

Written 2026-10-31 12:55 from cdm_ns_20261031.jsonl (882/882 jobs; run on the 5950X 01:11-05:25, 30 procs; `python3 cdm_ns_metrics_20261031.py`). Prereg: PREREG_CDM_20261031.md (stamped 00:58). Comparison arms U/B/A from marg_ns_20261030.jsonl (RESULT_NS_20261030).

## 1. Coverage
All 13 arm x sig_int cells complete at N=147: U/B/N/F at sig_int {0, 0.03, 0.06}; A at 0 (as pre-registered).

## 2. Preferred sig_int (sum logZ)
U -13065 / -10081 / -10248 -> 0.03.  B -13093 / -10191 / -10303 -> 0.03.
N -12152 / -9935 / -10275 -> 0.03.  F -11805 / -9818 / -10224 -> 0.03.
Every arm prefers 0.03 dex intrinsic scatter in v.

## 3. Metrics by the pre-registered rule
M1 (at preferred sig_int; A at 0 vs N at 0.03 as the prereg fixed):
  U - N  -145.8 +- 2.9  median +0.41  tail 45 > +2.5 / 52 < -2.5  bulk(|d|<=20) -129.1  -> N favoured
  A - N -1221.0 +- 3.2  median +1.28  tail 55 / 36                 bulk +233.7         -> N favoured (by rule; see sec 5)
  B - N  -255.6 +- 2.9  median +0.23  tail 26 / 45                 bulk -178.9         -> N favoured
M2 F - N +117.0 +- 3.0; median -0.19; bulk +74.3; 19 galaxies > +2.5, none < -2.5. The c-M prior costs 117 nats.
M3 distance shifts (rms z_D; counts beyond 2 / 3 sigma): U 1.39 (18/8); B 1.22 (15/3); A 1.26 (13/6); N 0.78 (4/0); F 0.63 (3/1).
M4 sum logL_max at sig 0: N -10458, F -10206 (U/B/A files lack logL_max; N - U not computed).
M5 posterior-median rms (dex, at preferred): U 0.056, B 0.054, A 0.039, N 0.043, F 0.040.
M6 concentration vs Dutton-Maccio (posterior medians, sig 0): N offset -0.03 dex, rms 0.18, 14/147 beyond 0.3 dex; F offset -0.16 dex, rms 0.39, 73/147 beyond 0.3 dex. Median log M200: N 11.43, F 11.62.

## 4. Decisions by the pre-registered rule
D1 relation vs CDM: N favoured over U (-146), over A (-1221), over B (-256). All three cross the -5 threshold.
D2 haloes leave the c-M relation: YES (F - N +117; F rms offset 0.39 dex > 0.2).
D3 N bought with distances: NO. N moves distances less than U (4 vs 18 beyond 2 sigma), not more.

## 5. Flags on the comparison itself
(a) A - N is NOT a matched comparison. The prereg fixed A at sig_int 0 and N at its preferred 0.03. The sig_int 0 -> 0.03 step is worth ~2200-3000 nats per arm (sec 2), so the -1221 is dominated by the intrinsic-scatter mismatch, not by the halo. POST HOC (not pre-registered, reported for honesty, carries no decision weight): at matched sig_int 0, A - N = +995 (median +0.99, bulk +150); A - F = +649 (median +0.56, bulk +129). At matched sig 0 the per-galaxy relation beats both NFW arms; at matched 0.03 the universal relation loses to N by 146. The fair statement needs an A arm at 0.03 and 0.06 (follow-up prereg).
(b) U - N at 0.03: median +0.41 (U wins more galaxies than it loses) but bulk -129 and sum -146. N's margin is carried by ~50 galaxies where NFW's second parameter fits a feature the one-parameter relation cannot (worst: NGC0247 -47, NGC5371 -36, F571-8 -32, NGC6015 -29, UGC00731 -29, NGC2841 -27). U's wins: IC2574 +105, NGC3109 +63, ESO563-G021 +45.
(c) B - N: the RAR fitting function loses to NFW by 256 at matched 0.03, with a worse bulk (-179) than U's. MOND's one-function fit is penalised the same way; this is the Li et al. 2020 "haloes fit rotation curves" result reproduced, not a GPRF-specific loss.
(d) M6: with the c-M prior the haloes stay on the line (rms 0.18 vs the 0.11 prior width); freed, half of them leave it by > 0.3 dex and the evidence rewards leaving (+117). D2 reproduces Li et al. 2020 qualitatively.
(e) M5: at the preferred scatter, rms ordering A 0.039 < F 0.040 < N 0.043 < B 0.054 < U 0.056. The per-galaxy relation (one free parameter) matches the free NFW halo (two) in residual rms.
(f) M4 incomplete: U/B/A jsonl lack logL_max. The MAP csv route in the prereg caveat was not run.

## 6. Against stated expectations
Prereg expected: D1 U favoured over N by a modest margin -- WRONG (N by 146). A favoured over N clearly -- WRONG by the rule, TRUE at matched sig_int (post hoc). B vs N indistinguishable -- WRONG (N by 256). D2 yes -- RIGHT. D3 no -- RIGHT. Occam did not pay for the two halo parameters: wide uniform priors on log M200 (5.5 dex) and the c-M prior still left N ahead of a one-parameter universal relation.

## 7. Reading (UNRATIFIED)
The honest B.2 sentence is no longer "CDM not scored". It is: on SPARC rotation curves alone, under identical nuisances, an NFW halo with the cosmological concentration prior is favoured over the universal-rho_c relation (DlnZ -146) and over the RAR fitting function (-256); a per-galaxy rho_c matches or beats it at equal intrinsic scatter (post hoc; matched-scatter run owed). The relation's case was never that it out-fits a two-parameter halo per galaxy; it is that one global number plus the baryons does as well as the RAR function (RESULT_NS D1) and that the halo parameters are then predicted, not fitted. This result does not change that framing but it removes any claim that the evidence favours the relation over CDM at the galaxy level.

## 8. Proposed note edits (awaiting JM)
- B.2: replace "CDM not scored" with the paragraph above (numbers D1 U-N, B-N, M5, M6; the A caveat).
- Summary table: CDM column entry "NFW + c-M prior: DlnZ +146 over universal rho_c, +256 over RAR fn; rms 0.043 vs 0.056".
- Sec 7: follow-up item -- A arm at sig_int 0.03/0.06; cored-halo arm (DC14/Einasto) per the prereg caveat.
- No change to abstract limits (iv)/(v).

## 9. Follow-up prereg candidates
PREREG_CDM2: A at sig_int {0.03, 0.06}; cored halo arm; logL_max for U/B/A (M4).
