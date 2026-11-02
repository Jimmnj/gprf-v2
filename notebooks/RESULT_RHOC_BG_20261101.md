# RESULT_RHOC_BG_20261101 -- background-normalisation rule as a zero-freedom evidence arm  [UNRATIFIED]

Written 2026-11-01 12:55. Prereg PREREG_RHOC_BG_20261031.md (stamped by launch 21:40 10/31). Run bg_ns_20261031.py, 294/294 jobs, 5950X, 247 min. Arm G: rho_nb = A exp[S(x_b) - dphi_nb/sigma^2] - rho_amb, log A = 6.75 fixed; implied rho_c = A e^{S(x_b0)}. Metrics computed inline (same code pattern as cdm3_metrics; to be frozen as bg_metrics if ratified).

## Preferred sig_int
G: 0.03 (-10228.8 vs -10302.6). All eight arms now prefer 0.03.

## J1 -- Delta logZ at 0.03 (sum; median; tail +/-; bulk |d|<=20)  [@0.06]
  G - U   -148.0   -0.08   23/38   -139.3   [-54.5]
  G - B    -38.2   -0.32   45/48    -33.7   [-0.1]
  G - N   -293.8   +0.41   50/54    -77.2   [-27.9]
  G - D   -683.6   -0.30   29/58   -203.7   [-213.7]
  A - G   +748.2   +0.71   58/7    +338.4   [+279.4]
  G - F   -410.7 ; G - K -602.3
G's losses vs U: UGC00731 -27, UGC05764 -23, UGC06787 -19, NGC2841 -16, UGC06667 -16. Wins: IC2574 +41, UGC08837 +15, IC4202 +14.
## J2 -- logL_max at 0.03: G - A -694; G - D -689; G - N -288. (U has no logL_max.)
## J3/J4 -- distances and rms: G zD rms 1.36, 18 / 6 beyond 2 / 3 sigma (U 18 / 8); median per-galaxy rms G 0.049 (U 0.046, A 0.033).
## J5 -- implied log rho_c under the rule: median 7.44 [16-84: 7.17, 7.62]; x_b0 median 2.15 (so e^{S(x_b0)} sits near its e^2 cap for most galaxies). Against the fitted per-galaxy values: Spearman -0.01, rms 0.61, median (A - G) +0.11.

## Decisions
K1 G vs U: U FAVOURED (-148 < -5).
K2 G vs B: B favoured at 0.03 (-38); indistinguishable at 0.06 (-0.1).
K3 G vs D: D favoured (-684).
K4 rule tracks the fit: NO (Spearman -0.01).

## Against stated expectations
K1 G over U by +30..+150: WRONG (U by 148). K2 G over B by > +110: WRONG (B by 38). K3 D favoured: RIGHT. K4 NO: RIGHT (worse than expected: zero correlation). J3 "between U and A": WRONG (same as U).

## Reading (UNRATIFIED)
1. The background-normalisation rule reproduces the universal core's amplitude (implied median 7.44 vs 7.45) with a 0.45-dex, x_b0-driven modulation on top. That modulation is uncorrelated with what the curves want (Spearman -0.01 against the per-galaxy fit) and costs 148 nats against the universal core and 38 against the RAR function. The CDR slope it carries at nominal nuisances (0.53, RESULT_RHOC_RULES) is bought with central structure the full curves reject once D, i and Upsilon are free.
2. Consequence for the note: the BG rule does NOT become the zero-freedom model. It stays a named candidate variant in sec 2.4 at most, with this result attached; the CDR concession stands. The zero-freedom relation remains universal rho_c: +110 over the RAR function, -146 vs NFW, -536 vs DC14.
3. Pattern across all rho_c rules now tested in the full machinery or at MAP (universal, relation 3, BG, CDR-borrowed, Milgrom): nothing built from M_b, Sigma_b0, Sigma_h, f_gas or x_b0 moves rho_c in the direction the individual curves want. The per-galaxy rho_c is a measurement the law does not yet predict (framing (be')). The CDR (central Sigma_dyn within R_d) and the curve evidence pull rho_c in different directions for the same galaxies; the only arm that satisfies both is the per-galaxy fit, which also trades rho_c against Upsilon and D.
4. Since x_b0 sits at the saturation cap (median 2.15) for most galaxies, the BG rule is nearly a constant 7.6 for high-x_b0 systems and lower only for the lowest-bump dwarfs -- which is where its wins are (IC2574, UGC08837, D631-7). A rule that lowers rho_c for low-x_b0 galaxies only is a weaker, one-sided version of the same idea; not pre-registered, noted for the open list, not to be tuned in-sample.

## Proposed note edits (add to the HELD list)
- Sec 2.4 candidate sentence: keep, append "scored against the universal core under marginalised nuisances it is disfavoured (DlnZ -148)".
- B.1: the rules paragraph gains the BG evidence number; the CDR/curve tension sentence per reading 3.
- Summary table: no new row. Sec 7: rho_c rule item unchanged (still first).

## Files
bg_ns_20261031.jsonl (294 rows; logLmax, lrc_implied, xb0), bg_ns_run.log, this file.
