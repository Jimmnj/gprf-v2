# RESULT_BURKMU_20261102 -- Burkert with the Donato mu_0 prior (arm M) vs the relation  [RATIFIED JM 2026-11-02 19:20]
# 2026-11-02 16:30 EST. Prereg PREREG_BURKMU_20261101 (stamped by launch; JM "done" 16:13). Run: burkmu_ns_20261101.jsonl, 294 jobs (147 galaxies x sig_int 0.03/0.06), 5950X, ~3 min. Metrics: burkmu_metrics_20261102.py -> burkmu_report_20261102.txt.

## 0. VERDICT
L1 PARITY. M - A = -3.0 +- 3.0 in summed ln Z at the preferred sig_int (0.03 for both). The relation's +146 over the flat-prior Burkert (K) was prior volume: M - K = +142.9 (L2, expected > +100). M - D (DC14) = +61.6, but the bulk (|d| <= 20) is -11 with 23 vs 20 tail galaxies: a tail effect, not a bulk preference (L3 sign: M over D, reported as such).
Expectation on record (M beats A by +50 to +100) WRONG by ~60-100: the fit gain over A (M ln L_max = K's, -6.5; K - A = +97) is paid back exactly by the volume of the two priors (r_c over 3 dex; mu_0 0.2 dex), leaving parity.

## 1. Numbers
  arm   sum ln Z (0.03)   vs M      ln L_max vs M (0.03)   |zD|>2   rms
  M     -9483.6           --        --                      0        0.040
  A     -9480.5           -3.0      (K-A +97 -> M-A ~ +90)  12       0.038
  D     -9545.2           +61.6     -90.9                   0        0.041
  K     -9626.5           +142.9    +6.5                    1        0.035
  U     -10080.7          +597.2    --                      18       0.056
Preferred sig_int: 0.03 for every arm (0.06 loses by 500-600 everywhere). Pairwise at 0.03 identical to "preferred".
Tails (M - A): 14 galaxies > +2.5, 18 < -2.5; median -0.57; bulk -62.7 over 145 -> the relation is favoured in the bulk by ~0.4/galaxy and M wins back the sum in two galaxies with |d| > 20. Honest reading: parity, with the relation slightly ahead galaxy by galaxy.

## 2. J4 (does the data move mu_0?)
Posterior percentiles only are stored, so the log mu_0 width is bounded assuming independence of log rho_0 and log r_c (an upper bound; they anticorrelate): median bound 0.24 dex vs the 0.2 prior; below 0.2 in 46/147. Posterior median log mu_0 = 2.17 (prior 2.15). In K (flat priors) the same bound is 0.31 and the median 2.20. Reading: r_c is compressed (0.11 from 0.87 half-width) and mu_0 sits on its prior: M is a ~1.2-parameter model in effect, which is why it ties a 1-parameter relation. A sample-level width needs a rerun storing samples (45 min in the VM, 3 min on the 5950X); not needed for L1-L3.

## 3. J3
M moves no galaxy beyond 2 sigma in distance; A moves 12 (known, Table B.3). Burkert's second parameter absorbs what the relation puts into distance.

## 4. What changes in the note (B.3; HELD for JM)
- Table B.3: new row after "this relation, rho_c free": "Burkert halo, core-surface-density prior & 2 & -9484 & +90 & 0 & 0.040" (ln L_max vs DC14: M - D = +90.9).
- Text "First, ..." : the sentence "a Burkert with the observed core-surface-density band as prior is the fair cored comparator and is pre-registered [...]" -> "a Burkert with the observed core-surface-density band as prior (Donato et al. 2009; log mu_0 Gaussian, 0.2 dex) is the fair cored comparator: it ties the relation in evidence (-9484 against -9481, within the sampling error of 3) with the same raw fit as the flat-prior Burkert (+90 over DC14) and no distance freedom; the relation's +146 over the flat-prior Burkert was the volume of those priors. The tie is between one fitted parameter and two with one of them pinned by an observed band: the posterior on mu_0 is no narrower than its prior in most galaxies." The clause "the relation's one parameter does what the feedback core's two do" stays (DC14 comparison unchanged).
- Sec 3.2 / summary table: wherever "best evidence among per-galaxy models" is claimed, it becomes "ties the primed Burkert, ahead of DC14 by 65". Needs a grep (held).
- Framing (bp): never "beats Burkert"; say "ties a Burkert whose core surface density is fixed by the Donato band; beats it only when that band is left free".

## 5. Ratification
  BM1 L1 PARITY (M - A -3.0 +- 3.0) (T2)
  BM2 L2 +142.9; L3 +61.6 tail-driven, bulk -11 (T2)
  BM3 [SUPERSEDED by BM3' sec 6] J4: mu_0 on its prior, r_c compressed; M ~ 1.2 effective parameters (T2 within the independence bound)
  BM4 Expectation (+50 to +100) wrong; reason stated (reading)
  BM5 B.3 row + text per sec 4 (JM rules)
# END [RATIFIED JM 2026-11-02 19:20]

## 6. CORRECTION 19:40 (JM 19:30 "rerun storing samples"; burkmu_ns_20261102.py/.jsonl, 294 jobs M+K at sig 0.03; burkmu_J4_report_20261102.txt, burkmu_J4_20261102.csv)
BM3 as written is WRONG. Sample-level log mu_0 posterior (M): median half-width 0.116 dex vs the 0.2 prior band; narrower than the prior in 144/147 (113 below 0.15, 58 below 0.1); median 2.17. corr(log rho_0, log r_c) = -0.96, which is why the independence bound (0.24) was useless. Effective parameters M = 0.62 (mu_0) + 0.91 (r_c) = 1.53; K = 0.97 + 0.87 = 1.84 (mu_0 width 0.147 even with flat priors). Evidence: rerun sum -9480.1 vs -9483.6 on 11/01 (sampler scatter 0.36/galaxy); M - A = +0.4 on the rerun, -3.0 on the first: PARITY either way (L1 unchanged).
  BM3' The data constrain mu_0 to ~0.12 dex at the Donato value and r_c to 0.1 dex; the primed Burkert is a 1.5-effective-parameter model tying a 0.9-effective-parameter relation (T2)
  Tex B.3 sentence "no narrower than its prior in most galaxies" must be replaced (HELD for JM "apply").
