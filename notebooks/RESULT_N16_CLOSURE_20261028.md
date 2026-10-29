# RESULT_N16_CLOSURE_20261028 [UNRATIFIED]

Run 2026-10-28 00:40 under PREREG_N16_20261028.md (STAMPED JM 2026-10-28 00:24 "accepted - go"; FROZEN).
Script n16_closure_20261028.py (parts S1, S3, S2, S2k). Reports n16_S1/S3/S2/S2k_report_20261028.txt. Data n16_S1_20261028.csv, n16_S3_dwarfs_20261028.csv, n16_closure_points_20261028.csv.
Harness check: the univ and rel3 cases reproduce pass3_relation3_20261022_report.txt exactly (rms 0.155 / 0.140; Spearman +0.25 / +0.16; quartile ratios identical). Eq. (2) sigma^2 equals TEST1 s2b to 2e-16.

## VERDICT: FAIL (on C4)

| ID | Metric | Closure | PASS line | FAIL line | Status |
|---|---|---|---|---|---|
| C0 | rms log rho_c, closure vs fitted (S1, 120) | 0.688 | <= 0.64 | > 0.80 | GAP |
| C1 | per-point rms (S2, 147 gal, 2696 pts) | 0.149 | <= 0.145 | > 0.155 | GAP |
| C2 | Spearman res vs log Sigma_h | +0.10 | <= 0.16 | > 0.25 | PASS |
| C3 | V ratio spread, low-high Sigma quartile | 0.056 | <= 0.104 | > 0.15 | PASS |
| C4 | dwarfs within bar | 6/16 | >= 11 | <= 10 | **FAIL** |
| C5 | kill rule (emb -0.105, iso -0.136) | not triggered | 0 | 1 | PASS |

Rule: any of C0, C1, C4, C5 fails -> FAIL. C4 fails. N16 as stated is closed.

## 1. The pre-stated expectation, checked

- OLS of the closure log rho_c on S1: -0.492 log M_b + 0.799 log Sigma_h + 0.030 f_gas.
  Pre-stated: -0.50 / +0.80 / 0. Confirmed to the second decimal.
  Relation (3): -0.365 / +0.268 / -0.929.
- The Sigma over-response did NOT dominate the SPARC curves. C2 and C3 pass; C3 is better than relation (3) (0.056 vs 0.104).
  Reason: the curves are weakly sensitive to rho_c (C0 sd 0.69 dex costs only 0.009 in C1).
- The missing gas term shows in the residuals: Spearman vs f_gas -0.28 (rel3 -0.13), vs log M_b +0.21 (rel3 +0.08).
- Level: mean closure - fitted = -0.03 dex. Best constant k = 1.03 (log k +0.015). The King-Burkert identification needs no order-unity factor. The k diagnostic changes nothing (C1 0.149, C4 7/16).

## 2. Where it fails: the dwarfs, and why

- Closure dwarf rho_c 6.06..8.53, against relation (3) 8.37..9.10. Too little core density.
- Predictions low: mean log(pred/obs) -0.112, rms 0.184 (relation (3): +0.028, 0.148).
- **The closure dispersions equal MOND's isolated dispersions, sigma^4 = (4/81) G a0 M_b, to 0.006 dex rms (max 0.024, AndXIX, r_h 3.1 kpc).**
  MOND iso on the same rows: rms 0.183, mean -0.110. The closure inherits MOND's dwarf prediction and its failure.
- Scaling reason (sketch, not derived): in the low-Sigma limit mu0_M ~ sqrt(Sigma_b Sigma_M), so rho_c ~ Sigma_b Sigma_M / sigma^2 ~ M_b a0 /(G r_h^2 sigma^2).
  In a core-dominated dwarf sigma_los^2 ~ G rho_c r_h^2 ~ a0 M_b / sigma^2 ~ sqrt(G a0 M_b). So sigma_los^4 ~ G a0 M_b, with r_h cancelling -- MOND's isolated law.
  The 1% coefficient match with 4/81 is not derived here. OPEN.
- Consequence: the closure removes the relation's dwarf advantage over MOND (0.148 vs 0.183). The data on dwarfs want rho_c above the MOND core surface density, by about 0.5-1 dex.
- Embedded vs isolated: both low (-0.105 / -0.136); no environment split, as before.

## 3. Diagnostics (no verdict)

| Variant | C0 rms | dwarfs | dwarf rms |
|---|---|---|---|
| baseline (Sigma_b 0.75, RAR mu) | 0.688 | 6/16 | 0.184 |
| Sigma_b 0.50 Sigma_h | 0.709 | 7/16 | 0.214 |
| Sigma_b 1.00 Sigma_h | 0.690 | 7/16 | 0.168 |
| simple mu | 0.689 | 6/16 | 0.184 |
| best k (1.03) | 0.687 | 7/16 | 0.180 |

No variant reaches 11/16. The failure does not ride on the inherited choices.
Note: 'lambda simple 9.21' in the report header is a grid artifact (for the simple mu, L = 1/(1+x) and int L dx diverges). The eq. 11 integral to X0 is finite and the simple-mu diagnostic is valid.

## 4. Reading points for JM

(a) N16 as pre-registered fails, on the dwarfs. Close it in sec 7, or reword B.1 ("first candidate") -- JM to rule.
(b) The scaling (-0.49, +0.80, 0) was predicted before the run and confirmed; the failure was predicted (FAIL), but its site was not: SPARC curves tolerate the closure (C2, C3 pass); dwarfs do not.
(c) New, not pre-registered: closure dwarfs = MOND isolated to 0.006 dex. A referee will read this as "a MOND core surface density in the relation reproduces MOND". Analytic coefficient OPEN.
(d) Reverse reading: the dwarfs need rho_c r_c above Milgrom's mu0 at low Sigma. That is the same fact as Fig 4's data mass term (+0.14) that MOND lacks.
(e) Nothing enters the tex until JM rules.
