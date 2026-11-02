# PREREG_RHOC_BG_20261031 -- background-normalisation rule as a zero-freedom evidence arm  [AWAITING STAMP]

Written 2026-10-31 21:30 before the run. Motivation: RESULT_RHOC_RULES (10/30) found the BG rule carries the CDR slope (0.53 vs 0.58) at no cost in RAR rms; RESULT_RHOC_SIGMAD (10/31) found no other rule does; the per-galaxy fit does both but partly by nuisance trade-off. The BG rule has never been scored in the evidence machinery.

## The rule
rho_nb = A exp[S(x_b) - dphi_nb/sigma^2] - rho_amb, log A = 6.75 (Msun/kpc^3), i.e. the exponential is normalised at the background (x_b = 0) instead of at the centre; the implied central value is rho_c = A e^{S(x_b0)}, between A and A e^2, set by the galaxy's own baryon bump height x_b0. One universal constant replaces the universal rho_c (10^7.45). No per-galaxy parameter. log A = 6.75 was chosen chi^2-optimal in-sample on 10/30 at nominal nuisances; it is FIXED here (not re-fitted), so the arm pays no Occam cost but the constant is not independent of SPARC -- same status as universal 7.45.

## Run (bg_ns_20261031.py; pilot 2 galaxies OK 21:25)
Arm G at sig_int {0.03, 0.06}; 147 galaxies, 294 jobs; nuisances, priors, likelihood, nlive 300, dlogz 0.1, rwalk as PREREG_NS. Output bg_ns_20261031.jsonl (with logLmax, implied log rho_c and x_b0 at posterior-median nuisances).
Command (5950X): cd /mnt/c/users/jm/source/repos/gprf/notebooks && nohup python3 bg_ns_20261031.py --sig 0.03,0.06 --procs 30 > bg_ns_run.log 2>&1 &

## Metrics (each arm at its preferred sig_int; also matched 0.03)
J1 sum Delta logZ (G - U), (G - B), (G - N), (G - D), (A - G); median, tails, bulk as before.
J2 logLmax sums (G - U) [U lacks logLmax at 0.03? -- U rows from the first NS run lack it; if absent, report G - A and G - D only].
J3 distance shifts for G (rms; >2 / >3 sigma) vs U (18/8) and A (12/1).
J4 posterior-median rms (G vs U 0.046 median, A 0.033).
J5 implied log rho_c distribution (median, 16-84) vs A's per-galaxy values: Spearman, rms.
J6 P5-type CDR slope is NOT recomputed here (already 0.53 at nominal, RESULT_RHOC_RULES).

## Decisions
K1 G vs U: G favoured if J1 (G-U) > +5; U favoured if < -5; else indistinguishable. (Both zero-freedom; this is the direct test of the rule against the universal core.)
K2 G vs B: same rule. K3 G vs D (cored CDM, two per-galaxy parameters): same rule, expected D.
K4 "rule tracks the fit": YES if J5 Spearman > 0.6 (same bar as the Milgrom rule, which gave 0.24-0.34).

## Expectations stated now
K1: G beats U modestly, +30 to +150 (the 10/30 MAP rms was 0.157 vs 0.182 pooled; the CDR-passing central structure should help the inner points); K2: G over B by more than U's +110; K3: D favoured by several hundred. K4: NO (Spearman 0.3-0.5; rho_c spans only e^2 = 0.87 dex under the rule, A spans 4 dex). J3: G's distance use between U and A.
If K1 is YES and K2 YES, the note's zero-freedom model becomes the BG rule and the CDR concession is withdrawn; if K1 is NO, the BG rule stays a candidate variant (sec 2.4) and the CDR remains a limit.

## Caveats fixed now
log A in-sample; not a derivation of rho_c. Spherical nb solve as in all arms. sig_int 0 omitted (never preferred).
