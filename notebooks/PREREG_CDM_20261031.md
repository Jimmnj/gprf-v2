# PREREG_CDM_20261031 -- LambdaCDM halo arm under the same evidence machinery  [STAMPED: JM 2026-10-31 00:58 "yes - so i can run it overnight?"]

Written 2026-10-31 01:10 before any run. Appendix B.2 says CDM is not scored in our harness; this fixes that for rotation curves only.

## Question
On SPARC, under identical nuisance priors, likelihood and sampler as PREREG_NS, how does an NFW halo compare in Bayesian evidence with (U) the relation at universal rho_c and (A) the relation at per-galaxy rho_c, and with (B) the RAR fitting function?

## Arms (new; same Gal.curves, same nuisances log Ud, log Ub, D, i and priors as PREREG_NS)
- N: NFW halo with the cosmological concentration-mass prior. V^2 = v_b^2 + V_NFW^2(r; M200, c). Priors: log M200/Msun ~ U[8.5, 14]; log c200 = 0.905 - 0.101 log(M200 / 10^12 h^-1 Msun) + delta, delta ~ N(0, 0.11) (Dutton & Maccio 2014, z=0, h=0.7). rho_crit from H0 = 70. Two halo parameters per galaxy.
- F: NFW halo with a flat concentration prior, log c200 ~ U[0, 2], log M200 ~ U[8.5, 14]. The "free halo" ceiling that Li et al. 2020 show the fits want.
- sig_int in {0, 0.03, 0.06} for both; nlive 300, dlogz 0.1, rwalk, as the NS run. Outputs cdm_ns_20261031.jsonl with the same fields plus logL_max and posterior lM200, lc.

## Metrics
M1 sum Delta logZ (U - N), (A - N) [A at sig 0 only], (B - N), each at the arms' preferred sig_int; per-galaxy tail audit as RESULT_NS.
M2 sum Delta logZ (F - N): the evidence cost of the c-M prior.
M3 N and F distance shifts z_D (rms, counts beyond 2 and 3 sigma) against U and B.
M4 sum Delta logL_max (N - U), (F - A): raw fit gain of the extra halo parameters before the Occam cost.
M5 posterior-median residual rms per arm.
M6 N posterior c vs M200 against the Dutton-Maccio line: median offset and rms (does the data pull the haloes off the relation?).

## Decisions
D1 "relation vs CDM": U favoured over N if sum Delta logZ (U-N) > +5; N favoured if < -5; else indistinguishable. Same for A-N and B-N.
D2 "haloes want to leave the c-M relation": yes if F - N > +5 AND M6 rms offset of F's log c from the line > 0.2 dex.
D3 distances: "N bought with distances" if N's count beyond 2 sigma > 1.5 x U's.
Expectations stated now: D1 U favoured over N by a modest margin (Occam on two parameters over wide priors); A favoured over N clearly; B vs N near indistinguishable. D2 yes. M4 N - U > 0 (two more parameters fit better) and F - A > 0. D3 no.

## Caveats fixed now
NFW only; no feedback-cored profile (DC14, Einasto) -- a cored halo arm is the natural follow-up if N loses. Priors on log M200 are wide and uniform; the evidence penalty for that width is part of what "two free parameters" costs and is reported separately via M4. Rotation curves only; nothing here bears on cosmological scales.
