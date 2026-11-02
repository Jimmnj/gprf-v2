# PREREG_RHOC_RULES_20261030 -- Replacing the per-galaxy rho_c fit by a rule  [STAMPED: JM 2026-10-30 20:53 "you're on a roll - so yes go"]

Written 2026-10-30 21:00 before any run. Follows RESULT_CDR_20261030 (universal rho_c and relation (3) fail the CDR slope).

## Question
Can the per-galaxy rho_c be replaced by (a) a physical rule with one universal constant, or (b) the CDR as a borrowed input, without losing the RAR scatter, and does each carry the CDR slope?

## Arms (nominal Upsilon_d=0.5, Upsilon_b=0.7, D, i; harness galaxy list; solver = marg_fit_20261030.Gal.model unchanged)
- A0: per-galaxy rho_c (1-D chi^2 fit; reference ceiling).
- U0: universal log rho_c = 7.45 (reference floor).
- R0: relation (3) (for the record).
- BG (physical): density normalised to the background instead of the galaxy centre, rho_nb = A exp[S(x_b) - delta phi/sigma^2] - rho_amb, i.e. log rho_c,i = log A + S(x_b0,i)/ln 10 with x_b0 the harness's own central baryonic potential in units of sigma^2. One universal constant log A, chosen to minimise the summed harness chi^2 over the sample on a grid 6.0..8.5 step 0.05 (the RAR-rms-optimal log A is reported alongside as a check, not used).
- CDR (borrowed): rho_c,i set so that the model's Toomre Sigma_dyn(0) equals the CDR value at the galaxy's Sigma_*(0), using OUR observed fit (beta 0.578, Sigma_0 1678, Sigma_crit 1263; same recipe as PREREG_CDR). 1-D root in 6<=log rho_c<=10; galaxies without a root are clipped and counted.

## Scores (per arm)
S1 RAR scatter: rms of log(g_obs/g_mod) over the harness-selected points; also over the low-acceleration subset g_b < 0.1 a_0.
S2 CDR slope beta (PREREG_CDR recipe, alpha=1) and C2 vs beta_OBS = 0.578+-0.048.
S3 median chi^2/N per galaxy.
S4 agreement with the free fit: Spearman and rms of (log rho_c,rule - log rho_c,A0).

## Decision rules
D1 "carries the CDR": passes C2.
D2 "holds the RAR": S1 <= S1(A0) + 0.010 dex AND S1 <= S1(U0) - 0.005 dex.
D3 an arm "replaces the fit" if D1 and D2 both pass.
Expectations stated in advance: BG fails D1 (beta ~0.3-0.4, range of e^S too small) and gives S1 between U0 and A0; CDR passes D1 by construction and lands near A0 on S1 but not within 0.010; S4 Spearman > 0.6 for CDR, < 0.5 for BG.

## Outputs
rhoc_rules_20261030.py, rhoc_rules_points_20261030.csv, rhoc_rules_summary_20261030.txt, RESULT_RHOC_RULES_20261030.md.
