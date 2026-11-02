# RESULT_RHOC_RULES_20261030  [UNRATIFIED]
Prereg PREREG_RHOC_RULES_20261030.md. Run 21:34 on marsGB WSL, 69 s, 147 galaxies, nominal parameters. Files: rhoc_rules_points_20261030.csv, rhoc_rules_summary_20261030.txt.

| arm | S1 RAR rms (all / g_b<0.1a0) | S2 beta | C2 | S3 med chi2/N | S4 rho_s vs A0 | sd log rho_c |
|---|---|---|---|---|---|---|
| A0 free per galaxy | 0.086 / 0.096 | 0.626+-0.035 | pass | 2.5 | 1 | 0.56 |
| U0 universal 7.45 | 0.155 / 0.180 | 0.276+-0.071 | fail | 14.4 | - | 0 |
| R0 relation (3) | 0.140 / 0.164 | 0.373+-0.048 | fail | 14.5 | 0.43 | 0.22 |
| BG background norm., log A = 6.75 | 0.157 / 0.182 | 0.534+-0.055 | pass | 15.7 | 0.16 | 0.21 |
| CDR borrowed rule | 0.200 / 0.186 | 0.596+-0.014 | pass | 19.0 | 0.32 | 0.92 (21 clipped) |
S1 here is the unbinned per-point rms over harness-selected points at nominal parameters; it is not the binned or marginalised figure quoted in the note and is for internal comparison only.

## Decisions
BG: D1 PASS (carries the CDR slope), D2 FAIL (RAR rms no better than universal), D3 FAIL.
CDR: D1 PASS (by construction), D2 FAIL (RAR rms worse than universal; imposed rho_c scatters 0.92 dex, 21 galaxies hit the bounds), D3 FAIL.
Neither rule replaces the per-galaxy fit. Relation (3) remains the best rule on S1 and the worst on S2.

## Against stated expectations
BG was expected to fail the CDR slope; it passes (0.53 vs 0.58). The factor e^{S(x_b0)} carries the population slope. It does not carry the per-galaxy values (rho_s 0.16), so it leaves the RAR scatter at the universal level. Prediction recorded as wrong.
CDR was expected to land near A0 on S1; it lands below U0. The central Toomre density is too noisy at the first radius to set rho_c galaxy by galaxy.

## Reading
1. The CDR slope is a weak population constraint: three different rho_c prescriptions pass it. It discriminates against a universal core and against relation (3), not among Sigma-dependent rules.
2. The RAR scatter is the strong constraint, and it is per galaxy. The fitted rho_c carries information that none of M_b, Sigma_h, f_gas, x_b0 or Sigma_*(0) supplies (sd 0.56 dex, of which relation (3) captures 0.22).
3. Two readings of (2), not separable at nominal parameters: rho_c tracks galaxy physics not in those variables, or rho_c absorbs distance, inclination and M/L errors (RESULT_MARG: universal rho_c fits at 2 sigma distance shifts). The NS evidence run (arm A vs U at marginalised nuisances) is the test that separates them.

## For the note (pending ratification)
Sec 3.1/B.1: one paragraph -- four rules tried (universal, relation 3, background normalisation, CDR-borrowed); all hold the CDR slope except the first two, none holds the RAR scatter; the per-galaxy rho_c stays, now as a measured quantity whose dependence on the baryons is only partly captured. Sec 7: the two readings and the NS test. Sec 2.4: background normalisation is a legitimate variant of the carrier picture and carries the CDR slope; recorded as a candidate, not adopted.
