# RESULT_OCCAM_20261101 -- posterior compression and effective parameter counts per arm (S2 of REVIEW_COLD_v2.0.1)  [UNRATIFIED]
Written 2026-11-01 19:45 (JM "go" 19:36). Script occam_20261101.py; inputs: the four NS jsonl files at sig_int 0.03, posterior 16/84 half-widths vs prior sigma (uniform: width/sqrt12; Gaussian: its sigma). No new runs.

## Numbers (147 galaxies; model parameters starred)
- Compression sigma_post/sigma_prior (median): A log rho_c 0.12 (constrained, rho<0.5, in 92%). D log M200 0.07 (93%), DC14 c-M offset 0.81 (10%). N log M200 0.11 (97%), offset 0.80 (8%). F log M200 0.22, log c 0.36. K log rho0 0.18, log r_c 0.15 (both ~90%).
- Effective constrained model parameters per galaxy, mean sum(1-rho^2): A 0.91 (nominal 1); D 1.31 (2); N 1.32 (2); F 1.64 (2); K 1.78 (2).
- Occam cost of model parameters, sum ln(sigma_prior/sigma_post): A 286; D 382; N 356; F 406; K 501. D - A = 96, against the evidence-minus-raw-fit gap of 60 (+65 evidence, +5 ln L_max); same order, crude estimator.
- Nuisance compression (median zD): A 0.73, D 0.83, N 0.83, K 0.96, F 0.95; zero-freedom arms U 0.55, B 0.60, G 0.57. The models with no per-galaxy parameter squeeze distance hardest; among the per-galaxy arms the relation squeezes it more than the haloes (12 vs 0 beyond 2 sigma, B.3).

## Reading
1. The MOND reviewer is right in substance: DC14 with the 0.11-dex c-M prior is 1.3 effective parameters, not 2; its concentration offset is constrained in one galaxy in ten. "The evidence margin +65 is the cost of the second parameter" should read "the Occam cost of a second parameter the data constrain in one galaxy in ten; by effective count the comparison is 0.9 against 1.3 parameters per galaxy".
2. The Burkert arm's evidence deficit is prior volume (5 dex x 3 dex, both parameters well constrained): Occam cost 501 vs A's 286. A Burkert with a mu_0 prior (S3) would recover most of that; expectation on record: Burkert-with-prior beats A in evidence by ~50-100 (raw fit +97 minus a smaller Occam gap).
3. Universal rho_c, the RAR function and the BG rule buy fit with distance (zD compression 0.55-0.60 vs 0.83 for the haloes). This is the B.3 "distance freedom" row seen from the other side and belongs next to it.

## Proposed B.3 edits (HELD)
- Replace "with the evidence margin (+65) equal to the cost of the second parameter and the raw fit equal" by the wording in reading 1. Add one sentence with the compression numbers and effective counts (0.9 / 1.3 / 1.3 / 1.6 / 1.8). Add a row to Table B.3: "effective parameters". Note the Burkert deficit as prior volume and point to the S3 arm.
