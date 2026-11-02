# PREREG_RHOC_SIGMAD_20261031 -- rho_c from Milgrom's central-surface-density bound  [AWAITING STAMP]

Written 2026-10-31 19:50 before any run. Prior art: Milgrom 2009, MNRAS 398, 1023 (arXiv:0909.5184): in MOND the phantom halo's central surface density is quasi-universal, Sigma_c* ~ gamma Sigma_M, Sigma_M = a0/(2 pi G) = 138 Msun/pc^2 (log 2.14), for systems with mean acceleration >~ a0 (gamma = 1 for mu = x/sqrt(1+x^2); 2/pi for the simple mu); for low-surface-density systems Sigma_c* ~ 2.4 (Sigma_b/Sigma_M)^(1/2) Sigma_M. Sigma_c is Donato+2009's rho_0 r_0; Milgrom defines it through the central column density Sigma(0) = int rho dz = (pi/2) Sigma_c (Burkert conversion). Matches Donato+2009 (2.15 +- 0.2) and Gentile+2009. Also: Kaplinghat, Ren & Yu 2020 (SIDM isothermal-Jeans: core density rises with baryon concentration). Not previously tried in RESULT_RHOC_RULES_20261030 (which tried universal, background-normalisation, CDR-borrowed, relation (3)).

## Question
Does the relation's fitted per-galaxy rho_c (arm A, RESULT_CDM2, posterior medians at sig_int 0.03) equal the value that makes the solution's central nb column density obey Milgrom's bound? If yes, rho_c is not a free lever: it is MOND's central-surface-density prediction, inherited.

## Rule (zero free parameters)
For each galaxy, at nominal nuisances (Ud 0.5, Ub 0.7, D_SPARC, i_SPARC), find rho_c such that
  (2/pi) * Sigma_nb(0) = Sigma_c*(Sigma_b0),   Sigma_nb(0) = 2 int_0^inf rho_nb(r) dr  (spherical solution, through the centre),
with Sigma_c*/Sigma_M = [ gamma^-2 + (2.4 y^(1/2))^-2 ]^(-1/2), y = Sigma_b0/Sigma_M  (smooth join of Milgrom's two limits; declared here),
Sigma_b0 = 0.5 SBdisk(r_1) + 0.7 SBbul(r_1) + Sigma_gas0, Sigma_gas0 = 1.33 M_HI / (2 pi (2 R_d)^2)  (declared approximation; no-gas variant as sensitivity).
Solver: the canon-form law (marg_fit_20261030 Gal.model: S = 2 tanh(x/2), sigma^2 = sqrt(G a0 M_b)/2, rho_amb 35.5), integrated beyond r_last with the point-mass baryon tail to r where rho_nb reaches the ambient floor (or 50 r_last). gamma = 1 primary; 2/pi reported.
Closed-form cross-check (not the rule): rho_c ~ 4 pi G Sigma_c^2 / sigma^2 -> 2 a0^(3/2) / (pi G^(3/2) M_b^(1/2)) at gamma = 1.

## Metrics (Stage 1, cheap: 147 root-finds)
P1 rms and median of (log rho_c,rule - log rho_c,A); P2 Spearman and slope of log rho_c,A vs log rho_c,rule; P3 the same against log M_b alone (does the rule beat M_b^-1/2 scaling?); P4 per-point RAR rms (dex in g) of the relation with rule rho_c vs universal 7.45 (0.057 at MAP) vs per-galaxy; P5 CDR slope with rule rho_c, by the RESULT_CDR recipe (obs 0.58 +- 0.05) -- if the recipe script is reusable; else Stage 2.
## Decisions (Stage 1)
Q1 "rule tracks the fit": YES if P2 Spearman > 0.6 AND P1 rms < 0.5 dex (A's own spread around universal is ~0.8 dex; to be stated from the data in the result). Q2 "no free amplitude": YES if |P1 median| < 0.3 dex at gamma = 1 or 2/pi. Q3 "rule fits curves": YES if P4 rms(rule) <= rms(universal). If Q1 and Q3 YES -> Stage 2: NS arm 'M' (rho_c from the rule inside the likelihood, nuisances marginalised) scored against A, U, D as in CDM3.
## Expectations stated now
Direction right (rho_c falls with M_b: relation (3) slope -0.365 vs rule -0.5); Spearman 0.4-0.7; P1 median offset +0.2..+0.5 dex (rule high: our core is not a Burkert core; the (2/pi) conversion is Burkert's); Q3 marginal. The low-Sigma branch should carry the CDR slope (P5 0.4-0.6). I expect Q1 borderline and Q2 NO at face value with a single global offset fixing it -- which would still be one global number, not one per galaxy.
