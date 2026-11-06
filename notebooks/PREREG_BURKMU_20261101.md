# PREREG_BURKMU_20261101 -- Burkert halo with the core-surface-density prior as a fair cored comparator (S3 of REVIEW_COLD_v2.0.1)  [UNRATIFIED]
Written 2026-11-01 19:58 (JM "all" 19:41). Stamp = JM launching the run. Motivation: RESULT_OCCAM_20261101 -- the free Burkert arm K loses to the relation by 146 in evidence while out-fitting it by 97 in ln L_max; its Occam cost (501 nats) is the 5 dex x 3 dex flat prior volume. The MOND reviewer called K "set up to lose". Arm M replaces the flat rho_0 prior with the observed core-surface-density band.

## Arm M
Burkert v^2 as in cdm3_ns_20261031.py (arm K code path). Priors: log r_c ~ U[-1, 2] kpc (as K); log mu_0 = log(rho_0 r_c) ~ N(2.15, 0.2) Msun pc^-2 (Donato et al. 2009, the band drawn in Fig. 4); rho_0 = mu_0 / r_c. Nuisances, likelihood, sampler identical to every other arm (log Upsilon_d, log Upsilon_b N(., 0.1); D, i Gaussian at SPARC errors; sig_int 0.03 and 0.06; dynesty nlive 300, dlogz 0.1, rwalk). 147 galaxies x 2 sig_int = 294 jobs. Script burkmu_ns_20261101.py (copy of cdm3_ns with the M transform); pilot OK in the VM (DDO064: M -9.95 vs K -10.87 at nlive 100, ln L_max equal).
Command (5950X, WSL): cd /mnt/c/Users/jm/source/repos/gprf/notebooks && python3 burkmu_ns_20261101.py --arms M --sig 0.03,0.06 --procs 30

## Metrics (at each arm's preferred sig_int; also at 0.03)
J1 sum Delta ln Z, median, tail counts, bulk (|d| <= 20) for M - A, M - K, M - D, M - U. J2 Delta ln L_max M - A, M - K. J3 distance shifts beyond 2 / 3 sigma for M. J4 posterior compression of log mu_0 (does the data move it off the prior?) and effective parameter count.
## Decisions
L1 M vs A: if M - A > +5, the relation's +146 over Burkert was prior volume and the fair cored two-parameter halo beats the one-parameter relation in evidence; B.3 and sec 3.2 say so. If |M - A| <= 5, parity. If M - A < -5, the relation's economy survives a physically-primed Burkert.
L2 M vs K: expected > +100 (the removed prior volume).
L3 M vs D: sign reported; both are cored two-parameter haloes with a physical prior.
## Expectations on record
L1: M beats A by +50 to +100 (ln L_max +97 for K, minus a smaller Occam gap; some of the 97 may go if the mu_0 prior fights the fit in the LSBs). L2: +150 to +250. L3: M over D by tens. J4: log mu_0 posterior narrower than 0.2 in most galaxies (the data do constrain it: K's two parameters were both compressed to ~0.15).
## Consequence either way
The B.3 Burkert sentence changes: "fitting slightly worse" is already gone (R12); after M it reads from this result. Table B.3 gains a row.
