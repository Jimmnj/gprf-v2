# PREREG_DIVERSITY_20261105 -- diversity at fixed baryons: matched SPARC pairs
Written 2026-11-05 18:15 EST BEFORE code. [UNRATIFIED]. Trigger: JM 18:04 "do both and verify angus".
## 1. Question. MOND and the empirical acceleration relation fix the curve once the baryons (and nuisances) are fixed. Eq. (1) allows two galaxies with the same baryon profile to differ in the transition region through rho_c. Do matched pairs with different inner curves exist, beyond what distance, inclination and M/L freedom can absorb?
## 2. Inputs. Marginalised MAP fits (marg_fit_params_20261030, arm A: rho_c with D, i, M/L free) for rho_c; per-galaxy nested-sampling evidences (marg_ns_20261030_pergal: logZ_A rho_c free, logZ_B empirical function, logZ_U universal; same nuisance priors; sig_int 0.03). Catalogue log M_b, log Sigma_h, f_gas.
## 3. Method. Pairs matched within |dlog M_b| <= 0.15, |dlog Sigma_h| <= 0.15, |df_gas| <= 0.15 (all three). For each pair: D = |dlog rho_c,A|. Random pairs: all unmatched pairs. Evidence check per galaxy: dZ = logZ_A - logZ_B (positive = the empirical function with free nuisances fits worse than the relation with rho_c free, after Occam).
## 4. Criteria (fixed now)
K1 EXISTENCE: >= 5 matched pairs with D >= 0.3 dex (and >= 2 with D >= 0.5).
K2 NOT NUISANCE: among K1 pairs, in at least half, one member has dZ >= +5 (the empirical function cannot absorb its curve with free D, i, M/L).
K3 BARYONS DON'T FIX rho_c: median D over matched pairs / median D over random pairs >= 0.7. (A baryon-determined rho_c would give << 1.)
Expectation (Claude): K1 pass (~10-20 pairs), K2 pass, K3 ~0.8. 
Decision: all pass -> sec 7 / B.1 sentence: "matched baryon pairs differ in rho_c by up to X dex and the empirical function cannot absorb the difference with free nuisances: diversity at fixed baryons is predicted by eq. (1) and forbidden by a baryon-only function." K2 fail -> the differences are nuisance-absorbable; no claim. K3 << 0.7 -> rho_c is more baryon-determined than B.1 says; report.
# END
