# RESULT_TWOCENTRE_20261102 -- Milky Way + LMC two-centre solve under ONE-FIELD  [RATIFIED JM 2026-11-02 19:20]
# gprf\notebooks\RESULT_TWOCENTRE_20261102.md. 2026-11-02 ~13:10 EST (Fable session). Prereg: PREREG_TWOCENTRE_20261102.md (stamped JM 12:25, unchanged).
# Scripts: twocentre_20261102.py (solver; args L= NMU= NR= lrc= MLMC= aLMC= D= gas= agas= tag=), twocentre_v1_shoot_20261102.py (V1), twocentre_twohalo_20261102.py (column ii).
# Outputs: twocentre_<tag>_20261102.json/.npz (primary, L24, L60, D500, rc788, twocomp), twocentre_twohalo_20261102.json, GPRF_TwoCentre_20261102.xlsx.

## 0. VERDICT (one line)
D2: REVERSE-IF FIRES. Under ONE-FIELD the relation places 1.2e9 Msun of nb excess within 30 kpc of the LMC; the LMC's model mass within 30 kpc is 4.2e9 against 4.3e10 for the two-halo reading by the same code (ratio 0.098 < 0.5) and 1-2e11 measured (25-50x short). Expectation (prereg sec 4: 1-3e9 within 30, D2) met.

## 1. Setup (as pre-registered; deviations in sec 4)
Law eq. (1)-(2), tanh2 cap, floor -rho_amb; one sigma 125.1 km/s (M_b,MW 6.15e10), one rho_c (10^7.45 primary), rho_amb 35.5. MW baryons spherical-equivalent (Hernquist a = 0.695 Rd for the three discs; bulge as given). LMC Plummer 3.0e9, a 1.5, at D 50 kpc on the symmetry axis. Multipoles l = 0..40 about the MW centre, 96 Gauss-Legendre nodes on [-1,1], 1200 log-radial points 0.05-3000 kpc; phi_b analytic; Picard with damping 0.3 (converged to 7.5e-5 in 60 iterations; the Newton-Krylov fallback was not needed for any run). Sphere integrals about the LMC on a dedicated (d, cos a) quadrature evaluating the law from the multipole coefficients (exact in angle).

## 2. Validity controls
- V1 (2D MW-alone vs 1D shooting, same spherical baryons): v_c at 8/12/16/20/27 kpc differs by +0.10/+0.06/-0.02/-0.03/-0.01 km/s (criterion 0.1: PASS at the edge); rho_nb at 50 kpc 2.08436e5 vs 2.08426e5 (5e-5); M_nb(<50 kpc about the MW centre) 3.001e11 vs 3.010e11 = 0.3% (criterion 0.1%: MISSED; cause: the 2D radial grid starts at 0.05 kpc and the inner-core mass below it is dropped; irrelevant at the LMC's position where the density agrees to 5e-5). Reported as a partial pass; the prereg's 0.1 km/s velocity criterion was written for the mw_decline disc model and is applied here to the same-baryon 1D solve, which is the meaningful control.
- V2 (LMC at 500 kpc): excess within 10/30/50 kpc of the LMC 1.3e6/1.2e7/3.4e7 (the response of a 2357 Msun/kpc^3 medium; nesting factor at centre 1.72); MW-alone run identical. PASS.
- V3 (L = 24/40/60, NMU 96/96/128): excess within 30 kpc 1.207194e9 / 1.207195e9 / 1.207200e9 (< 1e-5). PASS. Reason: the excess is dominated by the LOCAL response to the analytic LMC potential; the multipole-carried self-term is ~14% and smooth.

## 3. Results (primary: rho_c 10^7.45, Plummer 3e9)
Spheres about the LMC centre, Msun:
  R    M_b,LMC   MW-alone nb   two-centre nb   nb EXCESS   ONE-FIELD M_LMC(<R)   (i) linearised   (ii) two-halo M_LMC(<R)   ratio OF/2h
  10   2.90e9    8.80e8        9.97e8          1.17e8      3.02e9                1.11e8           1.44e10                   0.21
  30   2.99e9    2.56e10       2.68e10         1.21e9      4.20e9                1.06e9           4.29e10                   0.098
  50   3.00e9    1.52e11       1.55e11         3.60e9      6.60e9                3.58e9           4.69e10                   0.14
- Full-form nesting factor at the LMC centre: rho_nb(two-centre)/rho_nb(MW alone) = 1.67 (note sec 2.3 and sec 7 carry 1.36, the linearised NEST-1 number). MW nb density at 50 kpc 2.08e5 Msun/kpc^3.
- (i) linearised local integral (MW-alone nb x [exp(S(phi_LMC/sigma_MW^2)) - 1]) reproduces the full solve to 5-14%: the self-term adds 0.15e9 within 30 kpc. The full form is NOT larger than the linearised one by any factor that matters.
- (ii) two-halo column (pre-ONE-FIELD reading, same code): LMC own sigma 58.8, rho_amb 2.085e5, r_t (density match) 31.2 kpc, own nb 3.07e10 within r_t, plus 0.36 x MW nb inside r_t (9.2e9 at 30 kpc). The note's "~4e10" (sec 7) is this number (4.3e10 at 30 kpc). RESULT_RULEA's density-match 20.8 kpc / 2.2e10 used host_rho(200, 50) = 2.96e5 for the MW density; the self-consistent MW profile here gives 2.08e5 at 50 kpc, hence the larger r_t.
- MW circular speed in the plane perpendicular to the MW-LMC line, 8-27 kpc: raised by 0.2-0.6 km/s by the LMC (its own potential plus the response). The MW fit of RESULT_MW_DECLINE is not disturbed.
- Sensitivity: rho_c 10^7.88 (MW fitted): excess 0.89e8/0.84e9/2.6e9 (MW nb at 50 kpc 1.70e5); two-component LMC (stars 2.7e9 a 1.5 + gas 0.5e9 a 3.0): 1.06e8/1.09e9/3.25e9, nesting factor 1.58. Neither crosses the 0.5 boundary; D2 in every run.

## 4. Deviations from the prereg
- V1 mass criterion missed at 0.3% (sec 2); not a decision input. No other deviation. Picard converged; no Newton-Krylov. The thin-disc worst-case orientation run (prereg sec 2) is not triggered: the primary ratio 0.098 is a factor 5 from the boundary.

## 5. Reading
- ONE-FIELD removes the LMC's own sigma = 58.8 isothermal tail, and that tail WAS the two-halo number (3.1e10 of the 4.3e10). What replaces it is the MW medium's response to a 0.56 sigma_MW^2 potential dip: a factor 1.67 at the LMC centre, 1.1 by 10 kpc, 1.2e9 in all. The self-term (the only place a two-centre solve could have exceeded the linearised estimate) adds 14%. No Jeans runaway: the LMC's scale is well inside the local Jeans length (~31 kpc).
- So the canon REVERSE-IF condition ("fails the LMC mass worse than the two-halo reading") is met, by 10x at 30 kpc. The concession in sec 7 goes from "4-5x light" to "25-50x light" if ONE-FIELD stands.
- What the number does NOT say: it does not say the two-halo reading was right. The two-halo LMC was a second response with a second sigma, which ONE-FIELD identified as double counting on derivation grounds (resummation family). The choice is now between the derived rule and the LMC: (A) keep ONE-FIELD, widen the concession, and state in sec 7 that the LMC's perturber mass is the relation's largest single discrepancy, with the CDM/MOND routes as already worded; (B) reinstate a per-body response, which reopens S5 (which bodies, which sigma, double counting) and the cold-review objection it was adopted to answer. (C) a third reading -- the LMC carries its own sigma because it was an isolated galaxy before infall and the nb configuration has a memory (re-formed vs carried) -- is the momentum/history item (canon 2026-10-24) and is not available without an action. Claude's recommendation: (A); the relation is an equilibrium statement and the LMC is a first-infall body on which the field has not equilibrated; say so once, with the number.
- Ratio check (handoff item 2): under ONE-FIELD there is no (M_b)^{3/4} scaling for the LMC's nb mass at all, so the "1/8 vs 1/10" comparison has no footing and should be dropped from the owed list; the perturber mass is the excess above, full stop. The LIT item (what the reflex/stream papers actually constrain, and within what radius) stays owed before the wording is final.
- Referee flags: (a) the MW nb density at 50 kpc (2.1e5) is itself untested; a heavier MW disc (McGaugh 2019) raises sigma and the medium density, scaling the excess roughly linearly -- not by 10x. (b) The LMC's gas (0.5e9) and extended stellar envelope change the number by 10%. (c) An off-axis MW disc cannot change a monopole-dominated medium at 50 kpc by more than a few percent (Gate 2b q >= 0.9 beyond 5 Rd).

## 6. Ratification (for the ballot)
  TC1 V1-V3 controls as stated (V1 mass 0.3%, velocity pass) (T2)
  TC2 ONE-FIELD LMC: excess 1.17e8 / 1.21e9 / 3.60e9 within 10/30/50 kpc; model mass 4.2e9 at 30 kpc; nesting factor 1.67 at centre (T2)
  TC3 Two-halo same-code column 1.44e10 / 4.29e10 / 4.69e10; ratio 0.21 / 0.098 / 0.14 (T2)
  TC4 D2 REVERSE-IF FIRES by the stamped rule (T2 numbers; rule stamped 12:25)
  TC5 Linearised local integral within 14% of the full solve; the 1.36 in the note should read 1.67 (central) or be replaced by the integral (T2)
  TC6 Ratio-check item (1/8 vs 1/10) has no footing under ONE-FIELD: recommend DROP from the owed list (reading)
  TC7 [REJECTED JM 2026-11-02 19:27; withdrawn 13:34] Recommendation (A): keep ONE-FIELD, widen the LMC concession to 25-50x (superseded by TC10/TB4)

## 7. ADDENDUM 13:45 -- dwarf protocol under STRICT ONE-FIELD (JM 13:34 "go"). Script onefield_dwarfs_20261102.py, csv onefield_dwarfs_20261102.csv.
Control: the current protocol (per-body sigma from the dwarf's own M_b, host as ambient) reproduces the note: 9/16, rms 0.148 dex.
Strict ONE-FIELD (satellite has no response; host medium at the satellite's distance responds with host sigma^2 = v_c^2/2 to the satellite's phi_b): 0/16, rms 0.750 dex. Every dispersion falls to the stars-only Newtonian value (x_b0 = 1e-4 to 8e-3, so the host medium's local response is a 1e-4 to 1e-2 effect): Fornax 4.8 vs 11.7, Sculptor 2.4 vs 9.2, Leo I 3.9 vs 9.2, And XIX 0.6 vs 7.8, Crater II 0.4 vs 2.7. Only the diffuse objects (Crater II, And XIX, Sextans) pick up a host-medium excess comparable to their stellar mass, and it is far from enough.
Caveat on the three ISO objects (Cetus, Tucana, And XVIII at 450-880 kpc): there the ambient is near background and the object is an isolated body with its own sigma under any reading; the one-field column is not the right reading for them (the per-body column is). For the 12 satellites it is: 0/12.
READING: Appendix B's dwarf result exists only because satellites are treated as bodies with their own sigma. Strict ONE-FIELD removes it entirely. The LMC (30x) and the dwarfs (0/12, 0.75 dex) are the same question: is a self-bound satellite on its own orbit a body? Under (A') yes -> LMC two-halo 4e10 (3-5x short), dwarfs 9/16; under (A) no -> both lost.
  TC8 Current protocol control 9/16 rms 0.148 reproduced (T2)
  TC9 Strict ONE-FIELD dwarfs 0/16, rms 0.750 dex; satellites 0/12 (T2)
  TC10 Recommendation revised: (A') -- one-field within a body; self-bound satellites are bodies (stated assumption, S5 threshold); LMC = two-halo column; sec 2.3 "chosen by the analyst" becomes a stated rule (bound + own orbit + outside the host's baryons) (reading; JM rules)
# END -- [RATIFIED JM 2026-11-02 19:20]
