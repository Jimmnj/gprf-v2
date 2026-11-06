# RESULT_TWOCENTRE_B_20261102 -- the host's potential inside the satellite's own relation: NO EQUILIBRIUM  [RATIFIED JM 2026-11-02 19:20]
# 2026-11-02 15:40 EST. Prereg: PREREG_TWOCENTRE_B_20261102.md (JM 14:36 "let's try it"). Script twocentre_b_20261102.py (continuation in host strength hostf 0..1 from the spherical solve; Picard damped 0.2, Newton-Krylov check; mass-weighted residual; the exp clip ceiling counts as divergence). Outputs twocentre_b_<tag>_20261102.json.

## 0. VERDICT
BD3, and sharper than BD3: the satellite's own relation (sigma_sat) with the host's potential drop in phi_b has NO hydrostatic solution at the real host strength -- not for the LMC, not for Fornax. The equilibrium exists only up to ~20% of the Milky Way's potential gradient for the LMC (sigma 58.8 at 50 kpc) and below 5% for Fornax (sigma 17.3 at 147 kpc). Before it is lost, the host term REDUCES the LMC's mass (4.47e10 -> 2.91e10 on the grid at hostf 0.2): the far side empties, the near side piles up, the density peak slides off the LMC's centre toward the Milky Way, and at hostf ~0.25 the pile-up detaches and runs away. B2 (host nb in the uncapped self-term) tracks B1 and fails the same way.

## 1. Numbers (LMC, L 24 / 64 mu / 800 r; identical to the 16/48/500 run to 3 figures)
  hostf   converged   total excess   near hemisphere   far hemisphere   peak density location
  0.00    yes         4.47e10        6.6e9             3.81e10          centre (far-side node, 0.5 kpc)
  0.05    yes         4.02e10        9.7e9             3.05e10          centre
  0.10    yes         3.47e10        1.37e10           2.10e10          centre
  0.15    yes         3.00e10        1.90e10           1.10e10          0.3 kpc toward MW
  0.20    yes         2.91e10        2.38e10           5.3e9            0.9 kpc toward MW
  0.25    NO          --             (2.5e10 growing)  (3.0e9)          2.4 kpc toward MW and moving
Fornax (M_b 2.24e7, D 147): hostf 0.00 converges (2.2e9 on the grid); 0.05 diverges at once. Sculptor stage 0 under-resolved on the quick grid (not pursued; moot).

## 2. Controls
hostf 0 reproduces the spherical two-halo solve (3.07e10 within 30 kpc in the 1D shoot; 2D grid total 4.47e10 includes the tail to 45 kpc). The clip-at-zero convention and the position-dependent ambient change the stage-0 mass by < 2% vs uniform ambient. Resolution: L 16 vs 24, 500 vs 800 radial points: identical loss point (0.20 < f_crit < 0.25).

## 3. Reading (this is the result)
- The 9/22 picture read the nb mass as a compressible medium with pressure sigma^2 rho in the TOTAL potential. For a satellite that medium is a 59 km/s gas sitting in a 200 km/s potential gradient: it cannot be held. The relation's "drop from the centre" form with the host included is the hydrostatic condition for that gas, and the gas falls into the Milky Way. The scale is the Jacobi-like radius sigma_sat^2 D / v_c^2 ~ 4 kpc for the LMC, ~1 kpc for Fornax: beyond it the host's drop exceeds sigma_sat^2 and no equilibrium exists. The density-match truncation at 31 kpc never saw this because it never put the host's potential into the satellite's relation.
- So JM's hypothesis (something unique for a large satellite) has a definite answer in the gas reading: the host term is not a correction, it is fatal, and not only for the LMC -- for every satellite, more so the smaller its sigma. It cannot be the missing LMC mass.
- What survives is the OTHER reading, which is the one the protocol has used all along: eq. (1) as a RESPONSE RULE -- the body's nb mass is the response to the body's own phi_b with the body's own self-term, the host entering only as ambient density. In that reading the host's potential does not appear in the satellite's relation at all, there is no instability, and nothing is unique about a large satellite beyond the 1.67 nesting factor already measured. This is consistent with the canon stance that the nb mass is "re-formed in place rather than carried" (2026-10-24) and carries no momentum: a response has no hydrostatics to lose. It is INCONSISTENT with the 9/22 "compressible medium" language still in SNAP_0922 sec 0c and with any sentence in the note that calls the nb mass a pressure-supported medium in the total potential.
- Consequence for the ONE-FIELD ruling: (A') stands as the only working reading, and it should be stated as: a body's relation is solved in its OWN phi_b; the host is ambient density; nested components are in the body's phi_b. The phrase "phi_b = sum over all stars and gas" is to be scoped to the body's own baryons.
- The LMC stays where LIT_LMC_MASS left it: 1.6x short at 33 kpc, halo ends at 31 kpc vs detected to >= 53. Conceded.

## 4. Ratification
  TB1 Spherical control reproduced at hostf 0 (T2)
  TB2 No equilibrium for the LMC above hostf ~0.2-0.25, Fornax above 0.05; host term reduces the mass before the loss (T2)
  TB3 BD3; hypothesis "unique for a large satellite" NOT supported in the gas reading; nothing to add in the response reading (reading)
  TB4 Reading: eq. (1) is a response rule, not a medium in the total potential; "compressible medium" language retired; (A') scoped to the body's own phi_b (reading; JM rules; canon-level)
  TB5 LMC concession unchanged from LIT_LMC_MASS sec 4 (reading)
# END [RATIFIED JM 2026-11-02 19:20]
