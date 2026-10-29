# RESULT_CRATER2_TIDES_20261021 [UNRATIFIED]
Pre-registered 2026-10-21 17:35 (JM 'go'); crater2_tides_20261021.py; report crater2_tides_20261021_report.txt.
Inputs verified from Ji et al. 2021 (ApJ 921, 32), MW+LMC fiducial orbits: Crater II r_h 1.066 kpc, L_V 10^5.21, sigma 2.34 (+0.42/-0.30) [Caldwell 2017: 2.7+/-0.3], r_peri 21.7 (+5.4/-3.9); Antlia II r_h 2.54, L_V 10^5.88, sigma 5.98 (+0.37/-0.36), r_peri 38.6 (+7.2/-5.8). Ji et al.'s own tidal radius at pericentre: Crater II 386 (+120/-81) pc, Antlia II 1500 pc (their MW+LMC potential). Host = the law's own MW curve (RESULT_MW_DECLINE, universal rho_c): 169 km/s at 117 kpc, 187 at 21.7.

## Result (sigma_los within r_h, km/s; dex vs Ji value; protocol bar)
| case | Crater II, rel.(3) 10^8.70 | Crater II, univ. 10^7.45 | Antlia II, rel.(3) 10^8.43 | Antlia II, univ. 10^7.45 |
|---|---|---|---|---|
| A today, untruncated (formal law prediction) | 5.23 (+0.35) FAIL | 5.17 (+0.34) FAIL | 7.55 (+0.10) PASS | 7.90 (+0.12) FAIL |
| B stars+mass cut at r_J(peri), today's ambient | 3.83 (+0.21) FAIL | 3.10 (+0.12) FAIL (Caldwell: +0.06 PASS) | 5.91 (-0.01) PASS | 6.09 (+0.01) PASS |
| C as B, pericentre ambient (bound) | 3.61 (+0.19) FAIL | 2.78 (+0.07) PASS | 5.65 (-0.02) PASS | 5.80 (-0.01) PASS |
r_J(peri): Crater II 0.65-0.71 kpc (Ji: 0.39); Antlia II 1.8-2.1 kpc (Ji: 1.5). Enclosed mass positive throughout.

## Reading
- The law's equilibrium prediction for Crater II stays a 0.35 dex miss (concession unchanged as a prediction).
- Truncation at the pericentric Jacobi radius takes the universal-rho_c value from 5.2 to 3.1-2.8 km/s (Caldwell passes in B and C, Ji passes in C), and Antlia II, the control, passes in the same cases. So tides are a potential explanation inside the law, by the same mechanism CDM uses (Sanders+2018; Borukhovetskaya+2022; Fu+2019), not a law-specific one.
- relation (3) fails in every case (it predicts 10^8.70 for Crater II, too dense) -- the relation moves the wrong way at the diffuse end, as found before.
- Antlia II: the formal (case A) value is high by 0.10-0.12 dex; tides fix it too. Ji+2021 see a >5 sigma velocity gradient aligned with the orbit (Antlia II IS disrupting), so the control is consistent.

## Caveats (decisive)
1. r_J(peri) = 0.7 kpc is INSIDE today's half-light radius (1.07): a static cut at r_J is inconsistent with the observed photometry; the real system is post-pericentre and out of equilibrium. Cases B/C are bounds, not models.
2. The nb response is instantaneous, so the nb core regrows after pericentre while the stars relax on r_h/sigma ~ 400 Myr; a stripping simulation with the law's instantaneous response is required before any claim.
3. Ji's tidal radius (0.39 kpc) is smaller than ours (0.7): their MW+LMC potential is heavier at 22 kpc; with theirs the truncation would be stronger.
4. Pericentre uncertainty +5.4/-3.9 kpc -> r_J varies ~+/-15%.
## Label
Crater II: 'conceded (equilibrium); potential explanation by tides, shared with CDM, simulation owed'. Antlia II: added to the dwarf table as control (pass with tides; formal value +0.10/+0.12).
