# RESULT_CLUSTER_RELEASE_20261022 [UNRATIFIED -- order of magnitude]
cluster_release_20261022.py; report cluster_release_20261022_report.txt. JM 12:06 'sounds promising' (merger/stripping excess returned to the environment).
Setup: Coma-like cluster, 1357 member galaxies 1e8-3e11 Msun (Schechter alpha -1.1, M* 1e11 stellar, total stellar 1.5e13 VERIFY), isothermal cluster s = 1000 km/s,
R200 = 2 Mpc. Each galaxy: isolated nb M_iso = (4/3) sigma^2 r_t/G; in the cluster, tail cut where it meets the cluster density; released = difference, kept in the cluster.
| rho_amb | nb per 1e11 galaxy (isolated) | cluster sum isolated | retained on galaxies | RELEASED | released / needed (~8.5e14, VERIFY) |
|---|---|---|---|---|---|
| 35.5 (closure) | 2.8e13 | 5.5e15 | 1.7e14 | 5.4e15 | 6.3 |
| 150 (fit) | 1.4e13 | 2.7e15 | 1.7e14 | 2.5e15 | 3.0 |
Half of the budget comes from galaxies below ~4e10 Msun (M_nb ∝ M_b^3/4 favours small bodies).
READING: returning the excess to the environment turns the cluster regime from a deficit ('harder than MOND') into a SURPLUS of 3-6x. The surplus is the same
problem as the known Mpc-scale budget (3.7 rho_m at closure; sec 7): the isolated tails out to r_t (2-5 Mpc) are too massive in total. Matching Coma needs M ∝ rho_amb^-1/2 lower by
3x -> rho_amb ~ 1300, above the lensing 95% bound (630). Possible outs (untested): field tails overlap (L* r_t > typical separations), so summing isolated masses double-counts;
or a group/cluster body rule replaces member tails. Galaxies' Lagrangian region for 1e15 Msun is ~18 Mpc, so the whole of each tail is swept in -- no out from geometry.
VERDICT: promising in sign, too large in size; the binding constraint is the Mpc budget, not the cluster. Not for the note beyond one Open-list clause.

## Addendum 12:20 -- released nb as the cluster's OWN extended halo (JM: "could they extend the cluster's own halo")
Arrange the released mass as an isothermal halo of the cluster body, rho = s^2/(2 pi G r^2) (mass conserved by construction; dM/dr = 2 s^2/G):
| s (km/s) | source of s | M(<2 Mpc) | halo extends to (released 2.5e15 / 5.4e15) | rho at 5 / 10 Mpc |
|---|---|---|---|---|
| 1000 | observed Coma dispersion | 9.3e14 | 5.4 / 11.6 Mpc | 36 / 9 rho_m |
| 879 | BTFR sigma from cluster baryons incl. ICM (1.5e14) | 7.2e14 | 7.0 / 15.0 Mpc | 28 / 7 rho_m |
| 494 | BTFR sigma from stars only (1.5e13) | 2.3e14 | 22 / 48 Mpc | 9 / 2 rho_m |
READING: with the cluster's own BTFR sigma from its TOTAL baryons (879 km/s), the inner 2 Mpc holds 7.2e14 vs ~8.5e14 needed (-0.07 dex) with no free parameter, and the
surplus sits outside R200 out to 7-15 Mpc. That removes the 3-6x surplus INSIDE the cluster. Costs/flags: (1) using the ICM in the cluster's sigma conflicts with the source rule
(ICM does not source) unless 'arranging released nb' is distinguished from 'sourcing' -- a new rule, flag as ad hoc until justified; (2) the outskirts: 28 rho_m at 5 Mpc,
7 at 10 Mpc -- stacked cluster weak lensing and cluster-galaxy correlations measure outskirts to 10-30 Mpc (DES, HSC; VERIFY numbers): a direct test, likely denser than
observed at 5 Mpc (1-halo + 2-halo typically ~5-15 rho_m there, VERIFY); (3) why s(cluster) = BTFR is assumed, not derived; (4) the Mpc budget (3.7 rho_m) is moved, not removed.

## Addendum 13:05 -- overshoot vs undershoot reconciled (JM: "you say it overshoots - but you now say we still have the mond undershoot")
They are different radii. The law's own halo for the cluster body (sigma 879 km/s from total baryons 1.5e14) has total (4/3) sigma^2 r_t/G:
- rho_amb 35.5: r_t 28.4 Mpc, total 6.8e15; galaxies release 5.4e15 -> released/law = 0.79
- rho_amb 150: r_t 13.8 Mpc, total 3.3e15; galaxies release 2.5e15 -> released/law = 0.76
So the nb the member galaxies give up is ~3/4 of what the cluster body's halo holds out to its own r_t: conservation and the law are nearly self-consistent, at both ambient values. The 3-6x 'overshoot' was only against the mass needed INSIDE R200; the law puts most of the cluster halo at 2-28 Mpc.
Radial zones (isothermal at 879 km/s): < ~1 Mpc: 3.6e14 -- MOND's known inner deficit (~2x) remains, a diffuse core worsens it; 1-2 Mpc: total 7.2e14 inside 2 Mpc vs ~8.5e14 needed (parity); 2 Mpc to r_t: the released/extended halo, 28 rho_m at 5 Mpc -- the test against stacked cluster lensing.
