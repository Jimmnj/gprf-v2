# RESULT_LOCAL_NB_DENSITY_20261101 -- nb density at the Sun from the Milky Way solve (S6 of REVIEW_COLD_v2.0.1)  [UNRATIFIED]
Written 2026-11-01 19:52 (JM "go" 19:36). Same solve as mw_decline_20261021.py (Ou+2024 baryons, spherical response, rho_amb 35.5), density read at R0 = 8.2 kpc. Inline computation; to be frozen as a script if ratified.

## Result
- log rho_c 7.45 (universal): rho_nb(R0) = 0.0082 Msun/pc^3; spherical column |z| < 1.1 kpc ~ 18 Msun/pc^2.
- log rho_c 7.88 (MW fitted): rho_nb(R0) = 0.0100 Msun/pc^3; column ~ 22 Msun/pc^2.
- Reference (from memory, VERIFY): local dark-matter density 0.008-0.013 Msun/pc^3 (Bovy & Tremaine 2012 0.008; McKee+2015 0.013; de Salas & Widmark 2021 review 0.010-0.015); K_z total column within 1.1 kpc 68-74 Msun/pc^2 with baryons 47-55 -> non-baryonic column 13-27.
- Both values PASS the vertical-force range. The hostile reviewer's 0.028 at the fitted rho_c used rho_c/e; at 8 kpc the solution is already past the core (x_b has fallen), so the actual density is 2.8x lower than that estimate. The universal value sits at the low edge of the K_z range, the fitted value in the middle.
- Caveat: spherical response; a flattened nb distribution would raise the midplane density at fixed column. A disc-plane solve is owed (sec 3.1 (ii) midplane sentence already notes the geometry).

## Consequence for the note
One sentence in A.1 (Milky Way): "At the Sun the solve gives 0.008 (universal) to 0.010 (fitted) Msun/pc^3, inside the 0.008-0.013 of vertical-force determinations." Adds a passed check where a referee would otherwise compute a fail. HELD.
