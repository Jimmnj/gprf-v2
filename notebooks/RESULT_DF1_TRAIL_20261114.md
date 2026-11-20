# RESULT_DF1_TRAIL_20261114 [UNRATIFIED] -- NGC 1052-DF2, DF4 and DF9 at the June 2026 distances under ONE-FIELD SCOPED (A')
Queue item 3 (DF-1). JM 2026-11-14 14:43 "yes". Claude Opus 5.5, 14:50-15:05 EST.
Scripts: df1_trail_20261114.py (stages rep, main, alt, appB); df1_sat_20261114.py (saturation check). Report: df1_trail_20261114_report.txt.
Tables: df1_trail_20261114_main.csv, df1_trail_20261114_alt.csv, df1_trail_20261114_appB.csv. Canon, tex (md5 baecfbc1) and table of record untouched.

## 1. Protocol (rule 3a)
- Machinery verbatim: n9_dwarfs_20261003.py up to 'def run(' (solve, sigma_los, host_rho, rho_bg), exec'd as in dwarf_prereg_20261020.py,
  onefield_dwarfs_20261102.py and rhoc_sigma_b_oos_20261019.py. This is the note's dwarf protocol (Appendix B, sec 4.3): Plummer stars, sigma from eq. (2),
  isotropic Jeans, luminosity-weighted projected dispersion inside the aperture.
- A': each galaxy is solved in its own phi_b with its own sigma. The host enters only as ambient density. FIELD = rho_bg = 35.5 Msun/kpc^3.
  HOST bracket = NGC 1052 at 250 km/s and 80 kpc, plus rho_bg (the old protocol's value; an upper bound on the ambient).
- M_b = M* (gas-free; DF2 M_HI/M* < 0.016 at 3 sd). M* scales as D^2, R_e as D.
- rho_c: median 10^7.42 (V3 first minimum, fixed arm, 141 galaxies; RESULT_EP3_20261114). Rule (3) refit verbatim:
  log rho_c = 11.063 - 0.365 log M_b + 0.268 log Sigma_h - 0.929 f_gas (the note's equation to the digits).
  Trend-only (M_b, Sigma_h), as in rhoc_sigma_b_oos: 8.349 - 0.153 log M_b + 0.370 log Sigma_h.
- Inputs [T1; LIT_TRITONSTATION sec 8t; Danieli+2019 re-checked 2026-11-14 14:5x]:
  DF2: M* 2.0e8 and R_e,c 2.0 kpc at 20 Mpc (van Dokkum+2018; Cohen+2018 via Danieli+2019, whose r_1/2 = 4/3 R_e,c = 2.7 kpc).
       sigma 8.5 +2.3/-3.1 (Danieli+2019, ApJL 874, L12; arXiv v1 gave 8.4 +- 2.1). D 17.6 +- 0.6 (JWST TRGB, Tang+2026); alternative 21.7 (HST TRGB).
  DF4: M* 1.5e8 (2 x 7.5e7 inside r_1/2) and R_e 1.6 kpc at 20 Mpc (Shen+2023). sigma 8.0 +2.3/-1.9 (Shen+2023).
       D 18.8 +- 1.7 (SBF, Tang+2026); alternative 20.0 (TRGB, Danieli+2020).
  DF9: M* 1.42e8 (2 x 0.71e8 inside r_e) and R_e 1.1 kpc at 20.6 Mpc (Keim+2026). sigma 6.4 +4.0/-4.3 (Keim+2026).
       D 18.7 +- 1.5 (SBF, Tang+2026); alternative 20.6 (Keim's assumed value).
  Input check: the Wolf stars-only dispersion for DF9 at 20.6 Mpc is 8.3 km/s here and 8.3 in Keim+2026.

## 2. Reproduction before use
- DF2, old inputs (2e8, rh 2.7, host): 22.45 km/s at 10^7.45 (record 22.4). Required log rho_c 5.82 [nan, 6.22] (record identical).
- DF4, old inputs (1.5e8, rh 1.6, host): 17.31 km/s (record 17.3). Required 5.44 [nan, 6.52] (record identical).
- Appendix B, 16 rows: largest difference from dwarf_prereg_20261020.csv is 0.000 km/s.
- NOT reproduced: the 10-20 rerun's "stars-only Newton 6.6 / 7.8". That script was not saved. The value depends on how the ambient is subtracted
  when rho_c -> 0: 6.26 / 7.71 with the host ambient subtracted; 7.07 / 7.96 for clean Plummer Jeans. The note does not quote it. Record caveat only.

## 3. Results (convention C, see F1; FIELD ambient; primary distances)
| | D (Mpc) | M_b | R_e (kpc) | eq. (2) sigma | stars only | isolated MOND | at 10^7.42 | observed | required log rho_c [lower, +1 sd] | 7.42 - required |
|---|---|---|---|---|---|---|---|---|---|---|
| DF2 | 17.6 | 1.55e8 | 1.76 | 28.0 | 6.9 | 18.7 | 20.1 | 8.5 +2.3/-3.1 | 5.91 [<4.5, 6.42] | 1.51 (>= 1.00 at +1 sd) |
| DF4 | 18.8 | 1.33e8 | 1.50 | 26.9 | 6.9 | 18.0 | 18.3 | 8.0 +2.3/-1.9 | 5.86 [<4.5, 6.48] | 1.56 (>= 0.94) |
| DF9 | 18.7 | 1.17e8 | 1.00 | 26.1 | 8.0 | 17.4 | 14.9 | 6.4 +4.0/-4.3 | below stars only [<4.5, 6.74] | >= 0.68 at +1 sd |
(dispersions in km/s; lower bounds are open because each observed value minus 1 sd lies below the stars-only dispersion.)
- Rules: the trend-only rule predicts 7.42-7.58, i.e. the median's values (15-20 km/s). Rule (3) with f_gas = 0, the Appendix B protocol, predicts 8.31-8.46
  (24-27 km/s).
- SPARC context (ep3_val_20261114_part.csv, the 141): lowest fitted value 6.17 (UGC07125); 2.5th percentile 6.59. DF2 and DF4 need less than every one of
  the 141 (3.2 and 3.4 robust spreads below the median). At +1 sd (6.42, 6.48) one of 141 qualifies (UGC07125). DF9's bound (6.74) is met by 9 of 141.
- Ambient (A'): the HOST bracket raises the required value by 0.11-0.13 dex (DF2 6.02, DF4 5.99; +1 sd bounds 6.46, 6.52, 6.76). At fixed rho_c it LOWERS
  sigma by about 1%: truncation at r_t ~ 13 kpc instead of ~0.9 Mpc in the field. The environment is not a lever.
- Distance: DF2 at 21.7 Mpc gives 23.3 km/s at the median and needs 5.45 (1.97 dex below). The JWST distance narrows the gap by 0.46 dex.
  DF4 at 20.0 needs 5.71 (1.71 below). DF9 at 20.6: bound 6.57.
- Saturation: the n9 solver uses S(x) = x; the frozen R1 uses S = 2 tanh(x/2). With x_b0 = 0.48-0.74, the frozen form moves sigma at the median by
  +0.2 to +0.5% and the required rho_c by <= 0.003 dex (df1_sat_20261114.py; its identity-S variant reproduces n9 to 0.001 km/s). Immaterial here.
- Aperture: DF2's KCWI field covers R <= 0.7 R_e. With that aperture: 19.99 vs 20.07 km/s at the median; required 5.91 [<4.5, 6.44]. Immaterial.

## 4. Findings to flag before referees
(F1) RADIUS CONVENTION in the dwarf solver. n9 sets the Plummer scale a = rh/1.305, which reads rh as the 3-D half-light radius. The protocol feeds
  projected radii: Walker+2009 r_half; the M31 3-D values converted x3/4 to projected (dwarf_prereg header); DF4's R_e 1.6. The same scripts treat rh as
  projected elsewhere (Wolf r_1/2 = 4/3 rh; Sigma_h). Those stellar distributions are 1.305 times too compact. DF2's old input (Danieli's r_1/2 = 2.7 kpc,
  3-D) was consistent by accident. C = the note's stated protocol: a = R_e, aperture R <= R_e. P = as coded.
  Effect on Appendix B under C: every prediction rises 1-9%; pass 9/16 unchanged; rms 0.148 -> 0.152 dex; mean +0.028 -> +0.051 (embedded +0.056,
  isolated +0.038). Verdict and pattern unchanged. DF4's required value moves from 5.50 (P) to 5.86 (C), both at 18.8 Mpc.
  Options: (a) fix the one line (a = R_e) and carry the new Appendix B numbers into the v2.0.2 change list; (b) keep it and state the convention. Claude: (a).
(F2) Danieli+2019's journal version gives sigma = 8.5 +2.3/-3.1. The value 8.4 +- 2.1 is arXiv v1; the note and the 10-20 VERIFY LEDGER quote v1.
  Effect on the required rho_c: 0.03 dex. Cite the journal value.
(F3) The note's sec 4.3 says DF2 and DF4 need rho_c "about 1.5 and 1 dex" low. That paired DF2's central value with DF4's +1 sd bound.
  Like for like: 1.5 and 1.6 dex at the measured dispersions; at least 1.0 and 0.9 dex at +1 sd.
(F4) DF9's measured dispersion (6.4) is below its stars-only value (8.0 here; 8.3 in Keim+2026 at 20.6 Mpc). It bounds rho_c from above only.
(F5) CORRECTION to ledger (am') SCORECARD (Claude, 13:11): "the ambient term enhances" is wrong. Under A' the ambient is subtracted; it truncates.
  At 80 kpc from NGC 1052 the effect inside R_e is about 1% (0.1 dex in the rho_c needed). The clause's conclusion stands: a nearby normal galaxy does not
  help the relation.

## 5. READ [T3]
- At the SPARC median the relation gives 15-20 km/s, about what isolated MOND gives (17-19). The measured 6-9 km/s need rho_c about 1.5 dex below the
  median, below every SPARC galaxy's fitted value. With rho_c free the relation reproduces them; it does not predict them ((am') addendum 2).
- Under A' there is no environmental lever (0.1 dex) and no gas lever. MOND's lever is the external field, and the distances now disfavour the geometry it
  needs (8t). Neither model predicts the three. The bullet-dwarf account predicted DF9's low dispersion before it was measured (Keim+2026).
- The closer DF2 distance helps the relation by 0.46 dex. That is not enough.
- No new physics. This updates sec 4.3's numbers and adds DF9.

## 6. Owed
- JM: rule F1 (a or b) and F2. Then sec 4.3 in the tex pass (ledger (am') addendum 3 carries the numbers and the text).
- The apertures of the DF4 and DF9 measurements are not checked. DF2 shows the aperture is immaterial.

## 7. ADDENDUM 16:10 -- outer profile and globular clusters at the fitted rho_c (df1_gc_20261114.py) [UNRATIFIED]
> **JM 2026-11-14 16:03:** "even if the rho_c is a fitted number does the rest of the curve fit - are they eliptical"
> **JM 2026-11-14 16:04:** "df2 4 9"
- Morphology [T1]: round, gas-free, dispersion-supported ultra-diffuse spheroids. DF2: Sersic n 0.6, b/a 0.85 (van Dokkum+2018). DF4: n 0.79,
  b/a 0.89 (van Dokkum+2019). DF9: the closest analog of the two in luminosity, size and clusters (Keim+2026). No disc and no gas, so no rotation curve.
  DF2 shows a weak prolate-like rotation, about +-6 km/s (Emsellem+2019); Danieli+2019 saw none inside the KCWI field.
- Data beyond the stellar dispersion [T1]: DF2 GCs 7.8 +5.2/-2.2 (10 GCs; van Dokkum+2018 RNAAS) and 10.6 +3.9/-2.3 (7 GCs + 3 PNe; Emsellem+2019);
  GCs out to 7.6 kpc at 20 Mpc; GC half-number radius 1.4 R_e (Wasserman+2018). DF4 GCs 4.2 +4.4/-2.2 (7 GCs), out to >= 7 kpc. DF9: two clusters,
  7 +23/-6 (Keim+2026).
- At the fitted rho_c (convention C, field): core radius 12.6 kpc (DF2), 12.8 (DF4). nb/stars inside r: DF2 0.3 at R_e, 1.0 at 2 R_e, 5.0 at 4 R_e,
  11 at 10 kpc; DF4 0.2, 0.7, 3.4, 11. Field r_t about 0.9 Mpc; nb inside r_t 2.3e11 (DF2), 2.0e11 (DF4).
- Predicted GC dispersion inside R <= 6.7 kpc (DF2) and 6.6 kpc (DF4), Plummer tracer a_GC = 1.0 / 1.5 / 2.0 R_e:
  DF2 9.0 / 10.4 / 11.6 km/s (stars only 6.3-6.4); DF4 8.2 / 9.4 / 10.4 (stars only 6.2-6.3).
  At the +1 sd rho_c bounds: DF2 11.9 / 14.1 / 15.9; DF4 11.3 / 13.4 / 15.0.
  DF2: the fit lies inside 1 sd of both GC values. DF4: about 1.2 sd above 4.2 +4.4 at a_GC = 1.5 R_e; the +1 sd bound about 2 sd above.
  The GCs do not separate the fit from stars only.
- Circular speed at the fitted rho_c, DF2: 13 / 19 / 28 / 42 km/s at 2 / 5 / 10 / 30 kpc. The BTFR speed for M_b 1.55e8 is 40 km/s. Stars only: about
  5 km/s at 30 kpc.
- READ [T3]: the fit has zero degrees of freedom (one stellar number). Lowering rho_c does not remove the nb; it spreads it into a ~13 kpc core.
  The tail is set by eq. (2) and the ambient, not by rho_c, so beyond ~15 kpc the relation's DF2 is back on the BTFR. In the relation these galaxies
  are poor in nb inside ~2 R_e only. The outer mass is the discriminator: the relation gives the BTFR value; the bullet-dwarf account gives none.
  Current data stop near 4 R_e. Kumar+2024 (A&A) reach the same point for CDM: a cored halo with a 20 kpc scale length fits DF2's kinematics, and the
  mass models are largely unconstrained by the data.
- Not computed: the tidal-distortion constraint (Keim+2022: halo masses not much above the stellar masses, if the distortions are tides from NGC 1052).
  At the fitted rho_c the mass inside 7 kpc is about 6x the stars, which would need a closer pericentre. The new distances put DF2 1.8-3.8 Mpc in
  front of NGC 1052, which also weakens NGC 1052 as the tidal source.
- 17:05 (JM "can we say the gc's support the fit?"): likelihood ratio fit / stars-only, split-normal on each GC value, a_GC = 1.0 / 1.5 / 2.0 R_e [T2]:
  DF2 vD+2018 (10 GCs) 1.23 / 1.08 / 0.94; DF2 Emsellem+2019 (7 GCs + 3 PNe) 4.5 / 5.3 / 5.1 (stars only 1.8 sd low); DF4 vD+2019 (7 GCs) 0.73 / 0.56 / 0.42.
  The two DF2 sets share five clusters, so they cannot be multiplied. Net: no robust preference. ANSWER: no; "consistent with both", not "support".
  The same ~10 DF2 velocities gave 7.8 (likelihood) and 5.6 (ABC) in vD+2018 RNAAS, 14.3 +- 3.5 (Martin+2018) and 2.4-18.8 at 2 sd (Haghi+2019).
  A decisive test needs the ~4 km/s gap at 2-4 R_e measured to ~1.3 km/s: roughly 25-40 tracer velocities with 3-5 km/s errors [T2]. The gap grows
  outward (circular speed at 10 kpc: 28 km/s fitted vs 8 stars only), so the outermost tracers carry the most weight.
- 17:15: nb mass inside 50 / 100 / 200 kpc at the fitted rho_c (field): DF2 2.2e10 / 4.1e10 / 7.2e10 (265x M* inside 100 kpc); DF4 2.1e10 / 3.8e10 /
  6.7e10 (289x); DF9 at its +1 sd bound 1.7e10 / 3.1e10 / 5.9e10 (264x). Same order as a stellar mass-halo mass expectation, with a ~13 kpc core.

## 8. RULINGS
> **JM 2026-11-14 17:34:** "f1 accepted and start 4"  -> DF1-F1 (a).
> **JM 2026-11-14 18:38:** "accept all"  -> DF1-F2 (journal value), DF1-F3 (wording), ledger (am') with addenda 1-6. DF1-F5 was a logged correction.
This file as a record: not yet stamped.
