# RESULT TEST 5 -- globular-cluster nb fraction under the frozen law (2026-10-18 22:15, Claude, [UNRATIFIED]; JM 22:04 "ok")
Script test5_gc_nb_fraction_20261018.py; report test5_report_20261018.txt. GC parameters [mem, Harris-like; VERIFY]. MW isothermal v_c 200.
rho_c for a GC is NOT fixed by the law; bracketed by 10^7.45 (SPARC capped median) and 10^8.0 (dSph). Cap S <= 2 on the host ambient boost.
Bar (declared before run): nb/M* inside r_h < 0.1.
| GC | r_h pc | R kpc | r_J pc | nb/M* (<r_h) 7.45 / 8.0 | nb/M* (<r_J) 7.45 / 8.0 | stars=nb crossover pc |
| 47 Tuc | 4.2 | 7.4 | 111 | 1.1e-4 / 1.8e-4 | 1.05 / 1.58 | 36 / 28 |
| omega Cen | 7.6 | 6.4 | 173 | 1.7e-4 / 2.5e-4 | 0.99 / 1.39 | 63 / 49 |
| NGC 2419 | 24 | 90 | 639 | 3.7e-3 / 1.3e-2 | 22 / 41 | 74 / 57 |
| Pal 5 | 20 | 18 | 56 | 0.20 / 0.51 | 2.4 / 5.7 | 27 / 19 |
(sigma from GC BTFR ~3-11 km/s; r_c 110-840 pc >> r_J -> nb ~uniform at rho_c inside r_J.)
READING
1. CORES PASS by 3-4 dex for dense GCs: the cap removes the FRM-4 EXPONENT blow-up (omega Cen boost ~1e8 uncapped -> nb/M* ~2e-4).
2. GC vs UFD: same nb density (~0.03-0.1 Msun/pc^3) is negligible against GC stellar densities (~1e3 Msun/pc^3) and dominant in UFDs -> the
   GC/dwarf dichotomy follows from STELLAR density at a common rho_c (consistent with GAL-THRESH). Candidate framing, not yet a claim.
3. RISK -- outskirts and diffuse GCs: IF a GC carries a dSph-level rho_c out to its Jacobi radius, nb ~ M* inside r_J for 47 Tuc / omega Cen and
   >> M* for NGC 2419 and Pal 5. Pal 5 FAILS the core bar (0.20-0.51). Pal 5's thin, cold tidal tails and stars-only mass estimates [mem] would
   be hard to keep. Bound from Pal 5 (<r_h < 0.1): rho_c(GC) <~ 10^7.15.
4. => rho_c cannot be universal down to GCs at the dSph value. GCs need rho_c at least ~2-10x below dSph/SPARC values (or no own nb body at all).
   Back to the one open question: what sets rho_c. Crater II (needs 10^6.1) points the same way for diffuse, tidally affected systems.
PREDICTION (if any GC nb exists): flattening of the outer velocity-dispersion profile beyond the crossover radius (~20-75 pc). Prior art on flat
   GC outer dispersions (Scarpa et al.; tides / escapers / binaries alternatives) [mem, VERIFY].

## ADDENDUM 2026-10-19 00:30 -- R9 VERIFY (Claude Fable 5.1; JM 00:04 "do the verify for r9") [UNRATIFIED]
Source: Baumgardt GC catalogue, combined_table.txt (people.smp.uq.edu.au/HolgerBaumgardt/globular/; Baumgardt & Hilker 2018, Baumgardt+2019 arXiv:1908.02778), fetched 2026-10-19.
| GC | M used / cat (Msun) | r_h used / cat rh,l (pc) | R used / cat R_GC (kpc) | cat M/L_V | cat sig0 km/s | cat rho_h,m Msun/pc^3 |
| 47 Tuc | 7e5 / 8.53e5 +/- 0.05e5 | 4.2 / 4.03 | 7.4 / 7.52 | 1.87 +/- 0.09 | 11.9 | 10^2.61 |
| omega Cen | 3.5e6 / 3.94e6 +/- 0.17e6 | 7.6 / 7.56 | 6.4 / 6.50 | 2.80 +/- 0.10 | 18.0 | 10^2.62 |
| NGC 2419 | 9e5 / 7.83e5 +/- 1.2e5 | 24 / 19.76 | 90 / 95.93 | 1.61 +/- 0.24 | 5.1 | 10^0.71 |
| Pal 5 | 1.5e4 / 1.34e4 +/- 0.24e4 | 20 / 20.42 | 18 / 17.27 | 1.66 +/- 0.30 | 0.7 | 10^-1.12 |
All inputs within 20% of catalogue (r_h used = projected half-light, matching the Plummer a = r_h/1.305 convention). RERUN with catalogue values
(scratch t5v.py, cloud sandbox): nb/M*(<r_h) at 10^7.45 / 10^8.0: 47 Tuc 8.0e-5 / 1.3e-4; omega Cen 1.4e-4 / 2.1e-4; NGC 2419 2.4e-3 / 8.3e-3;
Pal 5 0.25 / 0.62. Crossover 37/64/67/27 pc (7.45). READINGS 1-2 UNCHANGED (cores pass by 3-4 dex). Pal 5 fail stands.
CORRECTION 1 (bound): the stated bound rho_c(GC) <~ 10^7.15 was a LINEAR extrapolation of the 0.20 and ignored the e^2 x rho_amb ambient-boost
floor in the script. Decomposition for Pal 5 (catalogue params): own-core 0.149 + ambient-boost 0.098 = 0.247 at 10^7.45; at 10^7.15 own 0.075 +
floor 0.098 = 0.17; at 10^6.5 total 0.11. WITH the e^2 boost Pal 5 fails the 0.1 bar at ANY rho_c (floor alone ~0.10, rho_amb(17 kpc) = 10^6.39).
WITHOUT the boost (host halo passes through the GC unboosted; floor 0.013) the bound is rho_c(GC) <~ 10^7.2. Which applies is exactly the NESTING
RULE owed under R7 (is a GC inside the MW a body whose baryon bump boosts the host's nb, or does the host nb pass through unchanged?). The Pal 5
bound is therefore CONDITIONAL on the nesting rule. Reading 4's '2-10x below dSph' survives either way; the specific 10^7.15 does not.
CORRECTION 2 (how hard the fail is): Pal 5 catalogue mass error 18% (M/L 1.66 +/- 0.30, N-body fit to kinematics + photometry = dynamical).
An nb excess of 0.25 is ~1.4 sigma against that error; 0.62 (dSph rho_c 10^8.0) is ~3 sigma. The 0.1 bar is stricter than the Pal 5 data.
HONEST STATEMENT: Pal 5 excludes a dSph-level rho_c (10^8) for GCs at ~3 sigma; the SPARC-median value (10^7.45) is in ~1.4 sigma tension only.
47 Tuc / omega Cen (M/L errors 5%) carry the clean core pass; Pal 5 carries only a soft bound. R10 (GAL-THRESH) rests on this soft bound -> PARK confirmed.
SCARPA VERIFY: Scarpa, Marconi & Gilmozzi 2003 (astro-ph/0305388): omega Cen dispersion flattens below a ~ 2.1 (+0.5/-0.4) e-8 cm/s^2; M15 reanalysis
~1.7e-8; authors call it not conclusive (tides, DM halo alternatives). Counter-literature (Lane+2010 Newtonian fits; Baumgardt+2005 tides) still [mem].
GPRF prediction (flattening beyond the crossover ~40-80 pc for omega Cen) is in the same radial range as Scarpa's claim -- prior art to cite, not evidence.
R9 RECOMMENDATION revised: RESTATE -- 'GC cores pass by 3-4 dex; Pal 5 excludes dSph-level rho_c for GCs at ~3 sigma; sharper bound awaits the nesting rule.'
