# RESULT_DWARF_PREREG_20261020 [UNRATIFIED] -- the pre-registered dwarf test, first pass on measured out-of-sample dwarfs. Script dwarf_prereg_20261020.py; report; csv.
RELATION (published here, from SPARC TEST1 tanh2, N=120): log rho_c = 11.063 - 0.365 log M_b + 0.268 log Sigma_h - 0.929 f_gas (Sigma_h = M_b/2pi R_h^2 in Msun/pc^2 [unit CORRECTED 23:15; script uses pc]; resid 0.64 dex).
PREDICTION: N9 solver (Plummer stars, M_b = 1.6 L_V, sigma^2 = sqrt(G a0 M_b)/2, host ambient), rho_c from the relation at f_gas = 0. NO per-object freedom.
DECLARED: per-object pass |dlog sigma| <= sqrt(0.10^2 + err^2); sample survives if >= 2/3 pass and embedded/isolated residuals agree in sign. Inputs [mem, VERIFY: McConnachie 2012 etc.].

## Result: 7/15 PASS -> FAIL by the declared rule. rms log(pred/obs): GPRF 0.179 (mean +0.02) | MOND iso (crude) 0.212 (-0.11) | MOND EFE (crude) 0.293 (-0.21) | Newton 0.651.
PASS: Sculptor +0.01, Carina 0.00, Sextans -0.06, Leo I +0.07, Leo II +0.04, And XIX 0.00, And XVIII -0.11.
FAIL: Fornax +0.15, And II +0.31 (massive, extended: over); Crater II +0.28, And XXI +0.25 (diffuse: over); Draco -0.17, UMi -0.20 (densest: under); Cetus +0.11; Tucana -0.34.
Embedded mean +0.06 (N=12), isolated -0.11 (N=3) -- the EFE-like direction, but carried entirely by Tucana at 15.8 km/s (Fraternali+2009); Taibi+2020 give 6.2 +/- 1.1
[mem, VERIFY] -> with that value the asymmetry vanishes. UNDECIDED pending verification.
## Reading
- A zero-parameter 4-dex extrapolation of a SPARC relation lands unbiased with 0.18 dex scatter; the relation's own 0.10 bar was too strict for that extrapolation.
  The residual pattern (over for extended/diffuse, under for the densest) says rho_c depends on something beyond (M_b, Sigma_h, f_gas) -- the 0.64 dex is real.
- The MOND columns are ONE-LINE estimates (isolated: 81 sigma^4/4Ga0 = M; EFE: quasi-Newtonian Wolf with G_eff = G a0/g_ext, capped at isolated). NOT literature
  predictions; MOND dSph fits in the literature (McGaugh & Milgrom 2013; Lelli+2017) do better than this column. The note may not say 'MOND fails dSphs'.
- For the note (6.6): report the relation, the protocol, this first pass as 'first pass: 7/15 within the bar, rms 0.18 dex, no bias', and list the objects to be
  predicted BEFORE measurement after a literature check (do not populate from memory).

## ADDENDUM 19:30 -- Tucana VERIFIED: Taibi et al. 2020 (A&A; VLT/FORS2, 39 members) sigma = 6.2 +1.6/-1.3 km/s, superseding Fraternali+2009's 15.8. Rerun:
8/15 PASS (still FAIL by the 2/3 rule); GPRF rms 0.157 dex (mean +0.05) | MOND iso 0.172 (-0.09) | MOND EFE 0.265 (-0.18) | Newton 0.620. Embedded +0.06 (N=12) vs
isolated +0.02 (N=3): NO EFE-like asymmetry; the first-pass direction was Tucana alone. Isolated dSphs (Cetus, Tucana, And XVIII) do not yet separate GPRF
from isolated MOND (both within errors on Tucana and And XVIII). Sign test: PASS (same sign); object count: FAIL. Report and csv regenerated.

## ADDENDUM 20:25 -- MW classical dSph inputs replaced by Walker+2009 Table 1 (VERIFIED). 8/15 pass; GPRF rms 0.158 (mean +0.035); MOND iso 0.186 (-0.10); MOND EFE 0.292;
Newton 0.651. Embedded +0.04 vs isolated +0.02. Same fails: Fornax +0.12, And II +0.31, Crater II +0.28, And XXI +0.25 (over); Draco -0.18, UMi -0.22 (under).
Remaining [mem]: Cetus, And II, And XXI, And XVIII. Report/csv regenerated.

## ADDENDUM 21:50 -- M31 dwarfs + Cetus VERIFIED (McGaugh & Milgrom 2013 I/II tables; Taibi+2018). And II stands (9.3e6; 7.3+/-0.8); And XXI r_h -> 0.77 kpc proj. (+0.19,
still over); And XVIII 9.7+/-2.3 (Tollerud) is DISPUTED (Collins+2013 <= 2.7 from 4 stars) -- kept, flagged; Cetus = 11.0 +1.6/-1.3 (Taibi+2018, 54 members) -> PASS (-0.01);
Kirby+2014 8.3 [mem] as alternate row (FAIL +0.11). Totals (15 objects + Cetus alt): 9/16 pass; GPRF rms 0.148 (mean +0.03); crude MOND iso 0.183, EFE 0.288, Newton 0.654.
Embedded +0.03 vs isolated +0.01. All inputs now verified except DF-class objects (not in this table) and Kirby's Cetus value. McGaugh & Milgrom's own MOND predictions
(a priori): And XXI iso 5.2 / EFE 3.7 (obs 4.5); And XVIII iso 5.6 (obs 9.7 or <= 2.7); And II iso 11.0 (obs 7.3) -- use THESE, not our crude column, for M31 dwarfs.
## ADDENDUM 23:40 -- McGaugh & Milgrom 2013 II Table 1 (fetched 10/20): And XIX L_V 4.1e5, r_1/2 2244 pc, sigma_iso 5.0 +1.0/-0.8, sigma_efe 2.6 +1.1/-0.8, obs 4.7 +1.6/-1.4 (Collins 2013)
[we use Collins+2020: 7.8 +1.7/-1.5, r_h 3.1 kpc]. Crater II dispersion saturates above rho_c ~10^7.5 near sigma_law 5.7 (RHOC scan 7.5:5.2, 7.75:5.3, 8.0:5.4) -- hence 5.1-5.4 at 10^8.7.
