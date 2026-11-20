# RESULT_FEATURE_TRANSMISSION_20261109 -- transmission of baryon features into rotation curves: data vs eq. (1) vs the RAR function (PREREG_FEATURE_TRANSMISSION_20261109, stamped 16:27) [RATIFIED 2026-11-09 22:03 -- JM "ratified"; BALLOT 22:03 line]
Written 2026-11-09 16:45 EST. Script: feature_transmission_20261109.py; table feature_transmission_20261109.csv (131 galaxies). ARM: fixed M/L (Yd 0.5, Yb 0.7), spherical pass2 solve_gal. T_x = signed feature amplitude of curve x divided by the baryon feature amplitude, window set from V_bar alone.

## 1. Sample
131 galaxies (Q <= 2, inc > 30, >= 8 points). Significant baryon feature (A_bar >= 2x V_obs noise, >= 3 points): 43 -- LOW g (< 0.3 a0) 21, MID 9, HIGH (> a0) 13. Scored set LOW+MID: 30.

## 2. Medians with bootstrap 16-84
                     T_obs             T_1 (univ 7.45)     T_1f (fitted)       T_R (RAR fn)
ALL sig (43)         0.21 [0.15,0.26]  0.67 [0.55,0.71]    0.69 [0.54,0.71]    0.64 [0.61,0.66]
SCORED LOW+MID (30)  0.18 [0.11,0.30]  0.47 [0.30,0.57]    0.41 [0.34,0.48]    0.56 [0.55,0.61]
LOW (21)             0.22 [0.15,0.38]  0.18 [0.17,0.32]    0.33 [0.26,0.39]    0.55 [0.54,0.55]
HIGH (13)            0.22 [0.15,0.61]  0.83 [0.81,0.84]    0.88 [0.77,1.09]    0.76 [0.70,0.94]
Paired differences, SCORED set: D_1 = T_obs - T_1 = -0.07 [-0.30,+0.02] covers 0; D_1f = -0.19 [-0.31,-0.05] excludes 0; D_R = -0.36 [-0.46,-0.30] excludes 0.
LOW only: D_1 +0.03 [-0.04,+0.05] covers 0; D_1f -0.02 [-0.19,+0.03] covers 0; D_R -0.33 [-0.39,-0.28] excludes 0.
HIGH only (control): D_1 -0.59 and D_R -0.43, both exclude 0.

## 3. Verdict against the pre-registered criteria
Universal arm (T_1): D_1 covers 0 AND D_R < 0 excluded -> criterion (ii) RELATION-LIKE on the scored set, and cleanly so in the LOW class (T_obs 0.22 vs T_1 0.18 vs T_R 0.55).
Fitted arm (T_1f): both D excluded -> (iii) for that arm; the data sit below both predictions.
The RAR function over-transmits in every class and every arm comparison: the observed curves carry about a third of the baryon-feature amplitude it would print at low g (0.18-0.22 vs 0.55).

## 4. The control that limits the reading [T3]
At HIGH g both laws must transmit nearly fully (T ~ 0.8-0.9: nu -> 1, nb share -> 0), yet the data give T_obs = 0.22 [0.15, 0.61] there too. The data under-transmit even where the physics leaves no room. Two readings, not separable here: (a) the V_bar features are partly artefacts (bulge/disc decomposition, M/L structure, photometric noise -- V_bar has no error column, and the significance cut used V_obs's errors); (b) beam smearing and binning smooth V_obs's features below V_bar's (the resolution limit 2508.03569 names). Either biases T_obs LOW for every model alike. So the firm statement is the ORDERING, not the agreement: the data transmit less than the RAR function by a wide margin in every class; eq. (1) sits closer to the data than the RAR function does at low g; "the data agree with eq. (1)" is NOT claimed, because the HIGH-g control shows T_obs is depressed by something other than dynamics.
This is the same conclusion as 2508.03569 (no clear evidence for Renzo's rule in SPARC at present noise), reached from the opposite direction: they count V_obs features without V_bar counterparts; here V_bar features fail to appear in V_obs at the amplitude either law predicts.

## 5. Named galaxies (T_obs / T_1 / T_1f / T_R; sig = significant baryon feature)
DDO 154 (LOW, sig): 0.77 / 0.32 / 0.39 / 0.53 -- the one named dwarf where the data transmit MORE than both laws (one galaxy; window 4-6 kpc).
IC 2574 (LOW, not sig at this cut; A_bar 0.08 vs noise 0.085): -0.10 / -0.06 / 0.17 / 0.54 -- the paper's "~40%" case; here the 2-4 kpc feature is below noise in V_obs.
NGC 2403 (MID, sig): 0.63 / 1.01 / 0.86 / 0.66. NGC 2841 (MID, sig): 0.67 / 0.52 / 0.72 / 0.67. NGC 3198 (MID, not sig): 0.31 / 0.90 / 0.91 / 0.61.
NGC 6946 (HIGH, sig): 0.62 / 0.81 / 0.74 / 0.76. NGC 7331 (HIGH, sig): 1.05 / 0.63 / 1.09 / 0.98.
UGC 2953, UGC 6787 (HIGH, not sig; 2508.03569's examples): -0.57 / 1.42 / 1.21 / 1.47 and 1.83 / 1.42 / 1.15 / 1.44 -- inner points, noise-dominated.
NGC 1560 (LOW; RESULT_FEATURE_AXISYM): 0.034/0.064 = 0.53 obs (below noise) / 0.20 / 0.27 / 0.56.

## 6. Consequences for v2.0.2 (App A.1)
Replace the 0.33-vs-0.24 correlation sentence with the transmission numbers: "On the 30 galaxies with a significant low- or mid-acceleration feature in V_bar, the observed curves carry a median 0.18 of the baryon feature's amplitude; eq. (1) at universal rho_c carries 0.47 (0.18 in the low-acceleration class), the acceleration relation 0.56. The data transmit less than either; the ordering, not the agreement, is the result, since at high acceleration -- where both must transmit fully -- the data still give 0.22, so part of the deficit is in the baryon features or the resolution." Keep "We do not claim Renzo's rule."
Exposure sentence: "Where a well-resolved curve shows a feature at the amplitude the acceleration relation prints, eq. (1) cannot: it transmits the baryon share only." (NGC 1560 number.)

## 7. Referee flags
F1 V_bar has no error column; the significance cut is one-sided. F2 Window choice (largest excursion, 2-8 points) is one of several; a per-galaxy fixed physical window was not tried. F3 Medians of ratios with small denominators are noisy; the 16-84 bands are bootstrap, not analytic. F4 Fixed arm only; the fitted arm transmits more (T_1f > T_1) because fitted rho_c is often above 7.45 in these galaxies. F5 No beam-smearing model. F6 2508.03569's DTW statistic not reproduced.

## 8. Ballot
FT-R1 sec 3 verdict: (ii) on the universal arm by the pre-registered criterion, REPORTED AS ORDERING ONLY because of the sec 4 control [KEEP with that qualification].
FT-R2 sec 4 reading (data depressed for all models; same conclusion as 2508.03569) [KEEP].
FT-R3 sec 6 A.1 replacement wording + exposure sentence [KEEP; ledger (ch)].
FT-R4 DDO 154 as the single over-transmitting named dwarf [KEEP as record; no claim].
# END -- secs 2-5 [UNRATIFIED]
POINTER 2026-11-18 (BF-79 (a); JM rulings 2026-11-18 20:43 "accepted" and 21:54 "1 go"; nothing above this line is changed). The numbers of secs 1-5 stand as ratified. The READING ratified 2026-11-09 22:03 as FT-R1 / FT-R2 ('ordering only') and the sec 6 wording ratified as FT-R3 ('the ordering, not the agreement, is the result'; the exposure sentence) are AMENDED for the paper. BF-79 (a): the observed curves carry less than either model in every class, the high-acceleration control included (0.22 where both models give about 0.8), so the statistic does not rank the models; the paper gives the low-class numbers (T_obs 0.22, eq. (1) 0.18, acceleration relation 0.55) and draws no ranking. Wording of record: NOTE_V2_FRAMING_20261014 ledger (dw) ADDENDA 7-8. The feature row of the Summary goes under Open (BF-82, 2026-11-18 18:45). Basis: LIT_RENZO_20261118 R4. Record: BALLOT RULINGS 2026-11-18 18:45, 19:19, 20:44, 21:57.
