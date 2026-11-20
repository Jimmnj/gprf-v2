# RESULT_DMS_20261105 -- DiskMass baryon fractions vs eq. (1), no free parameter (PREREG_DMS_20261105) [RATIFIED JM 2026-11-05 21:31 "accepted" (text applied 21:30)]
Run 2026-11-05 19:55 EST: dms_fb_20261105.py -> dms_fb_20261105.csv. 30 DMS galaxies; inputs parsed from the DMS VII/VI arXiv source tables (data/diskmass/). Universal rho_c = 10^7.45; baryons = exponential stellar disc (M_disk*, h_R) + point bulge + exponential atomic gas (2 h_R) + molecular (h_R); DMS gas masses already include He (x1.4; DMS VII l.499/523), no extra factor.

## 0. Two corrections to the pre-registration, made before scoring
(i) F_b in DMS is V_b/V_c at 2.2 h_R (a VELOCITY ratio; DMS VII sec "Disk and baryonic maximality"), not V_b^2/V_c^2 as written in PREREG sec 3. Scored as the velocity ratio. (ii) T2 (BTFR) is VOID: DMS inclinations are i_TF, derived from the K-band Tully-Fisher relation, so the deprojected speeds lie on the TF by construction. Pre-registration error; nothing scored for T2.

## 1. T1 -- F_b at 2.2 h_R
Predicted minus measured, 30 galaxies: median +0.091, rms 0.130, 16-84 [+0.02, +0.17]; excluding UGC 8196 (DMS F_b = 1.09 > 1, flagged by DMS) median +0.092, rms 0.125. Criteria: |median| <= 0.10 PASS (barely), rms <= 0.15 PASS. chi^2 against the DMS errors 85 for 30: the scatter exceeds the quoted errors (2.8 per galaxy), so this is a match at the 0.1 level with a systematic, not a fit.
Direction: the relation at the SPARC universal rho_c gives slightly too LITTLE nb mass in the inner disc of these galaxies; the DMS discs need 1.21x (median) the nb v^2 at 2.2 h_R that 10^7.45 supplies. Expectation (+0.05..+0.15) confirmed; SPARC's fitted rho_c for log M_b ~ 10.5 spirals sits at 7.5-7.7, which would close most of the gap (not fitted here).
Reference points: MOND, F_b = 1 by construction, is off by +0.435 (median), the Angus et al. 2015 problem in this variable. Maximum-disc Newtonian decompositions give F_b,max = 0.92 +- 0.08 (DMS); the sigma_z-based measurement is 0.57 (median).
T3 (report only): pISO Case-I asymptotic v^2 = 4 pi G rho_0 R_C^2 vs the relation's 2 sigma^2: median log ratio +0.12 (the DMS halos, fitted with the sub-maximal Upsilon, carry ~30% more v^2 than the BTFR tail; consistent with the +0.09 in F_b).

## 2. Reading [UNRATIFIED]
(a) A vertical-dynamics measurement that MOND misses by 0.4 in F_b, eq. (1) matches to 0.09 with nothing fitted: M_b and h_R in, F_b out. The nb mass in the inner disc is real mass in the right amount to within ~20% in v^2. This is the first test of the relation against a quantity not derived from rotation curves.
(b) The systematic (+0.09) has the sign and size expected from rho_c: these are massive spirals whose SPARC analogues want rho_c ~0.1-0.3 dex above universal. A per-galaxy rho_c would need DMS rotation curves (Atlas figures only; not tabulated) -- not done.
(c) Caveats: (1) i_TF-based V_c enters DMS's F_b (and all DMS dynamics); (2) Aniyan et al. 2016 argue DMS sigma_z are biased low by young stars, which would RAISE the DMS F_b toward our prediction (and worsen MOND less); unverified this session; (3) baryon model is schematic (pure exponential disc, point bulge, 2 h_R gas) -- DMS's own V_b uses their photometry; the difference is a few per cent in F_b; (4) the DMS error bars are formal; (5) rho_c universal, not fitted: the test is parameter-free but the 0.09 is partly the known rho_c-M_b trend.
(d) Text (HELD), new item in sec 6 or App B/E: "The DiskMass Survey measures the stellar vertical dispersions of 30 face-on spirals and from them the baryonic fraction of the rotation at 2.2 h_R, F_b = V_b/V_c = 0.57 (median), independently of the rotation curve. Solving eq. (1) at the universal rho_c for each galaxy's baryonic mass and scale length predicts F_b = 0.66 (median), 0.09 above the measurement with 0.13 rms; MOND, with no halo, predicts 1 (Angus et al. 2015). The offset has the sign expected if these spirals carry rho_c above the universal value, as their SPARC analogues do [RESULT_DMS_20261105]." Refs to add: Bershady et al. 2010, Martinsson et al. 2013a,b, Angus et al. 2015 (+ Milgrom 2015 reply, Aniyan et al. 2016 if read).
## 3. Ballot: DMS-1 T1 scoring (PASS with systematic); DMS-2 T2 void (prereg error); DMS-3 text (d). 
# END

## CORRECTION POINTER 2026-11-12 01:28 (DM-5; ruled JM 2026-11-11 23:15) -- no rewrite above
SUPERSEDED by CHECK_NP17_DISKMASS_20261111.md secs 1-2: (i) "MOND, F_b = 1 by construction, is off by +0.435" (sec 1) and reading (a) "MOND misses by 0.4" -- the "1" was Claude's PREREG_DMS sec 1 line, misattributed to Angus et al. 2015 (F1); in this variable MOND gives 0.66-0.85, offset +0.08 to +0.25 (F2); like for like MOND matches as well as eq. (1): +0.078 rms 0.106 (RAR function) against +0.091 rms 0.130 (F3). (ii) Text (d) "predicts F_b = 0.66 (median)": the median of the predictions is 0.695 (csv Fb_pred_v); 0.66 is 0.57 + 0.09 (F4). (iii) "independently of the rotation curve": only the disc mass is; F_b divides by the observed V_c (F5). (iv) Caveat (c)(2) on Aniyan understates: factor ~2 in surface density moves the datum to about 0.8 (F7). The T1 score (median +0.091, rms 0.130) STANDS.

## ADDENDUM 2026-11-12 01:28 (DM-6) -- MOND in the DMS variable, read from the csv [UNRATIFIED until JM stamps]
Script mond_fb_check_csv_20261112.py (= mond_fb_check_20261111.py reading dms_fb_20261105.csv directly; logic verbatim); report mond_fb_check_csv_20261112_report.txt; run 01:26:48 EST, device VM. Expectation (written first): identical to the 11/11 transcribed run. RESULT: IDENTICAL, every printed number (diff against the 11/11 footer empty).
| | median F_b | median offset from DMS | rms |
|---|---|---|---|
| DMS measured | 0.565 | - | - |
| eq. (1), universal rho_c | 0.695 | +0.091 | 0.130 |
| MOND like-for-like, RAR / simple / standard | 0.659 / 0.656 / 0.725 | +0.078 / +0.075 / +0.155 | 0.106 / 0.105 / 0.159 |
| MOND own masses, RAR / simple / standard | 0.727 / 0.719 / 0.854 | +0.127 / +0.120 / +0.250 | - |
Per galaxy, eq. (1) minus MOND (RAR) like for like: +0.005 (median). Excluding UGC 8196: +0.078 / +0.076 / +0.157, rms 0.102 / 0.099 / 0.160. Limits: algebraic MOND at 2.2 h_R, not a QUMOND disc; schematic baryon model inherited from dms_fb_20261105.py.
Reading [T3]: the DMS check does not separate eq. (1) from MOND. By the PREREG_DMS criteria (|median| <= 0.10, rms <= 0.15) MOND passes with the RAR and simple functions and fails the median with the standard function (+0.155). A check both pass from the same masses (DM-2 wording, ruled).
