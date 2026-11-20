# RESULT_FEATURE_AXISYM_20261109 -- NGC 1560 gas feature in disc geometry (PREREG_FEATURE_AXISYM_20261109, stamped 13:20; FA-4 amendment stamped 13:24) [RATIFIED 2026-11-09 22:03 -- JM "ratified"; BALLOT 22:03 line]
Written 2026-11-09 13:40 EST. Script: feature_axisym_20261109.py; outputs feature_axisym_sigma_20261109.csv (ring Sigma_gas, Sigma_disk), feature_axisym_curves_20261109.csv. ARM: fixed M/L_I 1.43, D 3.45, as PREREG_NGC1560.

## 1. Baryon model (FA-4): 120 log rings 0.1-25 kpc, NNLS, soft total-gas-mass constraint. Fidelity: gas rms 0.10 km/s, disc rms 0.35 km/s (gates 1.0), gas peak 3.39 -> 3.51 kpc (gate 0.3), ring gas mass 1.10e9 = 1.33 M_HI. PASS. The first (pre-amendment) model failed the gate (sec 10 of the prereg); nothing was scored from it.
## 2. Control: R1 re-solved = stored R1 (0.000 km/s). A0 (sphere on the ring-model in-plane V_bar) vs R1: max |dV| 0.39 km/s (gate 0.5). PASS -- the new harness is tied to the ratified result.
## 3. Runs (full, l <= 12; A1 converged in 9 iterations; nb isodensity q(2 Rd) = 0.74 against a thin disc)
curve              rms(dex)  chi2     A(3-6.5 kpc)  r_peak
V_obs              --        --       0.0344        4.64   (NULL-F: below 2x noise 0.036, as ratified)
R3 RAR function    0.0727    878.8    0.0356        3.01
R1 sphere 7.45     0.0673    1292.9   0.0124        3.01
A0 sphere/rings    0.0670    1289.1   0.0128        3.01
A1 disc 7.45       0.0612    1032.8   0.0127        3.01
A2 disc fit 7.15   0.0544    198.6    0.0174        3.01
A3 disc M* 1.0e9   0.0524    688.9    0.0103        3.01
A3 disc M* 1.6e9   0.0692    1408.6   0.0150        3.01
Pearson of detrended curves vs V_obs: A1 0.78, R1 0.78, R3 0.79.
## 4. Verdict against the pre-registered criteria
A(A1)/A(R3) = 0.36 ; A(A1)/A(R1 sphere) = 1.02 -> (iii) GEOMETRY IRRELEVANT. Disc geometry with l <= 12 changes the transmitted feature amplitude by 2%. RESULT_NGC1560 sec 5 (addendum 20:40) is CONFIRMED: the ~1/3 damping is the law's -- rho_nb responds to phi_b (one integral), M_nb is a second integral, S saturates, the core exceeds the feature width; the feature survives as the baryon share. Structural for any Poisson-obeying real mass. The spherical harness was not the problem here either (companion to RESULT_AXISYM_R1 sec 6).
Side results: fitted log rho_c 7.15 in disc geometry vs 7.10 spherical (+0.05, same sign and size as the SPARC shift +0.03); rms at universal 7.45 improves 0.067 -> 0.061 (ratio to R3 0.93 -> 0.84; PASS-U stands); chi2 1293 -> 1033 (still loses to R3 on chi2 -- say both, as ratified).
## 5. Consequences
- v2.0.2 NGC 1560 sentence: the damping clause is WRITTEN with this number (N1560-R3 wording HOLD released): "the relation transmits the gas feature at about one third of the RAR function's amplitude (0.013 vs 0.036 dex), in spherical and in disc geometry alike".
- Sec 7 feature-correlation sentence (0.33 vs 0.24): the NGC 1560 flag stands (0.78 vs 0.79); SPARC detrended check still owed before repeating it.
- Owed email for the ASCII mass model: no longer blocking for this test (reconstruction passed its gates); still a provenance item.
## 6. Referee flags
F1 Baryon surface densities are inverted from velocities (NNLS); the inversion is not unique, but the fidelity gates fix what matters for the potential (in-plane V to 0.1-0.35 km/s, total mass). F2 The M_HI soft constraint places unconstrained gas beyond the last data point; it affects phi_b(0) via the tail (small; the A3 bracket brackets a larger effect). F3 Thin-disc scale heights assumed (0.2 kpc gas; 0.19 kpc stars). F4 R3 remains the algebraic RAR function, not a QUMOND disc solve (few %, stated in the prereg).
## 7. Ballot
FX-1 sec 4 verdict (iii), sec 5 reading confirmed [KEEP].  FX-2 FA-4 amendment as the model of record [KEEP].  FX-3 v2.0.2 damping-clause wording (sec 5) [KEEP; ledger (cg)].  FX-4 side results (rho_c 7.15; rms 0.061) [KEEP as record].
# END -- secs 3-5 [UNRATIFIED]

## POINTER 2026-11-12 01:08 (DM-20) -- gas/stellar columns SWAPPED in this run [RULED OW-2 KEEP, JM 2026-11-12 01:25; no rewrite above]
feature_axisym_20261109.py reads the mislabelled csv as labelled (ERRATUM_NGC1560_COLUMNS_20261111): the "gas" rings are the stellar curve, the "disc" rings are the gas, M_b 2.40e9. Verdict (iii) (FX-1) survives as a like-for-like ratio (1.02). A3 rows INVALID as labelled (they scale the gas). FX-4 side results (10^7.15; rms 0.061; chi2 1033) computed at the swapped masses: superseded. "0.013" in sec 5 is A0/A1 at those masses; corrected spherical value 0.014. Detail: RESULT_NGC1560_OWED_20261112 sec 2.
