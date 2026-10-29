# RESULT_SFORM_FAMILY_20261028 [UNRATIFIED]

Run 2026-10-28 14:50 under PREREG_SFORM_20261028.md (STAMPED JM 14:44 "accepted - go"; FROZEN).
Script sform_family_20261028.py; reports sform_family_{univ,rel3}_report_20261028.txt; bins sform_family_bins_{univ,rel3}_20261028.csv.
Harness check: tanh/univ reproduces pass3 exactly (rms 0.155; all 13 binned medians identical).

## Verdict
- D1 'form immaterial': FAILS. hard cap M1 0.038 (rel3 0.035), rational M1 0.044 (rel3 0.045) > 0.03. erf 0.007, exp 0.021 pass. M3 <= 0.015 for all (passes).
- D2 'shape relieves the scatter tension': NOT MET. 0 of 6 window bins for any form (univ); 1 of 6 for every form including tanh (rel3).
- D3: tanh stays. No form adopted.

## Numbers (primary, universal rho_c 10^7.45; 147 gal, 2696 pts)
| form | rms | M1 max abs shift of bin median | M3 max abs shift of own scatter (window) | M4 bins |
|---|---|---|---|---|
| tanh | 0.1551 | 0 | 0 | 0/6 |
| hard | 0.1600 | 0.038 | 0.013 | 0/6 |
| erf | 0.1556 | 0.007 | 0.002 | 0/6 |
| exp | 0.1550 | 0.021 | 0.008 | 0/6 |
| rat | 0.1553 | 0.044 | 0.015 | 0/6 |
Secondary (relation (3) rho_c): rms 0.140 / 0.144 / 0.140 / 0.140 / 0.142; M1 0 / 0.035 / 0.006 / 0.024 / 0.045.

## Reading points
(a) The note's '0.02 dex' rests on the rational form on the old 120-galaxy sample. On the full 147-galaxy sample the same rational form moves the binned relation by 0.044 dex. The sentence is wrong as it stands.
(b) The shape TILTS the transition; it does not change the overall fit. rms moves by <= 0.005. A sharper cap raises the model in the bins below log g_bar/a0 = -0.75 by up to 0.03 and lowers it above -0.5 by up to 0.04; a softer cap does the reverse.
(c) The low-acceleration offsets the note quotes (0.04-0.08 dex, bins -2.0..-1.25) are form-dependent at the 0.01-0.03 level: hard cap -0.036 vs tanh -0.064 in -1.50..-1.25.
(d) The model's own scatter moves <= 0.015. Softer forms lower it in the two lowest window bins (rat 0.098 vs tanh 0.114), sharper raise it. No form brings model+errors below the data. The tension is not a shape artefact.
(e) For the note: replace 'immaterial (0.02 dex)' with the measured range, e.g. 'At fixed cap, the shape of S moves the binned relation by up to 0.04 dex (0.01-0.02 for erf and exponential forms, 0.04 for a hard cap or a rational form), tilting the transition without changing the per-point rms (0.155-0.160) or the scatter tension; tanh is fixed before any fit.' Places: Table 1 'form of S' row; sec 3 'its functional form is immaterial (0.02 dex, below)' and 'changing the FORM ... 0.02 dex'; Fig 1 caption 'form immaterial'. Counting argument (no fitted function) survives because tanh is fixed a priori and not chosen on this result; a referee can still call S a 0.04-dex systematic. JM to rule.
