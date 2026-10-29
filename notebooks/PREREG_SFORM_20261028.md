# PREREG_SFORM_20261028 [STAMPED JM 2026-10-28 14:44 "accepted - go"; FROZEN]

Drafted 2026-10-28 14:20 by Claude, answering JM "try 1" (14:12) = the shape-family test of S before any action work.
Nothing is run until JM stamps. After the stamp nothing below changes.

## 1. Question
Is the claim "the form of S is immaterial (0.02 dex)" (sec 3, Table 1) true across a family of shapes, and does the shape move the
model's own scatter in the tension window (log g_bar/a0 = -1.75 .. -0.25)? Previous evidence: one alternative (x/(1+x/2)) on the
older 120-galaxy sample (n12c_sform_sensitivity_20261020).

## 2. Forms (all: S(0)=0, S'(0)=1, cap 2 -- only the transition shape differs)
- tanh  : 2 tanh(x/2)                 [baseline, the note's form]
- hard  : min(x, 2)                   [sharpest]
- erf   : 2 erf(sqrt(pi) x / 4)
- exp   : 2 (1 - exp(-x/2))
- rat   : x / (1 + x/2)                [softest; the earlier alternative]
No other change: same solver, same sigma (eq. 2), same rho_amb.

## 3. Harness and sample (frozen)
pass2_sparc_full_20261022.py, unchanged except S. S2 = 147 galaxies, 2696 Lelli-style points (Q<=2, inc>=30, ev/v<0.1).
Primary: universal rho_c = 10^7.45 (the case the note's sensitivity sentence refers to). Secondary: rho_c from relation (3).
Error injection exactly as (15b): D, i, 0.1 dex M/L, velocity errors; 50 realisations, seed 7.

## 4. Metrics (each form against tanh; bins of 0.25 dex with N >= 10)
M1 max over bins |delta median log(g_obs/g_model)|
M2 per-point rms log(g_obs/g_model)
M3 model's own scatter ('law alone', sd of log g_model - log g_RARfn) per bin in the window: max |delta|
M4 model+errors scatter minus data scatter, per bin in the window

## 5. Decision rules (fixed now)
D1 'Form immaterial' STANDS if every form has M1 <= 0.03 AND M3 <= 0.02 (primary case). Otherwise the note's sentence is reworded to give the measured range; it does not say 'immaterial'.
D2 'Shape can relieve the scatter tension' if some form has model+errors <= data - 0.01 in at least 4 of the 6 window bins. Reported only.
D3 NO FORM IS ADOPTED ON THE RESULT. tanh stays the note's form whatever happens. Choosing a shape to fit would make S a fitted function
   (the counting argument, sec 3, would be lost); a referee-visible statement of the measured sensitivity is the only output.

## 6. Pre-stated expectation
M1 small for all smooth forms (<= 0.03); hard cap the largest. M3 shift small (<= 0.02); sign not predicted. D2 not met.
If D1 fails for 'hard' only, the note says so explicitly (smooth forms immaterial; a hard cap is not).

## 7. Outputs
sform_family_20261028.py, sform_family_20261028_report.txt, sform_family_bins_20261028.csv, RESULT_SFORM_FAMILY_20261028.md [UNRATIFIED].
Nothing enters the tex until JM rules on the result.
