# PREREG_FEATURE_AXISYM_20261109 -- does disc geometry restore the NGC 1560 gas feature in eq. (1)?
Written 2026-11-09 13:20 EST, BEFORE any disc-geometry solve of NGC 1560. [UNRATIFIED] until JM stamps. Trigger: RESULT_NGC1560_20261108 sec 5 (feature transmitted at ~1/3 of the RAR function's amplitude in the SPHERICAL harness; NULL-F on Gentile's curve) and SESSION_HANDOFF_20261109 queue item 4(ii). Ballot AX-7 (13:09): write, stamp before run.
Question: is the 1/3 damping a property of the law (two integrations + saturation) or of the spherical harness (a thin gas ring smeared into a shell)? RESULT_AXISYM_R1 sec 3(d) says the response is rounder than its source on SPARC; this test asks the same of a local feature.

## 1. Data (fixed; unchanged from PREREG_NGC1560 sec 1)
Gentile+2010 curve as extracted (data/ngc1560/ngc1560_gentile2010_rc_extracted_20261108.csv; 63 points); D 3.45 Mpc; M/L_I 1.43; M_HI 8.28e8 (x1.33); M_* primary 1.3e9 (bracket 1.0-1.6e9 as R4); M_b 2.40e9; sigma from eq. (2). In-plane model velocity = Gentile's V_bar^2 + nb part (same assumption as RESULT_AXISYM_R1 F2).

## 2. Baryon model in disc geometry (new)
Surface densities are NOT tabulated; they are reconstructed from Gentile's V_gas(r) and V_disk(r) by non-negative least squares on 40 logarithmic rings (0.1-15 kpc) using the same thin-ring potentials as the solver (ring_phi), gas scale height 0.2 kpc, stellar h_z = 0.196 Rd^0.633 with Rd = 1.3 kpc (Gentile/Broeils). FIDELITY CRITERION (fixed): the reconstructed in-plane V_gas and V_disk must reproduce the input to rms < 1.0 km/s over 1-9.5 kpc and keep the gas peak at 3.4 kpc within 0.3 kpc; otherwise the run is reported as "baryon model failed" and nothing is scored. The reconstruction is reported (Sigma_gas, Sigma_disk vs r) and is the test's weakest link (NNLS inversion of a rotation curve is ill-conditioned; regularisation = ring count and non-negativity only).
Control: the same ring model solved as a sphere (C0 on the reconstructed V_bar) must reproduce RESULT_NGC1560 R1 to < 0.5 km/s at every point; this ties the new harness to the ratified result.

## 3. Runs
A1 eq. (1), rho_c = 10^7.45, mode full (l <= 12), M_* primary.   A2 same, rho_c fitted (R2 analogue; bounded Brent [6,10]).
A3 A1 at M_* = 1.0e9 and 1.6e9.   A0 control: C0 on reconstructed V_bar (sec 2).   R3 RAR function (from RESULT_NGC1560, unchanged).

## 4. Criteria (fixed now; feature recipe identical to PREREG_NGC1560 (c))
Detrend log V by a per-curve quadratic in log r over 2 <= r <= 8 kpc; A = peak-to-trough of the detrended log V in 3 <= r <= 6.5 kpc. Reference values (ratified 01:06): A(V_obs) 0.034 (NULL-F: below 2x noise 0.036), A(R3) 0.036, A(R1 spherical) 0.012-0.018 (~1/3 of R3).
(i) GEOMETRY RESTORES: A(A1) >= 0.8 A(R3) -> the spherical harness was the problem; RESULT_NGC1560 sec 5's structural reading is RETRACTED and the v2.0.2 damping clause is dropped.
(ii) GEOMETRY HELPS, LAW STILL DAMPS: 1.5 A(R1 sph) <= A(A1) < 0.8 A(R3) -> report both; the sec 5 reading stands with the number corrected.
(iii) GEOMETRY IRRELEVANT: A(A1) < 1.5 A(R1 sph) -> sec 5 reading CONFIRMED; the 1/3 damping is the law's (two integrations + saturation), as the addendum argued.
Also report: fitted log rho_c in A2 vs 7.10 (spherical R2); rms and chi^2 of A1 vs R1; A(A3) bracket; Pearson of detrended A1 vs V_obs against R3 vs V_obs (reference 0.78 vs 0.79); the nb isodensity axis ratio at 2 Rd.

## 5. Expectation (Claude): (iii) or low (ii). A thin gas ring makes a sharper feature in phi_b than a shell, but eq. (1) still integrates twice and saturates; the SPARC sweep found l >= 2 terms worth < 0.01 dex. Expect A(A1) ~ 0.015-0.025 against R3's 0.036. If (i) occurs the structural claim goes, which is the point of pre-registering it.

## 6. Failure branches
(i) -> retraction line in RESULT_NGC1560 sec 5 + ledger; v2.0.2 NGC 1560 sentence loses the damping clause. (ii)/(iii) -> the v2.0.2 damping clause is written with the disc-geometry number (handoff item N1560-R3 wording HOLD released). Baryon model failed -> owed email for the ASCII mass model becomes blocking for this test.

## 7. Not done: QUMOND disc solve of NGC 1560 (R3 stays the algebraic function); nuisance marginalisation; the Broeils/Sanders curve.
## 8. Ballot: FA-1 stamp (keep / amend / kill). FA-2 NNLS ring reconstruction acceptable as the baryon model, with the sec 2 fidelity criterion? FA-3 thresholds 0.8 / 1.5 acceptable?
## 9. STAMP
> **JM 2026-11-09 13:20:** "stamped go"  -> FA-1 keep as written; FA-2 NNLS reconstruction with fidelity gate accepted; FA-3 thresholds 0.8/1.5 accepted.
# END

## 10. AMENDMENT FA-4 (2026-11-09 13:30, BEFORE any physics solve)
First run with the sec 2 model (40 rings, 0.1-15 kpc): FIDELITY GATE FAILED -- disc rms 2.68 km/s (gate 1.0), gas peak 3.39 -> 3.76 kpc (gate 0.3), ring gas mass 5.9e8 vs 1.1e9. Nothing scored. Diagnosis (fidelity only, no solve): the extracted V_disk has 1.3 km/s rms second differences; 40 log rings cannot follow it. Grid table: 80 rings/15 kpc gas 0.23 / disc 0.82 / peak 3.76; 120 rings/25 kpc + soft total-gas-mass constraint (1.33 M_HI): gas 0.10 / disc 0.35 / peak 3.51 / M_gas 1.10e9 -- passes every gate.
AMENDED MODEL: 120 log rings 0.1-25 kpc; soft constraint sum(Sigma_gas * ring area) = 1.33 M_HI (the curve does not constrain mass beyond the last point); everything else (thickness, fidelity gates, control, criteria, thresholds) UNCHANGED.
> **JM 2026-11-09 13:24:** "send an email and go"  -> FA-4 STAMPED; owed email (NGC 1560 ASCII mass model) to be sent now.
