# RESULT_N12_EMERGENT_NU_20261020 [UNRATIFIED] -- the interpolating function emerges from the frozen law with ONE universal constant. Pre-registered (script header).
Script n12_emergent_nu_20261020.py; per-point n12_emergent_nu_20261020.csv; bins n12_emergent_nu_bins_20261020.csv; raw n12_report_20261020.txt. Device VM, ~2.5 min.
Setup: TEST1 solver (tanh2 cap, BTFR sigma^2 = sqrt(G a0 M_b)/2, rho_amb 35.5), SPARC N=120 (Q<=2, inc>30, Yd 0.5, Yb 0.7), observed baryon profiles,
rho_c FIXED at 10^7.45 Msun/kpc^3 for every galaxy (no per-galaxy parameter); sensitivity 10^7.2 and 10^7.8; reference: fitted rho_c. nu_eff = g_model/g_bar at
the data radii, binned in log(g_bar/a0). Comparator: RAR function nu = 1/(1-exp(-sqrt(g_bar/a0))) (McGaugh, Lelli & Schombert 2016).
Declared PASS: in |log g_bar/a0| < 0.5, cross-galaxy scatter of log nu_eff <= 0.15 dex AND |median offset from RAR| <= 0.10 dex at rho_c 7.45.

## Result: PASS, by a wide margin.
| log g_bar/a0 bin | N | universal nu_eff median | scatter | offset from RAR | data offset from RAR | data scatter |
| -2.25..-2.00 | 13  | +1.076 | 0.183 | +0.034 | +0.056 | 0.322 |
| -2.00..-1.75 | 60  | +0.926 | 0.152 | -0.023 | +0.003 | 0.210 |
| -1.75..-1.50 | 194 | +0.853 | 0.146 | +0.008 | -0.058 | 0.180 |
| -1.50..-1.25 | 361 | +0.715 | 0.122 | -0.016 | -0.028 | 0.148 |
| -1.25..-1.00 | 390 | +0.584 | 0.098 | -0.036 | -0.003 | 0.164 |
| -1.00..-0.75 | 309 | +0.485 | 0.083 | -0.031 | -0.015 | 0.133 |
| -0.75..-0.50 | 286 | +0.396 | 0.113 | -0.023 | +0.010 | 0.149 |
| -0.50..-0.25 | 284 | +0.305 | 0.093 | -0.016 | -0.005 | 0.136 |
| -0.25..+0.00 | 253 | +0.225 | 0.078 | -0.016 | +0.014 | 0.172 |
| +0.00..+0.25 | 224 | +0.148 | 0.050 | -0.016 | +0.007 | 0.186 |
| +0.25..+0.50 | 163 | +0.088 | 0.029 | -0.012 | +0.017 | 0.128 |
| +0.50..+0.75 | 123 | +0.037 | 0.017 | -0.018 | -0.013 | 0.178 |
| +0.75..+1.00 | 82  | +0.015 | 0.007 | -0.015 | -0.002 | 0.139 |
Transition (|log g_bar/a0| < 0.5): scatter max 0.093 dex; |offset from RAR| max 0.016 dex. Over the full range -2.25..+1.0 (3.25 dex, N>=13 per bin):
|offset| <= 0.036 dex. Sensitivity: rho_c 7.2 -> 7.8 shifts nu_eff by -0.08..+0.09 dex, nearly uniformly across bins (a normalisation, not a shape change).
Fitted-rho_c model scatter 0.08-0.15; data scatter 0.13-0.19 (includes measurement error; McGaugh+2016 quote ~0.11 dex observed, ~0.05-0.08 intrinsic [mem, VERIFY]).

## Reading
1. The law reproduces McGaugh's interpolating function over three decades of g_bar with ONE global constant (rho_c = 10^7.45) plus a0 (through the BTFR).
   The two MOND limits are exact by construction (isothermal tail: g -> sqrt(g_bar a0); cap: nu -> 1); the TRANSITION SHAPE is an output, not an input.
2. The cross-galaxy scatter of the emergent nu at fixed g_bar (0.05-0.09 dex in the transition) is below the observed RAR scatter and consistent with the
   intrinsic scatter estimates: the law does not predict more RAR scatter than the data allow, even with rho_c universal.
3. rho_c acts on nu as a nearly uniform vertical shift; the shape of the transition comes from the cored-isothermal solve. The 0.65 dex per-galaxy spread
   of fitted rho_c therefore translates into ~0.1 dex of RAR scatter, not 0.65 -- which is why N15/N9b found the RAR reproduced.
## Caveats (referee)
- rho_c = 10^7.45 is the MEDIAN of the per-galaxy SPARC fits, so the normalisation is not independent of these data; the SHAPE is. An out-of-sample
  normalisation would use Donato/Milgrom's mu0 (which the fitted rho_c already matches, MU0-SPLIT).
- Comparator nu_RAR is itself a fit to SPARC; the honest comparison is to the binned DATA (last two columns): model-data offsets <= 0.034 dex in all bins
  with N >= 60 EXCEPT -1.75..-1.50 (0.066) [CORRECTED 22:40 after second cold review; earlier '<= 0.05' was wrong].
- Low-g bins (< -2) have few points and larger scatter in both model and data.
- Same M/L (Yd 0.5) and inclination cuts as TEST1; a0 = 1.2e-10 fixed.
## Claim for the note (framing (az)): 'With rho_c fixed to a single value for all galaxies, the law reproduces the radial acceleration relation --
including the shape of the transition between the Newtonian and deep regimes -- to within 0.04 dex over three decades in g_bar, with scatter below
that of the data. No interpolating function is assumed.' This is the headline candidate; it answers (aw)-style 'MOND with an extra parameter' directly:
the extra parameter is not needed for the RAR; it is needed for the residual per-galaxy curve shape (TEST1 chi 13.8 vs free-sigma 5.6) and for dwarfs.
- INTRINSIC SCATTER (added 12:58; JM 12:55 'so you're saying our btfr scatter is half of mond' -- no): the model's 0.05-0.09 dex is INTRINSIC
  (baryon-profile variety at universal rho_c); the data's 0.13-0.19 includes measurement error (McGaugh+2016: ~0.11 observed, intrinsic ~0.05-0.08
  [mem, VERIFY]). MOND predicts ZERO intrinsic RAR scatter. If intrinsic scatter is pinned below ~0.05 dex, universal rho_c is disfavoured (fitted
  rho_c survives, but the 'no parameter for the RAR' claim weakens). Not a BTFR statement: BTFR is an input in both frameworks.

## ADDENDUM 13:20 -- per-point residuals (JM 13:16 'can a fit per galaxy get closer'). rms log(g_obs/g_model), all / transition / x<-1:
universal 7.45: 0.172 / 0.154 / 0.187 | fitted rho_c: 0.107 / 0.099 / 0.104 | MOND RAR function: 0.159 / 0.156 / 0.165.
Per-galaxy rms: universal median 0.123, fitted 0.073; fitted better in 109/120, by >0.05 dex in 44. Fitted residual ~ measurement floor (10% v -> 0.09 dex).
READING: universal rho_c = MOND-level; per-galaxy rho_c halves the residual, which a one-function law cannot do (its scatter is irreducible). Success iff
rho_c is physical (aw). Note sentence drafted in framing (az').

## ADDENDUM 13:40 -- one-knob-each comparison (n12b_mond_freeML_20261020.py / n12b_report / csv; JM 'can you run it?'). Same 120 galaxies, 2837 points,
10% velocity weighting. Per-point rms: MOND Y fixed 0.159 | MOND Y_disk free 0.125 | MOND Y prior 0.1 dex 0.125 | GPRF rho_c free (Y 0.5) 0.108 | universal 0.173.
Per-galaxy median rms 0.125 / 0.081 / 0.082 / 0.073 / 0.123. GPRF better in 73/120 (vs free Y). MOND free Y median 0.42 (16-84% 0.16-0.69); 19/120 at bounds.
READING: (az') survives WEAKENED: ~14% lower rms with one knob each; 'halves' withdrawn. Li+2018 also free D and i -> matched D/i comparison OWED
before quoting against their 0.05. MOND's 19 bound hits = unphysical M/L demanded (say carefully). MOND's knob has external priors; ours none.

## ADDENDUM 13:50 -- gas split (n12b_gas_split_report_20261020.txt; JM 'what about gas'). Per-galaxy rms median, MOND Y fixed / Y free / GPRF universal / fitted:
f_gas<0.2: 0.115 / 0.061 / 0.112 / 0.071 (n=42, GPRF-fit better 22/42) | 0.2-0.5: 0.092 / 0.076 / 0.113 / 0.069 (35; 22/35) | >=0.5: 0.145 / 0.104 / 0.181 / 0.081 (43; 29/43).
READING: (1) the 73/120 one-knob advantage comes from gas-rich galaxies where MOND's M/L knob has no leverage -> NOT one-knob-each; (az''') DROPPED.
(2) Star-dominated, both knobs live: MOND wins 0.061 vs 0.071. (3) No knob either side: MOND beats universal rho_c in gas-rich galaxies (0.145 vs 0.181);
10^7.45 is too high there, consistent with RHOC-FGAS (gas-rich want lower rho_c). CAVEAT ON (az): universal-rho_c nu holds in the binned RAR but
underperforms the RAR function per galaxy in gas-rich systems by ~0.04 dex.

## ADDENDUM 16:20 -- S-form / cap sensitivity (n12c_sform_sensitivity_20261020.py; report; bins csv; cold review A1). Universal rho_c 7.45, 120 galaxies each.
max |median log nu_eff - tanh2| transition / full(-2..+1): rat2 0.020/0.051 (within 0.05/0.10 bar) | tanh1.5 0.038/0.043 (within) | tanh1 0.065/0.089 | tanh3 0.070/0.085 |
tanh4 0.134/0.153 | lin 0.220/0.252. FORM at fixed cap: immaterial (0.02 dex). CAP VALUE: ~0.07 dex per unit around 2, peaking at log g_bar/a0 -1..-0.5.
RAR offset of tanh2 is -0.016..-0.036; tanh1.5 adds +0.04 -> data prefer cap ~1.7-2. The cap at 2 is the identity x_b = 2 sqrt(g_b/a0) = 2 at g_b = a0 (FRM-4
GRADIENT 10/18), i.e. 'response stops growing once the baryonic acceleration exceeds a0' -- a0 entering through the same sigma, not a tuned number.
REVISED CLAIM (replaces 'no interpolating function assumed'): the transition shape follows from a potential-space law that saturates where g_b = a0; the
form of the saturation is immaterial; no function is fitted. Figure 1 gains a panel: nu_eff for cap 1/1.5/2/3/4 vs RAR.

## ADDENDUM 16:35 -- a-priori rho_c from Donato mu0 (n12d_apriori_rhoc_20261020.py; report; bins csv; cold review A2/#3). FAIL at the 0.05 bar.
Per-galaxy rho_c = (mu0/3 sigma)^2 4piG (median 7.23, range 6.56-8.29): shape WRONG (+0.16 at log g/a0 -1.6, -0.16 at -0.6; rms 0.198) -- it is the dead
a0-closure (mass slope -0.5). Global 7.23: uniformly -0.06..-0.09 below binned data (max 0.092); consistent with the sensitivity (0.22 dex in rho_c -> 0.07 in nu).
READING: shape derived; NORMALISATION is one global number from the data (SPARC median); an independent value (Donato) lands within ~0.1 dex. Half-sample
hold-out is not the issue (median SE ~0.08 dex); the single number is. Note wording: 'one global density, 10^7.45, the sample median; the Donato surface
density gives 10^7.2, within 0.1 dex in the relation.'

## ADDENDUM 19:50 -- MOND with Y_disk + D + i free under SPARC priors (n12e_mond_Di_free_20261020.py / report / csv; cold review A3, Li+2018 parity): per-point rms 0.112,
per-galaxy median 0.075, chi median 12.5 -- vs GPRF rho_c-only 0.108 / 0.073 / 13.8. A WASH. Inclination did not move (D carried the gain: D/D0 0.97, 0.80-1.07).
Li+2018's ~0.05 not reproduced (different error model/MCMC). CONCLUSION: 'comparable' is the only defensible word; one-knob claim withdrawn in every form; the matched
3-parameter GPRF run is moot. Note sec 3.0 table: add this row.
