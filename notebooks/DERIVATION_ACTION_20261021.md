# DERIVATION_ACTION_20261021 [UNRATIFIED] -- can the frozen law come from an action?
JM 2026-10-21 19:22 "can you derive a lagrangian"; 20:32 "ok" to the QUMOND-template attempt. Conventional-sign potentials (Phi <= 0); bump x = -Phi/sigma^2.

## 1. Target
psi: baryon Newtonian potential, lap psi = 4piG rho_b. chi: nb potential, lap chi = 4piG rho_nb. Matter feels Phi = psi + chi only.
Frozen law (ambient term dropped for clarity; it is linear and integrable): rho_nb = A exp[S(-psi/sigma^2) - chi/sigma^2], S(x) = 2 tanh(x/2).
Two nonlocalities: (i) psi is a functional of rho_b (the body's own baryon potential); (ii) sigma^2 = sqrt(G a0 M_b)/2 is a functional of the body's total M_b.

## 2. QUMOND template (Milgrom 2010) -- handles (i) without a reaction force
L = -(1/4piG) grad Phi . grad psi + (1/8piG)|grad psi|^2 - rho_b Phi + R.  delta Phi -> lap psi = 4piG rho_b (psi is the baryons' Newtonian field);
delta psi -> lap Phi = lap psi - 4piG dR/dpsi. Matter couples to Phi only, so no reaction force (contrast the naive action, RESULT_LAGRANGIAN_CHECKS (B)).
With R = P(psi): rho_nb = -dP/dpsi, any function of psi. WORKS -- but has no self-term, hence no isothermal 1/r^2 tail.

## 3. Adding the self-term: integrability
Add chi with canonical kinetic term and local couplings (and at most a linear multiplier field). The field equations close only if
   d rho_nb / d chi = 2 gamma d rho_nb / d psi  (gamma constant),  i.e.  rho_nb = F(psi + kappa chi)  for some function F and constant kappa.
For the frozen law (rho_chi / rho_psi) = cosh^2(x_b/2) = 1/S'(x_b): constant only if S is linear (verified symbolically, /tmp/sym.py output in ballot).
So: within this class THE FROZEN LAW HAS NO ACTION. The obstruction is exactly the saturation acting on the baryon potential alone.

## 4. The integrable alternative fails
kappa is fixed to 1 by the tail (v^2 = 2 sigma^2 needs the nb to respond to its own potential with the same sigma), so rho_nb = F(Phi_total):
saturation on the TOTAL potential -- the single-field perfect-fluid action L = -|grad Phi|^2/8piG - rho_b Phi + P(Phi).
Analytic: an isothermal sphere's own potential below its edge is ~2 sigma^2 ln(r_t/r), i.e. x_total > 2 for r < r_t/e (the self-term identity
delta phi_nb/sigma^2 ~ 12). With S saturating on x_total the response is flat (rho ~ rho_c) throughout r < r_t/e: no tail, mass grows as r^3.
Numerical (lawA_total_potential_20261021.py): the iteration for the nb central potential runs away (x_nb0 -> 1e6); MW rms 212 (universal) / 40 (best
rho_c 10^6.74) against 19.7 / 8.6 for the frozen law. DEAD.

## 5. Conclusion
Saturation on the baryon potential + an isothermal self-tail cannot coexist in an action with canonical fields and local couplings. Not proved for
non-canonical (AQUAL-type, gradient-dependent) couplings: the simplest gradient coupling (saturation through |grad psi|, i.e. through g_b, which is where the
cap physically sits: x_b = 2 sqrt(g_b/a0)) makes the gravitating nb density a total divergence (zero net mass) -- fails, but the class is not exhausted.
Nonlocality (ii), sigma(M_b), is untouched by any of this; MOND avoids it because sqrt(g_N a0) is local.
## 6. Routes
(a) Keep the law non-variational; state momentum conservation as open with the numbers (current draft sec 7).
(b) Give the nb mass its own equation of motion (a fluid with momentum, as in superfluid DM, Berezhiani & Khoury 2015): momentum conserved by construction,
    the response time (open item) becomes defined, saturation must come from a baryon-nb interaction term; loses 'no dynamics'. A different paper.
(c) Non-canonical couplings through g_b (AQUAL-like): unexplored; likely to bring an EFE, which would change sec 6.2.
