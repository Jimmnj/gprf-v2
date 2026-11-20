# ngc0289_diag_20261112.py -- why does NGC 0289 prefer log rho_c 10.37? POST-HOC diagnosis, not a test. [UNRATIFIED]
# Trigger: JM 2026-11-12 22:32 "is there anything different about 0289 that could explain it?"
# Solver/chi^2: the errprop_globalfit chi_fn (= pass2 solve_gal + pass2_table2 law_chi, verbatim; rule 3a followed 21:50).
# EXPECTATION (written 22:34 before the run): the two branches fit different radii -- the high branch (10.37) is bought by the INNER points
# (a steep inner rise above the baryons), the low branch (6.92) by the OUTER points (V_flat 163 sits below the BTFR value for its large M_b,
# which by the (V_f/V_BTFR)^6 lever wants a LOW rho_c). The galaxy is unusual in having a huge, extended HI disc and a bright compact centre.
import os, json, numpy as np, pandas as pd
os.environ['NOFIT']='1'; exec(open('errprop_globalfit_20261112.py').read()); kpc=ns['kpc']; Gc=ns['G']; a0k=ns['a0k']
g='NGC0289'; d=P0[P0.gal==g].sort_values('r'); t=T0.loc[g]; m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values
Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; Vbtfr=(Gc*a0k*Mb)**0.25
print('T %d | D %.1f +- %.1f (f_D Hubble flow) | i %.0f +- %.0f | L %.1fe9 | M_HI %.1fe9 | M_b %.2e | V_f %.0f vs V_BTFR %.0f (ratio %.2f) | Rd %.2f | Q %d'%(4,t.D,t.eD,t.inc,t.einc,t.L,t.MHI,Mb,t.Vf,Vbtfr,t.Vf/Vbtfr,t.Rd,t.Q))
f=chi_fn(g)
def curve(l):
    ns['P']=d; ns['T']=pd.DataFrame([t],index=[g]); o=solve_gal(g,l); return np.sqrt(o.vmod2.values)
vo=d.vobs.values; r=d.r.values; w=1/(0.1*vo)**2
print(' r(kpc)  Vobs   ev   Vgas  Vdisk  Vbul  Vbar | model 6.92  7.45  10.37 | chi2 contrib 6.92 / 10.37  (Lelli mask)')
C={l:curve(l) for l in (6.92,7.45,10.374)}
for k in range(len(r)):
    vb=np.sqrt(max(d.vbar2.values[k],0)); c1=w[k]*(C[6.92][k]-vo[k])**2; c3=w[k]*(C[10.374][k]-vo[k])**2
    print('%6.2f %6.1f %4.1f %6.1f %6.1f %5.1f %6.1f | %6.1f %6.1f %6.1f | %6.2f %6.2f %s'%(r[k],vo[k],d.ev.values[k],d.vg.values[k],d.vd.values[k],d.vb.values[k],vb,C[6.92][k],C[7.45][k],C[10.374][k],c1,c3,'' if m[k] else '(masked)'))
inner=r<2*t.Rd
for l in (6.92,7.45,10.374):
    c=w*(C[l]-vo)**2; print('log rho_c %.2f: chi2 inner (r<2Rd) %.1f | outer %.1f | total %.1f'%(l,c[m&inner].sum(),c[m&~inner].sum(),c[m].sum()))
# OUTPUT 2026-11-12 22:33 EST: see ngc0289_diag_20261112_report.txt. Marginalised arm (sig 0.03): U -86.1 (zD -1.24), A -85.9, B (RAR fn) -86.9, NFW -83.8, NFW-free -84.1, DC14 -84.8, Burkert -85.9.
# AGAINST EXPECTATION: WRONG on the mechanism -- the six inner points (r < 7 kpc, errors 25-32 km/s) are MASKED by the Lelli cut, so the fit sees only r > 7.5 kpc;
# the high branch is bought by the OUTER shape (a peak ~185 at 10-14 kpc, a dip to ~156 at 20-26 kpc, a rise to 182 at 40-44 kpc and a decline to ~158 at 62-71 kpc).
# RIGHT that the galaxy is unusual: HI to 71 kpc (R_HI 49.5 kpc, ~10 R_d), M_HI 2.7e10; V_f 163 below the BTFR value 184 (0.88); Hubble-flow distance 20.8 +- 5.2 Mpc (25%).
