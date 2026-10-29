# n12e_mond_Di_free_20261020.py [UNRATIFIED] PRE-REGISTERED 2026-10-20 19:45 (cold review A3; Li+2018 parity, MOND half). MOND RAR with Upsilon_disk, D, i FREE per
# galaxy under Gaussian priors from SPARC (e_D, e_Inc; log Y_disk 0.5 +/- 0.1 dex), same 120 galaxies and 10% velocity weighting. Scaling: V_bar components ~ sqrt(D)
# (M ~ D^2, r ~ D), radii ~ D; V_obs ~ 1/sin(i) relative to the catalogue inclination. Compare per-point rms log(g_obs/g_model) at the best fit (in the fitted frame)
# with GPRF rho_c-only 0.108 (Y, D, i fixed) and MOND Y-only 0.125. DECLARED: report; the matched GPRF 3-parameter run is NOT done (cost ~8e3 solves) and is OWED
# before the one-knob sentence is used in any form.
import numpy as np, pandas as pd
from scipy.optimize import minimize
from gate2c_robustness_20260925 import load, kpc, a0
rows=[]
for line in open('../sparc_galaxy_data/SPARC_Lelli2016c.mrt'):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: rows.append(dict(gal=p[0],D=float(p[2]),eD=float(p[3]),inc=float(p[5]),einc=float(p[6])))
        except ValueError: pass
E=pd.DataFrame(rows).set_index('gal'); P=load(); N12=pd.read_csv('n12_emergent_nu_20261020.csv'); gals=sorted(set(N12.gal))
nu=lambda gb: 1/(1-np.exp(-np.sqrt(gb/a0)))
out=[]
for g in gals:
    d=P[P.gal==g].sort_values('r'); r0=d.r.values; vobs0=d.vobs.values; w=1/(0.1*vobs0)**2
    base=np.sign(d.vg.values)*d.vg.values**2+0.7*d.vb.values**2; vd2=d.vd.values**2
    D0,eD,i0,ei=E.loc[g,['D','eD','inc','einc']]; i0r=np.radians(i0)
    def model(th):
        Y,fD,di=th; r=r0*fD; vb2=np.clip((base+Y*vd2)*fD,1e-3,None); gb=vb2/r*1e6/kpc
        vm=np.sqrt(gb*nu(gb)*r*kpc/1e6); vo=vobs0*np.sin(i0r)/np.sin(i0r+np.radians(di)); return vm,vo,gb
    def chi(th):
        Y,fD,di=th
        if Y<=0.05 or fD<=0.2 or not (5<i0+di<89): return 1e12
        vm,vo,_=model(th); return np.sum(w*(vm-vo)**2)+(np.log10(Y/0.5)/0.1)**2+((fD-1)*D0/eD)**2+(di/ei)**2
    res=minimize(chi,[0.5,1.0,0.0],method='Nelder-Mead',options={'xatol':1e-3,'fatol':1e-2,'maxiter':600})
    Y,fD,di=res.x; vm,vo,gb=model(res.x); lg=np.log10((vo**2)/(vm**2))
    out.append(dict(gal=g,Y=Y,fD=fD,di=di,rms=np.sqrt(np.mean(lg**2)),chi=np.sum(w*(vm-vo)**2),n=len(r0),res=lg))
R=pd.DataFrame(out).set_index('gal')
allrms=np.sqrt(np.mean(np.concatenate(R.res.values)**2))
print('MOND RAR, Y_disk + D + i free (SPARC priors), N=%d: per-point rms %.3f | per-galaxy median rms %.3f | chi median %.1f'%(len(R),allrms,R.rms.median(),R.chi.median()))
print('  fitted: Y median %.2f (16-84: %.2f-%.2f); D/D0 median %.2f (%.2f-%.2f); d_inc median %+.1f deg (%+.1f..%+.1f)'%(R.Y.median(),*np.percentile(R.Y,[16,84]),R.fD.median(),*np.percentile(R.fD,[16,84]),R.di.median(),*np.percentile(R.di,[16,84])))
print('  reference: MOND Y fixed 0.159 | MOND Y free 0.125 | GPRF rho_c free (Y,D,i fixed) 0.108 | GPRF universal 0.172')
R.drop(columns=['res']).to_csv('n12e_mond_Di_free_20261020.csv')
