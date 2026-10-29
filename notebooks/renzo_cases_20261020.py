# renzo_cases_20261020.py [UNRATIFIED] 2026-10-20 18:50 -- Renzo's rule on the canonical cases (Sancisi 2004: NGC 1560, NGC 6946, IC 2574; plus NGC 3198, NGC 2403
# as controls). From the N12 per-point table. Per galaxy: features = log V minus quadratic in log r; r(model, obs); and at the radius of the strongest
# baryonic feature, the feature amplitude in obs / GPRF fitted / GPRF universal / MOND / baryons-only. DECLARED: 'addressed' if, in the three canonical cases, the
# GPRF feature at the baryonic-feature radius has the observed sign and >= 50% of the observed amplitude, and r(GPRF) >= r(MOND) - 0.1.
import numpy as np, pandas as pd
from gate2c_robustness_20260925 import kpc, a0
D=pd.read_csv('n12_emergent_nu_20261020.csv')
nurar=lambda lgb: -np.log10(1-np.exp(-np.sqrt(10**lgb/a0)))
D['lV_obs']=0.5*(D.lgo+np.log10(D.r*kpc/1e6)); D['lV_mod']=0.5*(D.lgm+np.log10(D.r*kpc/1e6)); D['lV_bar']=0.5*(D.lgb+np.log10(D.r*kpc/1e6)); D['lV_mond']=0.5*(D.lgb+nurar(D.lgb)+np.log10(D.r*kpc/1e6))
def feat(x,y): return y-np.polyval(np.polyfit(x,y,2),x)
cases=['NGC1560','NGC6946','IC2574','NGC3198','NGC2403','NGC2841','DDO154','NGC7793']
present=[g for g in cases if g in set(D.gal)]; print('in sample:',present,'| missing:',[g for g in cases if g not in present])
for g in present:
    d=D[(D.gal==g)&(D.case=='fit')].sort_values('r'); u=D[(D.gal==g)&(D.case=='u7.45')].sort_values('r'); x=np.log10(d.r.values)
    fo=feat(x,d.lV_obs.values); ff=feat(x,d.lV_mod.values); fu=feat(x,u.lV_mod.values); fm=feat(x,d.lV_mond.values); fb=feat(x,d.lV_bar.values)
    i=np.argmax(np.abs(fb)); 
    print('%-8s n=%2d  r(fit,obs)=%+.2f r(univ,obs)=%+.2f r(MOND,obs)=%+.2f r(bar,obs)=%+.2f | strongest baryon feature at r=%.2f kpc: bar %+.3f obs %+.3f fit %+.3f univ %+.3f MOND %+.3f  (dex)'%(
        g,len(d),np.corrcoef(ff,fo)[0,1],np.corrcoef(fu,fo)[0,1],np.corrcoef(fm,fo)[0,1],np.corrcoef(fb,fo)[0,1],d.r.values[i],fb[i],fo[i],ff[i],fu[i],fm[i]))
    print('          point-by-point (r kpc: obs / fit / MOND / bar, mdex):',' '.join('%.1f:%+.0f/%+.0f/%+.0f/%+.0f'%(r_,1e3*a,1e3*b,1e3*c,1e3*e) for r_,a,b,c,e in zip(d.r.values,fo,ff,fm,fb)))
