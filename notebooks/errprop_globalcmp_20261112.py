# errprop_globalcmp_20261112.py -- compare the global fixed-arm fit with pass2_table2_part.csv (PREREG_ERRPROP sec 7). [UNRATIFIED]
import numpy as np, pandas as pd
G=pd.read_csv('errprop_globalfit_20261112_part.csv').drop_duplicates('gal').set_index('gal')
P2=pd.read_csv('pass2_table2_part.csv').drop_duplicates('gal').set_index('gal')
J=G.join(P2[['lrc','chi_law']]); J['d']=J.lrc_g-J.lrc; J['dchi']=J.chi_g-J.chi_law
bad=J[(J.d.abs()>0.1)&(J.dchi<-0.01)]
print('N %d | pass2 missed the global minimum (|d|>0.1 dex AND chi2 lower): %d'%(len(J),len(bad)))
for g,r in bad.sort_values('dchi').iterrows(): print('   %-10s pass2 %.2f (chi2 %.1f) -> global %.2f (chi2 %.1f)'%(g,r.lrc,r.chi_law,r.lrc_g,r.chi_g))
print('galaxies with >= 2 local minima on the grid: %d of %d'%((J.n_localmin>=2).sum(),len(J)))
hw=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
for lab,x in (('pass2 (Brent)',J.lrc),('global',J.lrc_g)):
    print('%-14s median %.2f | std %.2f | robust half-width %.2f | at upper bound (>=10.45) %d | at lower (<=4.55) %d'%(lab,x.median(),x.std(ddof=1),hw(x),(x>=10.45).sum(),(x<=4.55).sum()))
print('at bound (global):', ', '.join(J.index[(J.lrc_g>=10.45)|(J.lrc_g<=4.55)]))
