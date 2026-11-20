# errprop_analysis_20261112.py -- analysis of errprop_20261112_part.csv per PREREG_ERRPROP_20261112 sec 4-6. [UNRATIFIED]
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr
R=pd.read_csv('errprop_20261112_part.csv').drop_duplicates('gal').set_index('gal')
ref=pd.read_csv('pass2_table2_part.csv').set_index('gal').lrc
fD={}
for l in open('../sparc_galaxy_data/SPARC_Lelli2016c.mrt'):
    p=l.split()
    if len(p)>=18 and p[1].isdigit():
        try: fD[p[0]]=int(p[4])
        except: pass
R['fD']=[fD.get(g,0) for g in R.index]
dv=(R.lrc0-ref.reindex(R.index)).abs(); print('VALIDATION: %d of %d within 0.03 of pass2 (max |diff| %.3f)'%((dv<=0.03).sum(),len(R),dv.max()))
R['sD']=(R.lrcDp-R.lrcDm).abs()/2; R['sI']=(R.lrcI-R.lrc0).abs()
R['sU']=np.sqrt((R.lrcU-R.lrc0)**2+np.nan_to_num((R.lrcUb-R.lrc0).values)**2)
R['stot']=np.sqrt(R.sD**2+R.sI**2+R.sU**2)
cols=['lrc0','lrcDp','lrcDm','lrcI','lrcU','lrcUb']; R['flag']=(R[cols].ge(10.45)|R[cols].le(4.55)).any(axis=1)
print('flagged at a bound: %d -> %s'%(R.flag.sum(),', '.join(R.index[R.flag])))
U=R[~R.flag]
hw=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
print('median sigma (dex), unflagged N=%d: total %.2f | distance %.2f | inclination %.2f | M/L %.2f ; 16-84 of total %.2f-%.2f'%(len(U),U.stot.median(),U.sD.median(),U.sI.median(),U.sU.median(),np.percentile(U.stot,16),np.percentile(U.stot,84)))
print('robust spread lrc0: all %.2f | unflagged %.2f ; std unflagged %.2f'%(hw(R.lrc0),hw(U.lrc0),U.lrc0.std(ddof=1)))
def tau_mle(x,s):
    f=lambda p: 0.5*np.sum((x-p[0])**2/(s**2+p[1]**2)+np.log(s**2+p[1]**2))
    r=minimize(f,[np.median(x),0.3],method='Nelder-Mead'); return r.x
for lab,sub in (('unflagged',U),('direct distance (f_D != 1)',U[U.fD!=1]),('Hubble flow (f_D = 1)',U[U.fD==1])):
    mu,tau=tau_mle(sub.lrc0.values,sub.stot.values)
    print('%-28s N %3d  median sigma_tot %.2f  robust spread %.2f  MLE mu %.2f tau %.2f'%(lab,len(sub),sub.stot.median(),hw(sub.lrc0),mu,abs(tau)))
rho,p=spearmanr(U.stot,(U.lrc0-U.lrc0.median()).abs()); print('Spearman sigma_tot vs |lrc0 - median|: %.2f (p %.3f)'%(rho,p))
print('pinned three:', '; '.join('%s lrc0 %.2f D+/- %.2f/%.2f I %.2f U %.2f'%(g,*R.loc[g,['lrc0','lrcDp','lrcDm','lrcI','lrcU']]) for g in ('NGC4013','NGC5907','UGC01281')))
