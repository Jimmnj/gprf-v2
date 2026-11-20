# errprop3_analysis_20261113.py -- BD-3 analysis: ERRPROP stats with the 12 multi-minimum galaxies refit under V3 (first minimum). [UNRATIFIED]
# Stats code copied from errprop_analysis_20261112.py (tau_mle, hw, sD/sI/sU/stot, flag rule) so before/after differ only in the 12 rows.
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr
OLD=pd.read_csv('errprop2_20261112_part.csv').drop_duplicates('gal').set_index('gal')
N3=pd.read_csv('errprop3_first_20261113_part.csv').set_index('gal'); F=pd.read_csv('bound_check_first_20261112.csv').set_index('gal')
fD={}
for l in open('../sparc_galaxy_data/SPARC_Lelli2016c.mrt'):
    p=l.split()
    if len(p)>=18 and p[1].isdigit():
        try: fD[p[0]]=int(p[4])
        except: pass
cols=['lrc0','lrcDp','lrcDm','lrcI','lrcU','lrcUb']
def prep(R,edge=None):
    R=R.copy(); R['fD']=[fD.get(g,0) for g in R.index]
    R['sD']=(R.lrcDp-R.lrcDm).abs()/2; R['sI']=(R.lrcI-R.lrc0).abs()
    R['sU']=np.sqrt((R.lrcU-R.lrc0)**2+np.nan_to_num((R.lrcUb-R.lrc0).values)**2); R['stot']=np.sqrt(R.sD**2+R.sI**2+R.sU**2)
    R['flag']=(R[cols].ge(10.45)|R[cols].le(4.55)).any(axis=1)
    if edge is not None:  # V3 rows: grid runs to 13, so flag only a missing interior minimum or the lower bound
        for g in edge.index: R.loc[g,'flag']=bool(edge[g]) or bool((R.loc[g,cols].dropna()<=4.55).any())
    return R
hw=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
def tau_mle(x,s):
    f=lambda p: 0.5*np.sum((x-p[0])**2/(s**2+p[1]**2)+np.log(s**2+p[1]**2))
    r=minimize(f,[np.median(x),0.3],method='Nelder-Mead'); return r.x
def report(lab,R):
    U=R[~R.flag]; print('== %s: flagged %d -> %s'%(lab,R.flag.sum(),', '.join(R.index[R.flag])))
    print('   median sigma_tot %.2f (D %.2f I %.2f U %.2f); robust spread all %.2f unflagged %.2f; std unflagged %.2f; median lrc0 %.3f'%(U.stot.median(),U.sD.median(),U.sI.median(),U.sU.median(),hw(R.lrc0),hw(U.lrc0),U.lrc0.std(ddof=1),U.lrc0.median()))
    for l2,sub in (('unflagged',U),('direct (f_D!=1)',U[U.fD!=1]),('Hubble (f_D=1)',U[U.fD==1])):
        mu,tau=tau_mle(sub.lrc0.values,sub.stot.values); print('   %-16s N %3d  sigma_tot %.2f  spread %.2f  mu %.3f  tau %.3f'%(l2,len(sub),sub.stot.median(),hw(sub.lrc0),mu,abs(tau)))
    r,p=spearmanr(U.stot,(U.lrc0-U.lrc0.median()).abs()); print('   Spearman sigma_tot vs |dev|: %.2f (p %.3f)'%(r,p))
A=prep(OLD); NEW=OLD.copy()
for g in N3.index: NEW.loc[g,cols]=N3.loc[g,cols].values
B=prep(NEW,edge=N3.edge)
print('X1 baseline vs bound_check_first: max |diff| %.4f dex'%(N3.lrc0-F.lrc_first.reindex(N3.index)).abs().max())
print('\nThe 12 (old errprop2 -> V3):'); 
for g in N3.index: print('  %-9s lrc0 %5.2f->%5.2f  stot %5.2f->%5.2f  edge %d  V3 fits %s'%(g,A.loc[g,'lrc0'],B.loc[g,'lrc0'],A.loc[g,'stot'],B.loc[g,'stot'],N3.loc[g,'edge'],' '.join('%.2f'%x for x in N3.loc[g,cols].dropna())))
rest=OLD.drop(N3.index); jump=(rest[cols].sub(rest.lrc0,axis=0)).abs().max(axis=1)
print('\nX4 outside the V3 set: max perturbed jump %.2f dex (%s); galaxies with jump > 2.5: %s'%(jump.max(),jump.idxmax(),list(jump.index[jump>2.5])))
print(); report('BEFORE (errprop2, global min in [4.5,10.5])',A); report('AFTER BD-3 (17 rows refit under V3: the 12 + X4 screen DDO064 NGC3521 NGC5371 NGC6674 NGC6946)',B)
