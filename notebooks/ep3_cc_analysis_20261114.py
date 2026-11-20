# ep3_cc_analysis_20261114.py -- EP-3: clean_conc_20261109 summary block, old Brent fits vs V3 first-minimum fits. [UNRATIFIED]
import numpy as np, pandas as pd
from scipy.stats import spearmanr
C=pd.read_csv('clean_conc_part_20261109.csv').drop_duplicates('gal').set_index('gal'); E=pd.read_csv('ep3_cc_20261114_part.csv').set_index('key')
for d in (0,2):
    x=(E['fix%d_v3'%d]-E['fix%d_old'%d]); print('fix%d: N %d | |diff|>0.05: %d | edge %d'%(d,x.notna().sum(),(x.abs()>0.05).sum(),E['edge%d'%d].sum()))
    print(E.loc[x.abs()>0.05,['fix%d_old'%d,'fix%d_v3'%d]].round(2).to_string())
def resid(y,x): X=np.column_stack([np.ones(len(y)),x]); return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
def line(lab,D,y,x):
    D=D[[y,x,'lMb']].dropna(); n=len(D); lrho,lp=spearmanr(np.log10(D[x]),D[y])
    rp,pp=spearmanr(resid(D[y].values,D.lMb.values),resid(np.log10(D[x].values),D.lMb.values))
    return '  %-44s N %3d  rho %+.2f (p %.3f) | at fixed M_b %+.2f (p %.3f)'%(lab,n,lrho,lp,rp,pp)
for tag in ('OLD','V3'):
    R=C.copy()
    if tag=='V3': R['lrc_fix0']=E.fix0_v3; R['lrc_fix2']=E.fix2_v3
    Q1=R[R.Q==1]; print('==',tag)
    print(line('control r_peak/Rd all pts, fix0',R,'lrc_fix0','pk0'))
    print(line('A2 Q1 r_peak/Rd r[2:], fix2',Q1,'lrc_fix2','pk2')); print(line('A2 Q1 r_peak/Rd r[2:], fix0',Q1,'lrc_fix0','pk2'))
    print(line('A3 Q1 r_half/Rd r[2:], fix2',Q1,'lrc_fix2','hm2'))
    print(line('all Q r_peak/Rd r[2:], fix2',R,'lrc_fix2','pk2')); print(line('all Q r_half/Rd r[2:], fix2',R,'lrc_fix2','hm2'))
    print('  refit vs marg: Spearman %.2f; median |fix0-marg| %.2f; median |fix2-fix0| %.2f; fix0 median %.2f sd %.2f'%(spearmanr(R.lrc_fix0,R.lrc_marg)[0],(R.lrc_fix0-R.lrc_marg).abs().median(),(R.lrc_fix2-R.lrc_fix0).abs().median(),R.lrc_fix0.median(),R.lrc_fix0.std()))
    for lo,hi in ((0,1.5),(1.5,3),(3,99)): s=Q1[(Q1.pk2>=lo)&(Q1.pk2<hi)]; print('   bin %4.1f-%4.1f N %3d  median fix2 %.2f'%(lo,hi,len(s),s.lrc_fix2.median()))
