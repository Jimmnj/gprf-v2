# ep3_rule3_regress_20261114.py -- EP-3 [UNRATIFIED]: (a) the 120-galaxy regression behind rule (3) and B.1's slopes, with test1 fits (2026-10-18,
# bounds 4-11; two second-branch fits) vs the V3 table; (b) fig4_mu0_regress two-variable fit with V3 values. Read-only.
import numpy as np, pandas as pd
from gate2c_robustness_20260925 import load_table
G=4.30091e-6; rng=np.random.default_rng(1)
F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal'); T=load_table(); V3=pd.read_csv('pass2_table2_part.csv').set_index('gal').lrc
F['lMb']=np.log10(F.Mb); F['lSh']=np.log10(F.Mb/(2*np.pi*(1.678*T.loc[F.index,'Rd'].values*1e3)**2)); F['fg']=1.33*T.loc[F.index,'MHI'].values*1e9/F.Mb
F['v3']=F.index.map(V3); F['v3']=F.v3.fillna(F.l_tanh2)
def ols(y,cols,D):
    X=np.column_stack([np.ones(len(D))]+[D[c].values for c in cols]); b=np.linalg.lstsq(X,D[y].values,rcond=None)[0]
    bs=[]
    for _ in range(2000):
        i=rng.integers(0,len(D),len(D)); bs.append(np.linalg.lstsq(X[i],D[y].values[i],rcond=None)[0])
    lo,hi=np.percentile(bs,[2.5,97.5],axis=0); r=D[y].values-X@b; return b,lo,hi,r.std()
print('(a) 120-galaxy regressions of log rho_c  [rule (3) as frozen: 11.063 -0.365 lMb +0.268 lSh -0.929 fgas]')
for y in ('l_tanh2','v3'):
    for cols in (['lMb','lSh'],['lMb','lSh','fg']):
        b,lo,hi,s=ols(y,cols,F); print('  %-8s %-18s '%(y,'+'.join(cols))+' '.join('%s %+.3f [%+.2f,%+.2f]'%(c,b[k+1],lo[k+1],hi[k+1]) for k,c in enumerate(cols))+' | const %.3f | resid sd %.3f of %.3f'%(b[0],s,F[y].std()))
print('(b) two-variable fit of log mu0 (fig4_mu0_regress)')
B=pd.read_csv('test1c_burkert_data_20261018.csv',index_col=0); B=B[~B.r0_gt_rmax.astype(bool)]; B=B[B.index.isin(F.index)]
for y in ('l_tanh2','v3'):
    F['mu']=np.log10(3*np.sqrt(F.s2b)*np.sqrt(10**F[y]/(4*np.pi*G))/1e6)
    for lab,D in (('Burkert subset',F.loc[B.index]),('all 120',F)):
        b,lo,hi,s=ols('mu',['lMb','lSh'],D); print('  %-8s %-15s N=%3d dlogMb %+.2f [%+.2f,%+.2f] dlogSh %+.2f [%+.2f,%+.2f]'%(y,lab,len(D),b[1],lo[1],hi[1],b[2],lo[2],hi[2]))
F['mu_d']=np.nan; F.loc[B.index,'mu_d']=B.mu_data; b,lo,hi,s=ols('mu_d',['lMb','lSh'],F.loc[B.index]); print('  data     Burkert         N=%3d dlogMb %+.2f [%+.2f,%+.2f] dlogSh %+.2f [%+.2f,%+.2f]'%(len(B),b[1],lo[1],hi[1],b[2],lo[2],hi[2]))
