# pass8_bootstrap_20261025.py [UNRATIFIED] 2026-10-25 00:55. Sim-review item C2 (simulated Desmond; test ours): galaxy-level bootstrap.
# (a) rms differences between Table 2 models, resampling whole galaxies (points within a galaxy are correlated);
# (b) residual trends at galaxy level: per-galaxy mean residual vs log Sigma_h, f_gas, log M_b; within-galaxy slope vs log r.
import numpy as np, pandas as pd, ast
from scipy.stats import spearmanr
T2=pd.read_csv('pass2_table2_part.csv').set_index('gal')
cols={'law_rho_c_free':'law','RAR_fixed':'rfix','RAR_ML_free':'rfree','RAR_ML_D_i':'r3'}
R={k:{g:np.array(ast.literal_eval(T2.loc[g,c])) for g in T2.index} for k,c in cols.items()}
P=pd.read_csv('pass3_relation3_points_20261022.csv')
for case,lab in (('univ','law_universal'),('rel3','law_relation3')):
    d=P[P.case==case]; R[lab]={g:x.res.values for g,x in d.groupby('gal')}
gals=sorted(set.intersection(*[set(v) for v in R.values()]))
print('galaxies in common',len(gals))
rms=lambda k,gs: np.sqrt(np.mean(np.concatenate([R[k][g] for g in gs])**2))
for k in R: print('  %-16s rms %.4f'%(k,rms(k,gals)))
rng=np.random.default_rng(7); B=2000
pairs=[('law_rho_c_free','RAR_ML_free'),('law_universal','RAR_fixed'),('law_relation3','RAR_fixed'),('law_relation3','law_universal')]
boot={p:[] for p in pairs}
for _ in range(B):
    gs=[gals[i] for i in rng.integers(0,len(gals),len(gals))]
    for a,b in pairs: boot[(a,b)].append(rms(a,gs)-rms(b,gs))
print('(a) rms difference (A - B), galaxy bootstrap, 2000 resamples:')
for (a,b),v in boot.items():
    v=np.array(v); print('  %-16s - %-14s %+.4f  68%% [%+.4f,%+.4f]  95%% [%+.4f,%+.4f]  P(A<B)=%.2f'%(a,b,rms(a,gals)-rms(b,gals),*np.percentile(v,[16,84,2.5,97.5]),np.mean(v<0)))
# (b) galaxy-level trends for the universal law
d=P[P.case=='univ']; Gm=d.groupby('gal').agg(res=('res','mean'),Sh=('Sh','first'),fg=('fgas','first'),Mb=('Mb','first'))
slopes={}
for g,x in d.groupby('gal'):
    if len(x)>=4: slopes[g]=np.polyfit(np.log10(x.r),x.res,1)[0]
S=pd.Series(slopes)
print('(b) universal law, galaxy level (N=%d):'%len(Gm))
for lab,col in (('log Sigma_h',np.log10(Gm.Sh)),('f_gas',Gm.fg),('log M_b',np.log10(Gm.Mb))):
    r,p=spearmanr(Gm.res,col); bs=[]
    for _ in range(B):
        i=rng.integers(0,len(Gm),len(Gm)); bs.append(spearmanr(Gm.res.values[i],col.values[i])[0])
    print('  mean residual vs %-12s Spearman %+.2f (p %.0e), bootstrap 95%% [%+.2f,%+.2f]'%(lab,r,p,*np.percentile(bs,[2.5,97.5])))
bs=[np.median(S.values[rng.integers(0,len(S),len(S))]) for _ in range(B)]
print('  within-galaxy slope d(res)/d(log r): median %+.3f, bootstrap 95%% [%+.3f,%+.3f], N=%d galaxies (>=4 pts)'%(S.median(),*np.percentile(bs,[2.5,97.5]),len(S)))
