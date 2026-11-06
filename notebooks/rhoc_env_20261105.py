# rhoc_env_20261105.py [UNRATIFIED] -- PREREG_RHOC_ENV_20261105: fitted log rho_c (pass2_table2_part.csv, fixed M/L) vs Chae+2021 Table 3 log e_N,env.
import numpy as np, pandas as pd, io, contextlib, warnings; warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
rng=np.random.default_rng(20261105)
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
fit=pd.read_csv('pass2_table2_part.csv').set_index('gal')[['lrc','chi_law']]
env=pd.read_csv('../data/chae2021/chae2021_table3_eNenv.csv').set_index('galaxy')
G2=T.copy(); G2['Mb']=0.5*G2.L*1e9+1.33*G2.MHI*1e9; G2['lMb']=np.log10(G2.Mb); G2['fgas']=1.33*G2.MHI*1e9/G2.Mb
G2['lSig']=np.log10(G2.Mb/(2*np.pi*(1.678*G2.Rd*1e3)**2))
D=fit.join(G2[['lMb','lSig','fgas']],how='inner').join(env,how='inner').dropna(subset=['lrc','log_eN_env_maxclust'])
print('N galaxies with fitted rho_c and Chae e_N,env:',len(D))
def boot(fn,n=1000):
    v=[fn(D.iloc[rng.integers(0,len(D),len(D))]) for _ in range(n)]; return np.percentile(v,[2.5,97.5])
for col in ('log_eN_env_maxclust','log_eN_env_noclust'):
    r,p=spearmanr(D.lrc,D[col]); ci=boot(lambda d: spearmanr(d.lrc,d[col])[0])
    print('Q1 Spearman(log rho_c, %s) = %+.3f (p=%.3f), 95%% CI [%+.3f,%+.3f]'%(col,r,p,*ci))
    def ols(d,cols):
        X=np.column_stack([np.ones(len(d))]+[d[c].values for c in cols]); b=np.linalg.lstsq(X,d.lrc.values,rcond=None)[0]; return b[1]
    for cols,lab in (([col],'alone'),([col,'lMb','lSig','fgas'],'+ log Mb, log Sigma_h, f_gas')):
        b=ols(D,cols); ci=boot(lambda d: ols(d,cols)); print('Q2 slope d log rho_c / d %s (%s) = %+.3f dex/dex, 95%% CI [%+.3f,%+.3f]'%(col,lab,b,*ci))
    q=D[col].quantile([1/3,2/3]); lo=D[D[col]<=q.iloc[0]].lrc; hi=D[D[col]>=q.iloc[1]].lrc
    diff=hi.median()-lo.median()
    ci=boot(lambda d: d[d[col]>=d[col].quantile(2/3)].lrc.median()-d[d[col]<=d[col].quantile(1/3)].lrc.median())
    print('Q3 tercile contrast high-low median log rho_c = %+.3f (N %d/%d), 95%% CI [%+.3f,%+.3f]'%(diff,len(hi),len(lo),*ci))
print('\nmedian log rho_c %.2f, spread (16-84) %.2f dex; log e_N,env range %.2f..%.2f'%(D.lrc.median(),np.diff(D.lrc.quantile([.16,.84]))[0],D.log_eN_env_maxclust.min(),D.log_eN_env_maxclust.max()))
D.to_csv('rhoc_env_20261105.csv')
