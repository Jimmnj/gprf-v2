# pass5b_overshoot_20261022.py [UNRATIFIED] 22:15. Residual of the universal-rho_c law at R_last vs R_last/r0 (King radius r0 = 3 sigma/sqrt(4 pi G rho_c)).
# Isothermal v^2 overshoots its asymptote 2 sigma^2 near a few r0 (F11, 9/26); if LSB last points sit there, the LSB excess could be the overshoot.
import numpy as np, pandas as pd
from scipy.stats import spearmanr
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6
A=pd.read_csv('pass3_relation3_points_20261022.csv')
for case in ('univ','rel3'):
    D=A[A.case==case]; L=D.loc[D.groupby('gal').r.idxmax()].set_index('gal')
    s2=np.sqrt(a0k*G*L.Mb)/2; rc=10**L.lrc; r0=3*np.sqrt(s2/(4*np.pi*G*rc))
    x=np.log10(L.r/r0); y=np.log10(np.sqrt(L.vmod2)/L.vobs)
    r,p=spearmanr(y,x); rs,ps=spearmanr(x,np.log10(L.Sh))
    q=pd.qcut(x,4,labels=False)
    print('%s: log(R_last/r0) median %.2f (16-84 %.2f..%.2f); Spearman resid vs it %+.2f p=%.0e; vs Sigma_h it correlates %+.2f'%(case,x.median(),*x.quantile([.16,.84]),r,p,rs))
    print('   median V_model/V_obs by R_last/r0 quartile (inner->outer):',' '.join('%.3f'%10**y[q==k].median() for k in range(4)),' | quartile medians of R_last/r0:',' '.join('%.1f'%10**x[q==k].median() for k in range(4)))
    X=np.column_stack([np.ones(len(L)),x,np.log10(L.Sh)]); b=np.linalg.lstsq(X,y.values,rcond=None)[0]
    rng=np.random.default_rng(3); B=[np.linalg.lstsq(X[i],y.values[i],rcond=None)[0] for i in (rng.integers(0,len(L),len(L)) for _ in range(2000))]
    print('   OLS resid = %+.3f %+.3f(+-%.3f) log(R/r0) %+.3f(+-%.3f) log Sh'%(b[0],b[1],np.std(B,0)[1],b[2],np.std(B,0)[2]))
# quadratic control for the non-monotonic overshoot hump
print('--- quadratic in log(R/r0) + log Sh ---')
for case in ('univ','rel3'):
    D=A[A.case==case]; L=D.loc[D.groupby('gal').r.idxmax()].set_index('gal')
    s2=np.sqrt(a0k*G*L.Mb)/2; r0=3*np.sqrt(s2/(4*np.pi*G*10**L.lrc)); x=np.log10(L.r/r0); y=np.log10(np.sqrt(L.vmod2)/L.vobs).values
    X=np.column_stack([np.ones(len(L)),x,x**2,np.log10(L.Sh)]); b=np.linalg.lstsq(X,y,rcond=None)[0]
    rng=np.random.default_rng(4); B=np.array([np.linalg.lstsq(X[i],y[i],rcond=None)[0] for i in (rng.integers(0,len(L),len(L)) for _ in range(2000))])
    print('%s: x %+.3f(+-%.3f) x^2 %+.3f(+-%.3f) logSh %+.3f(+-%.3f); hump peak at R/r0 = %.1f'%(case,b[1],B[:,1].std(),b[2],B[:,2].std(),b[3],B[:,3].std(),10**(-b[1]/(2*b[2])) if b[2]<0 else np.nan))
