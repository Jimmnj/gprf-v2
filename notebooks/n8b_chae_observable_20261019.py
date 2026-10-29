# n8b_chae_observable_20261019.py [UNRATIFIED]  PRE-REGISTERED 2026-10-19 (ballot EFE-1 / SLATE)
# N8b: Chae's observable on SPARC vs Chae+2021 Table 3 environmental field e_N,env.
# Observables per galaxy: (A) outer RAR residual = mean of lgo - lgmond over points with g_bar < a0/3 (low-acceleration regime,
#     where MOND-EFE predicts a DIP growing with e_N; Chae+2020 Fig. style); (B) outer RC log-slope d log V / d log R over the
#     outer third of points. Regress each on log e_N (max clustering) controlling log M_b and log Sigma_b (M_b/2piRd^2).
# DECLARED: GPRF frozen law + NEST-1 C -> coefficient on log e_N = 0 for both. MOND EFE -> (A) negative (residual falls with e_N),
#     (B) negative (curves decline more in denser fields). Significance bar: 95% bootstrap CI excludes 0.
# Controls checked both ways (with/without Sigma_b). Also reported, NOT a test: Chae's RC-fitted e-tilde vs the same observables (circular).
import numpy as np, pandas as pd
from gate2c_robustness_20260925 import load, load_table, a0
P=load(); T=load_table(); rng=np.random.default_rng(11)
env=pd.read_csv('../data/chae2021/chae2021_table3_eNenv.csv').set_index('galaxy')
fit=pd.read_csv('../data/chae2021/chae2021_table2_fits.csv').set_index('galaxy')
rows=[]
for g,d in P.groupby('gal'):
    d=d.sort_values('r'); n=len(d)
    if n<6: continue
    low=d[d.gbar<a0/3]
    resA=low.res.mean() if len(low)>=3 else np.nan
    o=d.iloc[-max(3,n//3):]
    sl=np.polyfit(np.log10(o.r),np.log10(o.vobs),1)[0]
    Mb=0.5*d.L.iloc[0]*1e9+1.33*d.MHI.iloc[0]*1e9; Rd=d.Rd.iloc[0]
    rows.append(dict(gal=g,Q=d.Q.iloc[0],npts=n,nlow=len(low),resA=resA,slopeB=sl,lMb=np.log10(Mb),lSig=np.log10(Mb/(2*np.pi*(Rd*1e3)**2)) if Rd>0 else np.nan))
D=pd.DataFrame(rows).set_index('gal').join(env,how='inner').join(fit[['etilde']],how='left')
D.to_csv('n8b_chae_observable_20261019.csv')
def reg(y,cols,lab):
    m=np.isfinite(y)
    for c in cols: m&=np.isfinite(D[c].values)
    X=np.column_stack([np.ones(m.sum())]+[D[c].values[m] for c in cols]); yy=y[m]
    b=np.linalg.lstsq(X,yy,rcond=None)[0]
    bs=np.array([np.linalg.lstsq(X[i],yy[i],rcond=None)[0] for i in (rng.integers(0,len(yy),len(yy)) for _ in range(3000))])
    lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
    s=' | '.join('%s %+.3f [%+.3f,%+.3f]'%(c,b[k+1],lo[k+1],hi[k+1]) for k,c in enumerate(cols))
    print('  %-34s N=%3d  %s  resid %.3f'%(lab,m.sum(),s,(yy-X@b).std()))
print('N8b  galaxies with Chae e_N:',len(D),' (Q<=2, inc>30, >=6 pts); with >=3 low-acc points:',int(np.isfinite(D.resA).sum()))
for nm,y in (('(A) outer RAR residual [dex]',D.resA.values),('(B) outer RC log-slope',D.slopeB.values)):
    print(nm)
    reg(y,['log_eN_env_maxclust'],'e_N alone')
    reg(y,['log_eN_env_maxclust','lMb'],'e_N + log Mb')
    reg(y,['log_eN_env_maxclust','lMb','lSig'],'e_N + log Mb + log Sigma_b')
    reg(y,['log_eN_env_noclust','lMb','lSig'],'e_N(noclust) + Mb + Sigma_b')
    reg(y,['etilde','lMb','lSig'],'[circular] Chae e-tilde + Mb + Sigma_b')
    q1=(D.Q==1).values; print('  Q=1 only:'); 
    Dq=D[q1]; 
    m=np.isfinite(y[q1])&np.isfinite(Dq.lSig.values); X=np.column_stack([np.ones(m.sum()),Dq.log_eN_env_maxclust.values[m],Dq.lMb.values[m],Dq.lSig.values[m]]); yy=y[q1][m]
    b=np.linalg.lstsq(X,yy,rcond=None)[0]; bs=np.array([np.linalg.lstsq(X[i],yy[i],rcond=None)[0] for i in (rng.integers(0,len(yy),len(yy)) for _ in range(3000))]); lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
    print('    e_N + Mb + Sigma_b (Q=1)            N=%3d  e_N %+.3f [%+.3f,%+.3f]'%(m.sum(),b[1],lo[1],hi[1]))
lo,hi=np.percentile(D.log_eN_env_maxclust,[5,95]); print('e_N range 5-95%%: %.2f..%.2f dex (x%.1f); resA std %.3f; slopeB std %.3f'%(lo,hi,10**(hi-lo),np.nanstd(D.resA),D.slopeB.std()))
