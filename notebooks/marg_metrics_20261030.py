# marg_metrics_20261030.py -- metrics M1-M4 and decisions D1-D3 of PREREG_MARG_20261030 from marg_fit_{params,points}_20261030.csv. [UNRATIFIED]
import numpy as np, pandas as pd
from scipy.stats import spearmanr
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; la0=np.log10(a0)
P=pd.read_csv('marg_fit_params_20261030.csv'); Q=pd.read_csv('marg_fit_points_20261030.csv'); Q=Q[Q.sel&(Q.vbar2>0)].copy()
Q['lgb']=np.log10(Q.vbar2/Q.r*1e6/kpc); Q['lgo']=np.log10(Q.vobs**2/Q.r*1e6/kpc); Q['lgm']=np.log10(Q.vmod2/Q.r*1e6/kpc); Q['res']=Q.lgo-Q.lgm
edges=np.arange(-2.25,1.01,0.25)+la0; win=[(lo,hi) for lo,hi in zip(edges[:-1],edges[1:]) if lo-la0>=-1.751 and hi-la0<=-0.249]
out=[]
def say(s): print(s); out.append(s)
say('PREREG_MARG_20261030 metrics -- arms A (relation, rho_c free) and B (RAR function); nuisances Ud, Ub, D, i with Li+2018-style priors; N gal A=%d B=%d'%(P[P.arm=='A'].gal.nunique(),P[P.arm=='B'].gal.nunique()))
# M1
say('\nM1 per-point rms log(g_obs/g_model) after marginalisation (Lelli points): ' + ' | '.join('%s %.3f (N %d)'%(a,np.sqrt((Q[Q.arm==a].res**2).mean()),len(Q[Q.arm==a])) for a in 'AB'))
say('   baselines at fixed nuisances: universal rho_c 0.155, relation(3) 0.140, RAR function 0.133; per-galaxy rho_c (Table 2) 0.081 vs RAR M/L free 0.079')
say('   window bins (log g_bar/a0): residual sd  A | B')
# M2: measurement-only scatter by injection on the fitted curves (D, i, 0.1 dex M/L, velocity errors; 50 realisations, seed 7), per bin
rng=np.random.default_rng(7); meas={a:[] for a in 'AB'}
par={(r.gal,r.arm):r for r in P.itertuples()}
for a in 'AB':
    d=Q[Q.arm==a].copy(); eD=d.gal.map(lambda g: par[(g,a)].eD); D=d.gal.map(lambda g: par[(g,a)].D); inc=d.gal.map(lambda g: par[(g,a)].inc); ei=d.gal.map(lambda g: par[(g,a)].ei)
    sc=[]
    for k in range(50):
        gp={g:(rng.normal(),rng.normal(),rng.normal(0,0.1)) for g in d.gal.unique()}; eps=np.array([gp[g] for g in d.gal])
        dD=np.log10(1+eps[:,0]*eD/D).clip(-0.3,0.3); i2=np.radians((inc+eps[:,1]*ei).clip(20,90)); di=2*np.log10(np.sin(np.radians(inc))/np.sin(i2))
        dv=2*rng.normal(size=len(d))*d.ev/d.vobs/np.log(10); dML=eps[:,2]   # M/L moves g_bar; here it moves the model's log g by the fraction set by the stellar share (approximated as full, upper bound)
        pert=-dD+di+dv
        sc.append([pert[(d.lgb>=lo)&(d.lgb<hi)].std() for lo,hi in win])
    meas[a]=np.array(sc).mean(axis=0)
for j,(lo,hi) in enumerate(win):
    sA=Q[(Q.arm=='A')&(Q.lgb>=lo)&(Q.lgb<hi)].res; sB=Q[(Q.arm=='B')&(Q.lgb>=lo)&(Q.lgb<hi)].res
    say('   %+5.2f..%+5.2f  N %4d  resid sd %.3f | %.3f   meas-only %.3f | %.3f   intrinsic est. %.3f | %.3f'%(lo-la0,hi-la0,len(sA),sA.std(),sB.std(),meas['A'][j],meas['B'][j],
        np.sqrt(max(sA.std()**2-meas['A'][j]**2,0)),np.sqrt(max(sB.std()**2-meas['B'][j]**2,0))))
M1A=np.sqrt((Q[Q.arm=='A'].res**2).mean()); M1B=np.sqrt((Q[Q.arm=='B'].res**2).mean())
intA=[np.sqrt(max(Q[(Q.arm=='A')&(Q.lgb>=lo)&(Q.lgb<hi)].res.std()**2-meas['A'][j]**2,0)) for j,(lo,hi) in enumerate(win)]
say('M2 note: the measurement-only term here excludes the M/L term (it enters g_bar, i.e. the x axis, not the residual at fixed curve) and velocity errors are propagated fully; after marginalisation the fit has already absorbed part of these errors, so the subtraction can over-correct toward zero.')
# M3
A=P[P.arm=='A'].copy(); A['Dp']=A.D/A.D0; A['lMb']=np.log10(A.Mb); A['Sh']=A.Mb/(2*np.pi*(1.678*A.Rd*A.Dp*1e3)**2); A['lSh']=np.log10(A.Sh); A['fg']=1.33*A.MHI*1e9*A.Dp**2/A.Mb
A=A[np.isfinite(A.lSh)&(A.Rd>0)]
X=np.column_stack([np.ones(len(A)),A.lMb,A.lSh,A.fg]); y=A.lrc.values; b=np.linalg.lstsq(X,y,rcond=None)[0]; res=y-X@b
bs=[]
for k in range(2000):
    i=rng.integers(0,len(A),len(A)); bs.append(np.linalg.lstsq(X[i],y[i],rcond=None)[0])
bs=np.array(bs); lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
say('\nM3 fitted log rho_c (arm A, N=%d): sd %.3f (before, fixed nuisances: 0.68); median %.2f; sd after removing M_b and Sigma_h trends %.3f (before 0.65)'%(len(A),A.lrc.std(),A.lrc.median(),np.std(y-np.column_stack([np.ones(len(A)),A.lMb,A.lSh])@np.linalg.lstsq(np.column_stack([np.ones(len(A)),A.lMb,A.lSh]),y,rcond=None)[0])))
say('   OLS log rho_c = %.2f %+.3f log Mb %+.3f log Sigma_h %+.3f f_gas ; 95%%: Mb [%+.3f,%+.3f] Sigma_h [%+.3f,%+.3f] f_gas [%+.3f,%+.3f]; resid sd %.3f'%(b[0],b[1],b[2],b[3],lo[1],hi[1],lo[2],hi[2],lo[3],hi[3],res.std()))
say('   before (relation 3, fixed nuisances): -0.365 / +0.268 [95%% from B.1: Sigma_h +0.14..+0.62 with mass only] / -0.93 [-1.71,-0.14]')
# M4
for a in 'AB':
    d=P[P.arm==a]; zu=(d.lud-np.log10(0.5))/0.1
    say('M4 arm %s nuisance shifts: dlogUd/0.1dex median %+.2f rms %.2f | zD median %+.2f rms %.2f | zi median %+.2f rms %.2f | |z|>2: Ud %d, D %d, i %d of %d'%(a,zu.median(),np.sqrt((zu**2).mean()),d.zD.median(),np.sqrt((d.zD**2).mean()),d.zi.median(),np.sqrt((d.zi**2).mean()),(abs(zu)>2).sum(),(abs(d.zD)>2).sum(),(abs(d.zi)>2).sum(),len(d)))
rmsA=np.sqrt(np.mean(np.concatenate([((P[P.arm=='A'].lud-np.log10(0.5))/0.1)**2,P[P.arm=='A'].zD**2,P[P.arm=='A'].zi**2])))
rmsB=np.sqrt(np.mean(np.concatenate([((P[P.arm=='B'].lud-np.log10(0.5))/0.1)**2,P[P.arm=='B'].zD**2,P[P.arm=='B'].zi**2])))
# decisions
d1='RESOLVED' if (M1A<=M1B+0.01 and max(intA)<=0.04) else ('CONFIRMED' if (M1A>M1B+0.03 or sum(x>0.06 for x in intA)>=2) else 'PARTLY')
d3=rmsA>1.5*rmsB
if d3 and d1=='RESOLVED': d1='PARTLY (D3: nuisance guard)'
d2='SURVIVE' if (lo[2]>0 and hi[3]<0) else ('ERROR-DRIVEN (one or both intervals include zero)')
say('\nD1 scatter tension: %s   (M1 A %.3f vs B %.3f; max window intrinsic A %.3f)'%(d1,M1A,M1B,max(intA)))
say('D2 rho_c patterns: %s   (Sigma_h slope %+.3f [%+.3f,%+.3f]; f_gas %+.3f [%+.3f,%+.3f])'%(d2,b[2],lo[2],hi[2],b[3],lo[3],hi[3]))
say('D3 nuisance guard: pooled rms z  A %.2f vs B %.2f -> %s'%(rmsA,rmsB,'TRIGGERED' if d3 else 'ok'))
open('marg_fit_report_20261030.txt','w').write('\n'.join(out))
