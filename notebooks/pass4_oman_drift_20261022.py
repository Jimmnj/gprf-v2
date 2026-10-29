# pass4_oman_drift_20261022.py [UNRATIFIED] 2026-10-22 21:05. Handoff 10/22 sec 4 item 1 (check suggested in a SIMULATED Oman review; the check is ours).
# Question: is the LSB outer-speed over-prediction (law V_model > V_obs at R_last, trend with Sigma_h) an artefact of uncorrected
# pressure support (asymmetric drift) in the HI? Drift: V_c^2 = V_rot^2 + sigma^2 * n,  n = -dln(Sigma_gas sigma^2)/dlnR (sigma const).
# Outer HI disk exponential with h_g = R_HI/2  ->  n(R) = 2R/R_HI. Brackets: sigma 6/8/12 km/s; n x0.5/x1/x1.5.
# Inputs: pass3_relation3_points_20261022.csv (2696 pts x {univ, rel3}); SPARC MRT for R_HI. No refit: rho_c as in pass 2/3.
import numpy as np, pandas as pd
from scipy.stats import spearmanr
MRT='../sparc_galaxy_data/SPARC_Lelli2016c.mrt'
R={}
for line in open(MRT):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: R[p[0]]=float(p[14])
        except ValueError: pass
A=pd.read_csv('pass3_relation3_points_20261022.csv'); A['RHI']=A.gal.map(R)
print('galaxies',A.gal.nunique(),' missing/zero R_HI:',A[(A.RHI.isna())|(A.RHI<=0)].gal.nunique())
A.loc[A.RHI<=0,'RHI']=np.nan
lg=lambda v2,r: np.log10(v2/r)+6-19.48940
chk=A[A.case=='univ'].copy(); 
V=pd.read_csv('pass2_vlast_20261022.csv').set_index('gal')
def lastpts(D):
    return D.loc[D.groupby('gal').r.idxmax()].set_index('gal')
L0=lastpts(chk); rr=np.sqrt(L0.vmod2)/L0.vobs
print('check vs pass2 ratio (universal): max |diff| %.2e'%np.nanmax(np.abs(rr-V.ratio.reindex(rr.index))))
bins=[(-2.0,-1.75),(-1.75,-1.5),(-1.5,-1.25)]
def summary(D,sig,ns,lab):
    D=D.copy(); n=ns*2*D.r/D.RHI; n=n.fillna(ns*2*D.r/(3.0*D.Rd))   # fallback if no R_HI: h_g = 1.5 R_d
    D['vc2']=D.vobs**2+n*sig**2
    L=lastpts(D); lr=np.log10(np.sqrt(L.vmod2/L.vc2))
    q=pd.qcut(np.log10(L.Sh),4,labels=False)
    rho,p=spearmanr(lr,np.log10(L.Sh))
    D['lgoc']=lg(D.vc2,D.r); x=D.lgb+9.92082   # log g_bar/a0
    b=[(D.lgm-D.lgoc)[(x>=lo)&(x<hi)].mean() for lo,hi in bins]
    br=[(D.lrar-D.lgoc)[(x>=lo)&(x<hi)].mean() for lo,hi in bins]
    print('%-26s med %.3f | Spearman(Sh) %+.2f p=%.0e | Sh quartiles %s | low-g bins law-data %s | RAR-data %s'%(lab,10**lr.median(),rho,p,
          ' '.join('%.3f'%10**lr[q==k].median() for k in range(4)),' '.join('%+.3f'%v for v in b),' '.join('%+.3f'%v for v in br)))
    return L,lr
for case in ('univ','rel3'):
    D=A[A.case==case]
    print('\n=== case',case,'(ratio = V_model/V_c at R_last; quartiles low->high Sigma_h) ===')
    summary(D,0.0,1.0,'no drift')
    for sig in (6,8,12):
        for ns in (0.5,1.0,1.5):
            summary(D,sig,ns,'sigma %2d n x%.1f'%(sig,ns))
# regression on uncorrected ratio (Oman's ask): log ratio ~ (10/V_last)^2 * n_last + log Vmax + log Sigma_h
for case in ('univ','rel3'):
    D=A[A.case==case]; L=lastpts(D)
    vmax=D.groupby('gal').vobs.max().reindex(L.index)
    nl=(2*L.r/L.RHI).fillna(2*L.r/(3*L.Rd))
    y=np.log10(np.sqrt(L.vmod2)/L.vobs).values
    X=np.column_stack([np.ones(len(L)),(10/L.vobs)**2*nl,np.log10(vmax),np.log10(L.Sh)])
    names=['const','(10/V)^2 n','log Vmax','log Sh']
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    rng=np.random.default_rng(1); B=[]
    for _ in range(2000):
        i=rng.integers(0,len(y),len(y)); B.append(np.linalg.lstsq(X[i],y[i],rcond=None)[0])
    se=np.std(B,axis=0)
    print('\nOLS %s (N %d): '%(case,len(y))+'; '.join('%s %+.3f+-%.3f'%(n_,b_,s_) for n_,b_,s_ in zip(names,beta,se)))
    for k,(n_,xx) in enumerate(zip(names[1:],X[:,1:].T)):
        r_,p_=spearmanr(y,xx); print('   Spearman y vs %-12s %+.2f p=%.0e'%(n_,r_,p_))
    print('   corr between regressors: (10/V)^2n~logSh %+.2f, logVmax~logSh %+.2f'%(spearmanr(X[:,1],X[:,3])[0],spearmanr(X[:,2],X[:,3])[0]))
    print('   median predicted drift dV/V at R_last (sigma 8, n x1): all %.3f, low-Sh quartile %.3f, high %.3f'%(
        *[np.median(np.sqrt(1+64*nl[m]/L.vobs[m]**2)-1) for m in (slice(None),pd.qcut(np.log10(L.Sh),4,labels=False)==0,pd.qcut(np.log10(L.Sh),4,labels=False)==3)],))
