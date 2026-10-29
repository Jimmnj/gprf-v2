# pass3_relation3_20261022.py [UNRATIFIED] 2026-10-22 17:45 -- Lelli-2 'single most important item': the law with rho_c from relation (3)
# (log rho_c = 11.063 - 0.365 log M_b + 0.268 log Sigma_h - 0.929 f_gas; four global constants, zero per-galaxy freedom) on all 147 Lelli-style galaxies.
# Compare with universal 10^7.45 and with the per-galaxy fit: binned RAR offsets, per-point rms, V(R_last) ratio vs Sigma_h, residual Spearman vs R/Sigma_h/f_gas/M_b.
# Caveat: relation (3) was fitted to the TEST1 per-galaxy rho_c of the older 120-galaxy sample, so this is in-sample for those and out-of-sample for the 27 non-flat galaxies.
import numpy as np, pandas as pd, io, contextlib, warnings; warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
beta=np.array([11.063,-0.365,0.268,-0.929])
Lsel=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)&(P.vbar2>0)
gals=sorted(P[Lsel].gal.unique())
def lrc_rel(g):
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; Sh=Mb/(2*np.pi*(1.678*t.Rd*1e3)**2) if t.Rd>0 else np.nan; fg=1.33*t.MHI*1e9/Mb
    return beta[0]+beta[1]*np.log10(Mb)+beta[2]*np.log10(Sh)+beta[3]*fg if np.isfinite(Sh) else 7.45
out=[]
for g in gals:
    for lab,l in (('univ',7.45),('rel3',lrc_rel(g))):
        o=solve_gal(g,l)
        if o is None: continue
        o=o[(o.ev/o.vobs<0.1)&(o.vobs>0)&(o.vbar2>0)].copy(); o['case']=lab; o['lrc']=l; out.append(o)
A=pd.concat(out); A['lgb']=np.log10(A.vbar2/A.r*1e6/kpc); A['lgo']=np.log10(A.vobs**2/A.r*1e6/kpc); A['lgm']=np.log10(A.vmod2/A.r*1e6/kpc)
A['lrar']=np.log10(10**A.lgb/(1-np.exp(-np.sqrt(10**A.lgb/a0)))); A['res']=A.lgo-A.lgm; A['resR']=A.lgo-A.lrar
A.to_csv('pass3_relation3_points_20261022.csv',index=False)
la0=np.log10(a0); edges=np.arange(-2.25,1.01,0.25)+la0
print('binned median log(g_obs/g_model), Lelli-style points (N %d), 147 galaxies'%(len(A)//2))
print('  lg(gb/a0)     N   universal   relation(3)   RAR fn')
for lo,hi in zip(edges[:-1],edges[1:]):
    u=A[(A.case=='univ')&(A.lgb>=lo)&(A.lgb<hi)]; r3=A[(A.case=='rel3')&(A.lgb>=lo)&(A.lgb<hi)]
    if len(u)<10: continue
    print('  %+5.2f..%+5.2f %4d   %+.3f      %+.3f        %+.3f'%(lo-la0,hi-la0,len(u),u.res.median(),r3.res.median(),u.resR.median()))
for lab in ('univ','rel3'):
    d=A[A.case==lab]; print('%-10s per-point rms %.3f (RAR fn %.3f); Spearman res vs log r %+.2f, log Sigma_h %+.2f, f_gas %+.2f, log M_b %+.2f'%(lab,np.sqrt((d.res**2).mean()),np.sqrt((d.resR**2).mean()),
          spearmanr(np.log10(d.r),d.res)[0],spearmanr(np.log10(d.Sh),d.res)[0],spearmanr(d.fgas,d.res)[0],spearmanr(np.log10(d.Mb),d.res)[0]))
    G2=[]
    for g,dd in d.groupby('gal'):
        last=dd.sort_values('r').iloc[-1]; G2.append(dict(gal=g,Sh=last.Sh,Mb=last.Mb,ratio=np.sqrt(last.vmod2)/last.vobs))
    G2=pd.DataFrame(G2); q=pd.qcut(np.log10(G2.Sh),4,labels=['lowSigma','q2','q3','highSigma'])
    print('           V_model/V_obs at R_last: median %.3f; Spearman vs log Sigma_h %+.2f (p %.1e); by quartile %s'%(G2.ratio.median(),*spearmanr(np.log10(G2.Sh),np.log10(G2.ratio))[:2],G2.groupby(q).ratio.median().round(3).to_dict()))
d=A[A.case=='rel3']; print('relation(3) rho_c: median %.2f, range %.2f-%.2f'%(d.groupby('gal').lrc.first().median(),d.lrc.min(),d.lrc.max()))
