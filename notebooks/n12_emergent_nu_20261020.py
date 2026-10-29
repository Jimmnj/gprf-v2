# n12_emergent_nu_20261020.py [UNRATIFIED] PRE-REGISTERED 2026-10-20 12:40 (JM 12:30 'can we derive the mu/nu shape').
# Emergent interpolating function from the frozen law with UNIVERSAL rho_c (no per-galaxy parameter): solve every SPARC galaxy (TEST1 solver, tanh2 cap,
# BTFR sigma) at log rho_c = 7.45 (median) and 7.2 / 7.8 (sensitivity), record (g_bar, g_model) at the data radii; also at the FITTED rho_c (reference).
# nu_eff = g_model/g_bar binned in log g_bar. DECLARED: limits must be exact (nu->1 high g; g -> sqrt(g_bar a0) low g). TEST: in the transition bins
# (|log g_bar/a0| < 0.5) report (a) scatter of log nu_eff across galaxies at fixed g_bar (universal rho_c) vs the DATA scatter about the RAR (McGaugh+2016
# ~0.11 dex incl. errors); (b) offset of binned universal-rho_c nu_eff from the RAR function nu_RAR = 1/(1-exp(-sqrt(g_bar/a0))); (c) shift of the curve
# between rho_c 7.2 and 7.8. PASS for 'nu derived' = (a) <= 0.15 dex AND |(b)| <= 0.1 dex in all transition bins at rho_c 7.45.
import numpy as np, pandas as pd, sys
sys.argv=['x']; src=open('test1_capped_baryon_20261018.py').read().split('def fit(')[0]; exec(src)
from scipy.integrate import solve_ivp
from gate2c_robustness_20260925 import kpc, a0
F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal')
def curve(g,l):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None)
    chi,s2b,phib0,Mb,Q=setup(g); S=LAWS['tanh2']; rc=10**l; xb0=phib0/s2b; Sb0=S(xb0); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*(rc*np.exp(S(xb0-y[0]/s2b)-Sb0-y[1]/s2b)-rho_amb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    if not s.success or s.y.shape[1]!=len(r): return None
    gb=vb2/r*1e6/kpc; gm=(vb2+G*s.y[2]/r)/r*1e6/kpc; go=d.vobs.values**2/r*1e6/kpc
    return pd.DataFrame(dict(gal=g,r=r,lgb=np.log10(gb),lgm=np.log10(gm),lgo=np.log10(go),Q=Q))
out=[]
for g in F.index:
    for lab,l in (('u7.45',7.45),('u7.2',7.2),('u7.8',7.8),('fit',F.loc[g,'l_tanh2'])):
        c=curve(g,l)
        if c is not None: c['case']=lab; out.append(c)
D=pd.concat(out); D.to_csv('n12_emergent_nu_20261020.csv',index=False)
la0=np.log10(a0)
D['lnu']=D.lgm-D.lgb; D['lnu_obs']=D.lgo-D.lgb; D['lnu_rar']=-np.log10(1-np.exp(-np.sqrt(10**D.lgb/a0)))
edges=np.arange(-3.0,1.01,0.25)+la0  # log g_bar bins relative to a0
print('bins in log(g_bar/a0); nu_eff = g_model/g_bar; RAR nu = 1/(1-exp(-sqrt(g_bar/a0))); scatter = std over points of log nu_eff in bin')
print('%8s | %22s | %22s | %22s | %12s | %22s'%('lgb/a0','UNIV 7.45 med/scat/dRAR','UNIV 7.2 med/dRAR','UNIV 7.8 med/dRAR','FIT med/scat','DATA med/scat/dRAR'))
rows=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    s='%+5.2f..%+5.2f'%(lo-la0,hi-la0); vals=[]
    for case in ('u7.45','u7.2','u7.8','fit'):
        m=(D.case==case)&(D.lgb>=lo)&(D.lgb<hi); x=D.lnu[m]; dr=(D.lnu-D.lnu_rar)[m]
        vals.append((len(x),x.median() if len(x) else np.nan,x.std() if len(x)>1 else np.nan,dr.median() if len(x) else np.nan))
    m=(D.case=='u7.45')&(D.lgb>=lo)&(D.lgb<hi); xo=D.lnu_obs[m]; dro=(D.lnu_obs-D.lnu_rar)[m]
    print('%s | N%4d %+.3f %.3f %+.3f | %+.3f %+.3f | %+.3f %+.3f | %+.3f %.3f | %+.3f %.3f %+.3f'%(s,vals[0][0],vals[0][1],vals[0][2],vals[0][3],vals[1][1],vals[1][3],vals[2][1],vals[2][3],vals[3][1],vals[3][2],xo.median(),xo.std(),dro.median()))
    rows.append(dict(lo=lo-la0,hi=hi-la0,N=vals[0][0],u745_med=vals[0][1],u745_scat=vals[0][2],u745_dRAR=vals[0][3],u72_dRAR=vals[1][3],u78_dRAR=vals[2][3],fit_scat=vals[3][2],data_scat=xo.std(),data_dRAR=dro.median()))
R=pd.DataFrame(rows); R.to_csv('n12_emergent_nu_bins_20261020.csv',index=False)
t=R[(R.lo>=-0.5)&(R.hi<=0.5)]
print('TRANSITION |log g_bar/a0|<0.5: universal-7.45 scatter max %.3f, |dRAR| max %.3f; shift 7.2->7.8 max %.3f; data scatter %.3f; fitted-rho_c scatter %.3f'%(t.u745_scat.max(),t.u745_dRAR.abs().max(),(t.u78_dRAR-t.u72_dRAR).abs().max(),t.data_scat.max(),t.fit_scat.max()))
