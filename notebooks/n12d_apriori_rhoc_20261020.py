# n12d_apriori_rhoc_20261020.py [UNRATIFIED] PRE-REGISTERED 2026-10-20 16:25 (cold review A2/#3). rho_c set A PRIORI from Donato+2009 mu0 = 10^2.15 Msun/pc^2
# (independent of SPARC): mu0 = rho_c r_c with r_c = 3 sigma/sqrt(4 pi G rho_c) -> rho_c = (mu0/(3 sigma))^2 * 4 pi G, per galaxy from its BTFR sigma (case 'mu0'),
# and one global value at the sample-median sigma (case 'mu0g'). Compare emergent nu to the BINNED DATA (lgo - lgb), not the RAR function. Baseline tanh2.
# DECLARED: claim survives if |median model - median data| <= 0.05 dex in all bins with N >= 60 for case 'mu0'.
import numpy as np, pandas as pd, sys
sys.argv=['x']; src=open('test1_capped_baryon_20261018.py').read().split('def fit(')[0]; exec(src)
from scipy.integrate import solve_ivp
from gate2c_robustness_20260925 import kpc, a0
F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal'); S=LAWS['tanh2']
mu0=10**2.15*1e6   # Msun/kpc^2
lrc_mu0=np.log10((mu0/(3*np.sqrt(F.s2b)))**2*4*np.pi*G)   # per galaxy
sig_med=np.median(np.sqrt(F.s2b)); lrc_g=np.log10((mu0/(3*sig_med))**2*4*np.pi*G)
print('a-priori rho_c from mu0=10^2.15: per-galaxy log rho_c median %.2f (range %.2f..%.2f); global at sigma_med %.0f km/s: %.2f; SPARC fit median 7.45'%(lrc_mu0.median(),lrc_mu0.min(),lrc_mu0.max(),sig_med,lrc_g))
def curve(g,l):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None)
    chi,s2b,phib0,Mb,Q=setup(g); rc=10**l; xb0=phib0/s2b; Sb0=S(xb0); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*(rc*np.exp(S(xb0-y[0]/s2b)-Sb0-y[1]/s2b)-rho_amb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    if not s.success or s.y.shape[1]!=len(r): return None
    gb=vb2/r*1e6/kpc; gm=(vb2+G*s.y[2]/r)/r*1e6/kpc; go=d.vobs.values**2/r*1e6/kpc
    return pd.DataFrame(dict(gal=g,lgb=np.log10(gb),lnu=np.log10(gm/gb),lnu_obs=np.log10(go/gb)))
out=[]
for g in F.index:
    for case,l in (('mu0',lrc_mu0[g]),('mu0g',lrc_g)):
        c=curve(g,l)
        if c is not None: c['case']=case; out.append(c)
D=pd.concat(out); D.to_csv('n12d_apriori_rhoc_20261020.csv',index=False)
la0=np.log10(a0); edges=np.arange(-2.0,1.01,0.25)+la0
print('%12s %5s | %9s %9s | %9s %9s | %8s'%('lgb/a0','N','mu0-data','mu0g-data','mu0 scat','data scat','mu0 med'))
rows=[]
for lo,hi in zip(edges[:-1],edges[1:]):
    m=(D.case=='mu0')&(D.lgb>=lo)&(D.lgb<hi); mg=(D.case=='mu0g')&(D.lgb>=lo)&(D.lgb<hi)
    dm=D.lnu[m].median()-D.lnu_obs[m].median(); dg=D.lnu[mg].median()-D.lnu_obs[mg].median()
    print('%+5.2f..%+5.2f %5d | %+9.3f %+9.3f | %9.3f %9.3f | %+8.3f'%(lo-la0,hi-la0,m.sum(),dm,dg,D.lnu[m].std(),D.lnu_obs[m].std(),D.lnu[m].median()))
    rows.append(dict(lo=lo-la0,hi=hi-la0,N=int(m.sum()),mu0_minus_data=dm,mu0g_minus_data=dg))
R=pd.DataFrame(rows); R.to_csv('n12d_apriori_bins_20261020.csv',index=False)
big=R[R.N>=60]; print('bins N>=60: max |mu0 - data| %.3f ; max |mu0g - data| %.3f  (bar 0.05) -> %s'%(big.mu0_minus_data.abs().max(),big.mu0g_minus_data.abs().max(),'PASS' if big.mu0_minus_data.abs().max()<=0.05 else 'FAIL'))
d=D[D.case=='mu0']; print('per-point rms log(g_obs/g_model), a-priori mu0 rho_c: %.3f (universal 7.45 was 0.172; RAR fn 0.159)'%((d.lnu_obs-d.lnu).std()))
