# test1d_yd_robust_20261018.py [UNRATIFIED] Yd robustness of the TEST 1 'flat rho_c' finding (ballot GAL-THRESH / FR-AP gate).
# Variants DECLARED BEFORE RUN: Yd 0.3, 0.5 (reference; must reproduce TEST 1/1b/1c), 0.7, and 'mdep': Yd = 0.3 + 0.4*(log L - 8)/3 clipped 0.3..0.7
# (L = SPARC L[3.6] in Lsun). Yb 0.7. For each: M_b and BTFR sigma recomputed; capped law (tanh2) 1-param fit of rho_c; Burkert fit to the DATA;
# Burkert fit (V method) to the MODEL halo. Report slopes vs log M_b of model log rho_c, data log rho0, data log r0, data/model log mu0. Galaxy set and
# cuts = TEST 1 (load() with Yd 0.5 cuts; vbar2 recomputed per variant, clipped >= 0). Run from notebooks/: python3 test1d_yd_robust_20261018.py [variant]
import numpy as np, pandas as pd, sys, os
from scipy.optimize import minimize_scalar, least_squares
from scipy.integrate import solve_ivp
from multiprocessing import Pool
from gate2c_robustness_20260925 import load, G, a0k
from test1b_burkert_surface_20261018 import MB
rho_amb=35.5; S=lambda x:2*np.tanh(x/2)
N14=pd.read_csv('n14_sigma_free_20261001.csv',index_col=0); P0=load(); P0=P0[P0.gal.isin(N14.index)]
def ydof(var,L):
    if var=='mdep': return float(np.clip(0.3+0.4*(np.log10(L*1e9)-8)/3,0.3,0.7))
    return float(var)
def run(args):
    g,var=args; d=P0[P0.gal==g].sort_values('r'); L=float(d.L.iloc[0]); Yd=ydof(var,L)
    r=d.r.values; vo=d.vobs.values; vb2=np.clip(np.sign(d.vg)*d.vg**2+Yd*d.vd**2+0.7*d.vb**2,0,None).values; rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    Mb=Yd*L*1e9+1.33*float(d.MHI.iloc[0])*1e9; s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); w=1/(0.1*vo)**2
    def solve(l,teval):
        rc=10**l
        f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4),G*y[2]/max(x,1e-4)**2,4*np.pi*x*x*(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-rho_amb)]
        return solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=teval,rtol=1e-6,atol=1e-3,max_step=rl/200)
    def chi(l):
        s=solve(l,r)
        if not s.success or s.y.shape[1]!=len(r): return 1e12
        return float(np.sum(w*(np.sqrt(np.clip(vb2+G*s.y[2]/r,1,None))-vo)**2))
    res=minimize_scalar(chi,bounds=(4,11),method='bounded',options={'xatol':0.02}); l=res.x
    v2m=G*solve(l,r).y[2]/r; ok=v2m>0
    mV=least_squares(lambda p:(np.sqrt(G*MB(r[ok],10**p[0],10**p[1])/r[ok])-np.sqrt(v2m[ok]))/np.sqrt(v2m[ok]).max(),[np.log10(rl/3),l],bounds=([-2,3],[3,11])).x
    best=None
    for l0 in (-1,0,1):
        for p1 in (6,7.5,9):
            q=least_squares(lambda p:(np.sqrt(vb2+G*MB(r,10**p[0],10**p[1])/r)-vo)/(0.1*vo),[l0,p1],bounds=([-2,3],[3,11]))
            if best is None or q.cost<best.cost: best=q
    return dict(gal=g,var=var,Yd=Yd,Mb=Mb,Q=int(d.Q.iloc[0]),xb0=xb0,l_rhoc=l,chi=res.fun,lr0_m=mV[0],lrho0_m=mV[1],lr0_d=best.x[0],lrho0_d=best.x[1],
                chi_burk=2*best.cost,rmax=rl)
if __name__=='__main__':
    var=sys.argv[1]; out='test1d_yd_%s_20261018.csv'%var
    with Pool(2) as Pl: R=pd.DataFrame(Pl.map(run,[(g,var) for g in N14.index if g in set(P0.gal)]))
    R.to_csv(out,index=False); print(var,'done',len(R))
