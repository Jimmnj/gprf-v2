# test1b_burkert_surface_20261018.py [UNRATIFIED] Burkert-equivalent core surface density rho_0 r_0 of the TEST 1 model nb halos vs Donato+2009
# (MNRAS 397, 1169; log rho_0 r_0 = 2.15 +/- 0.2 Msun/pc^2, Burkert rho0 r0^3/((r+r0)(r^2+r0^2))). For each galaxy and law (lin, tanh2) at the
# TEST 1 best-fit rho_c: solve the nb profile, then fit Burkert (rho0, r0) two ways: (D) log density on 0..rmax; (V) halo v^2 = G M_h(<r)/r at the
# data radii (Donato's rotation-curve method). Run from notebooks/.
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from multiprocessing import Pool
import test1_capped_baryon_20261018 as T
from gate2c_robustness_20260925 import G
C=pd.read_csv('test1_capped_baryon_20261018.csv',index_col=0)
def MB(r,r0,rho0): x=r/r0; return np.pi*rho0*r0**3*(np.log(1+x*x)+2*np.log(1+x)-2*np.arctan(x))
def prof(g,law):
    d=T.P[T.P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    s2=C.loc[g,'s2b']; xb0=C.loc[g,'phib0']/s2; rc=10**C.loc[g,'l_'+law]; S=T.LAWS[law]; Sb0=S(xb0)
    rho=lambda x,y: rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4),G*y[2]/max(x,1e-4)**2,4*np.pi*x*x*(rho(x,y)-T.rho_amb)]
    xs=np.unique(np.concatenate([np.geomspace(max(r[0]/5,1e-3),rl,60),r]))
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=xs,rtol=1e-6,atol=1e-3,max_step=rl/200)
    dens=np.array([rho(x,s.y[:,i]) for i,x in enumerate(xs)])-T.rho_amb; M=s.y[2]
    return xs,dens,M,r
def one(g):
    out=dict(gal=g)
    for law in ('lin','tanh2'):
        xs,dens,M,r=prof(g,law); ok=dens>0
        rD=least_squares(lambda p: np.log10(dens[ok])-np.log10(10**p[1]*(10**p[0])**3/((xs[ok]+10**p[0])*(xs[ok]**2+10**(2*p[0])))),
                         [np.log10(xs[-1]/3),np.log10(dens[0])]).x
        Mr=np.interp(r,xs,M); v2=G*Mr/r; okv=v2>0
        rV=least_squares(lambda p: (np.sqrt(G*MB(r[okv],10**p[0],10**p[1])/r[okv])-np.sqrt(v2[okv]))/np.sqrt(v2[okv]).max(),
                         np.clip([rD[0],rD[1]],[-1.99,3.01],[2.99,10.99]),bounds=([-2,3],[3,11])).x
        for tag,p in (('D',rD),('V',rV)): out['lr0_%s_%s'%(tag,law)]=p[0]; out['mu_%s_%s'%(tag,law)]=p[0]+p[1]-6   # log(rho0 r0) Msun/pc^2: kpc->pc +3, kpc^-3->pc^-3 -9
        out['rmax']=r[-1]
    return out
if __name__=='__main__':
    with Pool(2) as Pl: R=pd.DataFrame(Pl.map(one,list(C.index))).set_index('gal')
    R=R.join(C[['Mb','Q']]); R.to_csv('test1b_burkert_surface_20261018.csv'); lM=np.log10(R.Mb)
    print('Donato+2009: log mu0D = 2.15 +/- 0.2 (Msun/pc^2), flat over 14 mag')
    for law in ('lin','tanh2'):
        for tag in ('D','V'):
            for sub,lab in ((R,'all'),(R[R.Q==1],'Q=1')):
                m=sub['mu_%s_%s'%(tag,law)]; A=np.polyfit(np.log10(sub.Mb),m,1)
                print('%-5s %s %-3s N=%3d  log mu0 median %.2f  16-84 %.2f..%.2f  std %.2f  slope vs log Mb %+.2f'%(law,tag,lab,len(sub),m.median(),*np.percentile(m,[16,84]),m.std(),A[0]))
    print('r0 (V,tanh2) / rmax median %.2f ; frac r0 > rmax %.2f'%((10**R.lr0_V_tanh2/R.rmax).median(),(10**R.lr0_V_tanh2>R.rmax).mean()))
