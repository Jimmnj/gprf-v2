# rhoc_sigma_b_oos_20261019.py [UNRATIFIED]  PRE-REGISTERED 2026-10-19 (ballot DF2-DETAIL / SLATE)
# Does rho_c track baryonic surface density, and does the SPARC trend PREDICT the low rho_c that Crater II and NGC1052-DF2/DF4 need?
# Step 1: SPARC frozen-law fits (TEST1 tanh2, N=120): log rho_c = a + b log M_b + c log Sigma_h, Sigma_h = M_b/(2 pi R_h^2), R_h = 1.678 R_d (proj. half-mass, disc proxy).
# Step 2: required rho_c for each dwarf from the N9 solver (Plummer stars, BTFR sigma, host ambient): scan log rho_c, interpolate to sigma_obs.
# Step 3: predicted log rho_c from step 1 at each dwarf's (M_b, Sigma_h); compare to required; 1-sigma band = SPARC residual + extrapolation error.
# DECLARED: PASS if predicted within required's error band (from sigma_obs +/- errors) for BOTH Crater II and DF2 (DF4 secondary, And XIX control must stay a pass).
# Dwarf inputs [mem, VERIFY]: Crater II M_b 2.6e5, r_h 1.07, MW vc 220 at 117 kpc, sigma 2.7+-0.3 (Caldwell+2017). And XIX M_b 1.6e6, r_h 3.1, M31 230 at 115, 7.8 +1.7/-1.5.
#   DF2 M_b 2.0e8, r_h 2.2, NGC1052 vc 250 at 80 kpc (projected), sigma 8.5 +2.3/-3.1 (Danieli+2019). DF4 M_b 1.5e8, r_h 1.6, same host, 4.2 +4.4/-2.2 (vD+2019).
import numpy as np, pandas as pd
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
from gate2c_robustness_20260925 import load_table
rng=np.random.default_rng(3)
T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd']]); F=F[F.Rd>0]
lM=np.log10(F.Mb.values); Rh=1.678*F.Rd.values*1e3; lS=np.log10(F.Mb.values/(2*np.pi*Rh**2)); y=F.l_tanh2.values
X=np.column_stack([np.ones(len(y)),lM,lS]); beta=np.linalg.lstsq(X,y,rcond=None)[0]; res=y-X@beta
bs=np.array([np.linalg.lstsq(X[i],y[i],rcond=None)[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(3000))])
lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
print('STEP 1  SPARC tanh2 rho_c, N=%d: log rho_c = %.2f + %+.2f log Mb + %+.2f log Sigma_h ; 95%%: b [%+.2f,%+.2f] c [%+.2f,%+.2f]; resid %.2f dex; median %.2f'%(len(y),beta[0],beta[1],beta[2],lo[1],hi[1],lo[2],hi[2],res.std(),np.median(y)))
A=np.polyfit(lS,y,1); print('        Sigma_h only: slope %+.2f resid %.2f | SPARC log Sigma_h range %.1f..%.1f (median %.1f)'%(A[0],(y-np.polyval(A,lS)).std(),lS.min(),lS.max(),np.median(lS)))
def required(name,Mb,rh,vc,D,so,ep,em,lrcs=np.arange(5.0,8.01,0.25)):
    a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb); amb=host_rho(vc,D)+rho_bg
    sig=[]
    for lrc in lrcs:
        rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb))); s=solve(Mb,a,s2,amb,lrc,rmax); sig.append(sigma_los(s,a,rh,rmax))
    sig=np.array(sig)
    f=lambda v: np.interp(v,sig,lrcs) if sig.min()<v<sig.max() else (np.nan)
    req=f(so); reqlo=f(so-em); reqhi=f(so+ep)
    lSh=np.log10(Mb/(2*np.pi*(rh*1e3)**2)); pred=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh
    pb=bs[:,0]+bs[:,1]*np.log10(Mb)+bs[:,2]*lSh; perr=np.sqrt(res.std()**2+pb.std()**2)
    print('%-10s Mb %.1e r_h %.2f log Sigma_h %5.2f | sigma(law) %4.1f | REQUIRED log rho_c %.2f [%.2f,%.2f] (sigma_obs %s +%s/-%s) | PREDICTED %.2f +/- %.2f (extrap %.2f) | at SPARC median %.2f -> sigma %.1f'%(
        name,Mb,rh,lSh,np.sqrt(s2),req,reqlo,reqhi,so,ep,em,pred,perr,pb.std(),np.median(y),np.interp(np.median(y),lrcs,sig)))
    print('           scan:',' '.join('%.2f:%.1f'%(l,s) for l,s in zip(lrcs,sig)))
print('STEP 2/3  (host vc, D as listed; ambient = host isothermal + rho_bg; truncated run)')
required('And XIX',1.6e6,3.1,230.,115.,7.8,1.7,1.5)
required('Crater II',2.6e5,1.07,220.,117.,2.7,0.3,0.3)
required('DF2',2.0e8,2.2,250.,80.,8.5,2.3,3.1)
required('DF4',1.5e8,1.6,250.,80.,4.2,4.4,2.2)
