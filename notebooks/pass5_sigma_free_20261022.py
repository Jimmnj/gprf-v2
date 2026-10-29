# pass5_sigma_free_20261022.py [UNRATIFIED] 2026-10-22 22:00. Handoff 10/22 sec 4 item 2 (diagnostic suggested in a SIMULATED Desmond review; test + numbers ours).
# Per galaxy: fit rho_c AND sigma^2 = q * sigma^2_BTFR (q free) to SPARC V_obs with catalogue errors (floor 3% of V). Same law, sample (Lelli-style, 147 gal) and
# solver as pass 2. Then: regress log q on log Sigma_h (+ f_gas, log M_b); compare chi2 vs rho_c-only fit. Resumable: appends to pass5_sigma_free_fits.csv.
# FIX 22:10: explicit initial_simplex (scipy NM default step from 0 is 0.00025 -> 61/147 stuck at log q = 0 in the first run; discarded).
# Usage: python3 pass5_sigma_free_20261022.py fit <seconds>   |   python3 pass5_sigma_free_20261022.py report
import numpy as np, pandas as pd, glob, os, sys, time, warnings; warnings.filterwarnings('ignore')
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from scipy.stats import spearmanr
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2)
MRT='../sparc_galaxy_data/SPARC_Lelli2016c.mrt'; SD='../sparc_galaxy_data/sparc_database'
T=[]
for line in open(MRT):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: T.append(dict(gal=p[0],inc=float(p[5]),L=float(p[7]),Rd=float(p[11]),MHI=float(p[13]),Q=int(p[17])))
        except ValueError: pass
T=pd.DataFrame(T).set_index('gal')
rows=[]
for f in glob.glob(os.path.join(SD,'*_rotmod.dat')):
    g=os.path.basename(f).replace('_rotmod.dat','')
    if g not in T.index: continue
    d=np.atleast_2d(np.loadtxt(f,comments='#'))
    rows.append(pd.DataFrame(dict(gal=g,r=d[:,0],vobs=d[:,1],ev=d[:,2],vg=d[:,3],vd=d[:,4],vb=d[:,5])))
P=pd.concat(rows).join(T,on='gal')
P['vbar2']=np.sign(P.vg)*P.vg**2+0.5*P.vd**2+0.7*P.vb**2
P=P[(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)]
gals=sorted(P.gal.unique())
def model(g,lrc,lq):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; s2=10**lq*np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-5,atol=1e-3,max_step=rl/100)
    if not s.success or s.y.shape[1]!=len(r): return None,d,Mb,t
    return np.sqrt(vb2+G*s.y[2]/r),d,Mb,t
def chi2(g,lrc,lq):
    v,d,_,_=model(g,lrc,lq)
    if v is None: return 1e9
    e=np.maximum(d.ev.values,0.03*d.vobs.values)
    return float(np.sum(((v-d.vobs.values)/e)**2))
if sys.argv[1]=='time':
    t0=time.time(); print(chi2(gals[0],7.45,0.0), '%.3fs per solve'%(time.time()-t0)); sys.exit()
out='pass5_sigma_free_fits.csv'
if sys.argv[1]=='fit':
    budget=float(sys.argv[2]); t0=time.time()
    done=set(pd.read_csv(out).gal) if os.path.exists(out) else set()
    for g in gals:
        if g in done: continue
        if time.time()-t0>budget: break
        # rho_c-only fit (q=1), then 2-param from grid start
        grid=[(chi2(g,lr,0.0),lr) for lr in np.arange(6.0,9.01,0.25)]
        c1,lr1=min(grid); r1=minimize(lambda x:chi2(g,x[0],0.0),[lr1],method='Nelder-Mead',options=dict(xatol=0.01,fatol=0.05,maxiter=40,initial_simplex=[[lr1],[lr1+0.15]]))
        best=(1e18,None)
        for lq in (-0.4,-0.2,0.0,0.2,0.4):
            for lr in (r1.x[0]-0.5,r1.x[0],r1.x[0]+0.5):
                c=chi2(g,lr,lq); 
                if c<best[0]: best=(c,(lr,lq))
        b0=list(best[1]); r2=minimize(lambda x:chi2(g,x[0],np.clip(x[1],-1.0,1.0)),b0,method='Nelder-Mead',options=dict(xatol=0.01,fatol=0.05,maxiter=160,initial_simplex=[b0,[b0[0]+0.15,b0[1]],[b0[0],b0[1]+0.15]]))
        d=P[P.gal==g]; t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9
        rec=dict(gal=g,N=len(d),Mb=Mb,Sh=Mb/(2*np.pi*(1.678*t.Rd*1e3)**2) if t.Rd>0 else np.nan,fgas=1.33*t.MHI*1e9/Mb,
                 lrc1=r1.x[0],chi1=r1.fun,lrc2=r2.x[0],lq=float(np.clip(r2.x[1],-1,1)),chi2=r2.fun)
        pd.DataFrame([rec]).to_csv(out,mode='a',header=not os.path.exists(out),index=False)
    n=len(pd.read_csv(out)) if os.path.exists(out) else 0; print('fitted %d/%d'%(n,len(gals)))
if sys.argv[1]=='report':
    F=pd.read_csv(out); F=F[np.isfinite(F.Sh)]
    print('N',len(F)); print('log q = log(sigma2_fit/sigma2_BTFR): median %+.3f, 16-84 %+.3f..%+.3f; at bound |lq|=1: %d'%(F.lq.median(),*F.lq.quantile([.16,.84]),(F.lq.abs()>0.99).sum()))
    dchi=F.chi1-F.chi2; print('chi2 gain from freeing sigma: median %.1f (dof 1); galaxies with gain>4: %d'%(dchi.median(),(dchi>4).sum()))
    for lab,x in (('log Sigma_h',np.log10(F.Sh)),('f_gas',F.fgas),('log M_b',np.log10(F.Mb))):
        r,p=spearmanr(F.lq,x); print('Spearman log q vs %-12s %+.2f p=%.0e'%(lab,r,p))
    X=np.column_stack([np.ones(len(F)),np.log10(F.Sh),F.fgas,np.log10(F.Mb)]); y=F.lq.values
    b=np.linalg.lstsq(X,y,rcond=None)[0]; rng=np.random.default_rng(2)
    B=[np.linalg.lstsq(X[i],y[i],rcond=None)[0] for i in (rng.integers(0,len(y),len(y)) for _ in range(2000))]
    print('OLS log q = %+.3f %+.3f(+-%.3f) logSh %+.3f(+-%.3f) fgas %+.3f(+-%.3f) logMb'%(b[0],b[1],np.std(B,0)[1],b[2],np.std(B,0)[2],b[3],np.std(B,0)[3]))
    q=pd.qcut(np.log10(F.Sh),4,labels=False); print('median log q by Sigma_h quartile (low->high):',' '.join('%+.3f'%F.lq[q==k].median() for k in range(4)))
    print('implied sigma ratio sqrt(q) by quartile:',' '.join('%.3f'%10**(F.lq[q==k].median()/2) for k in range(4)))
    # Sigma-modulated sigma: slope of log q on log Sigma_h alone
    b1=np.polyfit(np.log10(F.Sh),F.lq,1); print('log q = %.3f log Sigma_h %+.3f  -> sigma^2 ∝ Sigma_h^%.3f at fixed M_b'%(b1[0],b1[1],b1[0]))
