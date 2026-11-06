# h2_refit_20261103.py PREREG_H2_20261103 [results RATIFIED 2026-11-03 17:23]: add Leroy+2009 H2 (exponential disc) to the SPARC baryons for the 9 HERACLES x SPARC
# galaxies; refit at fixed M/L (law rho_c free; universal rho_c residuals; empirical function fixed/free Y). Resumable: rows -> h2_refit_part.csv.
import numpy as np, pandas as pd, os, sys, io, contextlib, time, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar, minimize
from scipy.special import i0,k0,i1,k1
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
H=pd.read_csv('../data/heracles2009/leroy2009_table3.csv',comment='#').set_index('gal')
GALS=['NGC2841','NGC2903','NGC2976','NGC3198','NGC3521','NGC4214','NGC5055','NGC6946','NGC7331']
def h2_model(g):
    h=H.loc[g]; Ds=T.loc[g].D; sc=Ds/h.D_leroy_Mpc
    M=5.5*h.LCO_1e7*1e7*sc**2
    l=h.lCO_kpc*sc if np.isfinite(h.lCO_kpc) else 0.2*h.r25_arcmin*np.pi/180/60*Ds*1e3
    return M,l,(not np.isfinite(h.lCO_kpc))
def vh2sq(r,M,l):  # Freeman 1970, razor-thin exponential disc; r,l kpc; returns (km/s)^2
    S0=M/(2*np.pi*l**2); y=np.clip(r/(2*l),1e-6,None)
    return 4*np.pi*G*S0*l*y**2*(i0(y)*k0(y)-i1(y)*k1(y))
def solve2(g,lrc,ext2,Mext,Ydisc=0.5):
    d=P[P.gal==g].sort_values('r'); r=d.r.values
    vb2=np.clip(np.sign(d.vg.values)*d.vg.values**2+Ydisc*d.vd.values**2+0.7*d.vb.values**2+ext2,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    t=T.loc[g]; Mb=Ydisc*t.L*1e9+1.33*t.MHI*1e9+Mext; s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    if not s.success or s.y.shape[1]!=len(r): return None
    return vb2+G*s.y[2]/r, vb2
def rarfn(gb): return gb/(1-np.exp(-np.sqrt(gb/a0)))
PART='h2_refit_part.csv'; done=set(pd.read_csv(PART).key) if os.path.exists(PART) else set(); t0=time.time()
for g in GALS:
    M,l,flag=h2_model(g); d=P[P.gal==g].sort_values('r'); r=d.r.values; vo=d.vobs.values; w=1/(0.1*vo)**2
    m=((d.ev/d.vobs<0.1)&(d.vobs>0)).values
    for case,ext2,Mext in (('noH2',0*r,0.0),('H2',vh2sq(r,M,l),M)):
        key=g+'_'+case
        if key in done: continue
        if time.time()-t0>140: print('TIME; rerun'); sys.exit()
        def chi_law(lrc,Y=0.5):
            o=solve2(g,lrc,ext2,Mext,Y); return 1e12 if o is None else float(np.sum((w*(np.sqrt(o[0])-vo)**2)[m]))
        res=minimize_scalar(lambda l_:chi_law(l_),bounds=(4.5,10.5),method='bounded',options={'xatol':0.02})
        o=solve2(g,7.45,ext2,Mext); vm2,vb2=o
        gb=vb2/r*1e6/kpc; go=vo**2/r*1e6/kpc; gm=vm2/r*1e6/kpc; inner=(gb>=a0)&m; outer=(gb<a0)&m
        res_m=np.log10(go/gm); res_r=np.log10(go/rarfn(gb))
        def chi_rar(Y):
            vb2y=np.clip(np.sign(d.vg.values)*d.vg.values**2+Y*d.vd.values**2+0.7*d.vb.values**2+ext2,0,None)
            vr=np.sqrt(rarfn(vb2y/r*1e6/kpc)*r*kpc/1e6); return float(np.sum((w*(vr-vo)**2)[m]))
        Yf=minimize_scalar(chi_rar,bounds=(0.05,2.0),method='bounded')
        # S1: law with rho_c and Y_disc free
        q=minimize(lambda p:chi_law(p[0],p[1]) if 0.05<p[1]<2 and 4.5<p[0]<10.5 else 1e13,[res.x,0.5],method='Nelder-Mead',options={'xatol':1e-2,'fatol':0.5,'maxiter':120})
        row=dict(key=key,gal=g,case=case,MH2=Mext,lH2=l,lflag=int(flag),Mb=0.5*T.loc[g].L*1e9+1.33*T.loc[g].MHI*1e9,
                 lrc=res.x,chi_law=res.fun,chi_law_u=chi_law(7.45),npts=int(m.sum()),nin=int(inner.sum()),
                 res_in_m=np.nanmedian(res_m[inner]) if inner.sum() else np.nan,res_out_m=np.nanmedian(res_m[outer]) if outer.sum() else np.nan,
                 res_in_r=np.nanmedian(res_r[inner]) if inner.sum() else np.nan,res_out_r=np.nanmedian(res_r[outer]) if outer.sum() else np.nan,
                 chi_rar_fix=chi_rar(0.5),Yfree=Yf.x,chi_rar_free=Yf.fun,lrc_Yfree=q.x[0],Y_law=q.x[1],chi_law_Yfree=q.fun)
        pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False); done.add(key)
        print('%-8s %-4s MH2 %.2e l %.2f lrc %.2f chi %.1f | u-res in %+.3f out %+.3f | rar fix %.1f free %.1f Y %.2f | S1 lrc %.2f Y %.2f chi %.1f  (%.0fs)'%(
            g,case,Mext,l,res.x,res.fun,row['res_in_m'],row['res_out_m'],row['chi_rar_fix'],Yf.fun,Yf.x,q.x[0],q.x[1],q.fun,time.time()-t0),flush=True)
print('DONE')
