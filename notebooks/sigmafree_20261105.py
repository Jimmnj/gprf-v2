# sigmafree_20261105.py [UNRATIFIED] -- PREREG_SIGMAFREE_20261105. Resumable: rows -> sigmafree_part.csv; rerun until DONE.
import numpy as np, pandas as pd, os, sys, io, contextlib, time, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize, minimize_scalar
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
W=0.06
def solve_s(g,lrc,s2):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    if not s.success or s.y.shape[1]!=len(r): return None
    return vb2+G*s.y[2]/r
PART='sigmafree_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set(); t0=time.time()
Lsel=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)&(P.vbar2>0)
for g in sorted(P[Lsel].gal.unique()):
    if g in done: continue
    if time.time()-t0>140: print('TIME; rerun'); sys.exit()
    d=P[P.gal==g].sort_values('r'); m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values; vo=d.vobs.values; w=1/(0.1*vo)**2
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; s2b=np.sqrt(a0k*G*Mb)/2; ls0=0.5*np.log10(s2b)
    def chi(lrc,ls):
        if not (4.5<lrc<10.5) or abs(ls-ls0)>0.6: return 1e13
        o=solve_s(g,lrc,10**(2*ls)); return 1e12 if o is None else float(np.sum((w*(np.sqrt(o)-vo)**2)[m]))
    r0=minimize_scalar(lambda l:chi(l,ls0),bounds=(4.5,10.5),method='bounded',options={'xatol':0.02})
    q=minimize(lambda p:chi(p[0],p[1])+((p[1]-ls0)/W)**2,[r0.x,ls0],method='Nelder-Mead',options={'xatol':1e-2,'fatol':0.05,'maxiter':150})
    qf=minimize(lambda p:chi(p[0],p[1]),[q.x[0],q.x[1]],method='Nelder-Mead',options={'xatol':1e-2,'fatol':0.05,'maxiter':150})
    row=dict(gal=g,Mb=Mb,Vf=t.Vf,lrc0=r0.x,chi0=r0.fun,lrc_p=q.x[0],dls_p=q.x[1]-ls0,chi_p=chi(q.x[0],q.x[1]),lrc_f=qf.x[0],dls_f=qf.x[1]-ls0,chi_f=qf.fun,npts=int(m.sum()))
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False); done.add(g)
    print('%-10s lrc0 %.2f chi %.1f | prior: lrc %.2f dls %+.3f chi %.1f | free: lrc %.2f dls %+.3f chi %.1f (%.0fs)'%(g,r0.x,r0.fun,q.x[0],q.x[1]-ls0,row['chi_p'],qf.x[0],qf.x[1]-ls0,qf.fun,time.time()-t0),flush=True)
print('DONE')
