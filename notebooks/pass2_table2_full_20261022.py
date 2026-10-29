# pass2_table2_full_20261022.py [UNRATIFIED] -- Table 2 on the full Lelli-style sample (147 galaxies with points, 2696 points). Resumable: per-galaxy rows
# appended to pass2_table2_part.csv; run repeatedly until 'DONE'. Per galaxy: law rho_c free (1D, 10% velocity weights as before); RAR function with
# Y_disk fixed 0.5, Y_disk free, and Y_disk+D+i free under SPARC priors (as n12e). Residuals log(g_obs/g_model) on Lelli-style points.
import numpy as np, pandas as pd, os, sys, json, time, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar, minimize
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()): exec(src)
PART='pass2_table2_part.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set()
Lsel=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)&(P.vbar2>0)
todo=[g for g in sorted(P[Lsel].gal.unique()) if g not in done]; t0=time.time()
def rarfn(gb): return gb/(1-np.exp(-np.sqrt(gb/a0)))
for g in todo:
    if time.time()-t0>150: break
    d=P[P.gal==g].sort_values('r'); m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values
    r=d.r.values; vo=d.vobs.values; w=1/(0.1*vo)**2; t=T.loc[g]
    def law_chi(l):
        o=solve_gal(g,l)
        if o is None: return 1e12,None
        return float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m])),o
    res=minimize_scalar(lambda l:law_chi(l)[0],bounds=(4.5,10.5),method='bounded',options={'xatol':0.02})
    c,o=law_chi(res.x); lgo=np.log10(vo**2/r*1e6/kpc)
    law_res=list((lgo-np.log10(o.vmod2.values/r*1e6/kpc))[m]) if o is not None else []
    def vb2(Y,D=1.0): return (np.sign(d.vg)*d.vg**2+Y*d.vd**2+0.7*d.vb**2).values*D   # V^2 ∝ D at fixed angle (M ∝ D^2, r ∝ D)
    def rar_res(Y,D=1.0,di=0.0):
        rr=r*D; gb=np.clip(vb2(Y,D),1e-9,None)/rr*1e6/kpc; inc=np.radians(t.inc); vobs=vo*np.sin(inc)/np.sin(np.clip(inc+np.radians(di),0.35,1.5708))
        go=vobs**2/rr*1e6/kpc; return np.log10(go)-np.log10(rarfn(gb)), vobs
    def rar_chi(Y,D=1.0,di=0.0):
        rs,vobs=rar_res(Y,D,di); gb=np.clip(vb2(Y,D),1e-9,None)/(r*D)*1e6/kpc; vm=np.sqrt(rarfn(gb)*r*D/1e6*kpc)
        return float(np.sum((((vm-vobs)/(0.1*vobs))**2)[m]))
    fx=list(rar_res(0.5)[0][m])
    Yf=minimize_scalar(lambda Y:rar_chi(Y),bounds=(0.05,2.0),method='bounded').x; fr=list(rar_res(Yf)[0][m])
    pri=lambda p: rar_chi(p[0],p[1],p[2])+(np.log10(max(p[0],1e-3)/0.5)/0.1)**2+((p[1]-1)/(t.eD/t.D))**2+(p[2]/max(t.einc,1.0))**2
    q=minimize(pri,[0.5,1.0,0.0],method='Nelder-Mead',options={'xatol':1e-3,'fatol':1e-2,'maxiter':600}).x; f3=list(rar_res(*q)[0][m])
    row=dict(gal=g,lrc=res.x,chi_law=c,Yfree=Yf,Y3=q[0],D3=q[1],di3=q[2],noVf=int(t.Vf<=0),
             law=json.dumps(law_res),rfix=json.dumps(fx),rfree=json.dumps(fr),r3=json.dumps(f3))
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
n=len(pd.read_csv(PART)); print('progress %d / %d'%(n,P[Lsel].gal.nunique()))
if n>=P[Lsel].gal.nunique():
    R=pd.read_csv(PART).drop_duplicates('gal')
    def rms(col,sub=R): a=np.concatenate([json.loads(x) for x in sub[col]]); return np.sqrt(np.mean(a**2)), len(a)
    print('DONE. per-point rms log(g_obs/g_model), Lelli-style sample, %d galaxies:'%len(R))
    for lab,col in (('RAR function, M/L fixed','rfix'),('RAR function, M/L_disc free','rfree'),('law, rho_c free','law'),('RAR function, M/L, D, i free (SPARC priors)','r3')):
        a,n2=rms(col); b,_=rms(col,R[R.noVf==0]); c2,_=rms(col,R[R.noVf==1]); print('  %-44s %.3f  (N pts %d; with V_f %.3f; without V_f %.3f)'%(lab,a,n2,b,c2))
    print('  fitted log rho_c: median %.2f, spread (std) %.2f; galaxies without V_f median %.2f'%(R.lrc.median(),R.lrc.std(),R[R.noVf==1].lrc.median()))
