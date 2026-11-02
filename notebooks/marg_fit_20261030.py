# marg_fit_20261030.py -- marginalised fit under PREREG_MARG_20261030.md (STAMPED JM 2026-10-30 12:20; FROZEN). [UNRATIFIED result]
# Arm A: the relation (eq. 1-2), log rho_c free flat [6,10]; Arm B: empirical RAR function (g_dagger = a0 fixed), no rho_c.
# Both arms: per galaxy log Ud ~ N(log 0.5, 0.1 dex), log Ub ~ N(log 0.7, 0.1 dex) (bulge galaxies), D ~ N(D_SPARC, e_D), i ~ N(i_SPARC, e_i) clipped [20,90].
# MAP: chi^2 with SPARC velocity errors + Gaussian prior penalties; Nelder-Mead, two starts, better kept.
# Harness choices carried from pass2_sparc_full_20261022.py: vbar^2 = vg|vg| + Ud vd^2 + Ub vb^2 (rotmod at Upsilon=1); M_b = (Ud L + 1.33 MHI) D'^2
# (the MRT has no separate bulge luminosity; the baseline harness uses 0.5 L the same way); sigma^2 = sqrt(G a0 M_b)/2 at the FITTED M_b; rho_amb 35.5.
# Distance: r -> r D', baryon v^2 -> v^2 D' (M ~ D^2, r ~ D); v_obs unchanged by D, rescaled by sin(i0)/sin(i); errors likewise.
# SUPPLEMENTARY (not in the prereg, added 2026-10-30 13:05 after the A/B run): arm U = relation at UNIVERSAL rho_c 10^7.45 with the same nuisances; arm R = rho_c from relation (3) at the fitted M_b, Sigma_h, f_gas.
# Usage: python3 marg_fit_20261030.py --arm both --procs 16 [--slice 0:50] [--pilot] [--out marg_fit]
import os, sys, glob, argparse, time, numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, minimize_scalar
from multiprocessing import Pool
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2)
LUD0,LUB0,SIGU=np.log10(0.5),np.log10(0.7),0.1
MRT='../sparc_galaxy_data/SPARC_Lelli2016c.mrt'; SD='../sparc_galaxy_data/sparc_database'
def load():
    T=[]
    for line in open(MRT):
        p=line.split()
        if len(p)>=18 and p[1].isdigit():
            try: T.append(dict(gal=p[0],D=float(p[2]),eD=float(p[3]),inc=float(p[5]),einc=float(p[6]),L=float(p[7]),Rd=float(p[11]),MHI=float(p[13]),Vf=float(p[15]),Q=int(p[17])))
            except ValueError: pass
    T=pd.DataFrame(T).set_index('gal'); rows=[]
    for f in glob.glob(os.path.join(SD,'*_rotmod.dat')):
        g=os.path.basename(f).replace('_rotmod.dat','')
        if g not in T.index: continue
        d=np.atleast_2d(np.loadtxt(f,comments='#'))
        rows.append(pd.DataFrame(dict(gal=g,r=d[:,0],vobs=d[:,1],ev=d[:,2],vg=d[:,3],vd=d[:,4],vb=d[:,5])))
    P=pd.concat(rows).join(T,on='gal')
    P['vbar2_0']=np.sign(P.vg)*P.vg**2+0.5*P.vd**2+0.7*P.vb**2
    sel=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)&(P.vbar2_0>0)
    gals=sorted(P[sel].gal.unique())
    return T,P,gals
T=P=None
def init(t,p):
    global T,P; T,P=t,p
class Gal:
    def __init__(s,g):
        d=P[(P.gal==g)].sort_values('r'); s.g=g; s.t=T.loc[g]
        s.r0=d.r.values; s.vobs0=d.vobs.values; s.ev0=d.ev.values; s.vg=d.vg.values; s.vd=d.vd.values; s.vb=d.vb.values
        s.bul=bool(np.any(s.vb>0)); s.D0=s.t.D; s.eD=max(s.t.eD,1e-3*s.D0); s.i0=s.t.inc; s.ei=max(s.t.einc,0.5)
        s.sel=(d.ev.values/d.vobs.values<0.1)&(d.vobs.values>0)   # Lelli point cut (at nominal values, fixed)
    def unpack(s,x,arm):
        k=0
        if arm=='A': lrc=x[k]; k+=1
        elif arm=='U': lrc=7.45
        elif arm=='R': lrc=None
        else: lrc=None
        lud=x[k]; k+=1; lub=x[k] if s.bul else LUB0; k+= (1 if s.bul else 0)
        zD=x[k]; zi=x[k+1]
        D=max(s.D0+zD*s.eD,0.1*s.D0); inc=float(np.clip(s.i0+zi*s.ei,20,90)); return lrc,lud,lub,D,inc
    def curves(s,lud,lub,D,inc):
        Dp=D/s.D0; r=s.r0*Dp
        vb2=(np.sign(s.vg)*s.vg**2+10**lud*s.vd**2+10**lub*s.vb**2)*Dp
        f=np.sin(np.radians(s.i0))/np.sin(np.radians(inc)); vobs=s.vobs0*f; ev=s.ev0*f
        Mb=(10**lud*s.t.L*1e9+1.33*s.t.MHI*1e9)*Dp**2
        return r,vb2,vobs,ev,Mb
    def model(s,x,arm):
        lrc,lud,lub,D,inc=s.unpack(x,arm); r,vb2,vobs,ev,Mb=s.curves(lud,lub,D,inc); vb2c=np.clip(vb2,0,None)
        if arm=='R':
            Dp=D/s.D0; Sh=Mb/(2*np.pi*(1.678*s.t.Rd*Dp*1e3)**2) if s.t.Rd>0 else np.nan; fg=1.33*s.t.MHI*1e9*Dp**2/Mb
            lrc=11.063-0.365*np.log10(Mb)+0.268*np.log10(Sh)-0.929*fg if np.isfinite(Sh) else 7.45
            if not (6<=lrc<=10): return None
        if arm=='B':
            gb=vb2c/r; gm=gb/(1-np.exp(-np.sqrt(np.clip(gb/a0k,1e-12,None)))); return r,vb2,vobs,ev,Mb,gm*r
        rl=r[-1]; rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c])
        vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2c[-1]*rl/x
        xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]
        s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
        fn=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
        so=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
        if not so.success or so.y.shape[1]!=len(r): return None
        return r,vb2,vobs,ev,Mb,vb2c+G*so.y[2]/r
    def obj(s,x,arm):
        m=s.model(x,arm)
        if m is None: return 1e9
        r,vb2,vobs,ev,Mb,vm2=m; vm=np.sqrt(np.clip(vm2,1,None))
        chi=np.sum(((vm-vobs)/ev)[s.sel]**2)
        lrc,lud,lub,D,inc=s.unpack(x,arm)
        pen=((lud-LUD0)/SIGU)**2+(((lub-LUB0)/SIGU)**2 if s.bul else 0)+((D-s.D0)/s.eD)**2+((inc-s.i0)/s.ei)**2
        if arm=='A': pen+= 0 if 6<=lrc<=10 else 1e3*(min(abs(lrc-6),abs(lrc-10)))**2
        return chi+pen
def fit_one(args):
    g,arm=args; t0=time.time(); s=Gal(g); nuis0=[LUD0]+([LUB0] if s.bul else [])+[0.0,0.0]
    if arm=='A':
        f1=lambda l: s.obj([l]+nuis0,arm); l0=minimize_scalar(f1,bounds=(6,10),method='bounded',options={'xatol':0.05}).x
        starts=[[l0]+nuis0,[l0+0.3,LUD0-0.1]+([LUB0] if s.bul else [])+[0.5,-0.5]]
    else:
        starts=[nuis0,[LUD0-0.1]+([LUB0] if s.bul else [])+[0.5,-0.5]]
    best=None
    for x0 in starts:
        r=minimize(lambda x:s.obj(x,arm),x0,method='Nelder-Mead',options={'xatol':0.01,'fatol':0.05,'maxiter':600,'maxfev':900})
        if best is None or r.fun<best.fun: best=r
    lrc,lud,lub,D,inc=s.unpack(best.x,arm); m=s.model(best.x,arm)
    if arm=='R':
        Dp=D/s.D0; Mb_=(10**lud*s.t.L*1e9+1.33*s.t.MHI*1e9)*Dp**2; Sh=Mb_/(2*np.pi*(1.678*s.t.Rd*Dp*1e3)**2) if s.t.Rd>0 else np.nan; lrc=11.063-0.365*np.log10(Mb_)+0.268*np.log10(Sh)-0.929*(1.33*s.t.MHI*1e9*Dp**2/Mb_) if np.isfinite(Sh) else 7.45
    if m is None: return None,None
    r,vb2,vobs,ev,Mb,vm2=m
    pts=pd.DataFrame(dict(gal=g,arm=arm,r=r,vobs=vobs,ev=ev,vbar2=vb2,vmod2=vm2,sel=s.sel,Mb=Mb))
    chi=np.sum(((np.sqrt(np.clip(vm2,1,None))-vobs)/ev)[s.sel]**2)
    par=dict(gal=g,arm=arm,lrc=lrc,lud=lud,lub=(lub if s.bul else np.nan),bul=s.bul,D=D,D0=s.D0,eD=s.eD,zD=(D-s.D0)/s.eD,inc=inc,i0=s.i0,ei=s.ei,zi=(inc-s.i0)/s.ei,
             Mb=Mb,L=s.t.L,MHI=s.t.MHI,Rd=s.t.Rd,chi2=chi,npts=int(s.sel.sum()),pen_obj=best.fun,nfev=best.nfev,sec=time.time()-t0)
    return par,pts
if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--arm',default='both'); ap.add_argument('--procs',type=int,default=2); ap.add_argument('--slice',default=None)
    ap.add_argument('--pilot',action='store_true'); ap.add_argument('--out',default='marg_fit'); a=ap.parse_args()
    T,P,gals=load()
    if a.pilot: gals=[g for g in gals if T.loc[g,'Q']==1][:10]
    if a.slice: lo,hi=map(int,a.slice.split(':')); gals=gals[lo:hi]
    arms=['A','B'] if a.arm=='both' else list(a.arm); jobs=[(g,arm) for arm in arms for g in gals]
    print('galaxies %d, jobs %d, procs %d'%(len(gals),len(jobs),a.procs),flush=True); t0=time.time()
    with Pool(a.procs,initializer=init,initargs=(T,P)) as pool:
        out=[]
        for k,(par,pts) in enumerate(pool.imap_unordered(fit_one,jobs),1):
            if par is not None: out.append((par,pts))
            if k%10==0 or k==len(jobs): print('  %d/%d done, %.0f s'%(k,len(jobs),time.time()-t0),flush=True)
    tag=('_pilot' if a.pilot else '')+(('_%s'%a.slice.replace(':','-')) if a.slice else '')
    pd.DataFrame([p for p,_ in out]).to_csv('%s_params_20261030%s.csv'%(a.out,tag),index=False)
    pd.concat([q for _,q in out]).to_csv('%s_points_20261030%s.csv'%(a.out,tag),index=False)
    print('wrote %s_params/points_20261030%s.csv; %.0f s total'%(a.out,tag,time.time()-t0))
