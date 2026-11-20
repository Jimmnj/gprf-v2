# axisym_r1_20261109.py [UNRATIFIED] -- axisymmetric solve of the FROZEN law (canon sec 6: R1 2026-10-18; 2026-10-21 clarification).
# Rule 3a read done 2026-11-09 (SNAP_20260922_001 sec 0c; superseded in form by R1 -- this file implements R1, not the 9/22 law).
# LAW:  rho_nb(R,z) = rho_c exp[ S(x_b) - S(x_b0) - dphi_nb/sigma^2 ] - rho_amb
#       x_b = phi_b/sigma^2 with phi_b > 0 the baryon bump potential (-> 0 at infinity), x_b0 = phi_b(0)/sigma^2,
#       S(x) = 2 tanh(x/2), dphi_nb = drop of the nb mass's OWN potential from the centre (>= 0),
#       sigma^2 = sqrt(G a0 M_b)/2 (BTFR input), rho_amb = 35.5 Msun/kpc^3 (pass2 value). Identical to pass2 solve_gal in 1D.
# GEOMETRY MODES: C0   spherical phi_b from SPARC V_bar (pass2 construction) -> MUST reproduce pass2 solve_gal (validation)
#                 C1   flattened phi_b (ring sums: disc Sigma=Yd*SBdisk, h_z=0.196 Rd^0.633; gas exp. disc M=1.33 MHI, Rg from Vgas,
#                      h_g=0.2 kpc; spherical bulge from Vbul), nb MONOPOLE only
#                 full flattened phi_b, nb multipoles l<=12 (Green's function, damped fixed point)
# ASSUMPTION (flagged): in-plane model v^2 = SPARC V_bar^2 + nb part; the ring-sum baryon model only SOURCES the nb response.
# OUTPUT: axisym_r1_<mode>_gals.csv (per galaxy: fitted log rho_c with pass2 mask/weights/bounds; chi2 at fit and at 7.45;
#         isodensity axis ratios q at 2/3.2/5 Rd), axisym_r1_<mode>_points.csv (vmod2 at fit and at 7.45). Resumable (part files).
# Usage: python axisym_r1_20261109.py C0|C1|full [--gal A,B] [--workers N] [--budget 140]
import sys, os, io, time, argparse, contextlib, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.special import eval_legendre
from scipy.optimize import minimize_scalar
from scipy.integrate import cumulative_trapezoid as ctz
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)          # T, P, solve_gal, G, a0k, ramb, S, gals
assert abs(T.loc['NGC2841'].Rd-3.64)<1e-6 and abs(T.loc['DDO154'].Rd-0.37)<1e-6, 'MRT column error'
import gate2b_axisym_20261004 as gb
from gate2b_axisym_20261004 import rotmod, freeman_v2, ring_sources, ring_phi, multipoles, NMU, MU, WMU, MU_ALL, NR, NRC, YD, YB
LMAX_FULL=12

def prep(g, mode, R_out):
    """Return rr, psib (drop of phi_b from centre, >=0, shape [NR, len(MU_ALL)]), phib0 (phi_b at centre), diag."""
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); rl=r[-1]
    rr=np.concatenate([[0],np.geomspace(0.005,R_out,NR-1)])
    diag=dict(mratio=np.nan,vbratio=np.nan,Rg=np.nan)
    if mode=='C0':
        # pass2 construction: v_bar^2 interpolated to rl, point-mass tail beyond; phi_b(0) = int_0^rl v^2/x dx + v^2(rl)
        xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,np.r_[0,r],np.r_[0,vb2])/xs,xs)+vb2[-1]
        v2=np.where(rr<=rl,np.interp(rr,np.r_[0,r],np.r_[0,vb2]),vb2[-1]*rl/np.maximum(rr,1e-9))
        integ=np.where(rr>0,v2/np.maximum(rr,1e-12),0); integ[0]=integ[1]
        psi=ctz(integ,rr,initial=0)
        return rr,np.repeat(psi[:,None],len(MU_ALL),axis=1),phib0,diag
    raw=rotmod(g); t=T.loc[g]; Rd=t.Rd; L=t.L; MHI=t.MHI
    R0,sb,vg,vbul=raw[:,0],raw[:,6],raw[:,3],raw[:,5]
    hz=0.196*Rd**0.633
    Rs=np.geomspace(0.005*Rd,max(10*Rd,1.5*R0[-1]),160)
    sbq=np.interp(Rs,R0,np.clip(sb,0,None)); pos=np.where(sb>0)[0]
    if len(pos):
        jl=pos[-1]; b=Rs>R0[jl]; sbq[b]=sb[jl]*np.exp(-(Rs[b]-R0[jl])/Rd)
    sbq[Rs<R0[0]]=max(sb[0],0)
    Rd_,zd_,md_,Mdisk=ring_sources(Rs,YD*sbq*1e6,hz)
    Mg=1.33*MHI*1e9; ok=vg>0
    if ok.sum()>=3 and Mg>0:
        Rg=np.exp(minimize_scalar(lambda lr: np.sum((np.sqrt(freeman_v2(R0[ok],Mg,np.exp(lr)))-vg[ok])**2),
                                  bounds=(np.log(0.2*Rd),np.log(30*Rd)),method='bounded').x)
    else: Rg=2*Rd
    Rgs=np.geomspace(0.005*Rg,max(10*Rg,1.5*R0[-1]),160)
    Rg_,zg_,mg_,_=ring_sources(Rgs,Mg/(2*np.pi*Rg**2)*np.exp(-Rgs/Rg),0.2)
    rc=np.concatenate([[0],np.geomspace(0.01,R_out,NRC)])
    Rf=(rc[:,None]*np.sqrt(1-MU_ALL**2)).ravel(); zf=(rc[:,None]*MU_ALL).ravel()
    phi=(ring_phi(Rf,zf,Rd_,zd_,md_,hz/2)+ring_phi(Rf,zf,Rg_,zg_,mg_,0.1)).reshape(len(rc),len(MU_ALL))   # Newtonian sign, ->0 at inf
    phib0=-phi[0,0]
    psic=phi-phi[0,0]
    psib=np.column_stack([np.interp(rr,rc,psic[:,j]) for j in range(len(MU_ALL))])
    if np.any(vbul>0):
        vb2b=YB*vbul**2
        v2=np.where(rr<=R0[-1],np.interp(rr,np.r_[0,R0],np.r_[0,vb2b]),vb2b[-1]*R0[-1]/np.maximum(rr,1e-9))
        integ=np.where(rr>0,v2/np.maximum(rr,1e-12),0); integ[0]=integ[1]
        psib=psib+ctz(integ,rr,initial=0)[:,None]
        xs=np.linspace(1e-4,R0[-1],20001); phib0+=np.trapezoid(np.interp(xs,np.r_[0,R0],np.r_[0,vb2b])/xs,xs)+vb2b[-1]
    v2m=rc*np.gradient(phi[:,NMU],rc)
    spc=YD*raw[:,4]**2+np.sign(vg)*vg**2; m=np.interp(R0,rc,v2m); sel=spc>0
    diag=dict(mratio=Mdisk/(YD*L*1e9) if L>0 else np.nan, vbratio=float(np.median(m[sel]/spc[sel])) if sel.any() else np.nan, Rg=Rg)
    return rr,psib,phib0,diag

def solve_halo(rr,psib,phib0,s2,lrc,amb,lmax,Phil=None,iters=25,damp=0.5):
    """Frozen law. Sb = S(x_b)-S(x_b0) from the baryon drop; nb part linear: -(p0 + sum_l Phi_l P_l)/s2."""
    rc=10**lrc; nr=len(rr); ls=np.arange(2,lmax+1,2)
    xb0=phib0/s2; Sb=S(xb0-psib/s2)-S(xb0)                       # [nr, len(MU_ALL)], <= 0
    Pl=np.array([eval_legendre(l,MU_ALL) for l in ls]) if len(ls) else None
    if Phil is None or not len(ls): Phil=np.zeros((len(ls),nr))
    dPhil=np.zeros_like(Phil); it=0
    for it in range(iters if len(ls) else 1):
        nbl=(Phil.T@Pl if len(ls) else np.zeros_like(psib))        # l>=2 nb potential (drop; Phi_l(0)=0)
        ex=Sb[:,:NMU]-nbl[:,:NMU]/s2
        f=lambda i,p: 4*np.pi*rr[i]**2*(rc*np.dot(WMU,np.exp(np.clip(ex[i]-p/s2,-700,50)))-amb)
        M=np.zeros(nr); p0=np.zeros(nr)
        for i in range(nr-1):
            h=rr[i+1]-rr[i]; fi=f(i,p0[i]); gi=G*M[i]/rr[i]**2 if rr[i]>0 else 0.0
            Mp=M[i]+h*fi; pp=p0[i]+h*gi
            M[i+1]=M[i]+0.5*h*(fi+f(i+1,pp)); p0[i+1]=p0[i]+0.5*h*(gi+G*Mp/rr[i+1]**2)
        psi_nb=p0[:,None]+nbl
        if not len(ls): break
        rho=rc*np.exp(np.clip(Sb[:,:NMU]-psi_nb[:,:NMU]/s2,-700,50))-amb
        new,dnew=multipoles(rr,rho,ls)
        delta=np.max(np.abs(new-Phil)); dPhil=dnew
        if delta<1e-3*s2: Phil=new; break
        Phil=Phil+damp*(new-Phil)
    return M,p0,Phil,dPhil,Sb,psi_nb,it

def fit_gal(args):
    g,mode=args; t0=time.time()
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vo=d.vobs.values; vb2=np.clip(d.vbar2.values,0,None)
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; s2=np.sqrt(a0k*G*Mb)/2; Rd=t.Rd
    rr,psib,phib0,diag=prep(g,mode,2.0*r[-1])
    lmax=LMAX_FULL if mode=='full' else 0
    ls=np.arange(2,lmax+1,2); P0=eval_legendre(ls,0.0) if len(ls) else np.zeros(0)
    st={'Phil':None}
    def model(lr):
        M,p0,Phil,dPhil,Sb,psi,it=solve_halo(rr,psib,phib0,s2,lr,ramb,lmax,st['Phil'])
        st['Phil']=Phil; st['last']=(M,p0,Phil,dPhil,Sb,psi,it)
        vnb2=G*np.interp(r,rr,M)/r
        if len(ls): vnb2=vnb2+r*np.interp(r,rr,(dPhil*P0[:,None]).sum(0))
        return vb2+vnb2
    m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values; w=1/(0.1*vo)**2          # pass2 refit mask/weights
    chi=lambda lr: float(np.sum((w*(np.sqrt(np.clip(model(lr),1e-6,None))-vo)**2)[m])) if m.sum()>=3 else np.nan
    v2u=model(7.45); chiu=chi(7.45)
    lr=minimize_scalar(chi,bounds=(4.5,10.5),method='bounded',options={'xatol':0.02}).x if m.sum()>=3 else np.nan
    v2f=model(lr) if np.isfinite(lr) else v2u*np.nan; M,p0,Phil,dPhil,Sb,psi,it=st['last']
    out=dict(gal=g,mode=mode,Q=int(t.Q),npts=int(m.sum()),lrc=lr,chi2_fit=chi(lr) if np.isfinite(lr) else np.nan,chi2_745=chiu,
             Mb=Mb,s2=s2,Rd=Rd,phib0_s2=phib0/s2,iters=it,secs=time.time()-t0,**diag)
    # isodensity axis ratio of the nb response at fit (log rho via exponent: Sb - psi_nb/s2) and of the baryon isopotential
    for c in (2.0,3.2,5.0):
        R0=c*Rd
        for tag,ps in (('q',-(Sb-psi/s2)),('qb',psib)):
            if mode!='C0' and R0<rr[-1]:
                tgt=np.interp(R0,rr,ps[:,NMU]); ax=np.maximum.accumulate(ps[:,NMU+1])
                out[f'{tag}{c:g}']=np.interp(tgt,ax,rr)/R0 if tgt<=ax[-1] else np.nan
            else: out[f'{tag}{c:g}']=np.nan
    dd=d[['gal','r','vobs','ev','vbar2']].copy(); dd['mode']=mode; dd['used']=m; dd['vmod2_745']=v2u; dd['vmod2_fit']=v2f
    return out,dd

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('mode'); ap.add_argument('--gal',default=None)
    ap.add_argument('--workers',type=int,default=1); ap.add_argument('--budget',type=float,default=140); a=ap.parse_args()
    gl=a.gal.split(',') if a.gal else gals
    PG=f'axisym_r1_{a.mode}_gals_part.csv'; PP=f'axisym_r1_{a.mode}_points_part.csv'
    done=set(pd.read_csv(PG).gal) if os.path.exists(PG) else set()
    todo=[g for g in gl if g not in done]; t0=time.time()
    print(f'axisym R1 | mode {a.mode} | todo {len(todo)} of {len(gl)} | workers {a.workers}',flush=True)
    def flush(res):
        for o,dd in res:
            pd.DataFrame([o]).to_csv(PG,mode='a',header=not os.path.exists(PG),index=False)
            dd.to_csv(PP,mode='a',header=not os.path.exists(PP),index=False)
            print(f"  {o['gal']:10s} lrc {o['lrc']:.2f} chi2 {o['chi2_fit']:.1f} (7.45: {o['chi2_745']:.1f}) it {o['iters']} {o['secs']:.1f}s",flush=True)
    if a.workers>1:
        with Pool(a.workers) as pool:
            for res in pool.imap_unordered(fit_gal,[(g,a.mode) for g in todo],chunksize=1):
                flush([res])
                if time.time()-t0>a.budget: print('budget reached; rerun to resume'); pool.terminate(); break
    else:
        for g in todo:
            flush([fit_gal((g,a.mode))])
            if time.time()-t0>a.budget: print('budget reached; rerun to resume'); break
    if not [g for g in gl if g not in (set(pd.read_csv(PG).gal) if os.path.exists(PG) else set())]: print('DONE')
