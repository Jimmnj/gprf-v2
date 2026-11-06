# twocentre_b_20261102.py [UNRATIFIED] PREREG_TWOCENTRE_B_20261102. Satellite's own relation with the host's potential drop in phi_b (B1 capped / B2 uncapped nb).
import numpy as np, sys, json, time
from scipy.special import eval_legendre
from scipy.integrate import cumulative_trapezoid as ctz, solve_ivp
from scipy.optimize import newton_krylov
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6; RHO_BG=35.5; S=lambda x:2*np.tanh(x/2)
args=dict(a.split('=') for a in sys.argv[1:])
L=int(args.get('L',40)); NMU=int(args.get('NMU',96)); NR=int(args.get('NR',1200)); LRC=float(args.get('lrc',7.45))
MSAT=float(args.get('Msat',3.0e9)); ASAT=float(args.get('asat',1.5)); D=float(args.get('D',50.0)); MODE=args.get('mode','B1'); TAG=args.get('tag','LMC_B1')
DAMP=float(args.get('damp',0.3)); ITERS=int(args.get('iters',80)); NOHOST=int(args.get('nohost',0)); NOSELF=int(args.get('noself',0)); HOSTF=float(args.get('hostf',1.0))
MW=[(1.55e10,0.70),(3.65e10,0.695*2.35),(8.2e9,0.695*18.24),(1.3e9,0.695*2.57)]; MB_MW=sum(m for m,_ in MW); S2MW=0.5*np.sqrt(G*a0k*MB_MW)
# --- MW unperturbed nb profile (1D shoot, universal rho_c), absolute potentials
def phi_mwb(r): return sum(-G*m/(r+a) for m,a in MW)
rcmw=10**7.45; xb0mw=-phi_mwb(1e-6)/S2MW
def fmw(r,y):
    rho=max(rcmw*np.exp(np.clip(S(-phi_mwb(r)/S2MW)-S(xb0mw)-y[0]/S2MW,-700,50))-RHO_BG,-RHO_BG); return [G*y[1]/r**2,4*np.pi*r*r*rho]
mwsol=solve_ivp(fmw,[1e-4,3000.],[0.,0.],rtol=1e-9,atol=1e-3,dense_output=True)
rg=np.logspace(-4,np.log10(3000.),3000); Mnb=np.array([mwsol.sol(r)[1] for r in rg]); rhonb=np.array([rcmw*np.exp(S(-phi_mwb(r)/S2MW)-S(xb0mw)-mwsol.sol(r)[0]/S2MW)-RHO_BG for r in rg])
gnb=G*Mnb/rg**2; phinb_abs=-(np.trapz(gnb,rg)-ctz(gnb,rg,initial=0))   # Phi_nb(r) = -int_r^Rmax g dr (zero at Rmax)
def phi_mwnb(d): return np.interp(d,rg,phinb_abs)
def rho_mwnb(d): return np.interp(d,rg,rhonb)
# --- satellite
S2=0.5*np.sqrt(G*a0k*MSAT); SIG=np.sqrt(S2); rc=10**LRC
rr=np.logspace(np.log10(0.05),np.log10(0.9*D),NR); mu,wmu=np.polynomial.legendre.leggauss(NMU); Pl=np.array([eval_legendre(l,mu) for l in range(L+1)])
R2,MU2=np.meshgrid(rr,mu,indexing='ij'); d2=np.sqrt(np.maximum(R2*R2+D*D-2*R2*D*MU2,1e-12))   # mu=+1 points TOWARD the host
phi_sat=-G*MSAT/np.sqrt(R2*R2+ASAT*ASAT)
if NOHOST: hostb=0*R2; hostnb=0*R2; amb=0*R2+max(rho_mwnb(D),0)+RHO_BG
else:
    hostb=HOSTF*(phi_mwb(d2)-phi_mwb(D)); hostnb=HOSTF*(phi_mwnb(d2)-phi_mwnb(D)); amb=np.maximum(rho_mwnb(d2),0)+RHO_BG     # drops from the satellite centre
def multipoles(rho):
    out=[]
    for l in range(L+1):
        rl=(2*l+1)/2*(rho@(wmu*Pl[l])); Iin=ctz(rl*rr**(l+2),rr,initial=0); b=rl*rr**(1.0-l); Iout=-ctz(b[::-1],rr[::-1],initial=0)[::-1]
        out.append(-4*np.pi*G/(2*l+1)*(Iin*rr**(-l-1.0)+rr**l*Iout))
    return np.array(out)
def law(phinb):
    if MODE=='B1': xb=-(phi_sat+hostb+hostnb)/S2; extra=0.0
    else: xb=-(phi_sat+hostb)/S2; extra=hostnb/S2
    xb0=xb[0].mean(); psi=(phinb-phinb[0].mean())*(0 if NOSELF else 1)
    return np.maximum(rc*np.exp(np.clip(S(xb)-S(xb0)-psi/S2-extra,-700,50))-amb,0.0)
VOL=2*np.pi*R2*R2*np.gradient(rr)[:,None]*wmu[None,:]
def hostterms(f):
    if NOHOST: return 0*R2,0*R2
    return f*(phi_mwb(d2)-phi_mwb(D)), f*(phi_mwnb(d2)-phi_mwnb(D))
def law_f(phinb,f):
    hb,hn=hostterms(f)
    if MODE=='B1': xb=-(phi_sat+hb+hn)/S2; extra=0.0
    else: xb=-(phi_sat+hb)/S2; extra=hn/S2
    xb0=xb[0].mean(); psi=(phinb-phinb[0].mean())*(0 if NOSELF else 1)
    return np.maximum(rc*np.exp(np.clip(S(xb)-S(xb0)-psi/S2-extra,-700,50))-amb,0.0)
def picard(rho,f,iters,damp):
    hist=[]; phinb=None
    for it in range(iters):
        phinb=(multipoles(rho).T@Pl); rn=law_f(phinb,f); res=np.sum(np.abs(rn-rho)*VOL)/np.sum(rho*VOL); hist.append(res); rho=rho+damp*(rn-rho)
        if res<1e-5: break
    return rho,hist
# stage 0: spherical start -- 1D isothermal solve in the LMC's own potential with ambient at D (truncated), then Picard with host OFF
amb0=max(rho_mwnb(D),0)+RHO_BG; xb0s=G*MSAT/ASAT/S2
def f1d(r,y):
    xb=G*MSAT/np.sqrt(r*r+ASAT*ASAT)/S2; return [G*y[1]/r**2,4*np.pi*r*r*max(rc*np.exp(np.clip(S(xb)-S(xb0s)-y[0]/S2,-700,50))-amb0,0.0)]
s1=solve_ivp(f1d,[1e-4,rr[-1]],[0.,0.],rtol=1e-8,atol=1e-3,dense_output=True)
psi1=np.array([s1.sol(r)[0] for r in rr]); rho=np.maximum(rc*np.exp(S(G*MSAT/np.sqrt(rr*rr+ASAT*ASAT)/S2)-S(xb0s)-psi1/S2)-amb0,0.0)[:,None]*np.ones((1,NMU))
t0=time.time(); hist=[]; conv=True; stages=[0.0] if NOHOST else [0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.5,0.7,0.85,1.0]
if HOSTF<1.0: stages=[x for x in stages if x<=HOSTF]+([HOSTF] if HOSTF not in stages else [])
for f in stages:
    rho,h=picard(rho,f,ITERS,DAMP); hist+=h; ok=h[-1]<1e-4
    if not ok:
        try:
            x=newton_krylov(lambda x: (law_f((multipoles(x.reshape(NR,NMU)).T@Pl),f)-x.reshape(NR,NMU)).ravel(),rho.ravel(),f_tol=1e-6*np.max(rho),maxiter=30,verbose=0)
            r2=x.reshape(NR,NMU); res2=np.sum(np.abs(law_f((multipoles(r2).T@Pl),f)-r2)*VOL)/np.sum(r2*VOL); 
            if res2<1e-4: rho=r2; ok=True; hist.append(res2)
        except Exception as e: pass
    if rho.max()>1e20: ok=False
    Mtot=float(np.sum(rho*VOL))
    im=np.unravel_index(np.argmax(rho),rho.shape); near=float(np.sum((rho*VOL)[MU2>0])); far=float(np.sum((rho*VOL)[MU2<0]))
    print(f'  stage hostf {f:.2f}: converged {ok} last res {h[-1]:.1e} (iters {len(h)})  total {Mtot:.3e} near {near:.3e} far {far:.3e}  max rho {rho.max():.3e} at r {rr[im[0]]:.2f} mu {mu[im[1]]:+.2f}')
    if not ok: conv=False; print('  STOP: no converged solution at this host strength'); break
phinb=(multipoles(rho).T@Pl); res=hist[-1]
print(f'[{TAG}] mode {MODE} sigma_sat {SIG:.1f} converged {conv} last stage hostf {f:.2f} {time.time()-t0:.1f}s  max rho {rho.max():.3e}')
# outputs: masses within spheres about the satellite (grid quadrature), hemispheres, extents on axis
vol=2*np.pi*R2*R2*np.gradient(rr)[:,None]*wmu[None,:]
def Msph(R,hemi=None):
    m=(R2<=R); 
    if hemi=='near': m&=(MU2>0)
    if hemi=='far': m&=(MU2<0)
    return float(np.sum(rho[m]*vol[m]))
Mb=lambda R: MSAT*R**3/(R*R+ASAT*ASAT)**1.5
near=rho[:,np.argmax(mu)]; far=rho[:,np.argmin(mu)]
def ext(prof):
    c=np.cumsum(prof*rr**2*np.gradient(rr)); return float(np.interp(0.999*c[-1],c,rr))
out=dict(tag=TAG,mode=MODE,sigma=float(SIG),Msat=MSAT,D=D,converged=bool(conv),hostf_reached=float(f),iters=len(hist),res=float(hist[-1]),
         M_exc={str(R):Msph(R) for R in (10.,30.,32.8,50.)},M_exc_near={str(R):Msph(R,'near') for R in (10.,30.,32.8,50.)},M_exc_far={str(R):Msph(R,'far') for R in (10.,30.,32.8,50.)},
         M_b={str(R):float(Mb(R)) for R in (10.,30.,32.8,50.)},extent_toward_host=ext(near),extent_away=ext(far),rho_amb_at_D=float(rho_mwnb(D)+RHO_BG))
# spherically averaged M(<r) for the dwarf Jeans control
Mavg=np.cumsum(np.sum(rho*vol,axis=1)); out['Mavg_r']=rr[::20].tolist(); out['Mavg']=Mavg[::20].tolist()
print(json.dumps({k:v for k,v in out.items() if k not in ('Mavg_r','Mavg')},indent=1)); json.dump(out,open(f'twocentre_b_{TAG}_20261102.json','w'))
