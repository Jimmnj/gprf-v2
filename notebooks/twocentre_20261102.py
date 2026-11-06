# twocentre_20261102.py -- PREREG_TWOCENTRE_20261102 (stamped JM 2026-11-02 12:25). Milky Way + LMC under ONE-FIELD.
# Law: rho_nb = rho_c exp[S(x_b)-S(x_b0) - psi_nb/sigma^2] - rho_amb, S(x)=2 tanh(x/2), x_b = -Phi_b/sigma^2, floor -rho_amb.
# One sigma (MW BTFR), one rho_c. Axisymmetric about the MW-LMC line; multipoles l=0..L about the MW centre, full mu range.
# Monopole shot radially (as gate2b_axisym); l>=1 iterated (Picard, damped) with Newton-Krylov fallback.
import numpy as np, sys, json, time
from scipy.special import eval_legendre
from scipy.integrate import cumulative_trapezoid as ctz
from scipy.optimize import newton_krylov
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; RHO_BG=35.5
S=lambda x: 2*np.tanh(x/2)
args=dict(a.split('=') for a in sys.argv[1:])
L=int(args.get('L',40)); NMU=int(args.get('NMU',96)); NR=int(args.get('NR',1200))
LRC=float(args.get('lrc',7.45)); MLMC=float(args.get('MLMC',3.0e9)); ALMC=float(args.get('aLMC',1.5)); D=float(args.get('D',50.0))
GAS=float(args.get('gas',0.0)); AGAS=float(args.get('agas',3.0))   # optional second LMC component
TAG=args.get('tag','primary'); DAMP=float(args.get('damp',0.3)); ITERS=int(args.get('iters',60))
# MW baryons: spherical equivalents (Ou+2024 Table 2 components; Hernquist a = 0.695 Rd for exponential discs)
MW=[(1.55e10,0.70),(3.65e10,0.695*2.35),(8.2e9,0.695*18.24),(1.3e9,0.695*2.57)]
MB_MW=sum(m for m,_ in MW); S2=0.5*np.sqrt(G*a0k*MB_MW); SIG=np.sqrt(S2)
rr=np.logspace(np.log10(0.05),np.log10(3000.),NR)
mu,wmu=np.polynomial.legendre.leggauss(NMU)            # full range [-1,1], weights sum to 2
Pl=np.array([eval_legendre(l,mu) for l in range(L+1)])   # (L+1,NMU)
R2,MU2=np.meshgrid(rr,mu,indexing='ij')
def phi_mw(r): return sum(-G*m/(r+a) for m,a in MW)
def dist_lmc(r,m): return np.sqrt(np.maximum(r*r+D*D-2*r*D*m,1e-12))
def phi_lmc(r,m,M=MLMC):
    d=dist_lmc(r,m); p=-G*M/np.sqrt(d*d+ALMC*ALMC)
    if GAS>0: p=p-G*GAS/np.sqrt(d*d+AGAS*AGAS)
    return p
def multipoles(rho):   # rho (NR,NMU) -> Phi_l(r) l=0..L, potential zero at infinity (outer region beyond rr[-1] ignored)
    out=[]
    for l in range(L+1):
        rl=(2*l+1)/2*(rho@(wmu*Pl[l]))
        Iin=ctz(rl*rr**(l+2),rr,initial=0)
        b=rl*rr**(1.0-l); Iout=-ctz(b[::-1],rr[::-1],initial=0)[::-1]
        out.append(-4*np.pi*G/(2*l+1)*(Iin*rr**(-l-1.0)+rr**l*Iout))
    return np.array(out)
def law(phib,phinb,rc):   # phib,phinb (NR,NMU) absolute potentials; centre values from l=0..L at r->0 (use first grid point)
    xb=-phib/S2; xb0=xb[0].mean(); psi=phinb-phinb[0].mean()
    return np.maximum(rc*np.exp(np.clip(S(xb)-S(xb0)-psi/S2,-700,50))-RHO_BG,-RHO_BG)
def solve(M_lmc,rc,label):
    phib=phi_mw(R2)+phi_lmc(R2,MU2,M_lmc)
    # start: spherical shoot with angle-averaged baryon potential (= monopole of phib), then Picard on all l
    phinb=np.zeros_like(phib); rho=law(phib,phinb,rc); hist=[]
    t0=time.time()
    for it in range(ITERS):
        Phil=multipoles(rho); new=(Phil.T@Pl)
        rho_new=law(phib,new,rc)
        res=np.max(np.abs(rho_new-rho)/(np.abs(rho)+1.0)); hist.append(res)
        rho=rho+DAMP*(rho_new-rho); phinb=new
        if res<1e-5: break
    conv=res<1e-4
    if not conv:
        print(f'[{label}] Picard not converged (res {res:.2e}); Newton-Krylov on rho')
        def F(x):
            r_=x.reshape(NR,NMU); return (law(phib,(multipoles(r_).T@Pl),rc)-r_).ravel()
        x=newton_krylov(F,rho.ravel(),f_tol=1e-6*np.max(np.abs(rho)),maxiter=40,verbose=0)
        rho=x.reshape(NR,NMU); phinb=(multipoles(rho).T@Pl); hist.append(np.max(np.abs(F(x))/(np.abs(rho).ravel()+1)))
        conv=True
    print(f'[{label}] iters {len(hist)} final res {hist[-1]:.2e}  {time.time()-t0:.1f}s')
    return dict(phib=phib,phinb=phinb,rho=rho,Phil=multipoles(rho))
def rho_at(sol,r,m):   # evaluate law at arbitrary points from multipole coefficients (exact in angle, log-interp in r)
    Phil=sol['Phil']; lr=np.log(rr)
    ph=np.array([np.interp(np.log(r),lr,Phil[l]) for l in range(L+1)])
    P=np.array([eval_legendre(l,m) for l in range(L+1)])
    phinb=np.sum(ph*P,axis=0); phib=phi_mw(r)+phi_lmc(r,m,sol['Mlmc'])
    rc=sol['rc']; xb=-phib/S2; xb0=sol['xb0']; psi=phinb-sol['phinb0']
    return np.maximum(rc*np.exp(np.clip(S(xb)-S(xb0)-psi/S2,-700,50))-RHO_BG,-RHO_BG)
def finish(sol,M_lmc,rc):
    sol['Mlmc']=M_lmc; sol['rc']=rc; sol['xb0']=(-sol['phib'][0]/S2).mean(); sol['phinb0']=sol['phinb'][0].mean(); return sol
def sphere_mass(f,Rs,nd=400,na=200):   # integrate f(r,mu) over spheres of radius Rs about the LMC centre (on axis at D)
    ca,wa=np.polynomial.legendre.leggauss(na); out={}
    for R in Rs:
        d=np.linspace(0,R,nd+1)[1:]-R/(2*nd); Dd,CA=np.meshgrid(d,ca,indexing='ij')
        r=np.sqrt(Dd*Dd+D*D+2*Dd*D*CA); m=(D+Dd*CA)/r
        val=f(r.ravel(),m.ravel()).reshape(nd,na)
        out[R]=float(np.sum(2*np.pi*Dd*Dd*val*(R/nd)*wa))
    return out
def vc_perp(sol,Rv):   # circular speed in the plane perpendicular to the MW-LMC line (mu=0): v^2 = r dPhi_tot/dr
    Phil=sol['Phil']; lr=np.log(rr); eps=1e-3
    def phit(r):
        ph=sum(np.interp(np.log(r),lr,Phil[l])*eval_legendre(l,0.0) for l in range(L+1))
        return ph+phi_mw(r)+phi_lmc(r,0.0,sol['Mlmc'])
    return np.sqrt(np.maximum(Rv*(phit(Rv*(1+eps))-phit(Rv*(1-eps)))/(2*eps*Rv),0))
def mb_lmc_sphere(R): # Plummer enclosed (+gas)
    m=MLMC*R**3/(R*R+ALMC*ALMC)**1.5
    if GAS>0: m+=GAS*R**3/(R*R+AGAS*AGAS)**1.5
    return m
if __name__=='__main__':
    rc=10**LRC; Rs=(10.,30.,50.)
    print(f'MW M_b {MB_MW:.3e}  sigma {SIG:.1f} km/s  rho_c 10^{LRC}  LMC M_b {MLMC+GAS:.2e} a {ALMC} at D {D}  L {L} NMU {NMU} NR {NR}  tag {TAG}')
    mw=finish(solve(0.0,rc,'MW alone'),0.0,rc)
    two=finish(solve(MLMC,rc,'MW+LMC'),MLMC,rc)
    f_mw=lambda r,m: rho_at(mw,r,m); f_two=lambda r,m: rho_at(two,r,m)
    M_mw=sphere_mass(f_mw,Rs); M_two=sphere_mass(f_two,Rs)
    # column (i): linearised local factor on the MW-alone nb density
    f_lin=lambda r,m: rho_at(mw,r,m)*(np.exp(S(-phi_lmc(r,m)/S2))-1.0)
    M_lin=sphere_mass(f_lin,Rs)
    Rv=np.array([8.,12.,16.,20.,27.])
    out=dict(tag=TAG,L=L,NMU=NMU,NR=NR,lrc=LRC,MLMC=MLMC+GAS,aLMC=ALMC,D=D,sigma=float(SIG),MB_MW=MB_MW,
             rho_mw_at_D=float(f_mw(np.array([D]),np.array([1.0]))[0]),rho_two_at_LMC=float(f_two(np.array([D-1e-3]),np.array([1.0]))[0]),
             M_mw={str(R):M_mw[R] for R in Rs},M_two={str(R):M_two[R] for R in Rs},M_exc={str(R):M_two[R]-M_mw[R] for R in Rs},
             M_lin={str(R):M_lin[R] for R in Rs},M_b_lmc={str(R):float(mb_lmc_sphere(R)) for R in Rs},
             vc_mw=vc_perp(mw,Rv).tolist(),vc_two=vc_perp(two,Rv).tolist(),Rv=Rv.tolist())
    out['nest_factor_centre']=out['rho_two_at_LMC']/out['rho_mw_at_D']
    print(json.dumps(out,indent=1))
    json.dump(out,open(f'twocentre_{TAG}_20261102.json','w'),indent=1)
    np.savez(f'twocentre_{TAG}_20261102.npz',rr=rr,mu=mu,rho_mw=mw['rho'],rho_two=two['rho'],Phil_two=two['Phil'])
