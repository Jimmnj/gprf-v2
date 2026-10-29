# n9_dwarfs_20261003.py  [UNRATIFIED]
# N9: satellite dSph velocity dispersion from the canon law with the HOST halo as ambient density.
# rho_nb = rho_amb[exp(phi_loc/s2) - 1], phi_loc = phi0 - psi, rho_c = rho_amb exp(phi0/s2) (central medium density).
# Baryons: Plummer sphere, M_b, scale a = r_half/1.305 (3D half-mass radius = 1.305 a).
# sigma (law) = BTFR input: s2 = 1/2 sqrt(G a0 M_b). rho_amb = host v_c^2/(4 pi G D^2) [+ rho_bg].
# Tracer: isotropic Jeans in the total potential, Plummer tracers, projected sigma_los within r_half.
# Compare: isolated (rho_amb = 35.5), truncated (host), MOND isolated (81 sigma^4/(4 G a0) = M -> sigma_iso),
# MOND EFE quasi-Newtonian (G_eff = G a0/g_ext), Wolf estimator check.
import numpy as np
from scipy.integrate import solve_ivp, quad
G=4.30091e-6; a0k=3.7e3; rho_bg=35.5
def host_rho(vc,D): return vc**2/(4*np.pi*G*D**2)
def solve(Mb,a,s2,rho_amb,lrc,rmax):
    rc=10**lrc
    rhob=lambda r: 3*Mb/(4*np.pi*a**3)*(1+r*r/(a*a))**-2.5
    def f(r,y):
        psi,Mtot=y; r=max(r,1e-5)
        return [G*Mtot/r**2, 4*np.pi*r*r*(rhob(r)+max(rc*np.exp(-psi/s2)-rho_amb,-rho_amb))]
    s=solve_ivp(f,[1e-5,rmax],[0,0],dense_output=True,rtol=1e-8,atol=1e-3,max_step=rmax/2000)
    return s
def sigma_los(s,a,rh,rmax):
    # isotropic Jeans: nu sigma_r^2 (r) = int_r^inf nu G M / r'^2 dr' ; nu Plummer (tracers follow light)
    nu=lambda r:(1+r*r/(a*a))**-2.5
    M=lambda r: s.sol(min(r,rmax))[1]
    rr=np.logspace(-2,np.log10(rmax),600)
    integrand=np.array([nu(r)*G*M(r)/r**2 for r in rr])
    cum=np.concatenate([[0],np.cumsum(0.5*(integrand[1:]+integrand[:-1])*np.diff(rr))])
    nus2=cum[-1]-cum            # int_r^rmax
    # project: sigma_los^2(R) = 2 int_R^inf nu sigma_r^2 r dr/sqrt(r^2-R^2) / Sigma(R); isotropic
    def Sig(R): return quad(lambda z:nu(np.sqrt(R*R+z*z)),0,rmax,limit=200)[0]*2
    def nus2_of(r): return np.interp(r,rr,nus2)
    def los2(R): return 2*quad(lambda z:nus2_of(np.sqrt(R*R+z*z)),0,rmax,limit=200)[0]/Sig(R)
    # luminosity-weighted within projected r_half
    Rs=np.linspace(0.02*rh,rh,25); w=np.array([Sig(R)*R for R in Rs]); l2=np.array([los2(R) for R in Rs])
    return np.sqrt(np.sum(w*l2)/np.sum(w))
def run(name,Mb,rh,vc_host,D,sig_obs,lrcs=(5.5,6.0,6.5,7.0,7.5)):
    a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb)
    rho_amb=host_rho(vc_host,D)+rho_bg
    Redge=np.sqrt(s2/(2*np.pi*G*rho_amb))
    sig_mond_iso=(4*G*a0k*Mb/81)**0.25
    gext=vc_host**2/D; Geff=G*a0k/gext
    r12=4/3*rh; sig_efe=np.sqrt(Geff*Mb/(3*r12))     # quasi-Newtonian Wolf
    sig_newt=np.sqrt(G*Mb/(3*r12))
    print(f'\n== {name}: Mb={Mb:.2e} r_half={rh} kpc host vc={vc_host} D={D} kpc  obs sigma={sig_obs}')
    print(f'   law sigma={np.sqrt(s2):.2f}  rho_amb={rho_amb:.2e}  R_edge={Redge:.2f} kpc  R_edge/r_half={Redge/rh:.2f}')
    print(f'   MOND isolated {sig_mond_iso:.2f} | MOND EFE quasi-Newtonian {sig_efe:.2f} (e_N={gext/a0k:.2f}) | Newton stars only {sig_newt:.2f}')
    print('   log rho_c   phi0/s2   sigma_los ISOLATED   sigma_los TRUNCATED   M_nb(<r_half)/M_b iso/trunc   M(<R_edge)')
    for lrc in lrcs:
        out=[]
        for amb in (rho_bg,rho_amb):
            rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb)))
            s=solve(Mb,a,s2,amb,lrc,rmax)
            sl=sigma_los(s,a,rh,rmax)
            Mh=s.sol(rh)[1]; Mbh=Mb*(rh/a)**3/(1+(rh/a)**2)**1.5
            out.append((sl,(Mh-Mbh)/Mb,s.sol(min(Redge,rmax))[1]))
        print(f'   {lrc:5.1f}      {np.log(10**lrc/rho_amb):5.2f}      {out[0][0]:6.2f}               {out[1][0]:6.2f}            {out[0][1]:6.1f} / {out[1][1]:6.1f}        {out[1][2]:.2e}')
# And XIX: Collins+2020 sigma 7.8 +1.7-1.5, r_half 3.1; M_V ~ -10.1 -> L_V ~ 1e6, M/L*=1.6 -> M_b 1.6e6 (check); D(M31) ~ 115 kpc (3D)
run('And XIX',1.6e6,3.1,230.,115.,'7.8 +1.7/-1.5')
run('And XIX low-mass / far',0.8e6,3.1,230.,130.,'7.8 +1.7/-1.5',lrcs=(6.0,7.0))
# Crater II: Caldwell+2017 sigma 2.7+-0.3, r_half 1.07 kpc, L_V 1.6e5 -> M_b 2.6e5; D(MW) 117 kpc
run('Crater II',2.6e5,1.07,220.,117.,'2.7 +- 0.3')
