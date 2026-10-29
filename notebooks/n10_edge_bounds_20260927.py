# n10_edge_bounds_20260927.py  [UNRATIFIED]
# N10: what the Mistele+2024 "flat to 1 Mpc" lensing result says about GPRF's R_edge, rho_amb,
# the C4 closure a0 = c sqrt(G rho_bg), and the halo mass budget.
# Law (canon 0c): rho_nb = rho_amb[exp(phi_loc/s2) - 1]; exterior solved as self-gravitating
# isothermal-with-floor: psi'' + 2 psi'/r = 4 pi G [rho_c exp(-psi/s2) - rho_amb], psi(0)=0.
# Units: kpc, km/s, Msun. G = 4.30091e-6 kpc (km/s)^2 / Msun.
import numpy as np
from scipy.integrate import solve_ivp
G=4.30091e-6; Mpc=1000.0
H0=70.0; rho_crit=3*(H0/1000.0)**2/(8*np.pi*G)   # Msun/kpc^3 (H0 in km/s/kpc)
Om=0.315; rho_m=Om*rho_crit
a0_SI=1.2e-10; c=2.998e5
# Mistele+2024 Table 2 (50-300 kpc): log Mb centers and Vflat (km/s)
bins=[(10.10,137.3,11.5),(10.66,182.1,10.8),(10.96,204.7,8.2),(11.29,260.2,6.2)]

def profile(s2, rho_amb, Mb, rmax_fac=1.5):
    """Solve PB with point-mass baryons Mb at center (softened) and rho_c chosen so that the
    isothermal asymptote v^2 -> 2 s2 (rho_c set from the flat-part normalization: we pick rho_c
    such that v^2 at 100 kpc ~ 2 s2; scan)."""
    def run(lrc):
        rc=10**lrc
        # exterior problem: baryons irrelevant at Mpc; pure self-gravitating isothermal with floor.
        f=lambda r,y:[ G*y[1]/r**2 , 4*np.pi*r*r*(rc*np.exp(-y[0]/s2)-rho_amb)]
        # note: y[0]=psi (potential rise), y[1]=M_nb excess enclosed (can go negative beyond edge)
        Redge_guess=np.sqrt(s2/(2*np.pi*G*rho_amb))
        r_end=rmax_fac*Redge_guess*3
        s=solve_ivp(f,[0.01,r_end],[0,0],dense_output=True,rtol=1e-8,atol=1e-2,max_step=r_end/1500)
        return s
    # choose rho_c so that v_c at 100 kpc equals sqrt(2 s2) (flat part), bisection in log rho_c
    def v100(lrc):
        s=run(lrc); r=100.0; return G*(Mb+s.sol(r)[1])/r
    # isothermal asymptote v^2 -> 2 s2 is independent of rho_c; rho_c only sets the core.
    # Use SPARC median log rho_c = 7.8 (RESULT_PB_SOLVE_CANON sec 2). No bisection.
    r0=5.0   # King core ~ 2 Rd (SPARC median r0/Rd ~ 1.7); sets rho_c, not v_flat
    lrc=np.log10(9*s2/(4*np.pi*G*r0**2))
    s=run(lrc)
    return s, lrc

def edge_and_decline(s,s2,Mb):
    r=np.logspace(1,np.log10(s.t[-1]),3000)
    M=s.sol(r)[1]; v2=G*(Mb+M)/r
    vf2=2*s2
    # R_edge: where excess density returns to zero = where M(<r) peaks
    iedge=np.argmax(M); Redge=r[iedge]
    def rwhere(frac):
        j=np.where((r>100)&(v2<frac*vf2))[0]
        return r[j[0]] if len(j) else np.nan
    return Redge, rwhere(0.9), rwhere(0.7), M[iedge]

print(f'rho_crit = {rho_crit:.1f} Msun/kpc^3 ; Omega_m rho_crit = {rho_m:.1f} ; 0.2 rho_crit = {0.2*rho_crit:.1f}')
# closure value of rho_bg from a0 = c sqrt(G rho): rho = a0^2/(G c^2) in SI -> Msun/kpc^3
Gsi=6.674e-11; Msun=1.989e30; kpc=3.0857e19
rho_closure=(a0_SI**2/(Gsi*(2.998e8)**2))*kpc**3/Msun
print(f'C4 closure: rho_bg = a0^2/(G c^2) = {rho_closure:.1f} Msun/kpc^3 = {rho_closure/rho_crit:.2f} rho_crit\n')

print('Per-bin analytic bound from "no decline to 1 Mpc" (R_edge > 1 Mpc, s2 = Vflat^2/2):')
for lm,vf,evf in bins:
    s2=vf**2/2; Mb=10**lm
    rho_max=s2/(2*np.pi*G*Mpc**2)
    print(f'  logMb={lm:5.2f} Vf={vf:6.1f}: rho_amb < {rho_max:7.1f} Msun/kpc^3 = {rho_max/rho_crit:5.2f} rho_crit ; '
          f'R_edge(rho_bg=closure)={np.sqrt(s2/(2*np.pi*G*rho_closure))/Mpc:4.2f} Mpc')

print('\nNumerical profile: decline onset and enclosed excess mass, for several rho_amb:')
print('  logMb  rho_amb  R_edge[Mpc]  R(v2=0.9vf2)  R(v2=0.7vf2)  M_nb(<R_edge)[Msun]  M_nb/M_b')
rows=[]
for lm,vf,evf in bins:
    s2=vf**2/2; Mb=10**lm
    for rho_amb in [rho_closure, 100.0, 300.0, 1000.0]:
        s,lrc=profile(s2,rho_amb,Mb)
        Re,r9,r7,Mn=edge_and_decline(s,s2,Mb)
        rows.append((lm,rho_amb,Re,Mn))
        print(f'  {lm:5.2f}  {rho_amb:7.1f}  {Re/Mpc:8.2f}     {r9/Mpc:8.2f}      {r7/Mpc:8.2f}     {Mn:12.3e}   {Mn/Mb:8.0f}')

# Mass budget: integrate n(Mb) M_nb(Mb) over the GSMF. Baldry+2012 double Schechter (h=0.7), per dex per Mpc^3.
# Convert stellar->baryonic crudely: Mb = 1.3 M*. M_nb(Mb) at fixed rho_amb ~ Mb^(3/4) (v^2 R_edge, s2 ~ Mb^(1/2)).
print('\nMass budget: sum over galaxies of M_nb(<R_edge) vs cosmic matter density (Omega_m rho_crit):')
lMs=np.linspace(8.0,11.8,77); Ms=10**lMs
Mstar=10**10.66; p1,a1,p2,a2=3.96e-3,-0.35,0.79e-3,-1.47
phi=np.log(10)*np.exp(-Ms/Mstar)*(p1*(Ms/Mstar)**(a1+1)+p2*(Ms/Mstar)**(a2+1))   # per dex per Mpc^3
Mb_grid=1.3*Ms
for rho_amb in [rho_closure,100.0,300.0,1000.0]:
    # M_nb(<R_edge) = v_flat^2 R_edge / G with v_flat^2 = sqrt(a0 G Mb) (BTFR, a0 in kpc units), R_edge = s2^(1/2)/sqrt(2 pi G rho_amb)
    a0k=a0_SI*kpc/1e6
    vf2=np.sqrt(a0k*G*Mb_grid); s2=vf2/2; Re=np.sqrt(s2/(2*np.pi*G*rho_amb))
    Mnb=vf2*Re/G
    rho_halos=np.trapezoid(phi*Mnb,lMs)/Mpc**3    # Msun/kpc^3
    print(f'  rho_amb={rho_amb:7.1f}: sum n*M_nb = {rho_halos:8.1f} Msun/kpc^3 = {rho_halos/rho_m:6.2f} x Omega_m rho_crit = {rho_halos/rho_amb:6.1f} x rho_amb')
