# cluster_release_20261022.py [UNRATIFIED] 2026-10-22 12:10 (JM 'sounds promising': nb returned to the environment).
# ORDER-OF-MAGNITUDE: a Coma-like cluster assembled from field galaxies. Each galaxy's isolated nb mass M_iso = (4/3) sigma^2 r_t/G (rho_amb = 35.5
# closure or 150 fit). Inside the cluster its tail is cut where it meets the cluster density rho_cl(R) = s_cl^2/(2 pi G R^2) (isothermal, s_cl = 1000 km/s);
# retained M_ret = (4/3) sigma^2 r_t(R)/G. Released = M_iso - M_ret, assumed to stay in the cluster. Galaxies: Schechter stellar mass function
# (alpha -1.1, M* 1e11 Msun stellar, 1e8-3e11), normalised to total stellar mass 1.5e13 (Coma, VERIFY), projected number density ~ isothermal
# (N(<R) ∝ R) out to R200 = 2 Mpc. Compare total released (and total retained) nb with the nb a Coma-like cluster needs: M200 ~ 1e15 minus
# baryons ~1.5e14 (ICM ~1.3e14 + stars) -> ~8.5e14 (VERIFY; lensing estimates for Coma range ~0.6-2e15).
import numpy as np
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6; s_cl=1000.; R200=2000.
lm=np.linspace(8,np.log10(3e11),400); M=10**lm; phi=(M/1e11)**(-1.1+1)*np.exp(-M/1e11)  # dN/dlogM ∝ (M/M*)^(alpha+1) e^-M/M*
N=phi/np.trapezoid(phi*M,lm)*1.5e13      # dN/dlogM normalised to total stellar 1.5e13
Ntot=np.trapezoid(N,lm)
s2=np.sqrt(G*a0k*M)/2; sig=np.sqrt(s2)
Rg=np.linspace(20,R200,200); wR=np.ones_like(Rg)/len(Rg)   # N(<R) ∝ R
print('galaxies 1e8-3e11: N = %.0f; total stellar %.2e'%(Ntot,np.trapezoid(N*M,lm)))
for ramb in (35.5,150.,355.):
    rt_iso=sig/np.sqrt(2*np.pi*G*ramb); Miso=4/3*s2*rt_iso/G
    rho_cl=s_cl**2/(2*np.pi*G*Rg**2)
    ret=np.array([np.sum(wR*4/3*s2_*np.minimum(sg/np.sqrt(2*np.pi*G*(rho_cl+ramb)),rt)/G) for s2_,sg,rt in zip(s2,sig,rt_iso)])
    tot_iso=np.trapezoid(N*Miso,lm); tot_ret=np.trapezoid(N*ret,lm); rel=tot_iso-tot_ret
    print('rho_amb %5.1f: field nb per L* (1e11) %.1e | cluster sum: isolated %.2e, retained %.2e, RELEASED %.2e Msun | released / needed(8.5e14) = %.1f'%(ramb,np.interp(11,lm,Miso),tot_iso,tot_ret,rel,rel/8.5e14))
    half=np.cumsum(N*Miso)/np.sum(N*Miso); print('      half of the released budget comes from galaxies below 10^%.1f Msun'%lm[np.searchsorted(half,0.5)])
