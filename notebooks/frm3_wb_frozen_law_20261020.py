# frm3_wb_frozen_law_20261020.py [UNRATIFIED] 2026-10-20 14:50 -- FRM-3 under the FROZEN LAW (R1 + NEST-1 C), replacing R2 (which was QUMOND+EFE floor-clipped, GE=1.2).
# Binary: total mass M, separation s. Two nb populations: (A) the binary's own (sub-system sigma^2 = sqrt(G a0 M)/2, own rho_c, baryon exponent capped at 2);
# (B) the MW's nb passing through (host sigma ~124 km/s, local rho ~0.01 Msun/pc^3, boost exp(phi_b/sigma_host^2) capped e^2).
# gamma-1 ~ M_nb(<s)/M for mass inside the orbit (shell theorem for the near-spherical bump). Report for rho_c = 10^7.45 (SPARC median), 10^8.0 (dSph), 10^8.5 (ceiling).
import numpy as np
G=4.30091e-6; a0=3.7e3; pc=1e-3; AU=4.848e-6*pc
M=2.0; s2=np.sqrt(G*a0*M)/2; sig=np.sqrt(s2)
print('binary M=%.0f Msun: sub-system sigma = %.3f km/s (sigma^2 = %.4f)'%(M,sig,s2))
rho_mw=0.01/pc**3   # Msun/kpc^3
sig_host=124.
for s_kAU in (7,14,28):
    s=s_kAU*1e3*AU; phib=G*M/s; xb=phib/s2
    print('\n s = %2d kAU = %.4f pc: phi_b = %.3f (km/s)^2, x_b = phi_b/sigma^2 = %.1f -> cap active (S=2) %s'%(s_kAU,s/pc,phib,xb,xb>2))
    for lrc in (7.45,8.0,8.5):
        rc=10**lrc; rho_max=np.e**2*rc            # ceiling of the binary's own population inside the capped region
        rcap=G*M/(2*s2)                            # radius inside which x_b > 2
        Min=4/3*np.pi*min(s,rcap)**3*rho_max + (4/3*np.pi*(s**3-min(s,rcap)**3)*rc if s>rcap else 0)
        print('   own nb, log rho_c %.2f: ceiling %.2f Msun/pc^3, cap radius %.3f pc, M_nb(<s) = %.1e Msun -> gamma-1 ~ %.1e'%(lrc,rho_max*pc**3,rcap/pc,Min,Min/M))
    boost=np.exp(phib/sig_host**2); Mmw=4/3*np.pi*s**3*rho_mw*boost
    print('   MW nb through-pass: boost %.6f, M(<s) = %.1e Msun -> gamma-1 ~ %.1e'%(boost,Mmw,Mmw/M))
print('\nFor gamma = 1.1 at 28 kAU one needs rho ~ %.0f Msun/pc^3 = 10^%.1f Msun/kpc^3 inside the orbit'%(0.1*M/(4/3*np.pi*(28e3*AU)**3)*pc**3, np.log10(0.1*M/(4/3*np.pi*(28e3*AU)**3))))
