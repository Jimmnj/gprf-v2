# test5_gc_nb_fraction_20261018.py [UNRATIFIED] Globular-cluster nb fraction under the FRM-4-FREEZE law (cap S<=2, sub-system BTFR sigma).
# rho_c for a GC is NOT fixed by the law (fitted in galaxies); two cases bracket it: 10^7.45 (SPARC capped median) and 10^8.0 (dSph, Gilmore/Strigari).
# nb density = max(cored isothermal(rho_c, sigma_GC), e^S_cap x rho_amb) is NOT used; instead: GC's own cored isothermal (r_c = 3 sigma/sqrt(4 pi G rho_c)),
# rho ~ rho_c/(1+(r/r_c)^2), truncated at the Jacobi radius r_J = R (m/(3 M_MW(<R)))^(1/3); plus host ambient boosted by at most e^2 (cap).
# Stars: Plummer, a = r_h/1.305. GC parameters [mem, Harris-catalogue-like; VERIFY]; MW isothermal v_c 200 km/s.
# Pass bar (declared): nb/M* inside r_h < 0.1 (GC M/L show no missing mass in cores).
import numpy as np
G=4.30091e-6; a0=3703.; vc=200.
GCS=[('47 Tuc',7e5,4.2e-3,7.4),('omega Cen',3.5e6,7.6e-3,6.4),('NGC 2419',9e5,2.4e-2,90.),('Pal 5',1.5e4,2.0e-2,18.)]   # M, r_h (kpc), R_gc (kpc)
def Menc_iso(r,rc0,rc):  # cored pseudo-isothermal enclosed mass
    return 4*np.pi*rc0*rc**3*(r/rc-np.arctan(r/rc))
for lrc in (7.45,8.0):
    rc0=10**lrc; print('--- rho_c = 10^%.2f Msun/kpc^3 (%.3f Msun/pc^3) ---'%(lrc,rc0/1e9))
    for nm,M,rh,R in GCS:
        s2=np.sqrt(G*a0*M)/2; rc=3*np.sqrt(s2/(4*np.pi*G*rc0)); MMW=vc**2*R/G; rJ=R*(M/(3*MMW))**(1/3)
        ramb=vc**2/(4*np.pi*G*R**2); a=rh/1.305
        Mstar=lambda r: M*r**3/(r**2+a**2)**1.5
        nb=lambda r: Menc_iso(min(r,rJ),rc0,rc)+np.exp(2)*ramb*4/3*np.pi*min(r,rJ)**3
        # crossover radius where Plummer density = rho_c/(1+(r/rc)^2)
        rr=np.geomspace(rh/10,rJ,4000); rs=3*M/(4*np.pi*a**3)*(1+rr**2/a**2)**-2.5; rn=rc0/(1+(rr/rc)**2)
        ix=np.argmax(rn>rs); rx=rr[ix] if rn[ix]>rs[ix] else np.nan
        print('%-10s M %.1e r_h %4.1f pc R %5.1f kpc | sigma_GC %.2f km/s r_c %4.0f pc r_J %4.0f pc rho_amb %.1e | nb/M* (<r_h) %.1e  (<r_J) %.2f | crossover %s pc'%(
            nm,M,rh*1e3,R,np.sqrt(s2),rc*1e3,rJ*1e3,ramb,nb(rh)/Mstar(rh),nb(rJ)/Mstar(rJ),'%.0f'%(rx*1e3) if np.isfinite(rx) else 'none<r_J'))
