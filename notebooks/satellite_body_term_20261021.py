# satellite_body_term_20261021.py [UNRATIFIED] 2026-10-21 01:40 -- the law-specific 'two-halo' term inside 1 Mpc (Mistele-persona review): under R7 each satellite of an isolated
# host is a body with its own nb halo, sigma_s^2 = sqrt(G a0 m)/2, truncated where its tail meets the HOST's tail density sigma_h^2/(2 pi G r^2) at the satellite's radius r:
# r_t,s = sigma_s r / sigma_h, M_nb,s = 2 sigma_s^2 r_t,s / G = 2 sigma_s^3 r / (G sigma_h). Sum over a satellite population below 10% of the host (Mistele isolation), uniformly
# distributed in a 1 Mpc sphere, mass function dN/dlog m ~ m^-0.3 (Schechter faint slope -1.3) from 1e7 to 0.1 M_host, normalised to N(>1e8) = 10 (rough; MW/M31-like within 1 Mpc)
# and 30 (upper). Output: satellite nb mass inside 1 Mpc / host nb mass inside 1 Mpc, and the resulting fractional boost of the host's enclosed mass at 0.5 and 1 Mpc. [counts mem, ROUGH]
import numpy as np
G=4.30091e-6; a0=3.7e3; rng=np.random.default_rng(2)
for Mh,Vf in ((1.26e10,137),(4.6e10,182),(9.1e10,205),(1.95e11,260)):
    sh=Vf/np.sqrt(2); Mhost=lambda R: 2*sh**2*R/G
    for N8 in (10,30):
        # sample masses from dN/dm ~ m^-1.3 on [1e7, 0.1 Mh]; N(>1e8) = N8 -> total N = N8 * (1e7^-0.3 - (0.1Mh)^-0.3)/(1e8^-0.3 - (0.1Mh)^-0.3)
        lo,hi=1e7,0.1*Mh; f=lambda m: m**-0.3; Ntot=int(round(N8*(f(lo)-f(hi))/(f(1e8)-f(hi))))
        res=[]
        for trial in range(200):
            u=rng.random(Ntot); m=(f(lo)-u*(f(lo)-f(hi)))**(-1/0.3)
            r=1e3*rng.random(Ntot)**(1/3)   # kpc, uniform in sphere of 1 Mpc
            ss=np.sqrt(np.sqrt(G*a0*m)/2); Ms=2*ss**3*r/(G*sh)
            rts=ss*r/sh
            res.append((Ms.sum(), Ms[r<500].sum()+ (Ms[(r>=500)]*0).sum(), (m.sum())))
        res=np.array(res); tot=res[:,0].mean(); in500=res[:,1].mean(); mb=res[:,2].mean()
        print('host Mb %.1e (V_f %3d): N(>1e8)=%2d (Ntot %3d, sat baryons %.2f of host): satellite nb inside 1 Mpc = %.2e = %.2f x host nb(<1 Mpc); inside 0.5 Mpc %.2f x host nb(<0.5 Mpc)'%(
            Mh,Vf,N8,Ntot,mb/Mh,tot,tot/Mhost(1e3),in500/Mhost(500)))
