# pass9_wb_rhoc_sweep_20261025.py [UNRATIFIED] 2026-10-25 01:05. Sim-review item C3 / pre-post item 11: gamma vs the binary's own rho_c.
# Same model and constants as frm3_wb_frozen_law_20261020.py (own population bounded by e^2 rho_c inside the cap radius, rho_c outside; gamma-1 = M_nb(<s)/M).
# M_nb is linear in rho_c, so the rho_c giving a threshold gamma-1 is closed-form. Also: relation (3) extrapolated to a 2 Msun binary.
import numpy as np
G=4.30091e-6; a0=3.7e3; pc=1e-3; AU=4.848e-6*pc; M=2.0; s2=np.sqrt(G*a0*M)/2; rcap=G*M/(2*s2)
def gm1_per_rhoc(s):
    v=lambda r:4/3*np.pi*r**3
    return (v(min(s,rcap))*np.e**2 + (v(s)-v(rcap) if s>rcap else 0))/M
print('binary 2 Msun; sub-system sigma %.3f km/s; cap radius %.1f kAU'%(np.sqrt(s2),rcap/AU/1e3))
print(' s(kAU) | gamma-1 at 10^7.45 | log rho_c for gamma-1 = 0.10 / 0.05 / 0.02 / 0.01')
for skau in (3,7,14,28):
    k=gm1_per_rhoc(skau*1e3*AU)
    print('  %4d  |   %.1e           |   %s'%(skau,k*10**7.45,'  '.join('%.2f'%np.log10(t/k) for t in (0.10,0.05,0.02,0.01))))
for lab,b in (('published (120-gal fit)',(11.063,-0.365,0.268,-0.929)),('full-sample refit (pass7)',(10.934,-0.381,0.374,-0.808))):
    for Sh in (30,130):
        l=b[0]+b[1]*np.log10(M)+b[2]*np.log10(Sh)
        print('relation (3) %-26s Sigma_h %3d Msun/pc^2, f_gas 0: log rho_c %.2f -> gamma at 28 kAU = %.2f'%(lab,Sh,l,1+gm1_per_rhoc(28e3*AU)*10**l))
print('Globular-cluster bound (TEST5, in note): log rho_c(GC) <~ 7.2')
