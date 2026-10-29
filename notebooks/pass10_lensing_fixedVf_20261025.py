# pass10_lensing_fixedVf_20261025.py [UNRATIFIED] 2026-10-25 01:20. Sim-review item C4 (simulated Mistele/Lelli; test ours): lensing fit with the
# amplitude FIXED by the law (V_f^2 = 2 sigma^2 = sqrt(G a0 M_b)) instead of free per bin. M_b = M_b,eff of Mistele et al. 2024b Table 2 ('All', R<300 kpc):
# log 10.10 / 10.66 / 10.96 / 11.29 (gas incl.: cold for LTGs, hot for ETGs). Bracket +-0.1 dex for stellar-mass/gas systematics.
# Same model, data (Table 1 'All', R >= 100 kpc), errors (stat (+) syst) and parser as mistele_fit_20261021.py. NOTE: V^4 = G a0 M_b is also MOND's asymptote.
import numpy as np
from mistele_table import load
G=4.30091e-6; a0=3.7e3; T=load('All')
lMb=[10.10,10.66,10.96,11.29]; Vflat_lens=[137.3,182.1,204.7,260.2]; eV=[11.5,10.8,8.2,6.2]
def vc_model(r,sig,lra):
    ra=10**lra; rt=sig/np.sqrt(2*np.pi*G*ra); rr=np.minimum(r,rt); M=2*sig**2*rr/G-4/3*np.pi*ra*rr**3
    return np.sqrt(np.maximum(G*M/r,0))
lras=np.arange(0.5,4.01,0.05)
print('Law amplitude vs lensing V_flat (Mistele+2024b Table 2):')
for b in range(4):
    Vl=(G*a0*10**lMb[b])**0.25; print('  bin %d: log M_b %.2f -> law V_f %.0f km/s (+-0.1 dex: %.0f-%.0f); lensing V_flat %.1f +- %.1f -> ratio %.3f (%+.1f sigma)'%(
        b+1,lMb[b],Vl,(G*a0*10**(lMb[b]-0.1))**0.25,(G*a0*10**(lMb[b]+0.1))**0.25,Vflat_lens[b],eV[b],Vl/Vflat_lens[b],(Vl-Vflat_lens[b])/eV[b]))
for dM in (-0.1,0.0,0.1):
    tot=np.zeros_like(lras); per=[]
    for b in range(1,5):
        d=T[(T.R>=100)&T['v%d'%b].notna()]; r=d.R.values; v=d['v%d'%b].values; e=np.sqrt(d['s%d'%b]**2+d['y%d'%b]**2).values
        sig=(G*a0*10**(lMb[b-1]+dM))**0.25/np.sqrt(2)
        chis=np.array([np.sum(((v-vc_model(r,sig,l))/e)**2) for l in lras]); tot+=chis; per.append((len(r),chis))
    N=sum(n for n,_ in per); i=np.argmin(tot); dchi=tot-tot.min()
    up=lambda lvl:10**lras[np.where(dchi<=lvl)[0].max()]
    at=lambda ra: tot[np.argmin(abs(lras-np.log10(ra)))]
    print('M_b shift %+.1f dex: joint chi2 min %.1f for %d points (no free amplitude) at rho_amb %.0f; 95%% upper %.0f; chi2 at 35.5 / 355 / 1000 = %.1f / %.1f / %.1f'%(
        dM,tot.min(),N,10**lras[i],up(4.0),at(35.5),at(355),at(1000)))
    print('   per bin chi2 at closure 35.5: '+' / '.join('%.1f (%d)'%(c[np.argmin(abs(lras-np.log10(35.5)))],n) for n,c in per))
