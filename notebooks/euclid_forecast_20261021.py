# euclid_forecast_20261021.py [UNRATIFIED] 2026-10-21 17:05 (JM 'can you do it?') -- forecast: with Mistele+2024-like stacks (Table 1 'All', four bins,
# R>=100 kpc) but errors scaled by f (1 = KiDS-1000; Euclid-wide ~15x area -> f~0.25) and radial coverage extended from the last bin to R_max with the
# same relative error as the last bin at fractional spacing 0.25 dex, what Delta chi2 separates (a) rho_amb=150 from 35.5, (b) 150 from no truncation,
# (c) 35.5 from no truncation? Fiducial truth = the law at the named rho_amb; V_f per bin fixed at its fitted value (optimistic: no normalisation freedom).
# Model: truncated isothermal tail as in mistele_fit. Delta chi2 of 9 ~ 3 sigma, 25 ~ 5 sigma, for one parameter.
import numpy as np, pandas as pd
G=4.30091e-6; from mistele_table import load as _mload
T=_mload('All')  # byte-position parser (fix 2026-10-22)
def vc(r,sig,ra):
    if ra<=0: return np.sqrt(2)*sig*np.ones_like(r)
    rt=sig/np.sqrt(2*np.pi*G*ra); rr=np.minimum(r,rt); M=2*sig**2*rr/G-4/3*np.pi*ra*rr**3; return np.sqrt(np.maximum(G*M/r,0))
VF={1:136,2:176,3:212,4:258}
def grid(b,Rmax):
    d=T[(T.R>=100)&T['v%d'%b].notna()]; r=d.R.values; e=np.sqrt(d['s%d'%b]**2+d['y%d'%b]**2).values
    rel=e[-1]/d['v%d'%b].values[-1]; ext=[]
    x=r[-1]
    while x*10**0.25<=Rmax: x*=10**0.25; ext.append(x)
    ext=np.array(ext); r2=np.concatenate([r,ext]); e2=np.concatenate([e,rel*VF[b]*np.ones_like(ext)]); return r2,e2
print('Delta chi2 (sum over 4 bins), V_f fixed. Rows: error scale f x radial reach. Columns: 355 vs 35.5 | 355 vs flat | 35.5 vs flat | 630 vs 355')
print('%5s %7s | %10s %10s %10s %10s'%('f','Rmax','355v35.5','355vFlat','35.5vFlat','630v355'))
for f in (1.0,0.5,0.25,0.1):
    for Rmax in (1400,2000,3000,5000):
        out=[]
        for A,B in ((355,35.5),(355,0),(35.5,0),(630,355)):
            dc=0
            for b in range(1,5):
                r,e=grid(b,Rmax); s=VF[b]/np.sqrt(2); dc+=np.sum(((vc(r,s,A)-vc(r,s,B))/(f*e))**2)
            out.append(dc)
        print('%5.2f %7.0f | %10.1f %10.1f %10.1f %10.1f'%(f,Rmax,*out))
