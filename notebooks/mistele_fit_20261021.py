# mistele_fit_20261021.py [UNRATIFIED] 2026-10-21 01:20 -- (A) likelihood fit of the frozen law's truncated tail to Mistele+2024 Table 1 ('All' sample, bins 1-4),
# errors stat (+) syst in quadrature; (B) Delta Sigma shape of the truncated tail at 1-5 Mpc vs an untruncated SIS. Model beyond ~50 kpc: excess density
# rho(r) = sigma^2/(2 pi G r^2) - rho_amb for r < r_t = sigma/sqrt(2 pi G rho_amb), 0 beyond; v_c^2 = G M(<r)/r with M from the excess. sigma^2 = V_f^2/2 with V_f free
# per bin (BTFR normalisation marginalised over) or fixed to the bin's flat value. Compare to the DEPROJECTED v_c as tabulated (caveat: their deprojection assumes a profile beyond
# the last bin; we fit v_c space, not Delta Sigma -- stated). Output: chi2 vs log rho_amb per bin; 68/95% upper bounds; joint.
import numpy as np, pandas as pd
G=4.30091e-6; from mistele_table import load as _mload
T=_mload('All')  # byte-position parser (fix 2026-10-22)
# bins: use rows with all four present (<=1351 kpc) plus 1836 where available
def vc_model(r,sig,lra):
    ra=10**lra; rt=sig/np.sqrt(2*np.pi*G*ra)
    # M(<r) of excess: integral 4 pi r^2 [sigma^2/(2piG r^2) - ra] = 2 sigma^2 r/G - 4/3 pi ra r^3, for r<rt; constant beyond
    rr=np.minimum(r,rt); M=2*sig**2*rr/G-4/3*np.pi*ra*rr**3
    return np.sqrt(np.maximum(G*M/r,0))
lras=np.arange(0.5,4.01,0.05); out={}
print('Fit of truncated isothermal tail to Mistele+2024 "All" bins (R >= 100 kpc), errors stat (+) syst; V_f fitted per bin at each rho_amb.')
print('%4s %6s | %22s | %28s | %s'%('bin','N','V_f(fit) at chi2min','rho_amb: best / 68% / 95% upper','chi2 at rho_amb = 35.5 / 150 / 1000 / 3000 (dof)'))
chi_tot=np.zeros_like(lras)
for b in range(1,5):
    d=T[(T.R>=100)&T['v%d'%b].notna()]; r=d.R.values; v=d['v%d'%b].values; e=np.sqrt(d['s%d'%b]**2+d['y%d'%b]**2).values
    chis=[];vfs=[]
    for lra in lras:
        best=(1e9,0)
        for vf in np.arange(100,320,1.0):
            m=vc_model(r,vf/np.sqrt(2),lra); c=np.sum(((v-m)/e)**2)
            if c<best[0]: best=(c,vf)
        chis.append(best[0]); vfs.append(best[1])
    chis=np.array(chis); vfs=np.array(vfs); chi_tot+=chis
    i=np.argmin(chis); dchi=chis-chis.min()
    def upper(lvl): 
        ok=np.where(dchi<=lvl)[0]; return 10**lras[ok.max()]
    def chi_at(ra): return chis[np.argmin(abs(lras-np.log10(ra)))]
    print('%4d %6d | %5.0f km/s          | %8.0f / %8.0f / %8.0f | %.1f / %.1f / %.1f / %.1f (%d)'%(b,len(r),vfs[i],10**lras[i],upper(1.0),upper(4.0),chi_at(35.5),chi_at(150),chi_at(1000),chi_at(3000),len(r)-2))
    out[b]=(chis,vfs)
d=chi_tot-chi_tot.min(); i=np.argmin(chi_tot)
print('JOINT: best rho_amb %.0f; 68%% upper %.0f; 95%% upper %.0f Msun/kpc^3; chi2min %.1f'%(10**lras[i],10**lras[np.where(d<=1)[0].max()],10**lras[np.where(d<=4)[0].max()],chi_tot.min()))
# (B) Delta Sigma shape
def DeltaSigma(R,sig,lra):
    ra=10**lra; rt=sig/np.sqrt(2*np.pi*G*ra)
    def Sigma(Rp):
        if Rp>=rt: return 0.0
        z=np.linspace(0,np.sqrt(rt**2-Rp**2),400); r=np.sqrt(Rp**2+z**2); rho=sig**2/(2*np.pi*G*r**2)-ra
        return 2*np.trapezoid(rho,z)
    Rs=np.linspace(1,R,400); S=np.array([Sigma(x) for x in Rs]); Sbar=2*np.trapezoid(S*Rs,Rs)/R**2
    return Sbar-Sigma(R)
sig=182/np.sqrt(2)
print('\n(B) Delta Sigma of the truncated tail / untruncated SIS, bin 2 (V_f 182):')
for lra in (np.log10(35.5),np.log10(150),np.log10(1000)):
    rt=sig/np.sqrt(2*np.pi*G*10**lra)
    vals=[DeltaSigma(R*1e3,sig,lra)/DeltaSigma(R*1e3,sig,-3) for R in (0.5,1,2,3,5)]
    print('  rho_amb %5.0f (r_t %.1f Mpc): ratio at R = 0.5/1/2/3/5 Mpc = '%(10**lra,rt/1e3)+' '.join('%.2f'%x for x in vals))
