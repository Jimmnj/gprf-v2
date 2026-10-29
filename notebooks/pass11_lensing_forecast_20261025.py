# pass11_lensing_forecast_20261025.py [FIX 01:05: SIS comparison now analytic; log-spaced integration -- first run's Delta Sigma column was a grid artefact]
# [UNRATIFIED] 2026-10-25 01:00. Forecast beyond Mistele+2024b: v_c(R) and Delta Sigma relative to an untruncated SIS
# (= MOND-like isolated asymptote, same amplitude) for the four bins, amplitude fixed by the law (V^4 = G a0 M_b, M_b,eff of their Table 2 +0.1 dex, the value that fits <1 Mpc),
# at rho_amb = 35.5 (closure) and 160 / 355 (present 95% bounds, amplitude-fixed fits). Also the MOND external-field radius scaling for comparison (r_EFE ∝ M_b^1/2 at fixed e_N).
import numpy as np
G=4.30091e-6; a0=3.7e3
lMb=[10.10,10.66,10.96,11.29]
def DS(R,sig,ra):
    # Delta Sigma = mean Sigma inside R minus Sigma(R) for the truncated excess profile; log-spaced line-of-sight grid
    rt=sig/np.sqrt(2*np.pi*G*ra)
    def Sig(Rp):
        if Rp>=rt: return 0.0
        zmax=np.sqrt(rt**2-Rp**2); z=np.concatenate([[0],np.geomspace(1e-3*Rp,zmax,800)]); r=np.sqrt(Rp**2+z**2)
        return 2*np.trapezoid(np.maximum(sig**2/(2*np.pi*G*r**2)-ra,0),z)
    Rs=np.geomspace(1.0,R,400); S=np.array([Sig(x) for x in Rs])
    inner=2*np.trapezoid(S*Rs,Rs)/R**2 + (Sig(1.0)*1.0**2)/R**2   # small-R core contribution ~ Sigma(1 kpc)*1 kpc^2
    return inner-Sig(R)
DS_SIS=lambda R,sig: sig**2/(2*G*R)   # exact: for an SIS, Delta Sigma(R) = Sigma(R) = sigma^2/(2 G R)
def vc(r,sig,ra):
    rt=sig/np.sqrt(2*np.pi*G*ra); rr=min(r,rt); return np.sqrt(max(G*(2*sig**2*rr/G-4/3*np.pi*ra*rr**3)/r,0))
Rs=[1,2,3,5]
for ra in (35.5,160,355):
    print('\nrho_amb = %.1f Msun/kpc^3'%ra)
    print(' bin log Mb |  V_f  | r_t (Mpc) | v_c/V_f at R = 1/2/3/5 Mpc | Delta Sigma / SIS at R = 1/2/3/5 Mpc')
    for l in lMb:
        Mb=10**(l+0.1); Vf=(G*a0*Mb)**0.25; sig=Vf/np.sqrt(2); rt=sig/np.sqrt(2*np.pi*G*ra)
        v=[vc(R*1e3,sig,ra)/Vf for R in Rs]; d=[DS(R*1e3,sig,ra)/DS_SIS(R*1e3,sig) for R in Rs]
        print('  %5.2f   | %4.0f  |   %4.2f    | %s | %s'%(l,Vf,rt/1e3,' '.join('%.2f'%x for x in v),' '.join('%.2f'%x for x in d)))
print('\nDownturn-radius ratio, most/least massive bin (log Mb 11.29 vs 10.10): law (M_b^1/4) %.2f ; MOND EFE at fixed e_N (M_b^1/2) %.2f'%(10**((11.29-10.10)/4),10**((11.29-10.10)/2)))
