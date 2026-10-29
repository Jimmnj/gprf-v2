# crater2_tides_20261021.py [UNRATIFIED] PRE-REGISTERED 2026-10-21 17:35 (JM 'go': Crater II tidal run with Antlia II control).
# Inputs VERIFIED 2026-10-21 from Ji et al. 2021 (ApJ 921, 32; arXiv:2106.12656), MW+LMC fiducial orbits:
#   Crater II: D_h 117.5, r_h 1.066 kpc, L_V 10^5.21, sigma 2.34 (+0.42/-0.30) [2.7+/-0.3 Caldwell 2017 also shown], r_peri 21.7 (+5.4/-3.9), r_apo 136.7
#   Antlia II: D_h 131.8, r_h 2.541 kpc, L_V 10^5.88, sigma 5.98 (+0.37/-0.36) (gradient removed), r_peri 38.6 (+7.2/-5.8), r_apo 136.2
# Solver: frozen dwarf protocol (n9_dwarfs solve/sigma_los; Plummer stars; BTFR sigma with M_b = 1.6 L_V; rho_c from relation (3) AND universal 10^7.45).
# Host: the LAW's own Milky Way (mw_decline_20261021.csv, universal rho_c): v_c and ambient rho = v^2/(4 pi G r^2) at the dwarf's radius.
# CASES (declared): A today (ambient at present radius, untruncated) = the protocol value; B tracers + mass cut at the Jacobi radius at pericentre
#   r_J = (G m(<r_J) r_p^2 / 2 v_p^2)^(1/3) (flat host), ambient still today's -- the 'stars stripped, nb regrown' reading (nb response instantaneous);
#   C as B with the ambient at pericentre (upper bound on stripping, 'nb frozen at pericentre'). Formal law prediction = A (instantaneous response);
#   B is the tidal reading; C a bound. Verdict per case: |log pred/obs| <= sqrt(0.10^2 + err^2) as in the protocol. Antlia II must pass in the
#   SAME case that rescues Crater II, or the rescue is not a tidal explanation. Sanity: report min M(<r) inside the cut (must be > 0).
import numpy as np, pandas as pd
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
from gate2c_robustness_20260925 import load_table
T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
lM=np.log10(F.Mb.values); lS=np.log10(F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2)); fg=1.33*F.MHI.values*1e9/F.Mb.values; y=F.l_tanh2.values
X=np.column_stack([np.ones(len(y)),lM,lS,fg]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
MW=pd.read_csv('mw_decline_20261021.csv'); vMW=lambda r: np.interp(r,MW.R.values,MW.law745.values)
def jacobi(s,rp,rmax):
    vp=vMW(rp); rr=np.logspace(-2,np.log10(rmax),400); m=np.array([s.sol(r)[1] for r in rr])
    rj=(G*np.maximum(m,0)*rp**2/(2*vp**2))**(1/3); k=np.where(rj<=rr)[0]; return rr[k[0]] if len(k) else rmax
dw=[('CraterII',10**5.21,1.066,2.34,0.42,0.30,117.5,21.7),('CraterII_C17',10**5.21,1.066,2.7,0.3,0.3,117.5,21.7),('AntliaII',10**5.88,2.541,5.98,0.37,0.36,131.8,38.6)]
print('relation (3) beta:',np.round(beta,3))
for name,LV,rh,so,ep,em,Dh,rp in dw:
    Mb=1.6*LV; a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb); Rgc=Dh  # galactocentric ~ heliocentric at >100 kpc (within 8%)
    lSh=np.log10(Mb/(2*np.pi*(rh*1e3)**2)); lrel=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh
    err=np.sqrt(0.10**2+(0.5*(ep+em)/so/np.log(10))**2)
    print('\n%s  M_b %.2e  r_h %.2f  sigma_BTFR %.2f  obs %.2f (+%.2f/-%.2f)  bar %.2f dex  r_peri %.1f  v_host(today/peri) %.0f/%.0f'%(name,Mb,rh,np.sqrt(s2),so,ep,em,err,rp,vMW(Rgc),vMW(rp)))
    for lab,lrc in (('relation(3) %.2f'%lrel,lrel),('universal 7.45',7.45)):
        amb_t=host_rho(vMW(Rgc),Rgc)+rho_bg; amb_p=host_rho(vMW(rp),rp)+rho_bg
        rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb_t)))
        sA=solve(Mb,a,s2,amb_t,lrc,rmax); gA=sigma_los(sA,a,rh,rmax)
        rJ=jacobi(sA,rp,rmax); sB=solve(Mb,a,s2,amb_t,lrc,rJ); gB=sigma_los(sB,a,rh,rJ)
        sC0=solve(Mb,a,s2,amb_p,lrc,rmax); rJc=jacobi(sC0,rp,rmax); sC=solve(Mb,a,s2,amb_p,lrc,rJc); gC=sigma_los(sC,a,rh,rJc)
        mminB=min(sB.sol(r)[1] for r in np.linspace(0.05,rJ,50)); mminC=min(sC.sol(r)[1] for r in np.linspace(0.05,rJc,50))
        for cl,g,rj,mm in (('A today',gA,np.nan,np.nan),('B r_J(peri)',gB,rJ,mminB),('C peri ambient',gC,rJc,mminC)):
            d=np.log10(g/so); print('  %-18s %-15s sigma %5.2f  (%+.2f dex) %s   r_J %s  minM %s'%(lab,cl,g,d,'PASS' if abs(d)<=err else 'FAIL',
                  '%.2f'%rj if rj==rj else '  - ', '%.1e'%mm if mm==mm else '  -'))
