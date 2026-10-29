# mistele_dsigma_fit_20261021.py [UNRATIFIED] 2026-10-21 02:45 -- Delta Sigma-space fit to Mistele+2024 Table 1 with the satellite-body term (Mistele-persona 1a/1b; owed).
# Forward model per mass bin: 3D excess density = host tail [sigma^2/(2piG r^2) - rho_amb, to r_t] + satellite term [rho_sat(r) ~ r for r < 1 Mpc, normalised to f_sat x host M(<1 Mpc);
# f_sat from satellite_body_term (N(>1e8)=10): 0.46/0.28/0.25/0.17 for bins 1-4; scaled by a factor s_sat in {0, 1, 3}]. Project to Sigma(R), form Delta Sigma(R) = Sbar(<R) - Sigma(R).
# Deproject AS MISTELE DOES (their systematic column): assume Delta Sigma ~ 1/R beyond the last bin (SIS), recover Sbar via dSbar/dR = -2 Delta Sigma/R, Sigma = Sbar - Delta Sigma,
# Abel-invert Sigma -> rho -> M(<r) -> v_c, sampled at their radii; chi2 against their v_c with STATISTICAL errors only (the systematic is the same assumption applied to both).
# V_f free per bin; grid in log rho_amb. Output: 95% upper bound on rho_amb per bin and joint, for each s_sat. Caveat: our deprojection is a generic spherical inversion; theirs may differ in detail.
import numpy as np, pandas as pd
G=4.30091e-6; from mistele_table import load as _mload
T=_mload('All')  # byte-position parser (fix 2026-10-22)
FSAT={1:0.46,2:0.28,3:0.25,4:0.17}; FSAT={1:0.46,2:0.28,3:0.25,4:0.17}
rgrid=np.logspace(0,4.3,480)   # 1 kpc .. 20 Mpc
def rho3d(r,sig,ra,fsat):
    rt=sig/np.sqrt(2*np.pi*G*ra); host=np.where(r<rt,sig**2/(2*np.pi*G*r**2)-ra,0.0)
    Mh1=2*sig**2*1e3/G-4/3*np.pi*ra*min(1e3,rt)**3
    # satellite term: rho_sat = k r for r<1 Mpc; M(<1Mpc)= k pi 1e12 -> k = fsat*Mh1/(pi*1e12)
    k=fsat*Mh1/(np.pi*1e12); sat=np.where(r<1e3,k*r,0.0); return host+sat
def project(rho_r,Rs):
    # Sigma(R) = 2 int_R^inf rho(r) r dr / sqrt(r^2-R^2), on rgrid
    out=[]
    for R in Rs:
        m=rgrid>R; r=rgrid[m]; integrand=rho_r[m]*r/np.sqrt(r**2-R**2+1e-9); out.append(2*np.trapezoid(integrand,r))
    return np.array(out)
Rs=np.logspace(0,np.log10(2000),160)
def observe(sig,ra,fsat,Rdata):
    rho=rho3d(rgrid,sig,ra,fsat); Sig=project(rho,Rs)
    Sbar=np.array([2*np.trapezoid(Sig[:i+1]*Rs[:i+1],Rs[:i+1])/Rs[i]**2 if i>0 else Sig[0] for i in range(len(Rs))]); DS=Sbar-Sig
    # 'observed' Delta Sigma sampled at data radii; deproject assuming DS ~ 1/R beyond Rmax
    Rmax=Rdata.max(); DSd=np.interp(Rdata,Rs,DS)
    Rmin=Rdata.min(); Rext=np.logspace(0,np.log10(50*Rmax),400); DSext=np.where(Rext<Rmin,DSd[0]*Rmin/Rext,np.where(Rext<=Rmax,np.interp(Rext,Rdata,DSd),DSd[-1]*Rmax/Rext))
    # Sbar(R) = 2 int_R^inf DS/R' dR' (Sbar -> 0 at infinity for 1/R tail? Sbar ~ 2 DS for SIS -> use Sbar(Rfar) = 2 DS(Rfar))
    Sb=np.array([2*np.trapezoid(DSext[i:]/Rext[i:],Rext[i:])+2*DSext[-1] for i in range(len(Rext))]); Sg=Sb-DSext
    # Abel: rho(r) = -1/pi int_r^inf dSigma/dR / sqrt(R^2-r^2) dR
    dS=np.gradient(Sg,Rext); rr=Rext[:-40]; rho=[]
    for r in rr:
        t=np.concatenate([[0],np.logspace(np.log10(0.01*r),np.log10(Rext[-1]),300)]); Rq=np.sqrt(r**2+t**2); f=np.interp(Rq,Rext,dS)
        rho.append(-1/np.pi*np.trapezoid(f/Rq,t))
    rho=np.array(rho)
    M=np.concatenate([[0],np.cumsum(0.5*(4*np.pi*rho[1:]*rr[1:]**2+4*np.pi*rho[:-1]*rr[:-1]**2)*np.diff(rr))]); M+=4*np.pi*rho[0]*rr[0]**3   # inner mass for rho ~ 1/r^2 inside the first grid point
    vc=np.sqrt(np.maximum(G*M/rr,0)); return np.interp(Rdata,rr,vc)
import sys
lras=np.arange(1.0,4.01,0.2)
for ssat in (float(sys.argv[1]),):
    chi_tot=np.zeros_like(lras); line=[]
    for b in range(1,5):
        d=T[(T.R>=100)&T['v%d'%b].notna()]; r=d.R.values; v=d['v%d'%b].values; e=d['s%d'%b].values
        chis=[]
        for lra in lras:
            best=1e9
            for vf in np.arange(110,300,8.0):
                m=observe(vf/np.sqrt(2),10**lra,ssat*FSAT[b],r); c=np.sum(((v-m)/e)**2); best=min(best,c)
            chis.append(best)
        chis=np.array(chis); chi_tot+=chis; dch=chis-chis.min(); up95=10**lras[np.where(dch<=4)[0].max()]; print('  s_sat %.0f bin %d: chi2min %.1f at rho_amb %.0f; 95%% upper %.0f'%(ssat,b,chis.min(),10**lras[np.argmin(chis)],up95),flush=True)
    d=chi_tot-chi_tot.min(); print('s_sat=%.0f JOINT: best rho_amb %.0f, 68%% upper %.0f, 95%% upper %.0f, chi2min %.1f (41 pts)'%(ssat,10**lras[np.argmin(chi_tot)],10**lras[np.where(d<=1)[0].max()],10**lras[np.where(d<=4)[0].max()],chi_tot.min()),flush=True)
