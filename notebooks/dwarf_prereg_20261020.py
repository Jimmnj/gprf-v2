# dwarf_prereg_20261020.py [UNRATIFIED] PRE-REGISTERED 2026-10-20 19:10 (outline 6.6 / cold review B8). THE DWARF TEST.
# STEP 1: rho_c relation from SPARC frozen-law fits (TEST1 tanh2): log rho_c = a + b log M_b + c log Sigma_h + d f_gas (Sigma_h = M_b/2pi R_h^2, R_h = 1.678 R_d).
# STEP 2: for each dwarf, rho_c from the relation (f_gas = 0, Sigma_h from M_b, r_h), sigma_los within r_h from the N9 solver (Plummer stars, BTFR sigma, host ambient).
#         NO per-object freedom. Beside: Newton stars-only, MOND isolated, MOND EFE (quasi-Newtonian Wolf), all from N9 formulas.
# STEP 3: DECLARED PASS/FAIL per object: |log sigma_pred - log sigma_obs| <= sqrt(0.10^2 + err_obs^2) where 0.10 dex = relation residual (0.65 dex in rho_c -> ~0.1 in sigma,
#         from the N9 scan slope d log sigma/d log rho_c ~ 0.2 at rho_c ~7). SAMPLE VERDICT: the law survives if >= 2/3 of objects pass AND the embedded and
#         isolated subsets do not differ in sign of residual (GPRF predicts no g_ext dependence; MOND-EFE predicts embedded suppressed).
# Inputs: MW classical dSphs VERIFIED 2026-10-20 20:15 against Walker+2009 Table 1 (L_V, r_half, sigma); M31 dwarfs VERIFIED vs McGaugh & Milgrom 2013 I/II tables (McConnachie 2012 photometry; r_1/2 3-D -> projected x3/4): And II 9.3e6/1.66 kpc(3D)/7.3+/-0.8 (Ho+2012); And XXI 4.6e5/1023 pc(3D)/4.5 +1.2/-1.0 (Collins+2013, 32 stars); And XVIII 6.3e5/417 pc(3D)/9.7+/-2.3 (Tollerud+2012, 22) DISPUTED (Collins+2013 <= 2.7, 4 stars). Cetus VERIFIED Taibi+2018: 11.0 +1.6/-1.3 (54 members); Kirby+2014 8.3+/-1.0 [mem] kept as alternate row. others [mem, VERIFY: McConnachie 2012; Kirby+2014; Collins+2013/2020; Fraternali+2009 (Tucana) -> REPLACED 19:25 by Taibi+2020 A&A: 6.2 +1.6/-1.3, 39 members, VERIFIED]: L_V, r_h (kpc), sigma_obs +/- , host, D_host (kpc). M_b = 1.6 L_V.
import numpy as np, pandas as pd
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
from gate2c_robustness_20260925 import load_table
T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
lM=np.log10(F.Mb.values); lS=np.log10(F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2)); fg=1.33*F.MHI.values*1e9/F.Mb.values; y=F.l_tanh2.values
X=np.column_stack([np.ones(len(y)),lM,lS,fg]); beta=np.linalg.lstsq(X,y,rcond=None)[0]; res=(y-X@beta).std()
print('STEP 1 relation (N=%d): log rho_c = %.3f %+.3f log Mb %+.3f log Sigma_h %+.3f f_gas ; resid %.2f dex'%(len(y),*beta,res))
dw=[ # name, L_V, r_h kpc, sigma, +err, -err, host vc, D kpc, class
('Fornax',1.4e7,0.668,11.7,0.9,0.9,220.,147.,'MW'),('Sculptor',1.4e6,0.260,9.2,1.1,1.1,220.,86.,'MW'),('Carina',2.4e5,0.241,6.6,1.2,1.2,220.,105.,'MW'),
('Sextans',4.1e5,0.682,7.9,1.3,1.3,220.,86.,'MW'),('Draco',2.7e5,0.196,9.1,1.2,1.2,220.,76.,'MW'),('UrsaMinor',2.0e5,0.280,9.5,1.2,1.2,220.,76.,'MW'),
('LeoI',3.4e6,0.246,9.2,1.4,1.4,220.,254.,'MW'),('LeoII',5.9e5,0.151,6.6,0.7,0.7,220.,233.,'MW'),('CraterII',1.6e5,1.07,2.7,0.3,0.3,220.,117.,'MW'),
('AndXIX',1.0e6,3.1,7.8,1.7,1.5,230.,115.,'M31'),('AndII',9.3e6,1.18,7.3,0.8,0.8,230.,184.,'M31'),('AndXXI',4.6e5,0.77,4.5,1.2,1.0,230.,134.,'M31'),
('Cetus',2.6e6,0.70,11.0,1.6,1.3,220.,755.,'ISO'),('Cetus_K14',2.6e6,0.70,8.3,1.0,1.0,220.,755.,'ISO'),('Tucana',5.6e5,0.28,6.2,1.6,1.3,220.,880.,'ISO'),('AndXVIII',6.3e5,0.31,9.7,2.3,2.3,230.,450.,'ISO')]
rows=[]
for name,LV,rh,so,ep,em,vc,Dh,cls in dw:
    Mb=1.6*LV; lSh=np.log10(Mb/(2*np.pi*(rh*1e3)**2)); lrc=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh+beta[3]*0.0
    a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb); amb=host_rho(vc,Dh)+rho_bg
    rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb))); s=solve(Mb,a,s2,amb,lrc,rmax); sg=sigma_los(s,a,rh,rmax)
    gext=vc**2/Dh; r12=4/3*rh; sN=np.sqrt(G*Mb/(3*r12)); sMi=(4*G*a0k*Mb/81)**0.25; sMe=np.sqrt(min(G*a0k/gext,1e9)*Mb/(3*r12)) if gext>0 else np.nan
    sMe=min(sMe,sMi)  # EFE cannot exceed isolated
    err=np.sqrt(0.10**2+(0.5*(ep+em)/so/np.log(10))**2); dl=np.log10(sg/so); ok=abs(dl)<=err
    rows.append(dict(name=name,cls=cls,Mb=Mb,rh=rh,lSig=lSh,lrc=lrc,gext_a0=gext/a0k,s_obs=so,s_GPRF=sg,s_Newton=sN,s_MONDiso=sMi,s_MONDefe=sMe,dlog=dl,bar=err,PASS=ok))
    print('%-10s %-3s Mb %.1e r_h %.2f gext/a0 %.2f | rho_c %.2f | obs %5.1f +%.1f/-%.1f | GPRF %5.1f (%+.2f dex, bar %.2f) %s | Newton %4.1f MOND iso %4.1f EFE %4.1f'%(name,cls,Mb,rh,gext/a0k,lrc,so,ep,em,sg,dl,err,'PASS' if ok else 'FAIL',sN,sMi,sMe))
R=pd.DataFrame(rows); R.to_csv('dwarf_prereg_20261020.csv',index=False)
print('\nSAMPLE: %d/%d pass. Residual mean: embedded (MW+M31) %+.2f dex (N=%d), isolated %+.2f (N=%d). GPRF predicts equal; MOND-EFE predicts embedded lower than isolated-MOND.'%(
    R.PASS.sum(),len(R),R[R.cls!='ISO'].dlog.mean(),(R.cls!='ISO').sum(),R[R.cls=='ISO'].dlog.mean(),(R.cls=='ISO').sum()))
for lab,col in (('GPRF',R.s_GPRF),('MOND iso',R.s_MONDiso),('MOND EFE',R.s_MONDefe),('Newton',R.s_Newton)):
    d=np.log10(col/R.s_obs); print('  %-9s rms log(pred/obs) %.3f  mean %+.3f'%(lab,np.sqrt((d**2).mean()),d.mean()))
