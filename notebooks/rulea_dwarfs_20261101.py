# rulea_dwarfs_20261101.py [UNRATIFIED] PREREG_RULEA_20261101. Rule A (host-field truncation) on the 16-dwarf protocol + LMC.
import numpy as np, pandas as pd
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
from gate2c_robustness_20260925 import load_table
T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
lM=np.log10(F.Mb.values); lS=np.log10(F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2)); fg=1.33*F.MHI.values*1e9/F.Mb.values; y=F.l_tanh2.values
X=np.column_stack([np.ones(len(y)),lM,lS,fg]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
exec(open('dwarf_prereg_20261020.py').read().split('rows=[]')[0].split('dw=[')[1].join(['dw=[','']) ) if False else None
dw=[('Fornax',1.4e7,0.668,11.7,0.9,0.9,220.,147.,'MW'),('Sculptor',1.4e6,0.260,9.2,1.1,1.1,220.,86.,'MW'),('Carina',2.4e5,0.241,6.6,1.2,1.2,220.,105.,'MW'),
('Sextans',4.1e5,0.682,7.9,1.3,1.3,220.,86.,'MW'),('Draco',2.7e5,0.196,9.1,1.2,1.2,220.,76.,'MW'),('UrsaMinor',2.0e5,0.280,9.5,1.2,1.2,220.,76.,'MW'),
('LeoI',3.4e6,0.246,9.2,1.4,1.4,220.,254.,'MW'),('LeoII',5.9e5,0.151,6.6,0.7,0.7,220.,233.,'MW'),('CraterII',1.6e5,1.07,2.7,0.3,0.3,220.,117.,'MW'),
('AndXIX',1.0e6,3.1,7.8,1.7,1.5,230.,115.,'M31'),('AndII',9.3e6,1.18,7.3,0.8,0.8,230.,184.,'M31'),('AndXXI',4.6e5,0.77,4.5,1.2,1.0,230.,134.,'M31'),
('Cetus',2.6e6,0.70,11.0,1.6,1.3,220.,755.,'ISO'),('Cetus_K14',2.6e6,0.70,8.3,1.0,1.0,220.,755.,'ISO'),('Tucana',5.6e5,0.28,6.2,1.6,1.3,220.,880.,'ISO'),('AndXVIII',6.3e5,0.31,9.7,2.3,2.3,230.,450.,'ISO')]
class Frozen:  # wrap a solution with M(r) frozen beyond rA
    def __init__(s,sol,rA,rmax,Mb,a): s.s=sol; s.rA=rA; s.rmax=rmax; s.Mb=Mb; s.a=a
    def Mbar(s,r): return s.Mb*(r/s.a)**3/(1+(r/s.a)**2)**1.5
    def sol(s,r):
        if r<=s.rA: return s.s.sol(r)
        MA=s.s.sol(s.rA)[1]; return np.array([0.0,MA+s.Mbar(r)-s.Mbar(s.rA)])  # nb frozen beyond rA, baryons continue
rows=[]
for name,LV,rh,so,ep,em,vc,Dh,cls in dw:
    Mb=1.6*LV; lSh=np.log10(Mb/(2*np.pi*(rh*1e3)**2)); lrc=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh
    a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb); amb=host_rho(vc,Dh)+rho_bg; gh=vc**2/Dh
    rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb))); s=solve(Mb,a,s2,amb,lrc,rmax)
    rr=np.logspace(-3,np.log10(rmax),1500); g=np.array([G*s.sol(r)[1]/r**2 for r in rr]); gmax=g.max()
    sg0=sigma_los(s,a,rh,rmax)
    if gmax<gh:
        rA=0.0; sgA=np.sqrt(G*Mb/(3*4/3*rh))  # Newton stars only (Wolf)
        # better: Jeans with baryons only
        sb=solve(Mb,a,s2,amb,-99,rmax); sgA=sigma_los(sb,a,rh,rmax)
    else:
        rA=rr[np.where(g>=gh)[0][-1]]
        if rA>=rmax*0.99: sgA=sg0
        else:
            # baryon mass beyond rA: Plummer fraction
            fb_out=1-(rA/a)**3/(1+(rA/a)**2)**1.5
            sgA=sigma_los(Frozen(s,rA,rmax,Mb,a),a,rh,rmax)
    err=np.sqrt(0.10**2+(0.5*(ep+em)/so/np.log(10))**2); dl0=np.log10(sg0/so); dlA=np.log10(sgA/so)
    rows.append(dict(name=name,cls=cls,gh=gh,gmax=gmax,rA=rA,rh=rh,s_obs=so,s_cur=sg0,s_A=sgA,dl_cur=dl0,dl_A=dlA,bar=err,PASS_cur=abs(dl0)<=err,PASS_A=abs(dlA)<=err))
    print('%-10s g_host %6.0f g_own,max %6.0f  r_A %5.2f kpc (r_h %.2f) | obs %4.1f | current %5.1f %s | rule A %5.1f %s'%(name,gh,gmax,rA,rh,so,sg0,'P' if abs(dl0)<=err else 'F',sgA,'P' if abs(dlA)<=err else 'F'))
R=pd.DataFrame(rows); R.to_csv('rulea_dwarfs_20261101.csv',index=False)
print('\nCURRENT: %d/16 pass, rms %.3f  | RULE A: %d/16 pass, rms %.3f, mean %+.3f'%(R.PASS_cur.sum(),np.sqrt((R.dl_cur**2).mean()),R.PASS_A.sum(),np.sqrt((R.dl_A**2).mean()),R.dl_A.mean()))
# LMC
Mb=3e9; s2=0.5*np.sqrt(G*a0k*Mb); gh=200.**2/50.; amb=host_rho(200.,50.)+rho_bg
for lrc in (7.45,):
    a=1.5; rmax=60.; s=solve(Mb,a,s2,amb,lrc,rmax); rr=np.logspace(-2,np.log10(rmax),2000); g=np.array([G*s.sol(r)[1]/r**2 for r in rr])
    rA=rr[np.where(g>=gh)[0][-1]]; rt=np.sqrt(s2/(2*np.pi*G*amb))
    print('LMC: sigma %.1f  g_host %.0f  g_own,max %.0f  r_A %.1f kpc  M_tot(<r_A) %.2e  | density-match r_t %.1f kpc M_tot(<r_t) %.2e | need 1-2e11'%(np.sqrt(s2),gh,g.max(),rA,s.sol(rA)[1],rt,s.sol(min(rt,rmax))[1]))
