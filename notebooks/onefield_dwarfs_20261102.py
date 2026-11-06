# onefield_dwarfs_20261102.py [UNRATIFIED] JM 2026-11-02 13:34 "go". Dwarf protocol under STRICT ONE-FIELD: the satellite has no response of its own;
# the host medium (density rho_host at the satellite's distance, host sigma) responds to the satellite's phi_b. Local response (linear in the
# small x_b regime: rho_exc = rho_host[exp(x_b - psi_nb/sigma_host^2) - 1]); same n9 machinery with s2 -> host sigma^2, rho_c -> rho_host exp(x_b0).
# Compare: observed, current protocol (per-body sigma, as in the note), strict one-field, Newton stars-only (Wolf). Host sigma^2 = v_c,host^2/2.
import numpy as np, pandas as pd
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
from gate2c_robustness_20260925 import load_table
T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
lM=np.log10(F.Mb.values); lS=np.log10(F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2)); fg=1.33*F.MHI.values*1e9/F.Mb.values; y=F.l_tanh2.values
X=np.column_stack([np.ones(len(y)),lM,lS,fg]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
exec('dw=['+open('rulea_dwarfs_20261101.py').read().split("\ndw=[")[1].split('class Frozen')[0])
rows=[]
for name,LV,rh,so,ep,em,vc,Dh,cls in dw:
    Mb=1.6*LV; lSh=np.log10(Mb/(2*np.pi*(rh*1e3)**2)); lrc=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh
    a=rh/1.305; s2=0.5*np.sqrt(G*a0k*Mb); amb=host_rho(vc,Dh)+rho_bg
    rmax=max(8*rh,3*np.sqrt(s2/(2*np.pi*G*amb)))
    s_cur=solve(Mb,a,s2,amb,lrc,rmax); sg_cur=sigma_los(s_cur,a,rh,rmax)
    s2h=vc**2/2; xb0=G*Mb/a/s2h; lrc_of=np.log10(amb*np.exp(xb0))      # one-field: central medium density = ambient x exp(x_b0)
    s_of=solve(Mb,a,s2h,amb,lrc_of,8*rh); sg_of=sigma_los(s_of,a,rh,8*rh)
    s_N=solve(Mb,a,s2h,amb,-99,8*rh); sg_N=sigma_los(s_N,a,rh,8*rh)      # baryons only (rho_c -> 0 => excess = -amb floor ~ 0)
    Mexc=s_of.sol(rh)[1]-s_N.sol(rh)[1]
    err=lambda d: ep if d>0 else em
    dl_cur=np.log10(sg_cur/so); dl_of=np.log10(sg_of/so)
    rows.append(dict(name=name,cls=cls,obs=so,cur=sg_cur,onefield=sg_of,newton=sg_N,xb0=xb0,Mexc_rh_over_Mb=Mexc/Mb,
                     dl_cur=dl_cur,dl_of=dl_of,PASS_cur=abs(sg_cur-so)<=err(sg_cur-so),PASS_of=abs(sg_of-so)<=err(sg_of-so)))
    print('%-10s obs %4.1f | current %5.1f %s | ONE-FIELD %5.2f %s | Newton %5.2f | x_b0 %.1e  M_exc(<rh)/M_b %.1e'%(name,so,sg_cur,'P' if rows[-1]['PASS_cur'] else 'F',sg_of,'P' if rows[-1]['PASS_of'] else 'F',sg_N,xb0,Mexc/Mb))
R=pd.DataFrame(rows); R.to_csv('onefield_dwarfs_20261102.csv',index=False)
print('PASS current %d/%d  rms dex %.3f | PASS one-field %d/%d  rms dex %.3f'%(R.PASS_cur.sum(),len(R),np.sqrt((R.dl_cur**2).mean()),R.PASS_of.sum(),len(R),np.sqrt((R.dl_of**2).mean())))
