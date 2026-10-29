# n9c_thindisc_phib_20261017.py  [UNRATIFIED]  N9c: N9b with the EXACT thin-disc central potential for the stellar disc.
# For any axisymmetric razor-thin disc, Phi(0) = -2 pi G int_0^inf Sigma(R) dR. SBdisk(R) is in the rotmod files (L/pc^2); Sigma = Yd*SB.
# Inner 0..r1: SB(r1) constant. Outer: exponential tail from the last two positive points. Bulge and gas: spherical proxy (int v^2/r dr + Keplerian),
# as in N9b. Compared against the N14 free fits (Yd 0.5), same as N9b sec 1.
import numpy as np, pandas as pd, os
Gpc=4.30091e-3; rho_amb=35.5; Yd=0.5; SPARC_DIR=os.path.join('..','sparc_galaxy_data','sparc_database')
N14=pd.read_csv('n14_sigma_free_20261001.csv',index_col=0); rows=[]
def sph(v2,r): return v2[0]/2+np.trapezoid(v2/r,r)+v2[-1]
for g,row in N14.iterrows():
    f_=os.path.join(SPARC_DIR,g+'_rotmod.dat')
    if not os.path.exists(f_): continue
    d=np.atleast_2d(np.loadtxt(f_,comments='#')); r=d[:,0]; vg,vd,vb,SBd,SBb=d[:,3],d[:,4],d[:,5],d[:,6],d[:,7]
    Rpc=r*1e3; Sig=Yd*SBd
    inner=Sig[0]*Rpc[0]; mid=np.trapezoid(Sig,Rpc)
    pos=np.where(Sig>0)[0]; outer=0.0
    if len(pos)>=2:
        i1,i2=pos[-2],pos[-1]; h=(Rpc[i2]-Rpc[i1])/np.log(Sig[i1]/Sig[i2]) if Sig[i1]>Sig[i2] else 0
        outer=Sig[i2]*h if h>0 else 0.0
    Phi_disc_exact=2*np.pi*Gpc*(inner+mid+outer)
    Phi_disc_sph=sph(Yd*vd**2,r); Phi_bul=sph(0.7*vb**2,r); Phi_gas=sph(np.sign(vg)*vg**2,r)
    Phi_new=Phi_disc_exact+Phi_bul+Phi_gas; Phi_old=Phi_disc_sph+Phi_bul+Phi_gas
    rows.append(dict(gal=g,lrhoc=row.lrhoc_free,s2=row.s2_free,lPhi_old=np.log10(Phi_old),lPhi_new=np.log10(Phi_new),
                     disc_ratio=Phi_disc_exact/max(Phi_disc_sph,1e-9),disc_share=Phi_disc_exact/Phi_new,tail_share=2*np.pi*Gpc*outer/max(Phi_disc_exact,1e-9),
                     pred_old=np.log10(rho_amb)+Phi_old/(row.s2_free*np.log(10)),pred_new=np.log10(rho_amb)+Phi_new/(row.s2_free*np.log(10))))
C=pd.DataFrame(rows).set_index('gal'); C.to_csv('n9c_thindisc_phib_20261017.csv')
print('N =',len(C)); print('disc exact / spherical proxy: median %.2f  16-84 %.2f..%.2f ; disc share of Phi_b,0 median %.2f ; exp-tail share median %.2f'%(C.disc_ratio.median(),*np.percentile(C.disc_ratio,[16,84]),C.disc_share.median(),C.tail_share.median()))
def reg(x,y,lab):
    m=np.isfinite(x)&np.isfinite(y); A=np.polyfit(x[m],y[m],1); res=y[m]-np.polyval(A,x[m]); print('%-36s slope %+.3f  scatter %.2f dex  r=%.2f'%(lab,A[0],res.std(),np.corrcoef(x[m],y[m])[0,1]))
for tag in ('old','new'):
    reg(C['lPhi_'+tag].values,C.lrhoc.values,'log rho_c vs log Phi_b,0 (%s)'%tag)
    reg(C['pred_'+tag].values,C.lrhoc.values,'fit vs canon pred (%s)'%tag)
    dd=C['pred_'+tag]-C.lrhoc; print('   pred-fit median %+.2f  std %.2f'%(dd.median(),dd.std()))
