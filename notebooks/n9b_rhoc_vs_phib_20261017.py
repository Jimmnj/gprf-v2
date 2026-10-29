# n9b_rhoc_vs_phib_20261017.py  [UNRATIFIED]  N9b: is rho_c predictable from the baryons alone?
# Canon definition: rho_c = rho_amb exp(phi_0/sigma^2), phi_0 = total central potential depth. Test: replace phi_0 by the BARYON central
# depth Phi_b,0 = int_0^inf vbar^2/r dr (spherical-equivalent; inner solid-body, outer Keplerian) and see whether
# log rho_c,pred = log rho_amb + Phi_b,0/(sigma^2 ln10) reproduces the N14 free-fit log rho_c. Also plain regressions on log Phi_b,0, log SBd, log Mb.
import numpy as np, pandas as pd, os
G=4.30091e-6; rho_amb=35.5
SPARC_DIR=os.path.join('..','sparc_galaxy_data','sparc_database'); MRT=os.path.join('..','sparc_galaxy_data','SPARC_Lelli2016c.mrt')
tab={}
for line in open(MRT):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: tab[p[0]]=dict(Rd=float(p[11]),SBd=float(p[12]),L=float(p[7]),MHI=float(p[13]),Q=int(p[17]))
        except ValueError: pass
N14=pd.read_csv('n14_sigma_free_20261001.csv',index_col=0); rows=[]
for g,row in N14.iterrows():
    f_=os.path.join(SPARC_DIR,g+'_rotmod.dat')
    if g not in tab or not os.path.exists(f_): continue
    d=np.atleast_2d(np.loadtxt(f_,comments='#')); r=d[:,0]; vg,vd,vb=d[:,3],d[:,4],d[:,5]
    vb2=np.clip(np.sign(vg)*vg**2+0.5*vd**2+0.7*vb**2,0,None)
    inner=vb2[0]/2                                  # solid body 0..r1: int (v1^2 r/r1^2)/r dr = v1^2/2
    mid=np.trapezoid(vb2/r,r)
    outer=vb2[-1]                                   # Keplerian beyond rmax
    Phi=inner+mid+outer                             # km^2/s^2, depth (positive)
    res=dict(gal=g,SBd=tab[g]['SBd'],Mb=row.Mb,lrhoc_free=row.lrhoc_free,lrhoc_btfr=row.lrhoc_btfr,s2_free=row.s2_free,s2_btfr=row.s2_btfr,
             lPhi=np.log10(Phi),Q=tab[g]['Q'],rmax=r[-1],vb2_last_frac=outer/Phi)
    for tag,s2 in (('free',row.s2_free),('btfr',row.s2_btfr)):
        res['lrc_pred_'+tag]=np.log10(rho_amb)+Phi/(s2*np.log(10))
    rows.append(res)
C=pd.DataFrame(rows).set_index('gal'); C.to_csv('n9b_rhoc_vs_phib_20261017.csv')
print('N =',len(C))
def reg(x,y,lab):
    m=np.isfinite(x)&np.isfinite(y); A=np.polyfit(x[m],y[m],1); res=y[m]-np.polyval(A,x[m])
    print('%-42s slope %+.3f  intercept %+.2f  scatter %.2f dex  r=%.2f'%(lab,A[0],A[1],res.std(),np.corrcoef(x[m],y[m])[0,1]))
print('--- regressions of log rho_c (free) ---')
reg(C.lPhi.values,C.lrhoc_free.values,'vs log Phi_b,0')
reg(np.log10(C.SBd).values,C.lrhoc_free.values,'vs log SBd')
reg(np.log10(C.Mb).values,C.lrhoc_free.values,'vs log Mb')
reg(np.log10(C.s2_free).values,C.lrhoc_free.values,'vs log sigma^2 (free)')
print('--- regressions of log rho_c (BTFR sigma) ---')
reg(C.lPhi.values,C.lrhoc_btfr.values,'vs log Phi_b,0'); reg(np.log10(C.SBd).values,C.lrhoc_btfr.values,'vs log SBd')
print('--- direct prediction log rho_c = log rho_amb + Phi_b,0/(sigma^2 ln10) ---')
for tag in ('free','btfr'):
    d=C['lrc_pred_'+tag]-C['lrhoc_'+tag]
    print('%s: pred-fit median %+.2f  16-84 %+.2f..%+.2f  std %.2f  | pred med %.2f, fit med %.2f'%(tag,d.median(),*np.percentile(d,[16,84]),d.std(),C['lrc_pred_'+tag].median(),C['lrhoc_'+tag].median()))
    reg(C['lrc_pred_'+tag].values,C['lrhoc_'+tag].values,'  fit vs pred (%s)'%tag)
print('raw spread log rho_c free: std %.2f, btfr: std %.2f'%(C.lrhoc_free.std(),C.lrhoc_btfr.std()))
print('Phi_b,0 range: log %.2f..%.2f ; outer Keplerian share median %.2f'%(C.lPhi.min(),C.lPhi.max(),C.vb2_last_frac.median()))
# Q=1 only
q=C[C.Q==1]; print('--- Q=1 only, N=%d ---'%len(q)); reg(q.lPhi.values,q.lrhoc_free.values,'vs log Phi_b,0 (free)'); reg(q['lrc_pred_free'].values,q.lrhoc_free.values,'fit vs pred (free)')
