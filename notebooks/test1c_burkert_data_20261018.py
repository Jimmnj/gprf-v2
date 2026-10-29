# test1c_burkert_data_20261018.py [UNRATIFIED] CONTROL for test1b: Burkert halo fit to the SPARC DATA itself (same 120 galaxies, Yd 0.5/Yb 0.7, same
# cuts, 10% velocity weights as N14): v_obs^2 = v_bar^2 + G M_B(<r)/r. Gives the observed log(rho0 r0) and its mass slope in OUR sample/method.
import numpy as np, pandas as pd
from scipy.optimize import least_squares
import test1_capped_baryon_20261018 as T
from test1b_burkert_surface_20261018 import MB
from gate2c_robustness_20260925 import G
C=pd.read_csv('test1_capped_baryon_20261018.csv',index_col=0); B=pd.read_csv('test1b_burkert_surface_20261018.csv',index_col=0); rows=[]
for g in C.index:
    d=T.P[T.P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); vo=d.vobs.values
    best=None
    for l0 in (-1,0,1):
        for p1 in (6,7.5,9):
            res=least_squares(lambda p:(np.sqrt(vb2+G*MB(r,10**p[0],10**p[1])/r)-vo)/(0.1*vo),[l0,p1],bounds=([-2,3],[3,11]))
            if best is None or res.cost<best.cost: best=res
    rows.append(dict(gal=g,mu_data=best.x[0]+best.x[1]-6,lr0_data=best.x[0],chi_burk=2*best.cost,rmax=r[-1],r0_gt_rmax=10**best.x[0]>r[-1]))
R=pd.DataFrame(rows).set_index('gal').join(C[['Mb','Q','chi_tanh2']]).join(B[['mu_V_tanh2']]); R.to_csv('test1c_burkert_data_20261018.csv')
ok=~R.r0_gt_rmax
for sub,lab in ((R,'all'),(R[ok],'r0<rmax'),(R[ok&(R.Q==1)],'r0<rmax Q=1')):
    for col,nm in (('mu_data','DATA'),('mu_V_tanh2','MODEL tanh2')):
        A=np.polyfit(np.log10(sub.Mb),sub[col],1); res=sub[col]-np.polyval(A,np.log10(sub.Mb))
        print('%-12s %-12s N=%3d median %.2f std %.2f slope vs log Mb %+.2f (scatter about trend %.2f)'%(lab,nm,len(sub),sub[col].median(),sub[col].std(),A[0],res.std()))
    dd=sub.mu_V_tanh2-sub.mu_data; print('%-12s MODEL-DATA per galaxy: median %+.2f std %.2f r=%.2f'%(lab,dd.median(),dd.std(),np.corrcoef(sub.mu_V_tanh2,sub.mu_data)[0,1]))
print('median chi: Burkert(2 par) %.1f vs GPRF tanh2 BTFR (1 par) %.1f ; log Mb range %.1f..%.1f'%(R.chi_burk.median(),R.chi_tanh2.median(),np.log10(R.Mb).min(),np.log10(R.Mb).max()))
