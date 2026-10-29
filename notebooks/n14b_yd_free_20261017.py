# n14b_yd_free_20261017.py  [UNRATIFIED]  N14b: N14 free fit with Yd also free (0.3-0.8), Yb fixed 0.7. For N9b residual test.
# Same chi as n14 (10% velocity weights, canon PB ODE). Start from the N14 free (log s2, log rho_c) point and Yd 0.5; one Nelder-Mead.
# Run from notebooks/: python3 n14b_yd_free_20261017.py  (Pool of 2; appends rows to n14b_part.csv as they finish; final n14b_yd_free_20261017.csv)
import numpy as np, pandas as pd, os, sys, glob
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from multiprocessing import Pool
G=4.30091e-6; rho_amb=35.5; Yb=0.7
SPARC_DIR=os.path.join('..','sparc_galaxy_data','sparc_database')
N14=pd.read_csv('n14_sigma_free_20261001.csv',index_col=0)
def fit(g):
    d=np.atleast_2d(np.loadtxt(os.path.join(SPARC_DIR,g+'_rotmod.dat'),comments='#'))
    r=d[:,0]; vobs=d[:,1]; vg,vd,vb=d[:,3],d[:,4],d[:,5]
    w=1/(0.1*vobs)**2; row=N14.loc[g]
    def chi(p):
        ls,l,Yd=p; Yd=np.clip(Yd,0.3,0.8); l=np.clip(l,4,11); s2=10**ls; rc=10**l
        vb2=np.clip(np.sign(vg)*vg**2+Yd*vd**2+Yb*vb**2,0,None)
        rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2])
        vbar=lambda x: np.interp(x,rr,vv) if x<=r[-1] else vb2[-1]*r[-1]/x
        f=lambda x,y:[(vbar(max(x,1e-4))+G*y[1]/max(x,1e-4))/max(x,1e-4),4*np.pi*x*x*(rc*np.exp(-y[0]/s2)-rho_amb)]
        s=solve_ivp(f,[1e-4,r[-1]],[0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=r[-1]/200)
        m=vb2+G*s.y[1]/r
        return np.sum(w*(np.sqrt(np.clip(m,1,None))-vobs)**2)
    p0=[np.log10(row.s2_free),row.lrhoc_free,0.5]
    res=minimize(chi,p0,method='Nelder-Mead',options={'xatol':0.01,'fatol':0.1,'maxfev':400})
    ls,l,Yd=res.x; out=dict(gal=g,s2_yd=10**ls,lrhoc_yd=np.clip(l,4,11),Yd=float(np.clip(Yd,0.3,0.8)),chi_yd=res.fun,chi_free=row.chi_free,nfev=res.nfev)
    with open('n14b_part.csv','a') as fh: fh.write(','.join(str(out[k]) for k in ('gal','s2_yd','lrhoc_yd','Yd','chi_yd','chi_free','nfev'))+'\n')
    return out
if __name__=='__main__':
    done=set()
    if os.path.exists('n14b_part.csv'): done=set(pd.read_csv('n14b_part.csv',header=None)[0])
    gals=[g for g in N14.index if g not in done and os.path.exists(os.path.join(SPARC_DIR,g+'_rotmod.dat'))]
    with Pool(2) as P:
        for o in P.imap_unordered(fit,gals): print(o['gal'],'%.2f %.2f %.1f %.1f'%(o['lrhoc_yd'],o['Yd'],o['chi_yd'],o['chi_free']),flush=True)
    C=pd.read_csv('n14b_part.csv',header=None,names=['gal','s2_yd','lrhoc_yd','Yd','chi_yd','chi_free','nfev']).set_index('gal').sort_index()
    C.to_csv('n14b_yd_free_20261017.csv'); print('DONE',len(C))
