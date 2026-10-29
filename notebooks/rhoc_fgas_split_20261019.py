# rhoc_fgas_split_20261019.py [UNRATIFIED] -- follow-up to RESULT_RHOC_SIGMA_B_OOS: does the SPARC frozen-law rho_c residual (after M_b, Sigma_h) depend on gas fraction?
# DECLARED: if gas-poor SPARC galaxies have LOWER rho_c at fixed (M_b, Sigma_h), a gas lever could move DF2/DF4 (f_gas=0) in the right direction; needs ~-1.6 dex at f_gas 0.
import numpy as np, pandas as pd
from gate2c_robustness_20260925 import load_table
rng=np.random.default_rng(7)
T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
lM=np.log10(F.Mb.values); lS=np.log10(F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2)); y=F.l_tanh2.values
fg=1.33*F.MHI.values*1e9/F.Mb.values; F['fgas']=fg
X=np.column_stack([np.ones(len(y)),lM,lS]); b=np.linalg.lstsq(X,y,rcond=None)[0]; r=y-X@b
print('N=%d  f_gas range %.2f..%.2f median %.2f'%(len(y),fg.min(),fg.max(),np.median(fg)))
for lab,m in (('f_gas<0.2',fg<0.2),('0.2-0.5',(fg>=0.2)&(fg<0.5)),('>=0.5',fg>=0.5)):
    print('  %-10s N=%3d  mean resid %+.2f +/- %.2f dex  (mean log rho_c %.2f)'%(lab,m.sum(),r[m].mean(),r[m].std()/np.sqrt(m.sum()),y[m].mean()))
for cols,lab in ((['fgas'],'resid ~ f_gas'),(['fgas','lM','lS'],'log rho_c ~ f_gas + log Mb + log Sigma_h')):
    Z=np.column_stack([np.ones(len(y))]+[F[c].values if c=='fgas' else {'lM':lM,'lS':lS}[c] for c in cols]); yy=r if cols==['fgas'] else y
    bb=np.linalg.lstsq(Z,yy,rcond=None)[0]; bs=np.array([np.linalg.lstsq(Z[i],yy[i],rcond=None)[0] for i in (rng.integers(0,len(yy),len(yy)) for _ in range(3000))]); lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
    print('  %-42s f_gas coeff %+.2f [%+.2f,%+.2f] dex per unit f_gas; resid %.2f'%(lab,bb[1],lo[1],hi[1],(yy-Z@bb).std()))
print('  Implied shift from SPARC median f_gas to f_gas=0: see coeff x median. DF2 needs about -1.6 dex.')
