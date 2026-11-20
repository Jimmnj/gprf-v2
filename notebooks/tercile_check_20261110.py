# tercile_check_20261110.py [UNRATIFIED] -- consistency check for B.1 (ledger (ck), BALLOT 2026-11-10 17:15)
# Sigma_h-tercile medians of the fitted log rho_c, raw and at fixed mass, both arms.
# Input: rhoc_corr_20261105_table.csv (fixed arm lrc = pass2_table2_part; marg arm lrc_marg = marg_fit_params arm A, MAP).
import numpy as np, pandas as pd
rng=np.random.default_rng(10)
D=pd.read_csv('rhoc_corr_20261105_table.csv').set_index('gal')
assert abs(D.loc['NGC2841','Rd']-3.64)<1e-6 and abs(D.loc['DDO154','Rd']-0.37)<1e-6
def ols(y,X):
    A=np.column_stack([np.ones(len(y))]+list(X)); b=np.linalg.lstsq(A,y,rcond=None)[0]; return b,y-A@b
def boot(d,cols,k,n=2000):
    out=[]
    for _ in range(n):
        s=d.iloc[rng.integers(0,len(d),len(d))]; out.append(ols(s[cols[0]].values,[s[c].values for c in cols[1:]])[0][k])
    return np.percentile(out,[2.5,97.5])
def terc(x,y,label):
    q=np.quantile(x,[1/3,2/3]); g=np.digitize(x,q); m=[np.median(y[g==i]) for i in range(3)]; n=[(g==i).sum() for i in range(3)]
    xm=[np.median(x[g==i]) for i in range(3)]
    bs=[]
    for _ in range(4000):
        i=rng.integers(0,len(x),len(x)); xx,yy=x[i],y[i]; gg=np.digitize(xx,np.quantile(xx,[1/3,2/3])); bs.append(np.median(yy[gg==2])-np.median(yy[gg==0]))
    lo,hi=np.percentile(bs,[2.5,97.5])
    print('  %-34s N %s  x-med %+.2f %+.2f %+.2f  y-med %+.3f %+.3f %+.3f  hi-lo %+.3f [%+.2f,%+.2f]'%(label,n,*xm,*m,m[2]-m[0],lo,hi))
for arm,col in [('FIXED','lrc'),('MARG (MAP)','lrc_marg')]:
    d=D[[col,'lMb','lSig','fgas']].dropna(); y=d[col].values
    print('\n=== %s  N=%d  median %.3f  sd %.3f'%(arm,len(d),np.median(y),y.std()))
    b2,r2=ols(y,[d.lMb.values,d.lSig.values]); c2=boot(d,[col,'lMb','lSig'],2)
    b3,r3=ols(y,[d.lMb.values,d.lSig.values,d.fgas.values]); c3=boot(d,[col,'lMb','lSig','fgas'],2)
    print('  OLS 2-var: lMb %+.3f lSig %+.3f [%+.2f,%+.2f]   resid sd %.3f'%(b2[1],b2[2],*c2,r2.std()))
    print('  OLS 3-var: lMb %+.3f lSig %+.3f [%+.2f,%+.2f] fgas %+.3f   resid sd %.3f'%(b3[1],b3[2],*c3,b3[3],r3.std()))
    print('  variance removed by lSig beyond lMb: %.1f%%'%(100*(1-r2.var()/ols(y,[d.lMb.values])[1].var())))
    terc(d.lSig.values,y,'raw Sigma_h terciles')
    ry=ols(y,[d.lMb.values])[1]; rx=ols(d.lSig.values,[d.lMb.values])[1]
    terc(rx,ry+np.median(y),'at fixed M_b (residual terciles)')
    ry=ols(y,[d.lMb.values,d.fgas.values])[1]; rx=ols(d.lSig.values,[d.lMb.values,d.fgas.values])[1]
    terc(rx,ry+np.median(y),'at fixed M_b and f_gas')
    print('  range of lSig: 16-84 %.2f dex; of lSig|M_b residual: %.2f dex'%(np.diff(np.percentile(d.lSig,[16,84]))[0],np.diff(np.percentile(ols(d.lSig.values,[d.lMb.values])[1],[16,84]))[0]))
