# two-variable fit log mu0 = a + b logMb + c logSh, bootstrap errors. [UNRATIFIED]
import numpy as np, pandas as pd
D=pd.read_csv('fig4_mu0_mond_20261027.csv').set_index('gal')
B=pd.read_csv('test1c_burkert_data_20261018.csv',index_col=0); B=B[~B.r0_gt_rmax.astype(bool)]
B=B[B.index.isin(D.index)]
rng=np.random.default_rng(1)
def fit(x1,x2,y):
    X=np.c_[np.ones_like(x1),x1,x2]; c=np.linalg.lstsq(X,y,rcond=None)[0]
    bs=[]
    for _ in range(2000):
        i=rng.integers(0,len(y),len(y)); bs.append(np.linalg.lstsq(X[i],y[i],rcond=None)[0])
    return c,np.std(bs,axis=0)
for k,idx,y in (('data (Burkert)',B.index,B.mu_data.values),('this work',B.index,D.loc[B.index,'mu_thiswork'].values),('MOND',B.index,D.loc[B.index,'mu_mond'].values),('this work all120',D.index,D.mu_thiswork.values),('MOND all120',D.index,D.mu_mond.values)):
    c,e=fit(D.loc[idx,'lMb'].values,D.loc[idx,'lSh'].values,y)
    print('%-17s N=%3d  dlogMb %+.2f+-%.2f  dlogSh %+.2f+-%.2f'%(k,len(y),c[1],e[1],c[2],e[2]))
print('corr(lMb,lSh) =',np.corrcoef(D.lMb,D.lSh)[0,1].round(2))
