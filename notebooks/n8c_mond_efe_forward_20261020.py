# n8c_mond_efe_forward_20261020.py [UNRATIFIED] 2026-10-20 16:45 -- MOND-EFE FORWARD MODEL of N8b observable (A) (cold review A6; N8b 'owed').
# Chae+2021 Eq. (qumond1d) verbatim from data/chae2021/RCLSSr1.tex: g = nu(y+yNe) gN + [nu(y+yNe) - nu(yNe)] gNe, nu = nu0 (simple IF), y = gN/a0, yNe = e_N.
# For each N8b galaxy, at its observed points with g_bar < a0/3, compute the MOND-EFE prediction at its Chae e_N (max-clust) and the no-EFE prediction;
# observable (A)_model = mean[ log g_EFE - log g_RAR(McGaugh fn) ] -- same reference as the data residual. Regress (A)_model on log e_N (+ log Mb, log Sigma_b)
# with the SAME code path as N8b -> MOND's predicted slope, to compare with the observed -0.12 [-0.28,+0.04] and GPRF's 0.
# Note: simple-IF nu0 differs from the RAR (McGaugh) fn; the residual relative to the RAR fn is taken for both EFE and no-EFE and the DIFFERENCE is the EFE effect,
# so the IF mismatch cancels to first order. Also report the pure EFE dip: mean[log g_EFE - log g_noEFE].
import numpy as np, pandas as pd
from gate2c_robustness_20260925 import load, a0
P=load(); rng=np.random.default_rng(11)
env=pd.read_csv('../data/chae2021/chae2021_table3_eNenv.csv').set_index('galaxy')
N8=pd.read_csv('n8b_chae_observable_20261019.csv',index_col=0)
nu0=lambda y: 0.5+np.sqrt(0.25+1/y)
nurar=lambda y: 1/(1-np.exp(-np.sqrt(y)))
rows=[]
for g,d in P.groupby('gal'):
    if g not in N8.index or g not in env.index: continue
    low=d[d.gbar<a0/3]
    if len(low)<3: continue
    y=low.gbar.values/a0; eN=10**env.loc[g,'log_eN_env_maxclust']
    g_efe=nu0(y+eN)*y+(nu0(y+eN)-nu0(eN))*eN       # in units of a0
    g_no=nu0(y)*y; g_rar=nurar(y)*y
    rows.append(dict(gal=g,dip=np.mean(np.log10(g_efe/g_no)),resA_model=np.mean(np.log10(g_efe/g_rar)),resA_noefe=np.mean(np.log10(g_no/g_rar)),
                     log_eN=np.log10(eN),lMb=N8.loc[g,'lMb'],lSig=N8.loc[g,'lSig'],resA_data=N8.loc[g,'resA']))
D=pd.DataFrame(rows).set_index('gal'); D.to_csv('n8c_mond_efe_forward_20261020.csv')
def reg(y,cols,lab):
    m=np.isfinite(y)
    for c in cols: m&=np.isfinite(D[c].values)
    X=np.column_stack([np.ones(m.sum())]+[D[c].values[m] for c in cols]); yy=y[m]
    b=np.linalg.lstsq(X,yy,rcond=None)[0]
    bs=np.array([np.linalg.lstsq(X[i],yy[i],rcond=None)[0] for i in (rng.integers(0,len(yy),len(yy)) for _ in range(3000))]); lo,hi=np.percentile(bs,[2.5,97.5],axis=0)
    print('  %-44s N=%3d  log_eN coeff %+.3f [%+.3f,%+.3f]'%(lab,m.sum(),b[1],lo[1],hi[1]))
print('N=%d galaxies. EFE dip (log g_EFE - log g_noEFE at g_bar<a0/3): median %.3f, range %.3f..%.3f; e_N median %.4f'%(len(D),D.dip.median(),D.dip.min(),D.dip.max(),10**D.log_eN.median()))
print('MOND-EFE PREDICTED slope of observable (A) vs log e_N:')
reg(D.resA_model.values,['log_eN'],'EFE model, e_N alone')
reg(D.resA_model.values,['log_eN','lMb','lSig'],'EFE model, e_N + log Mb + log Sigma_b')
reg(D.dip.values,['log_eN','lMb','lSig'],'pure dip (EFE - noEFE), + Mb + Sigma_b')
print('Reference: no-EFE model (should be ~0 slope; nonzero = IF-mismatch leakage):'); reg(D.resA_noefe.values,['log_eN','lMb','lSig'],'no-EFE simple IF vs RAR fn')
print('DATA (N8b): e_N + log Mb + log Sigma_b slope -0.120 [-0.279,+0.044]'); reg(D.resA_data.values,['log_eN','lMb','lSig'],'data, this subsample')
print('Separation: (MOND slope - 0) vs data CI half-width 0.16.')
