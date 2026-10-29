# make_fig1_full_20261022.py -- Fig 1 remade on the full Lelli-style sample (147 galaxies, 2696 points): (a) universal 10^7.45, relation (3), fitted rho_c,
# vs binned data and the McGaugh+16 function; (b) residuals log(g_obs/g_model) vs log Sigma_h and log r per point (universal and relation 3), binned medians.
import numpy as np, pandas as pd, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
a0=1.2e-10; kpc=3.0857e19; la0=np.log10(a0)
A=pd.read_csv('pass3_relation3_points_20261022.csv')
R2=pd.read_csv('pass2_table2_part.csv').drop_duplicates('gal')
# per-galaxy fitted residuals, aligned to univ points by galaxy order
U=A[A.case=='univ'].copy(); fitres=[]
for g,d in U.groupby('gal',sort=False):
    r=json.loads(R2[R2.gal==g].law.iloc[0]); assert len(r)==len(d),(g,len(r),len(d)); fitres+=list(r)
U['res_fit']=fitres
edges=np.arange(-2.25,1.01,0.25)+la0; x=0.5*(edges[:-1]+edges[1:])-la0
def binned(D,col):
    return np.array([D[col][(D.lgb>=lo)&(D.lgb<hi)].median() if ((D.lgb>=lo)&(D.lgb<hi)).sum()>=10 else np.nan for lo,hi in zip(edges[:-1],edges[1:])])
nurar=lambda l: -np.log10(1-np.exp(-np.sqrt(10**l)))
base=nurar(x)
R3=A[A.case=='rel3']
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
dat=base+binned(U,'resR'); ax[0].plot(x,dat,'ks',ms=6,label='SPARC data (binned median), 147 galaxies')
ax[0].plot(x,base,'k--',lw=1,label='empirical function (McGaugh+16)')
ax[0].plot(x,dat-binned(U,'res'),'-',color='tab:red',lw=2,label=r'this work, universal $\rho_c=10^{7.45}$')
ax[0].plot(x,dat-binned(R3,'res'),'-',color='tab:blue',lw=1.6,label=r'this work, $\rho_c$ from relation (3)')
ax[0].plot(x,dat-binned(U,'res_fit'),'-',color='tab:green',lw=1.2,label=r'this work, $\rho_c$ fitted per galaxy')
ax[0].set_xlabel(r'$\log(g_{\rm bar}/a_0)$'); ax[0].set_ylabel(r'$\log(g_{\rm obs}/g_{\rm bar})$'); ax[0].legend(fontsize=8); ax[0].set_title('(a) the relation, full sample')
# (b) residual vs Sigma_h
for D,c,lab in ((U,'tab:red','universal'),(R3,'tab:blue','relation (3)')):
    q=pd.qcut(np.log10(D.Sh),8); m=D.groupby(q).res.median(); xs=[iv.mid for iv in m.index]
    ax[1].plot(xs,m.values,'o-',color=c,label='this work, '+lab)
q=pd.qcut(np.log10(U.Sh),8); m=U.groupby(q).resR.median(); ax[1].plot([iv.mid for iv in m.index],m.values,'ks--',ms=5,label='empirical function')
ax[1].axhline(0,color='0.7',lw=0.8); ax[1].set_xlabel(r'$\log\,\Sigma_h$ ($M_\odot$ pc$^{-2}$)'); ax[1].set_ylabel(r'median $\log(g_{\rm obs}/g_{\rm model})$'); ax[1].legend(fontsize=8); ax[1].set_title('(b) residual against surface density')
fig.tight_layout(); fig.savefig('figs_v2/fig1_emergent_nu.png',dpi=160); print('fig1 saved')
# Fig 2 on the full sample: chi2 per galaxy universal vs fitted (10% weights) -- from pass2 csv
U['w']=1/(0.1*U.vobs)**2; U['vm']=np.sqrt(U.vmod2)
chi_u=U.groupby('gal').apply(lambda d: np.sum(d.w*(d.vm-d.vobs)**2))
chi_f=R2.set_index('gal').chi_law.reindex(chi_u.index)
fig,ax=plt.subplots(figsize=(5.5,4)); bins=np.linspace(0,2.5,26)
ax.hist(np.log10(chi_u.clip(0.1)),bins,alpha=.5,label=r'universal $\rho_c$ (median %.1f)'%chi_u.median()); ax.hist(np.log10(chi_f.clip(0.1)),bins,alpha=.6,label=r'$\rho_c$ fitted (median %.1f)'%chi_f.median())
ax.set_xlabel(r'$\log\,\chi^2$ per galaxy (10% velocity weights), 147 galaxies'); ax.set_ylabel('galaxies'); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig('figs_v2/fig2_chi.png',dpi=160); print('fig2 saved; chi2 medians univ %.1f fitted %.1f'%(chi_u.median(),chi_f.median()))
