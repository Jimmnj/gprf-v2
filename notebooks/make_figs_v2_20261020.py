# make_figs_v2_20261020.py -- six figures for NOTE_V2_DRAFT from existing CSVs. [UNRATIFIED]
import sys, numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from gate2c_robustness_20260925 import load_table, kpc, a0
G=4.30091e-6; out='figs_v2/'
ONLY=set(sys.argv[1:])  # 2026-10-27: e.g. fig3 fig5; fig1/fig2 now from make_fig1_full_20261022.py, fig4 from fig4_mu0_mond_20261027.py
def _save(fig,path,**k):
    if any(t in path for t in ONLY): fig.savefig(path,**k); print('wrote',path)
nurar=lambda x: -np.log10(1-np.exp(-np.sqrt(10**x)))
# Fig 1
R=pd.read_csv('n12_emergent_nu_bins_20261020.csv'); R=R[R.N>=10].reset_index(drop=True); C=pd.read_csv('n12c_sform_bins_20261020.csv')
x=0.5*(R.lo+R.hi); base=nurar(x)
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
ax[0].plot(x,base+R.data_dRAR,'ks',ms=6,label='SPARC data (binned median)'); ax[0].plot(x,base,'k--',lw=1,label='empirical RAR function (McGaugh+16)')
for col,lab,c in (('u72_dRAR',r'this work, $\rho_c=10^{7.2}$','tab:blue'),('u745_dRAR',r'this work, $\rho_c=10^{7.45}$ (median)','tab:red'),('u78_dRAR',r'this work, $\rho_c=10^{7.8}$','tab:green')):
    ax[0].plot(x,base+R[col],'-',color=c,lw=2 if '745' in col else 1.2,label=lab)
ax[0].set_xlabel(r'$\log(g_{\rm bar}/a_0)$'); ax[0].set_ylabel(r'$\log(g_{\rm obs}/g_{\rm bar})$'); ax[0].legend(fontsize=8); ax[0].set_title('(a) emergent relation, no per-galaxy parameter')
xc=0.5*(C.lo+C.hi); b2=np.interp(xc,x,base+R.u745_dRAR)
for col,lab in (('tanh1','saturation 1'),('tanh1.5','1.5'),('tanh3','3'),('tanh4','4'),('rat2','rational, 2'),('lin','no saturation')):
    ax[1].plot(xc,b2+C[col],lw=1.2,label=lab)
ax[1].plot(xc,b2,'tab:red',lw=2,label='tanh, 2 (adopted)'); ax[1].plot(x,base+R.data_dRAR,'ks',ms=5,label='data')
ax[1].set_xlabel(r'$\log(g_{\rm bar}/a_0)$'); ax[1].set_ylabel(r'$\log(g_{\rm obs}/g_{\rm bar})$'); ax[1].legend(fontsize=7,ncol=2); ax[1].set_title(r'(b) saturation value and form, $\rho_c=10^{7.45}$')
fig.tight_layout(); _save(fig,out+'fig1_emergent_nu.png',dpi=160); plt.close()
# Fig 2
F=pd.read_csv('test1_capped_baryon_20261018.csv')
fig,ax=plt.subplots(figsize=(5.5,4)); bins=np.linspace(0,2.5,26)
ax.hist(np.log10(F.chi_lin),bins,alpha=.5,label=r'no saturation (median %.1f)'%F.chi_lin.median()); ax.hist(np.log10(F.chi_tanh2),bins,alpha=.6,label=r'saturated (median %.1f)'%F.chi_tanh2.median())
ax.hist(np.log10(F.chi_tf),bins,histtype='step',lw=1.5,label=r'saturated, $\sigma$ free (median %.1f)'%F.chi_tf.median())
ax.set_xlabel(r'$\log\,\chi^2$ per galaxy (10% velocity weights)'); ax.set_ylabel('galaxies'); ax.legend(fontsize=8); fig.tight_layout(); _save(fig,out+'fig2_chi.png',dpi=160); plt.close()
# Fig 3
D=pd.read_csv('n12_emergent_nu_20261020.csv'); D['V_obs']=np.sqrt(10**D.lgo*D.r*kpc/1e6); D['V_mod']=np.sqrt(10**D.lgm*D.r*kpc/1e6); D['V_bar']=np.sqrt(10**D.lgb*D.r*kpc/1e6)
T=load_table(); fg=(1.33*T.MHI*1e9/(0.5*T.L*1e9+1.33*T.MHI*1e9))
picks=[('DDO154','gas-rich dwarf'),('NGC2403','intermediate'),('NGC2841','massive, bulge')]
fig,ax=plt.subplots(1,3,figsize=(12,3.8))
for a,(g,lab) in zip(ax,picks):
    d=D[(D.gal==g)&(D.case=='fit')].sort_values('r'); u=D[(D.gal==g)&(D.case=='u7.45')].sort_values('r')
    a.plot(d.r,d.V_obs,'ko',ms=3,label='observed'); a.plot(d.r,d.V_bar,'c:',label='baryons'); a.plot(u.r,u.V_mod,'tab:blue',lw=1.2,label=r'this work, $\rho_c=10^{7.45}$'); a.plot(d.r,d.V_mod,'tab:red',lw=2,label=r'this work, $\rho_c$ fitted')
    a.set_title('%s (%s, $f_{gas}$=%.2f)'%(g,lab,fg[g])); a.set_xlabel('r (kpc)'); a.set_ylabel('V (km/s)'); a.set_ylim(0,None)
ax[0].legend(fontsize=8); fig.tight_layout(); _save(fig,out+'fig3_curves.png',dpi=160); plt.close()
# Fig 4
B=pd.read_csv('test1c_burkert_data_20261018.csv',index_col=0); B=B[~B.r0_gt_rmax.astype(bool)]
mu_law=np.log10(3*np.sqrt(F.s2b)*np.sqrt(10**F.l_tanh2/(4*np.pi*G))/1e6); lM=np.log10(F.Mb)
fig,ax=plt.subplots(figsize=(5.5,4))
ax.plot(np.log10(B.Mb),B.mu_data,'ks',ms=4,alpha=.6,label='SPARC Burkert fits (data)'); ax.plot(lM,mu_law,'r.',ms=5,alpha=.7,label=r'this work, from fitted $\rho_c$')
A=np.polyfit(lM,mu_law,1); xx=np.linspace(7,11.5,10); ax.plot(xx,np.polyval(A,xx),'r-',lw=1,label='this work trend (slope %.2f)'%A[0])
ax.axhline(2.14,color='b',ls='--',label=r'MOND $a_0/2\pi G$ (Milgrom 2009)'); ax.axhspan(1.95,2.35,color='b',alpha=.08,label='Donato+09 $2.15\pm0.2$')
ax.set_xlabel(r'$\log M_b$ ($M_\odot$)'); ax.set_ylabel(r'$\log\,\mu_0$ ($M_\odot$ pc$^{-2}$)'); ax.legend(fontsize=7); fig.tight_layout(); _save(fig,out+'fig4_mu0.png',dpi=160); plt.close()
# Fig 5
fig,ax=plt.subplots(figsize=(5.5,4)); s=np.array([7,14,28])
ax.plot([1,40],[1,1],'r-',lw=2,label=r'this work: $\gamma=1.000\,(+<0.003)$'); ax.plot(s,[1.281,1.354,1.361],'b-o',label='MOND + Galactic field (QUMOND, 1 $M_\odot$)')
ax.axhspan(0.97,1.03,color='g',alpha=.12,label='Cookson+26 / Banik+24 / Pittordis+25: Newtonian'); ax.errorbar([20],[1.43],yerr=0.06,fmt='m^',label='Chae 2023 ($g_N\lesssim10^{-10}$)')
ax.set_xscale('log'); ax.set_xlabel('separation (kAU)'); ax.set_ylabel(r'$\gamma$'); ax.set_ylim(0.9,1.6); ax.legend(fontsize=7); fig.tight_layout(); _save(fig,out+'fig5_wb.png',dpi=160); plt.close()
# Fig 6
W=pd.read_csv('dwarf_prereg_20261020.csv'); W=W[W.name!='Cetus_K14']
fig,ax=plt.subplots(figsize=(5.5,4.5))
for cls,c,lab in (('MW','tab:blue','MW satellites'),('M31','tab:green','M31 satellites'),('ISO','tab:red','isolated')):
    w=W[W.cls==cls]; ax.errorbar(w.s_obs,w.s_GPRF,xerr=[w.s_obs-w.s_obs*10**(-np.sqrt(np.clip(w.bar**2-0.01,0,None))),w.s_obs*10**(np.sqrt(np.clip(w.bar**2-0.01,0,None)))-w.s_obs],fmt='o',color=c,label=lab)
    for _,r in w.iterrows(): ax.annotate(r['name'],(r.s_obs,r.s_GPRF),fontsize=6,xytext=(3,3),textcoords='offset points')
xx=np.linspace(1,20,5); ax.plot(xx,xx,'k-',lw=1); ax.fill_between(xx,xx*10**-0.1,xx*10**0.1,color='k',alpha=.08,label=r'$\pm0.1$ dex')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'observed $\sigma$ (km/s)'); ax.set_ylabel(r'predicted $\sigma$, relation (3), no free parameter'); ax.legend(fontsize=8); fig.tight_layout(); _save(fig,out+'fig6_dwarfs.png',dpi=160); plt.close()
print('figs written')
