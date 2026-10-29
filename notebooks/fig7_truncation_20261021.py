# fig7_truncation_20261021.py [UNRATIFIED] -- Fig. 7: Mistele+2024 stacked circular velocities (Table 1, 'All', four mass bins) vs the law's truncated
# isothermal tail at rho_amb = 35.5 (closure), 355 (fit, after parser fix 2026-10-22), 630 (95% bound). LINEAR axes (JM 16:54). v normalised by the bin's fitted V_f at rho_amb=35.5.
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
G=4.30091e-6; from mistele_table import load as _mload
T=_mload('All')  # byte-position parser (fix 2026-10-22)
def vc_model(r,sig,ra):
    rt=sig/np.sqrt(2*np.pi*G*ra); rr=np.minimum(r,rt); M=2*sig**2*rr/G-4/3*np.pi*ra*rr**3; return np.sqrt(np.maximum(G*M/r,0)),rt
fig,ax=plt.subplots(2,2,figsize=(9,6.6),sharex=True,sharey=True); ax=ax.ravel()
rr=np.linspace(20,3500,2000)
for b in range(1,5):
    d=T[T['v%d'%b].notna()]; r=d.R.values; v=d['v%d'%b].values; e=np.sqrt(d['s%d'%b]**2+d['y%d'%b]**2).values
    # V_f: chi2 fit at rho_amb 35.5 over R>=100 kpc
    m=r>=100; best=(1e9,0)
    for vf in np.arange(80,350,0.5):
        c=np.sum(((v[m]-vc_model(r[m],vf/np.sqrt(2),35.5)[0])/e[m])**2)
        if c<best[0]: best=(c,vf)
    vf=best[1]; a=ax[b-1]
    a.errorbar(r/1000,v/vf,yerr=e/vf,fmt='o',ms=4,color='k',capsize=2,label='Mistele et al. 2024, bin %d'%b)
    for ra,ls,lab in ((35.5,'-','35.5 (closure)'),(355,'--','355 (fit)'),(630,':','630 (95% bound)')):
        vm,rt=vc_model(rr,vf/np.sqrt(2),ra); a.plot(rr/1000,vm/vf,ls,color='C0',label=r'$\rho_{\rm amb}$ = %s'%lab)
    a.axhline(1,color='0.7',lw=0.8); a.set_ylim(0,1.35); a.set_xlim(0,3.5)
    a.text(0.03,0.06,r'bin %d, $V_f$ = %.0f km/s'%(b,vf),transform=a.transAxes,fontsize=9)
    if b==1: a.legend(fontsize=8,loc='upper right',frameon=False)
for i in (2,3): ax[i].set_xlabel('R (Mpc)')
for i in (0,2): ax[i].set_ylabel(r'$v_c\,/\,V_f$')
fig.suptitle('Isothermal tail with ambient truncation vs stacked lensing circular velocities (linear axes)',fontsize=10)
fig.tight_layout(); fig.savefig('figs_v2/fig7_truncation.png',dpi=150); print('saved')
