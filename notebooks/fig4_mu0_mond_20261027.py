# fig4_mu0_mond_20261027.py -- Fig 4 redo: halo core surface density vs M_b, with MOND's per-galaxy value
# (Milgrom 2009, MNRAS 398, 1023, eqs 10-11, verified from arXiv:0909.5184 2026-10-27). [UNRATIFIED]
# MOND: Sigma_c* = Sigma_M {(6/pi) X0 [1-mu(X0)] + int_0^X0 L dx},  X0 mu(X0) = (2/3) Sigma_b/Sigma_M,
#   Sigma_M = a0/2piG, L = x mu'/mu, Sigma_b = rho R (mean density x extent).
# Our choice: R = R_h = 1.678 R_d (SPARC), rho = (M_b/2)/(4pi/3 R_h^3) -> Sigma_b = (3/4) Sigma_h, Sigma_h = M_b/2piR_h^2.
# mu from the RAR function (McGaugh+16), the same nu used elsewhere in the note. Sensitivity: Sigma_b x0.5, x2.
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from gate2c_robustness_20260925 import load_table
G=4.30091e-6
a0SI=1.2e-10; SigM=a0SI/(2*np.pi*6.674e-11)*(3.0857e16**2)/1.989e30   # Msun/pc^2
nu=lambda y: 1/(1-np.exp(-np.sqrt(y)))
yg=np.logspace(-8,4,200001); Xg=yg*nu(yg); mug=1/nu(yg)
lnmu=np.log(mug); lnX=np.log(Xg); Lg=np.gradient(lnmu,lnX)
cumL=np.concatenate([[Xg[0]*Lg[0]],Xg[0]*Lg[0]+np.cumsum(0.5*(Lg[1:]+Lg[:-1])*np.diff(Xg))])  # int_0^X L dx (L->1 as x->0)
def sigc(Sb):
    y=(2/3)*Sb/SigM; X0=y*nu(y); mu0=1/nu(y)
    return SigM*((6/np.pi)*X0*(1-mu0)+np.interp(X0,Xg,cumL))
lam=cumL[-1]; print('Sigma_M=%.1f (log %.3f); lambda(RAR mu)=%.3f gamma=%.3f'%(SigM,np.log10(SigM),lam,2*lam/np.pi))
F=pd.read_csv('test1_capped_baryon_20261018.csv'); T=load_table()
F['Rh']=1.678*T.loc[F.gal,'Rd'].values*1e3; F['Sh']=F.Mb/(2*np.pi*F.Rh**2); F['Sb']=0.75*F.Sh
lM=np.log10(F.Mb.values)
mu_law=np.log10(3*np.sqrt(F.s2b)*np.sqrt(10**F.l_tanh2/(4*np.pi*G))/1e6).values
mu_mond=np.log10(sigc(F.Sb.values))
B=pd.read_csv('test1c_burkert_data_20261018.csv',index_col=0); B=B[~B.r0_gt_rmax.astype(bool)]
def slope(x,y): A=np.polyfit(x,y,1); return A
out=[]
for k,v in (('this work',mu_law),('MOND per galaxy',mu_mond),('MOND Sb x0.5',np.log10(sigc(0.5*F.Sb.values))),('MOND Sb x2',np.log10(sigc(2*F.Sb.values)))):
    A=slope(lM,v); out.append((k,A[0],np.median(v),np.std(v)))
Ab=slope(np.log10(B.Mb),B.mu_data); out.append(('Burkert data',Ab[0],np.median(B.mu_data),np.std(B.mu_data)))
for o in out: print('%-18s slope %+.3f  median %.2f  sd %.2f'%o)
# per-galaxy agreement with the Burkert data (same galaxies)
m=F.gal.isin(B.index).values; d=B.loc[F.gal[m],'mu_data'].values
for k,v in (('this work',mu_law[m]),('MOND per galaxy',mu_mond[m])):
    r=v-d; print('%-16s vs data N=%d: median offset %+.2f, rms %.2f, Spearman(res,Sh)=%+.2f'%(k,m.sum(),np.median(r),np.sqrt(np.mean(r**2)),spearmanr(r,np.log10(F.Sh[m]))[0]))
# Sigma dependence at fixed mass: residual about each set's own mass trend vs log Sh
for k,x,v in (('this work',lM,mu_law),('MOND',lM,mu_mond),('data',np.log10(B.Mb.values),B.mu_data.values)):
    Sh=np.log10(F.Sh.values) if k!='data' else np.log10(F.set_index('gal').loc[B.index,'Sh'].values)
    r=v-np.polyval(slope(x,v),x); s=spearmanr(r,Sh); print('%-10s resid vs log Sh at fixed M_b: Spearman %+.2f p=%.1e'%(k,s[0],s[1]))
fig,ax=plt.subplots(figsize=(5.8,4.3)); xx=np.linspace(7,11.5,10)
ax.plot(np.log10(B.Mb),B.mu_data,'ks',ms=4,alpha=.55,label='SPARC Burkert fits (Li et al. 2020)')
ax.plot(lM,mu_law,'r.',ms=6,alpha=.75,label=r'this work, fitted $\rho_c$ (slope %.2f)'%out[0][1])
ax.plot(lM,mu_mond,'b^',ms=3.5,alpha=.6,mfc='none',label='MOND per galaxy, Milgrom 2009 (slope %.2f)'%out[1][1])
ax.plot(xx,np.polyval(slope(lM,mu_law),xx),'r-',lw=1); ax.plot(xx,np.polyval(slope(lM,mu_mond),xx),'b-',lw=1)
ax.plot(xx,np.polyval(Ab,xx),'k:',lw=1,label='data trend (slope %.2f)'%Ab[0])
ax.axhline(np.log10(SigM),color='b',ls='--',lw=.8,label=r'$\Sigma_M=a_0/2\pi G$'); ax.axhspan(1.95,2.35,color='gray',alpha=.12,label=r'Donato et al. 2009, $2.15\pm0.2$')
ax.set_xlabel(r'$\log M_b$ ($M_\odot$)'); ax.set_ylabel(r'$\log\,\mu_0$ ($M_\odot$ pc$^{-2}$)'); ax.set_ylim(1.4,4.3); ax.legend(fontsize=6.5,loc='upper left',framealpha=.9); fig.tight_layout(); fig.savefig('figs_v2/fig4_mu0.png',dpi=160); plt.close()
pd.DataFrame(dict(gal=F.gal,lMb=lM,lSh=np.log10(F.Sh),mu_thiswork=mu_law,mu_mond=mu_mond)).to_csv('fig4_mu0_mond_20261027.csv',index=False)
print('fig4 written')
