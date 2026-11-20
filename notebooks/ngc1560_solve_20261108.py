# ngc1560_solve_20261108.py -- PREREG_NGC1560_20261108 (STAMPED JM 2026-11-08 19:11). Solver = marg_fit_20261030.py model() verbatim
# (harness of v2.0.1), galaxy input replaced by data/ngc1560/ngc1560_gentile2010_rc_extracted_20261108.csv. Rule 3a: SNAP_20260922_001 sec 0c read 19:20.
# Runs R1 (rho_c 7.45), R2 (rho_c fitted), R3 (RAR function), R4 (M_* bracket). Criteria (a)-(d) computed here. [UNRATIFIED result]
import numpy as np, sys
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2)
UD=1.43; MHI=2.36e5*3.45**2*294.6
d=np.loadtxt('../data/ngc1560/ngc1560_gentile2010_rc_extracted_20261108.csv',delimiter=',')
r,vobs,lo,hi,vg,vd=d[:,0],d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]; ev=(hi-lo)/2
vb2=np.sign(vg)*vg**2+np.sign(vd)*vd**2            # already at M/L 1.43
vb2c=np.clip(vb2,0,None)
def relation(lrc,Mb):
    rl=r[-1]; rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2c[-1]*rl/x
    xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]
    s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    fn=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    so=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    assert so.success and so.y.shape[1]==len(r)
    return np.sqrt(np.clip(vb2c+G*so.y[2]/r,1,None)), np.sqrt(s2)
def rarfn():
    gb=vb2c/r; gm=gb/(1-np.exp(-np.sqrt(np.clip(gb/a0k,1e-12,None)))); return np.sqrt(gm*r)
def rms(v): return np.sqrt(np.mean(np.log10(vobs/v)**2))
def chi2(v): return np.sum(((v-vobs)/ev)**2)
def detrend(v,w):
    c=np.polyfit(np.log10(r[w]),np.log10(v[w]),2); return np.log10(v[w])-np.polyval(c,np.log10(r[w]))
w=(r>=2)&(r<=8); fw=(r[w]>=3)&(r[w]<=6.5)
def feature(v):
    f=detrend(v,w); return f[fw].max()-f[fw].min(), r[w][fw][f[fw].argmax()]
v3=rarfn()
print('R3 RAR fn : rms %.4f dex  chi2 %.1f'%(rms(v3),chi2(v3)))
out={}
for Ms in [1.3e9,1.0e9,1.6e9]:
    Mb=Ms+1.33*MHI
    v1,sig=relation(7.45,Mb)
    f=minimize_scalar(lambda l:chi2(relation(l,Mb)[0]),bounds=(6,10),method='bounded',options={'xatol':1e-3})
    v2,_=relation(f.x,Mb)
    tag='PRIMARY' if Ms==1.3e9 else 'bracket'
    print('M_*=%.1e M_b=%.2e sigma=%.1f km/s | R1 (7.45): rms %.4f chi2 %.1f ratio %.2f | R2: log rho_c %.2f rms %.4f chi2 %.1f ratio %.2f  [%s]'%(Ms,Mb,sig,rms(v1),chi2(v1),rms(v1)/rms(v3),f.x,rms(v2),chi2(v2),rms(v2)/rms(v3),tag))
    out[Ms]=(v1,v2,f.x)
v1,v2,lrc=out[1.3e9]
print('--- (c) feature 3-6.5 kpc, detrended 2-8 kpc: amplitude (dex), r_peak')
for nm,v in [('V_obs',vobs),('V_bar',np.sqrt(vb2c)),('R1',v1),('R2',v2),('R3',v3)]:
    A,rp=feature(v); print('  %-6s A=%.4f  r_peak=%.2f'%(nm,A,rp))
noise=np.median(ev[w][fw]/vobs[w][fw])/np.log(10); print('  noise (median sigma_V/V in window, dex) %.4f ; 2x = %.4f'%(noise,2*noise))
fo=detrend(vobs,w); print('  Pearson detrended R1 vs obs %.2f ; R3 vs obs %.2f ; V_bar vs obs %.2f'%(np.corrcoef(detrend(v1,w),fo)[0,1],np.corrcoef(detrend(v3,w),fo)[0,1],np.corrcoef(detrend(np.sqrt(vb2c),w),fo)[0,1]))
print('--- (d) outer slope dlogV/dlogr, last 5 pts: obs %.3f  R1 %.3f  R2 %.3f  R3 %.3f'%tuple(np.polyfit(np.log10(r[-5:]),np.log10(v[-5:]),1)[0] for v in [vobs,v1,v2,v3]))
np.savetxt('ngc1560_curves_20261108.csv',np.c_[r,vobs,ev,np.sqrt(vb2c),v1,v2,v3],fmt='%.3f',delimiter=',',header='r,vobs,ev,vbar,R1_745,R2_fit_lrc%.2f,R3_RAR'%lrc,comments='')
