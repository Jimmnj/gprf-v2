# ngc1560_curves_supp_20261111.py -- PREREG_NGC1560_20261108 sec 13 (declared 2026-11-11 15:00 before the run). [UNRATIFIED; supplementary; display only]
# = ngc1560_dist_supp_20261110.py with ONE line added (np.savetxt of the curves, marked ADDED). Solver text untouched. Generated mechanically from that file.
# ngc1560_dist_supp_20261110.py -- SUPPLEMENTARY distance sensitivity (PREREG_NGC1560_20261108 sec 11, declared 2026-11-10 19:18 before the run).
# Solver = ngc1560_solve_20261108.py relation()/rarfn() VERBATIM (marg_fit_20261030 model), with the inputs rescaled for distance.
# 19:35 extension (PREREG sec 12): argument D:ML sets the I-band M/L (default 1.43); stellar V^2 and M_* scale with ML/1.43.
# Scaling from D0 = 3.45: r -> r f; Vgas^2, Vdisk^2 -> x f; V_obs, errors unchanged; M_HI, M_* -> x f^2 (f = D/D0). Rule 3a: SNAP_20260922_001 sec 0c read 19:18.
import numpy as np, sys, time
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2)
D0=3.45; MHI0=2.36e5*D0**2*294.6; MS0=1.3e9
d=np.loadtxt('../data/ngc1560/ngc1560_gentile2010_rc_extracted_20261108.csv',delimiter=',')
r0,vobs,lo,hi,vg,vd=d[:,0],d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]; ev=(hi-lo)/2
vb20=np.sign(vg)*vg**2+np.sign(vd)*vd**2
def run(D,ML=1.43):
    f=D/D0; m=ML/1.43; r=r0*f; vb2c=np.clip((np.sign(vg)*vg**2+m*np.sign(vd)*vd**2)*f,0,None); Mb=(MS0*m+1.33*MHI0)*f*f
    def relation(lrc):
        rl=r[-1]; rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2c])
        vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2c[-1]*rl/x
        xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2c[-1]
        s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
        fn=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
        so=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
        assert so.success and so.y.shape[1]==len(r)
        return np.sqrt(np.clip(vb2c+G*so.y[2]/r,1,None)), np.sqrt(s2)
    gb=vb2c/r; v3=np.sqrt(gb/(1-np.exp(-np.sqrt(np.clip(gb/a0k,1e-12,None))))*r)
    rms=lambda v: np.sqrt(np.mean(np.log10(vobs/v)**2)); chi2=lambda v: np.sum(((v-vobs)/ev)**2)
    v1,sig=relation(7.45)
    ft=minimize_scalar(lambda l:chi2(relation(l)[0]),bounds=(6,10),method='bounded',options={'xatol':1e-3}); v2,_=relation(ft.x)
    o=r>2*f   # the points beyond 2 kpc at D0
    med=lambda v: np.median((v[o]-vobs[o]))
    sl=lambda v: np.polyfit(np.log10(r[-5:]),np.log10(v[-5:]),1)[0]
    np.savetxt('ngc1560_curves_supp_D%.2f_ML%.2f_20261111.csv'%(D,ML),np.c_[r,vobs,ev,np.sqrt(vb2c),v1,v2,v3],fmt='%.3f',delimiter=',',header='r,vobs,ev,vbar,R1_745,R2_fit_lrc%.2f,R3_RAR'%ft.x,comments='')   # ADDED 2026-11-11: the only change
    print('M/L %.2f '%ML,end=''); print('D %.2f  M_b %.2e sigma %.1f | R3 rms %.4f chi2 %.0f | R1 rms %.4f chi2 %.0f ratio %.2f | R2 lrc %.2f rms %.4f chi2 %.0f | median(model-obs) outer: R1 %+.1f R3 %+.1f km/s | outer slope obs %.2f R1 %.2f R3 %.2f | Vbar/Vobs last %.2f'%(D,Mb,sig,rms(v3),chi2(v3),rms(v1),chi2(v1),rms(v1)/rms(v3),ft.x,rms(v2),chi2(v2),med(v1),med(v3),sl(vobs),sl(v1),sl(v3),np.sqrt(vb2c[-1])/vobs[-1]),flush=True)
for a in sys.argv[1:]:
    D,ML=(a.split(':')+['1.43'])[:2]; t=time.time(); run(float(D),float(ML)); print('   (%.0f s)'%(time.time()-t),flush=True)
