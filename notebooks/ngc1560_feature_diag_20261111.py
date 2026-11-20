# NGC 1560 feature DIAGNOSTIC, 2026-11-11 (Claude). POST-HOC, NOT pre-registered, NOT scored. [UNRATIFIED]
# Input: notebooks/ngc1560_curves_20261108.csv (ratified run outputs). No solver is run; nothing touches the constitutive law.
# Run 2026-11-11 between 14:23 and 14:26 EST (clock readings either side) in the cloud workspace on a transcribed copy of the csv (validated by step 1, which reproduces
# the pre-registered numbers). The only change from the run copy is the input file name below.
# EXPECTATIONS written before the run (in the script header, same window; the "14:26" first typed here was a guess):
# E1 prereg (c) numbers reproduce: A obs 0.034 @4.64; Vbar 0.064 @3.01; R1 0.012; R2 0.018; R3 0.036; Pearson 0.78/0.79/0.79.
# E2 3.011 kpc is the FIRST point of the 3-6.5 window: the model/baryon "peak" is a window edge.
# E3 regression transmission (detrended log V on detrended log Vbar, 2-8 kpc): R3 0.5-0.6; R1 0.15-0.30; R2 0.20-0.35;
#    observed 0.3-0.6, error 0.10-0.15 (nearer R3 than R1, not decisively).
# E4 local dip 4.64 -> 5.0-5.4 kpc: observed larger than both models; R1 shows none.
import numpy as np
d=np.genfromtxt('ngc1560_curves_20261108.csv',delimiter=',',names=True)
r=d['r']; cols={'obs':d['vobs'],'bar':d['vbar'],'R1':d['R1_745'],'R2':d['R2_fit_lrc710'],'R3':d['R3_RAR']}
ev=d['ev']
def detr(v,lo=2.0,hi=8.0,deg=2):
    m=(r>=lo)&(r<=hi); x=np.log10(r[m]); y=np.log10(v[m]); p=np.polyfit(x,y,deg); return m, y-np.polyval(p,x)
print("1. prereg (c) recipe: quadratic in log r over 2-8 kpc; amplitude = peak-to-trough in 3-6.5 kpc")
D={}
for k,v in cols.items():
    m,res=detr(v); rr=r[m]; w=(rr>=3.0)&(rr<=6.5); D[k]=(rr,res)
    ip=np.argmax(res[w]); it=np.argmin(res[w])
    print("   %-3s A=%.4f dex  peak at %.3f kpc (window index %d of %d)  trough at %.3f kpc"%(k,res[w].max()-res[w].min(),rr[w][ip],ip,w.sum(),rr[w][it]))
rr,ro=D['obs']
for k in ('R1','R3','bar'):
    print("   Pearson detrended %s vs obs (2-8 kpc): %.2f"%(k,np.corrcoef(D[k][1],ro)[0,1]))
m=(r>=3)&(r<=6.5); print("   noise 2 x median(ev/vobs) in window = %.4f (fraction) = %.4f dex"%(2*np.median(ev[m]/cols['obs'][m]), 2*np.median(ev[m]/cols['obs'][m])/np.log(10)))
print("2. transmission = OLS slope of detrended log V on detrended log Vbar, 2-8 kpc (same detrend)")
xb=D['bar'][1]
def slope(y,x): return np.sum(x*y)/np.sum(x*x)
for k in ('obs','R1','R2','R3'):
    y=D[k][1]; s=slope(y,xb)
    if k=='obs':
        mm=(r>=2)&(r<=8); sig=ev[mm]/cols['obs'][mm]/np.log(10)          # per-point error in dex
        rng=np.random.default_rng(1); S=[]
        for _ in range(4000):                                            # perturb V_obs by its errors, re-detrend
            vv=cols['obs'].copy(); vv[mm]=vv[mm]*10**(rng.normal(0,sig)); _,rs=detr(vv); S.append(slope(rs,xb))
        jk=[slope(np.delete(y,i),np.delete(xb,i)) for i in range(len(y))]; n=len(y)
        print("   %-3s %.3f  (error-perturbation sd %.3f; jackknife se %.3f; n=%d)"%(k,s,np.std(S),np.sqrt((n-1)/n*np.sum((np.array(jk)-np.mean(jk))**2)),n))
    else: print("   %-3s %.3f"%(k,s))
print("3. raw read-offs (no detrend), km/s")
def at(v,x): return v[np.argmin(abs(r-x))]
def mean(v,a,b): mm=(r>=a)&(r<=b); return v[mm].mean()
for k,v in cols.items():
    print("   %-3s V(4.64)=%.2f  mean V(5.0-5.4)=%.2f  dip=%.2f | slope 4.0-6.0: %.2f km/s/kpc  slope 6.0-6.9: %.2f  ratio %.1f"%(
        k,at(v,4.641),mean(v,5.0,5.4),mean(v,5.0,5.4)-at(v,4.641),
        np.polyfit(r[(r>=4.0)&(r<=6.03)],v[(r>=4.0)&(r<=6.03)],1)[0], np.polyfit(r[(r>=6.0)&(r<=6.91)],v[(r>=6.0)&(r<=6.91)],1)[0],
        np.polyfit(r[(r>=6.0)&(r<=6.91)],v[(r>=6.0)&(r<=6.91)],1)[0]/np.polyfit(r[(r>=4.0)&(r<=6.03)],v[(r>=4.0)&(r<=6.03)],1)[0]))
mm=(r>=4.0)&(r<=6.03); p=np.polyfit(r[mm],cols['obs'][mm],1,w=1/ev[mm],cov=True); print("   obs slope 4.0-6.0 weighted: %.2f +- %.2f km/s/kpc"%(p[0][0],np.sqrt(p[1][0,0])))
mm=(r>=6.0)&(r<=6.91); p=np.polyfit(r[mm],cols['obs'][mm],1,w=1/ev[mm],cov=True); print("   obs slope 6.0-6.9 weighted: %.2f +- %.2f km/s/kpc"%(p[0][0],np.sqrt(p[1][0,0])))
print("4. window sensitivity of step 1 (peak-to-trough), lower edge moved")
for lo in (2.5,3.0,3.5,4.0):
    out=[]
    for k in ('obs','bar','R1','R3'):
        rr,res=D[k]; w=(rr>=lo)&(rr<=6.5); out.append("%s %.3f@%.2f"%(k,res[w].max()-res[w].min(),rr[w][np.argmax(res[w])]))
    print("   window %.1f-6.5: "%lo+"  ".join(out))
# OUTPUT 2026-11-11, run between 14:23 and 14:26 EST:
# 1. prereg (c) recipe: quadratic in log r over 2-8 kpc; amplitude = peak-to-trough in 3-6.5 kpc
#    obs A=0.0344 dex  peak at 4.641 kpc (window index 13 of 28)  trough at 5.018 kpc
#    bar A=0.0635 dex  peak at 3.011 kpc (window index 0 of 28)  trough at 5.394 kpc
#    R1  A=0.0124 dex  peak at 3.011 kpc (window index 0 of 28)  trough at 5.394 kpc
#    R2  A=0.0175 dex  peak at 3.011 kpc (window index 0 of 28)  trough at 5.269 kpc
#    R3  A=0.0356 dex  peak at 3.011 kpc (window index 0 of 28)  trough at 5.394 kpc
#    Pearson detrended R1 vs obs (2-8 kpc): 0.78
#    Pearson detrended R3 vs obs (2-8 kpc): 0.79
#    Pearson detrended bar vs obs (2-8 kpc): 0.79
#    noise 2 x median(ev/vobs) in window = 0.0830 (fraction) = 0.0361 dex
# 2. transmission = OLS slope of detrended log V on detrended log Vbar, 2-8 kpc (same detrend)
#    obs 0.354  (error-perturbation sd 0.139; jackknife se 0.036; n=45)
#    R1  0.188
#    R2  0.256
#    R3  0.560
# 3. raw read-offs (no detrend), km/s
#    obs V(4.64)=63.51  mean V(5.0-5.4)=61.45  dip=-2.05 | slope 4.0-6.0: 1.96 km/s/kpc  slope 6.0-6.9: 6.05  ratio 3.1
#    bar V(4.64)=32.85  mean V(5.0-5.4)=30.71  dip=-2.14 | slope 4.0-6.0: -0.64 km/s/kpc  slope 6.0-6.9: 5.38  ratio -8.4
#    R1  V(4.64)=73.56  mean V(5.0-5.4)=75.30  dip=1.75 | slope 4.0-6.0: 4.36 km/s/kpc  slope 6.0-6.9: 4.96  ratio 1.1
#    R2  V(4.64)=61.16  mean V(5.0-5.4)=63.06  dip=1.90 | slope 4.0-6.0: 4.90 km/s/kpc  slope 6.0-6.9: 6.44  ratio 1.3
#    R3  V(4.64)=69.78  mean V(5.0-5.4)=68.94  dip=-0.84 | slope 4.0-6.0: 2.34 km/s/kpc  slope 6.0-6.9: 9.32  ratio 4.0
#    obs slope 4.0-6.0 weighted: 1.89 +- 0.48 km/s/kpc
#    obs slope 6.0-6.9 weighted: 5.96 +- 0.84 km/s/kpc
# 4. window sensitivity of step 1 (peak-to-trough), lower edge moved
#    window 2.5-6.5: obs 0.034@4.64  bar 0.064@2.88  R1 0.013@2.88  R3 0.036@2.88
#    window 3.0-6.5: obs 0.034@4.64  bar 0.063@3.01  R1 0.012@3.01  R3 0.036@3.01
#    window 3.5-6.5: obs 0.034@4.64  bar 0.060@3.76  R1 0.011@3.76  R3 0.033@3.76
#    window 4.0-6.5: obs 0.034@4.64  bar 0.054@4.01  R1 0.010@4.01  R3 0.030@4.01
