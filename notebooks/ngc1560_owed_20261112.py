# ngc1560_owed_20261112.py -- owed items (a) and (c) of SESSION_HANDOFF_20261111 sec 3.4 / SESSION_HANDOFF_20261112 queue item 2. [UNRATIFIED]
# POST-PROCESSING ONLY. Inputs: the saved curve csvs. No solver is run; nothing touches the constitutive law (rule 3a not triggered).
# Recipe = ngc1560_feature_diag_20261111.py (H4 transmission, H5 shelf/step slopes, peak-to-trough amplitude), applied to the CORRECTED curves
# (ngc1560_corr_*_20261111.csv, ERRATUM_NGC1560_COLUMNS_20261111). Window and segments chosen on the 3.45-Mpc radius scale r0 (as in
# ngc1560_corrected_20261111.py); slopes in km/s per kpc of the run's own radius scale (as in the plot report).
# EXPECTATIONS, written 2026-11-12 01:06 EST (clock read 01:06:18) BEFORE the run:
# X1 (validation) old csv ngc1560_curves_20261108.csv reproduces H4: obs 0.354 (sd 0.139, jk 0.036), R1 0.188, R2 0.256, R3 0.560; H5 obs 1.89 +- 0.48, R3 2.34, R1 4.36, R2 4.90.
# X2 (validation) corrected (3.45, 1.43, M_* 7.0e8): obs and R3 transmissions EQUAL the old values (V_bar unchanged by the swap at M/L 1.43);
#    shelf slopes reproduce the plot report: R1 3.60, R2 4.31, R3 2.34, obs 1.89 +- 0.48.
# X3 corrected R1 universal transmission 0.19-0.24 (old 0.188; amplitude rose 0.0124 -> 0.014); R2 fitted 0.22-0.28.
# X4 (3.45, 1.00) and (2.94, 0.98): lighter stars -> larger V_bar feature -> obs transmission FALLS to 0.25-0.33; R3 0.50-0.60; R1 0.15-0.21.
#    Obs error: jackknife 0.03-0.04, perturbation 0.12-0.17.
# X5 (item c) R1 universal peak-to-trough at the corrected primary point is 0.014 (0.0135-0.0145), not 0.013; "0.013" in the tex (l.106, A.3)
#    comes from the DISC-geometry run at the pre-registered (swapped) masses (A0 0.0128 / A1 0.0127), not from any corrected run.
import numpy as np, glob
def load(fn):
    d=np.genfromtxt(fn,delimiter=',',skip_header=1); return d[:,0],d[:,1],d[:,2],{'bar':d[:,3],'R1':d[:,4],'R2':d[:,5],'R3':d[:,6]}
runs=[('ngc1560_curves_20261108.csv',3.45,'OLD pre-registered (swapped M_* 1.3e9)')]+[(f,float(f.split('_D')[1].split('_')[0]),'corrected') for f in sorted(glob.glob('ngc1560_corr_D*_20261111.csv'))]
for fn,D,lab in runs:
    r,vobs,ev,C=load(fn); f=D/3.45; r0=r/f; C['obs']=vobs
    w=(r0>=2)&(r0<=8)
    def detr(v): x=np.log10(r[w]); y=np.log10(v[w]); return y-np.polyval(np.polyfit(x,y,2),x)
    fw=(r0[w]>=3)&(r0[w]<=6.5)
    xb=detr(C['bar']); sl=lambda y,x: np.sum(x*y)/np.sum(x*x)
    print('\n== %s  [%s]  D %.2f'%(fn,lab,D))
    print('   amplitude (peak-to-trough, 3-6.5 kpc at 3.45 Mpc): '+'  '.join('%s %.4f'%(k,detr(C[k])[fw].max()-detr(C[k])[fw].min()) for k in ('obs','bar','R1','R2','R3')))
    out=[]
    for k in ('obs','R1','R2','R3'):
        y=detr(C[k]); s=sl(y,xb)
        if k=='obs':
            sig=ev[w]/vobs[w]/np.log(10); rng=np.random.default_rng(1); S=[]
            for _ in range(4000):
                vv=vobs[w]*10**rng.normal(0,sig); x=np.log10(r[w]); yy=np.log10(vv); S.append(sl(yy-np.polyval(np.polyfit(x,yy,2),x),xb))
            n=len(y); jk=np.array([sl(np.delete(y,i),np.delete(xb,i)) for i in range(n)])
            out.append('obs %.3f (pert sd %.3f, jk %.3f, n %d)'%(s,np.std(S),np.sqrt((n-1)/n*np.sum((jk-jk.mean())**2)),n))
        else: out.append('%s %.3f'%(k,s))
    print('   transmission (OLS slope, detrended log V on detrended log V_bar, 2-8 kpc): '+' | '.join(out))
    seg=lambda a,b: (r0>=a)&(r0<=b)
    s1=seg(4.0,6.03); s2=seg(6.0,6.91)
    p1=np.polyfit(r[s1],vobs[s1],1,w=1/ev[s1],cov=True); p2=np.polyfit(r[s2],vobs[s2],1,w=1/ev[s2],cov=True)
    print('   shelf 4.0-6.0 kpc (km/s/kpc): obs %.2f +- %.2f | '%(p1[0][0],np.sqrt(p1[1][0,0]))+' '.join('%s %.2f'%(k,np.polyfit(r[s1],C[k][s1],1)[0]) for k in ('R3','R1','R2')))
    print('   step  6.0-6.9 kpc (km/s/kpc): obs %.2f +- %.2f | '%(p2[0][0],np.sqrt(p2[1][0,0]))+' '.join('%s %.2f'%(k,np.polyfit(r[s2],C[k][s2],1)[0]) for k in ('R3','R1','R2')))
    for k in ('R3','R1','R2'):
        print('   %s shelf pull vs obs (weighted obs error): %+.1f sigma'%(k,(np.polyfit(r[s1],C[k][s1],1)[0]-p1[0][0])/np.sqrt(p1[1][0,0])),end='')
    print()
