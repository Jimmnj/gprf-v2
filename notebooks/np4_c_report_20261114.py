# np4_c_report_20261114.py [UNRATIFIED] test1d_report_20261018.py statistics (sl() verbatim, rng seed 4) on the old csvs (reproduction) and the V3 part files.
import numpy as np, pandas as pd
six=['CamB','D512-2','NGC4068','NGC6789','UGC05750','UGC07232']
def sl(x,y,seed=4):
    rng=np.random.default_rng(seed); A=np.polyfit(x,y,1)[0]; bs=[np.polyfit(x[i],y[i],1)[0] for i in (rng.integers(0,len(x),len(x)) for _ in range(1000))]; lo,hi=np.percentile(bs,[2.5,97.5])
    return '%+.2f [%+.2f,%+.2f]'%(A,lo,hi)
f=lambda x:(np.percentile(x,84)-np.percentile(x,16))/2
for v in ('0.7','mdep'):
    O=pd.read_csv('test1d_yd_%s_20261018.csv'%v).set_index('gal'); N=pd.read_csv('np4_test1d_v3_%s_20261114_part.csv'%v).set_index('gal').reindex(O.index); N['Mb']=O.Mb; N['rmax']=O.rmax; N['lr0_d']=O.lr0_d.where(N.lr0_d.isna(),N.lr0_d); N['lrho0_d']=O.lrho0_d.where(N.lrho0_d.isna(),N.lrho0_d)
    d=(N.brent-O.l_rhoc).abs(); print('== Yd %s: reproduction: max |Brent now - record| %.4f (N %d) | moved by V3 (|v3 - brent| > 0.05): %s | edge %d | failed Burkert-of-model %d | six in set: %s'%(v,d.max(),len(O),{g:(round(N.brent[g],2),round(N.l_rhoc[g],2)) for g in N.index[(N.l_rhoc-N.brent).abs()>0.05]},int(N.edge.sum()),int(N.lr0_m.isna().sum()),[g for g in six if g in O.index]))
    bad=list(N.index[(N.l_rhoc>11)|(N.brent<4.05)]); NX=N.drop(bad); NL=N.copy(); NL.loc[bad,'l_rhoc']=4.5
    for tag,R in (('OLD',O),('OLD without %s'%bad,O.drop(bad)),('V3 as defined (first interior minimum)',N),('V3 without %s'%bad,NX),('V3, %s at the scan floor 4.5'%bad,NL)):
        ok=(10**R.lr0_d<R.rmax)&R.lr0_m.notna(); S=R[ok]; x=np.log10(S.Mb.values); xa=np.log10(R.Mb.values); mu_d=S.lr0_d+S.lrho0_d-6; mu_m=S.lr0_m+S.lrho0_m-6
        print('   %-52s N=%3d/%3d | MODEL log rho_c slope %s median %.2f rob %.2f sd %.2f | mu0 data %s | mu0 model %s median %.2f'%(tag,len(S),len(R),sl(xa,R.l_rhoc.values),R.l_rhoc.median(),f(R.l_rhoc),R.l_rhoc.std(),sl(x,mu_d.values),sl(x,mu_m.values),mu_m.median()))
print('record: 0.7: -0.26 [-0.38,-0.15], mu0 data +0.13, model +0.11, N 97; mdep: -0.30 [-0.40,-0.18], +0.12, +0.09, N 100; medians 7.21 / 7.28')
