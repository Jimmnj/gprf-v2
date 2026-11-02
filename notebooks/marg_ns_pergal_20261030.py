"""Per-galaxy audit of the NS evidence comparison (PREREG_NS_20261030, diagnostic; not a scored metric).
Joins NS logZ (U,B,A at sig_int 0) with MAP chi2 from marg_fit_params_20261030.csv.
Usage: python3 marg_ns_pergal_20261030.py [marg_ns_20261030.jsonl]"""
import json,sys,pandas as pd,numpy as np
f=sys.argv[1] if len(sys.argv)>1 else 'marg_ns_20261030.jsonl'
R=[json.loads(l) for l in open(f) if l.strip()]
D=pd.DataFrame([dict(gal=d['gal'],arm=d['arm'],sig=d['sig_int'],logZ=d['logZ'],rms=d['rms_med'],zD=d['post']['zD'][1],
   wD=d['post']['zD'][2]-d['post']['zD'][0],wlud=d['post']['lud'][2]-d['post']['lud'][0],npts=d['npts']) for d in R])
D=D[D.sig==0]
M=pd.concat([pd.read_csv('marg_fit_params_20261030.csv'),pd.read_csv('marg_supp_params_20261030.csv')])[['gal','arm','chi2','pen_obj','zD']]
W=D.pivot(index='gal',columns='arm',values=['logZ','rms','zD','wD','wlud']); W.columns=['%s_%s'%(a,b) for a,b in W.columns]
C=M.pivot(index='gal',columns='arm',values=['chi2','pen_obj']); C.columns=['%s_%s'%(a,b) for a,b in C.columns]
T=W.join(C).join(D.groupby('gal').npts.first())
T['dlogZ_UB']=T.logZ_U-T.logZ_B; T['dchi2_UB']=T.chi2_U-T.chi2_B; T['dpen_UB']=T.pen_obj_U-T.pen_obj_B
if 'logZ_A' in T: T['dlogZ_AU']=T.logZ_A-T.logZ_U; T['dchi2_AU']=T.chi2_A-T.chi2_U
T=T.sort_values('dlogZ_UB',ascending=False)
out=f.replace('.jsonl','_pergal.csv'); T.to_csv(out,float_format='%.3f')
g=T.dlogZ_UB.dropna(); print('U-B: N=%d sum %+.1f  top10 sum %+.1f (%.0f%%)  bottom10 sum %+.1f'%(len(g),g.sum(),g.head(10).sum(),100*g.head(10).sum()/g.sum() if g.sum() else 0,g.tail(10).sum()))
print('rank corr dlogZ_UB vs -dpen_UB/2 (MAP objective incl. priors): %.2f'%g.corr(-T.loc[g.index,'dpen_UB']/2,method='spearman'))
print('rank corr dlogZ_UB vs -dchi2_UB/2: %.2f'%g.corr(-T.loc[g.index,'dchi2_UB']/2,method='spearman'))
cols=['npts','dlogZ_UB','dchi2_UB','dpen_UB','rms_U','rms_B','zD_U','zD_B','wD_U','wD_B']
pd.set_option('display.width',200); print('\nTop 12 (U favoured):'); print(T[cols].head(12).round(2)); print('\nBottom 8 (B favoured):'); print(T[cols].tail(8).round(2))
if 'dlogZ_AU' in T:
    a=T.dlogZ_AU.dropna(); print('\nA-U: N=%d sum %+.1f  top5 sum %+.1f; rank corr with -dchi2_AU/2: %.2f'%(len(a),a.sum(),a.sort_values(ascending=False).head(5).sum(),a.corr(-T.loc[a.index,'dchi2_AU']/2,method='spearman')))
print('wrote',out)
