"""PREREG_BURKMU_20261101 metrics J1-J4 / L1-L3. Reads marg_ns (A,U,B), cdm3_ns (D,K), burkmu_ns (M)."""
import json,numpy as np,pandas as pd
rows=[]
for f in ('marg_ns_20261030.jsonl','cdm3_ns_20261031.jsonl','cdm_ns_20261031.jsonl','burkmu_ns_20261101.jsonl'):
    for l in open(f):
        if not l.strip(): continue
        d=json.loads(l); p=d['post']
        rows.append(dict(gal=d['gal'],arm=d['arm'],sig=d['sig_int'],logZ=d['logZ'],err=d['logZerr'],rms=d['rms_med'],zD=p['zD'][1],lLmax=d.get('logLmax',np.nan),
                         wrc=(p['lrc'][2]-p['lrc'][0])/2 if 'lrc' in p else np.nan, wrho=(p['lrho0'][2]-p['lrho0'][0])/2 if 'lrho0' in p else np.nan,
                         lmu=(p['lrho0'][1]+p['lrc'][1]-6) if ('lrc' in p and 'lrho0' in p) else np.nan))
D=pd.DataFrame(rows); out=[]; say=lambda t:(print(t),out.append(t))
say('PREREG_BURKMU metrics  (294 M jobs; coverage '+', '.join('%s@%g:%d'%(a,s,n) for (a,s),n in D.groupby(['arm','sig']).gal.nunique().items() if a in 'MKDAU')+')')
pref={}
for a in 'MKDAU':
    sub=D[D.arm==a]; sigs=[s for s in sorted(sub.sig.unique()) if s>0]
    tot={s:sub[sub.sig==s].logZ.sum() for s in sigs}; pref[a]=max(tot,key=tot.get)
    say('arm %s: '%a+'  '.join('sig %g %.1f'%(s,tot[s]) for s in sigs)+' -> preferred %g'%pref[a])
def pair(a1,a2,s1=None,s2=None):
    s1=pref[a1] if s1 is None else s1; s2=pref[a2] if s2 is None else s2
    x=D[(D.arm==a1)&(D.sig==s1)].set_index('gal'); y=D[(D.arm==a2)&(D.sig==s2)].set_index('gal'); g=x.index.intersection(y.index)
    d=(x.loc[g].logZ-y.loc[g].logZ); e=np.sqrt((x.loc[g].err**2+y.loc[g].err**2).sum()); return d,e
say('\nJ1 Delta ln Z (M minus arm) at preferred sig_int, and at 0.03:')
res={}
for a in 'AKDU':
    for lab,s in (('pref',None),('0.03',0.03)):
        d,e=pair('M',a,None if lab=='pref' else 0.03, None if lab=='pref' else 0.03); bulk=d[abs(d)<=20]; res[(a,lab)]=d.sum()
        say('  M - %s [%s]: N=%d sum %+.1f +- %.1f | median %+.2f | >+2.5: %d, <-2.5: %d | bulk(|d|<=20) %+.1f (N=%d)'%(a,lab,len(d),d.sum(),e,d.median(),(d>2.5).sum(),(d<-2.5).sum(),bulk.sum(),len(bulk)))
say('\nJ2 Delta ln L_max (where stored):')
for a in 'KD':
    x=D[(D.arm=='M')&(D.sig==0.03)].set_index('gal'); y=D[(D.arm==a)&(D.sig==0.03)].set_index('gal'); g=x.index.intersection(y.index)
    say('  M - %s at 0.03: sum %+.1f over %d (A lacks logLmax; RESULT_OCCAM: K - A = +97)'%(a,(x.loc[g].lLmax-y.loc[g].lLmax).sum(),len(g)))
say('\nJ3 distance shifts at preferred sig_int:')
for a in 'MKDAU':
    sub=D[(D.arm==a)&(D.sig==pref[a])]; say('  %s: rms zD %.2f  >2sig %d  >3sig %d  rms(resid) %.3f'%(a,np.sqrt(np.mean(sub.zD**2)),(abs(sub.zD)>2).sum(),(abs(sub.zD)>3).sum(),sub.rms.mean()))
say('\nJ4 posterior widths (half 16-84, dex) at preferred sig_int; prior: M log mu0 N(.,0.2), log rc U[-1,2] (half-width 0.87); K both U (1.5 / 0.87):')
for a in 'MK':
    sub=D[(D.arm==a)&(D.sig==pref[a])]; w=np.sqrt(sub.wrc**2+sub.wrho**2)
    say('  %s: median w(log rc) %.2f, w(log rho0) %.2f; independence bound on w(log mu0) median %.2f, <0.2 in %d/%d; median log mu0 %.2f (prior 2.15)'%(a,sub.wrc.median(),sub.wrho.median(),w.median(),(w<0.2).sum(),len(w),sub.lmu.median()))
say('\nDECISIONS (stamped rule: |M-A|<=5 parity; M-A>+5 M wins; M-A<-5 relation economy survives)')
mA=res[('A','pref')]; say('L1 M vs A: %+.1f -> %s'%(mA,'M BEATS A' if mA>5 else ('PARITY' if abs(mA)<=5 else 'RELATION ECONOMY SURVIVES')))
say('L2 M vs K: %+.1f (expected > +100)'%res[('K','pref')]); say('L3 M vs D: %+.1f'%res[('D','pref')])
open('burkmu_report_20261102.txt','w').write('\n'.join(out))
