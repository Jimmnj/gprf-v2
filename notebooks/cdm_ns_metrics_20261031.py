"""PREREG_CDM_20261031 metrics. Reads cdm_ns_20261031.jsonl (arms N,F) and marg_ns_20261030.jsonl (arms U,B,A)."""
import json,sys,numpy as np,pandas as pd
rows=[]
for f in ('marg_ns_20261030.jsonl','cdm_ns_20261031.jsonl'):
    for l in open(f):
        if not l.strip(): continue
        d=json.loads(l); rows.append(dict(gal=d['gal'],arm=d['arm'],sig=d['sig_int'],logZ=d['logZ'],err=d['logZerr'],rms=d['rms_med'],zD=d['post']['zD'][1],
            lLmax=d.get('logLmax',np.nan),lM=d['post'].get('lM200',[None,np.nan])[1],lc=d.get('lc_med',np.nan),line=d.get('cM_line',np.nan)))
D=pd.DataFrame(rows); out=[]; say=lambda t:(print(t),out.append(t))
say('PREREG_CDM_20261031 metrics'); say('coverage: '+', '.join('%s@%g:%d'%(a,s,n) for (a,s),n in D.groupby(['arm','sig']).gal.nunique().items()))
pref={}
for a in 'UBNF':
    sub=D[D.arm==a]; sigs=sorted(sub.sig.unique())
    if not len(sigs): continue
    gs=set.intersection(*[set(sub[sub.sig==s].gal) for s in sigs]); tot={s:sub[(sub.sig==s)&(sub.gal.isin(gs))].logZ.sum() for s in sigs}; pref[a]=max(tot,key=tot.get)
    say('arm %s N=%d: '%(a,len(gs))+'  '.join('sig %g %.1f'%(s,tot[s]) for s in sigs)+' -> preferred %g'%pref[a])
pref['A']=0.0
def pair(a1,a2):
    x=D[(D.arm==a1)&(D.sig==pref[a1])].set_index('gal'); y=D[(D.arm==a2)&(D.sig==pref[a2])].set_index('gal'); g=x.index.intersection(y.index)
    d=(x.loc[g].logZ-y.loc[g].logZ); e=np.sqrt((x.loc[g].err**2+y.loc[g].err**2).sum()); return d,e
say('\nM1/M2 Delta logZ at preferred sig_int:')
for a1,a2 in (('U','N'),('A','N'),('B','N'),('F','N'),('U','F'),('A','F')):
    if a1 not in pref or a2 not in pref: continue
    d,e=pair(a1,a2); bulk=d[abs(d)<=20]
    verdict='%s favoured'%a1 if d.sum()>5 else ('%s favoured'%a2 if d.sum()<-5 else 'indistinguishable')
    say('  %s - %s: N=%d sum %+.1f +- %.1f | median %+.2f | >+2.5: %d, <-2.5: %d | bulk(|d|<=20) %+.1f (N=%d) -> %s'%(a1,a2,len(d),d.sum(),e,d.median(),(d>2.5).sum(),(d<-2.5).sum(),bulk.sum(),len(bulk),verdict))
say('\nM3 distance shifts at preferred sig_int:'); cnt={}
for a in 'UBANF':
    if a not in pref: continue
    sub=D[(D.arm==a)&(D.sig==pref[a])]; cnt[a]=int((abs(sub.zD)>2).sum()); say('  %s: rms %.2f  >2sig %d  >3sig %d'%(a,np.sqrt(np.mean(sub.zD**2)),cnt[a],(abs(sub.zD)>3).sum()))
say('\nM4 Delta logL_max (raw fit) at sig 0, N vs F only (U/B/A files lack logLmax; use chi2 from MAP csv if needed):')
for a in 'NF':
    sub=D[(D.arm==a)&(D.sig==0)]; say('  %s: sum logLmax %.1f over %d'%(a,sub.lLmax.sum(),len(sub)))
say('\nM5 rms at preferred: '+'  '.join('%s %.3f'%(a,D[(D.arm==a)&(D.sig==pref[a])].rms.mean()) for a in 'UBANF' if a in pref))
say('\nM6 concentration vs Dutton-Maccio line (posterior medians, sig 0):')
for a in 'NF':
    sub=D[(D.arm==a)&(D.sig==0)]; off=sub.lc-sub.line; say('  %s: median offset %+.2f dex, rms %.2f, |off|>0.3: %d/%d; median log M200 %.2f'%(a,off.median(),np.sqrt(np.mean(off**2)),(abs(off)>0.3).sum(),len(off),sub.lM.median()))
say('\nDECISIONS')
for a1 in 'UAB':
    d,e=pair(a1,'N'); say('D1 %s vs N: %s (%+.1f)'%(a1,'%s favoured'%a1 if d.sum()>5 else ('N favoured' if d.sum()<-5 else 'indistinguishable'),d.sum()))
d,e=pair('F','N'); offF=(D[(D.arm=='F')&(D.sig==0)].lc-D[(D.arm=='F')&(D.sig==0)].line); say('D2 haloes leave c-M: %s (F-N %+.1f; F rms offset %.2f dex)'%('YES' if d.sum()>5 and np.sqrt(np.mean(offF**2))>0.2 else 'NO',d.sum(),np.sqrt(np.mean(offF**2))))
say('D3 N bought with distances: %s (N %d vs U %d beyond 2 sigma)'%('YES' if cnt.get('N',0)>1.5*max(cnt.get('U',0),1) else 'NO',cnt.get('N',0),cnt.get('U',0)))
open('cdm_ns_report_20261031.txt','w').write('\n'.join(out))
