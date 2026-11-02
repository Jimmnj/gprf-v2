# marg_ns_metrics_20261030.py -- metrics M1-M5 and decisions D1-D4 of PREREG_NS_20261030 from marg_ns_20261030.jsonl. Usage: python3 marg_ns_metrics_20261030.py [file]
import sys, json, numpy as np, pandas as pd
f=sys.argv[1] if len(sys.argv)>1 else 'marg_ns_20261030.jsonl'
R=[json.loads(l) for l in open(f) if l.strip()]
D=pd.DataFrame([dict(gal=d['gal'],arm=d['arm'],sig=d['sig_int'],logZ=d['logZ'],err=d['logZerr'],rms=d['rms_med'],zD=d['post']['zD'][1],zDlo=d['post']['zD'][0],zDhi=d['post']['zD'][2],zi=d['post']['zi'][1],lud=d['post']['lud'][1],lrc=(d['post'].get('lrc',[None,None,None])[1]),sec=d['sec']) for d in R])
out=[]; say=lambda s:(print(s),out.append(s))
say('PREREG_NS_20261030 metrics from %s: %d results'%(f,len(D)))
say('coverage: '+', '.join('%s@%g: %d gal'%(a,s,n) for (a,s),n in D.groupby(['arm','sig']).gal.nunique().items()))
# common galaxy set per comparison
def common(a1,s1,a2,s2):
    x=D[(D.arm==a1)&(D.sig==s1)].set_index('gal'); y=D[(D.arm==a2)&(D.sig==s2)].set_index('gal'); g=x.index.intersection(y.index); return x.loc[g],y.loc[g]
sigs=sorted(D[D.arm=='U'].sig.unique())
say('\nM2 summed logZ by sig_int (common galaxies within arm):')
pref={}
for a in 'UB':
    sub=D[D.arm==a]; gs=set.intersection(*[set(sub[sub.sig==s].gal) for s in sigs]) if len(sigs)>1 else set(sub.gal)
    tot={s:sub[(sub.sig==s)&(sub.gal.isin(gs))].logZ.sum() for s in sigs}; best=max(tot,key=tot.get); pref[a]=best
    say('  arm %s (N=%d): '%(a,len(gs))+'  '.join('sig %g: %.1f'%(s,tot[s]) for s in sigs)+'  -> preferred sig_int %g'%best)
say('\nM1 Delta logZ (U - B) summed, per sig_int:')
for s in sigs:
    x,y=common('U',s,'B',s); d=x.logZ-y.logZ; e=np.sqrt((x.err**2+y.err**2).sum())
    say('  sig %g: N=%d  sum %+.1f +- %.1f  | galaxies with Delta>+2.5: %d, < -2.5: %d  | median per-galaxy %+.2f'%(s,len(d),d.sum(),e,(d>2.5).sum(),(d<-2.5).sum(),d.median()))
x,y=common('U',pref['U'],'B',pref['B']); d1=(x.logZ-y.logZ).sum() if len(x) else np.nan
say('  at each arm\'s preferred sig_int (U@%g vs B@%g): sum %+.1f (N=%d)'%(pref['U'],pref['B'],d1,len(x)))
say('\nM3 Delta logZ (A - U) at sig_int 0:')
if (D.arm=='A').any():
    x,y=common('A',0.0,'U',0.0); d3=(x.logZ-y.logZ); say('  N=%d  sum %+.1f +- %.1f | >+2.5: %d, <-2.5: %d'%(len(d3),d3.sum(),np.sqrt((x.err**2+y.err**2).sum()),(d3>2.5).sum(),(d3<-2.5).sum()))
else: d3=pd.Series(dtype=float); say('  arm A not run')
say('\nM4 posterior-median distance shifts z_D (prior sigmas), at each arm\'s preferred sig_int:')
cnt={}
for a in 'UBA':
    s=pref.get(a,0.0) if a!='A' else 0.0; sub=D[(D.arm==a)&(D.sig==s)]
    if not len(sub): continue
    cnt[a]=(sub.zD.abs()>2).sum(); say('  arm %s@%g: N=%d  rms %.2f  |zD|>2: %d  |zD|>3: %d  median 16-84 width %.2f'%(a,s,len(sub),np.sqrt((sub.zD**2).mean()),cnt[a],(sub.zD.abs()>3).sum(),(sub.zDhi-sub.zDlo).median()))
say('\nM5 posterior-median residual rms: '+'  '.join('%s@%g %.3f'%(a,s,D[(D.arm==a)&(D.sig==s)].rms.mean()) for (a,s) in sorted(set(zip(D.arm,D.sig)))))
say('\nDECISIONS'+('  [PREVIEW: coverage incomplete, %d/147 U, %d/147 B, %d/147 A]'%tuple(D[(D.arm==a)&(D.sig==0)].gal.nunique() for a in 'UBA') if min(D[(D.arm==a)&(D.sig==0)].gal.nunique() for a in 'UBA')<147 or len(sigs)<3 else ''))
D1='FAVOURS universal-rho_c relation' if d1>5 else ('FAVOURS empirical function' if d1<-5 else 'INDISTINGUISHABLE')
D4=('bought with distances' if cnt.get('U',0)>1.5*max(cnt.get('B',0),1) else 'distance shifts comparable')
say('D1 (Q1): %s (sum Delta logZ %+.1f)  -- D4: %s (U %d vs B %d beyond 2 sigma)'%(D1,d1,D4,cnt.get('U',0),cnt.get('B',0)))
if len(sigs)<3: D2='INSUFFICIENT COVERAGE (sig_int values run: %s)'%sigs
else: D2='CONFIRMED' if (pref.get('U',0)>=0.06 and pref.get('B',0)<=0.03) else ('NOT CONFIRMED' if pref.get('U')==pref.get('B') else 'MIXED (U %g, B %g)'%(pref.get('U'),pref.get('B')))
say('D2 (Q2): scatter tension %s (preferred sig_int U %s, B %s)'%(D2,pref.get('U'),pref.get('B')))
if len(d3): say('D3 (Q3): per-galaxy rho_c %s (sum Delta logZ A-U %+.1f)'%('SUPPORTED' if d3.sum()>5 else ('NOT SUPPORTED' if d3.sum()<-5 else 'INDIFFERENT'),d3.sum()))
say('\nrun time: total %.1f core-hours; median per job U %.0f s, B %.0f s, A %.0f s'%(D.sec.sum()/3600,*(D[D.arm==a].sec.median() if (D.arm==a).any() else 0 for a in 'UBA')))
open(f.replace('.jsonl','_report.txt'),'w').write('\n'.join(out))
