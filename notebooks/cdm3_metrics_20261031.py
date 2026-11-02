# cdm3_metrics_20261031.py -- PREREG_CDM3 metrics G1-G6, decisions H1-H4. Reads marg_ns_20261030.jsonl, cdm_ns_20261031.jsonl, cdm3_ns_20261031.jsonl.
import json,collections,statistics as st,numpy as np
Z=collections.defaultdict(dict); LL=collections.defaultdict(dict); ZD=collections.defaultdict(dict); RMS=collections.defaultdict(dict); XD={}; GAM={}
for f in ['marg_ns_20261030.jsonl','cdm_ns_20261031.jsonl','cdm3_ns_20261031.jsonl']:
    for l in open(f):
        r=json.loads(l); k=(r['arm'],r['sig_int']); Z[k][r['gal']]=r['logZ']
        if 'logLmax' in r: LL[k][r['gal']]=r['logLmax']
        ZD[k][r['gal']]=r['post']['zD'][1]; RMS[k][r['gal']]=r.get('rms_med')
        if r['arm']=='D' and r['sig_int']==0.03: XD[r['gal']]=r['X']; GAM[r['gal']]=r['abg'][2]
print('coverage:',{'%s@%g'%k:len(v) for k,v in sorted(Z.items()) if k[0] in 'DK'})
pref={}
for arm in 'AUBNFDK':
    cand={s:sum(Z[(arm,s)].values()) for s in (0,0.03,0.06) if len(Z[(arm,s)])==147}
    pref[arm]=max(cand,key=cand.get); print('arm %s sum logZ'%arm,{s:round(v,1) for s,v in cand.items()},'-> preferred',pref[arm])
def S(a,b):
    g=set(Z[a])&set(Z[b]); d={x:Z[a][x]-Z[b][x] for x in g}; v=sorted(d.values())
    return sum(v),st.median(v),sum(x>2.5 for x in v),sum(x<-2.5 for x in v),sum(x for x in v if abs(x)<=20),len(g),d
print('\nG1 Delta logZ at preferred sig_int (and matched 0.03 / 0.06):')
for a,b in [('A','D'),('A','K'),('U','D'),('B','D'),('D','N'),('K','D'),('A','N'),('F','D'),('K','N')]:
    s,m,p,n,bulk,ng,_=S((a,pref[a]),(b,pref[b])); s3=S((a,0.03),(b,0.03))[0]; s6=S((a,0.06),(b,0.06))[0]
    print('  %s - %s: %+8.1f  median %+.2f  tail +%d/-%d  bulk %+.1f   | @0.03 %+.1f  @0.06 %+.1f'%(a,b,s,m,p,n,bulk,s3,s6))
print('\nG2 logLmax at 0.03:')
for a,b in [('A','D'),('A','K'),('D','N'),('K','D')]:
    g=set(LL[(a,0.03)])&set(LL[(b,0.03)]); print('  %s - %s: %+.1f (n %d)'%(a,b,sum(LL[(a,0.03)][x]-LL[(b,0.03)][x] for x in g),len(g)))
print('\nG3 distance shifts @0.03 (rms; >2 / >3 sigma) and G4 median per-galaxy rms:')
for arm in 'AUBNFDK':
    z=np.array(list(ZD[(arm,0.03)].values())); print('  %s: zD rms %.2f  >2sig %d  >3sig %d   rms_med %.3f'%(arm,np.sqrt(np.mean(z*z)),(abs(z)>2).sum(),(abs(z)>3).sum(),np.median([v for v in RMS[(arm,0.03)].values() if v is not None])))
X=np.array(list(XD.values())); gm=np.array(list(GAM.values()))
print('\nG5 DC14 @0.03: X median %.2f [16-84: %.2f, %.2f]; clamped (X<=-4.1 or >=-1.3): %d/147; gamma median %.2f, gamma<0.5: %d'%(np.median(X),*np.percentile(X,[16,84]),((X<=-4.1)|(X>=-1.3)).sum(),np.median(gm),(gm<0.5).sum()))
print('\nG6 per-galaxy logZ @0.03 (A, N, D, K):')
for g in ['NGC0247','F571-8','NGC2403','UGC12732','UGC00128','IC2574','NGC3109','IC4202','ESO563-G021']:
    print('  %-12s'%g,'  '.join('%s %8.1f'%(a,Z[(a,0.03)][g]) for a in 'ANDK'))
dA_D=S(('A',pref['A']),('D',pref['D']))[0]; dA_K=S(('A',pref['A']),('K',pref['K']))[0]; dD_N=S(('D',pref['D']),('N',pref['N']))[0]; dU_D=S(('U',pref['U']),('D',pref['D']))[0]; dA_N=S(('A',pref['A']),('N',pref['N']))[0]
v=lambda d,a,b: a+' favoured' if d>5 else (b+' favoured' if d<-5 else 'indistinguishable')
print('\nDECISIONS'); print('H1 A vs D: %s (%+.1f); A vs K: %s (%+.1f)'%(v(dA_D,'A','D'),dA_D,v(dA_K,'A','K'),dA_K))
print('H2 feedback core closes the gap: %s (D-N %+.1f; A-D %+.1f vs 0.5 x A-N = %+.1f)'%('YES' if dD_N>5 and dA_D<0.5*dA_N else 'NO',dD_N,dA_D,0.5*dA_N))
print('H3 U vs D: %s (%+.1f)'%(v(dU_D,'U','D'),dU_D))
for arm in 'DK':
    z=np.array(list(ZD[(arm,0.03)].values())); c=(abs(z)>2).sum(); print('H4 %s bought with distances: %s (%d vs 1.5 x N = 6)'%(arm,'YES' if c>6 else 'NO',c))
