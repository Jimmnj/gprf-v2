# occam_20261101.py [UNRATIFIED] S2 of REVIEW_COLD_v2.0.1: posterior widths, compression and effective parameter counts per arm at sig_int 0.03, from the NS jsonl files.
import json, numpy as np, collections, pandas as pd
rows=[]
for f in ['marg_ns_20261030.jsonl','cdm_ns_20261031.jsonl','cdm3_ns_20261031.jsonl','bg_ns_20261031.jsonl']:
    for l in open(f):
        try: r=json.loads(l)
        except: continue
        if abs(r.get('sig_int',-1)-0.03)<1e-6: rows.append(r)
U=lambda w: w/np.sqrt(12)
# prior sigma per parameter; Burkert's 'lrc' key is the core radius (U[-1,2]), the relation's 'lrc' is log rho_c (U[6,10])
def prior_sig(arm,k):
    if k=='lrc': return U(4.0) if arm=='A' else U(3.0)
    return {'lM200':U(5.5),'dlc':0.11,'lc':U(2.0),'lrho0':U(5.0),'lud':0.1,'lub':0.1,'zD':1.0,'zi':1.0}[k]
model={'A':['lrc'],'D':['lM200','dlc'],'N':['lM200','dlc'],'F':['lM200','lc'],'K':['lrho0','lrc'],'U':[],'B':[],'G':[]}
out=collections.defaultdict(list)
for r in rows:
    arm=r['arm']; g=r['gal']
    for k,(p16,p50,p84) in r['post'].items():
        sp=0.5*(p84-p16); rho=sp/prior_sig(arm,k)
        out[(arm,k)].append((g,rho,np.log(max(rho,1e-6))))
print('Arm  param   N   median sigma_post/sigma_prior   frac constrained (rho<0.5)   sum ln(prior/post)  [model params only -> Occam cost per arm]')
tot={}
for arm in 'ABDKFNUG':
    s=0.0
    for k in model[arm]+['lud','lub','zD','zi']:
        v=out.get((arm,k))
        if not v: continue
        rh=np.array([x[1] for x in v]); ln=np.array([x[2] for x in v])
        flag='*' if k in model[arm] else ' '
        print('%s  %-6s%s %3d   %5.2f                         %4.2f                      %7.1f'%(arm,k,flag,len(v),np.median(rh),(rh<0.5).mean(),-ln.sum()))
        if k in model[arm]: s+=-ln.sum()
    tot[arm]=s
print('\nOccam cost of MODEL parameters, summed over galaxies (nats = sum ln(sigma_prior/sigma_post)):',{a:round(v,1) for a,v in tot.items()})
# effective constrained-parameter count per galaxy: sum(1-rho^2) over model params (Bayesian dimensionality approx)
print('\nEffective number of constrained model parameters per galaxy (mean of sum(1-rho^2)):')
for arm in 'ADKFN':
    per=collections.defaultdict(float)
    for k in model[arm]:
        for g,rho,ln in out[(arm,k)]: per[g]+=max(0,1-rho**2)
    print('  %s: %.2f  (nominal %d)'%(arm,np.mean(list(per.values())),len(model[arm])))
# nuisance compression: do arms squeeze D/i/Upsilon differently?
print('\nNuisance posterior/prior width (median) per arm:')
for arm in 'ABDKFNUG':
    print('  %s: '%arm+'  '.join('%s %.2f'%(k,np.median([x[1] for x in out[(arm,k)]])) for k in ['lud','zD','zi'] if out.get((arm,k))))
