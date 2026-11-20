# clean_conc_20261109.py [UNRATIFIED] -- clean-data check of RESULT_RHOC_CONSTRAINED_20261108 sec 7 (concentration vs rho_c).
# PRECONDITION (JM 2026-11-09 00:00 "it was the ratio"): 9/23 r/Rd signal died in Q=1 with the innermost ~2 points dropped.
# Arms: (A0) reproduce sec 7: constrained set (marg width<0.5), marg rho_c median, r_peak,bar/Rd from all points;
#       (A1) Q=1 only, same quantities;  (A2) Q=1, fixed-M/L 1D refit of rho_c (pass2 method) with the innermost TWO points
#       dropped from the fit mask, index r_peak from r[2:];  (A3) smoother index: half-baryonic-mass radius / Rd (r[2:]);
#       controls: fixed-M/L refit WITHOUT dropping (separates method from drop); full-sample versions of each index.
# Resumable: per-galaxy refits appended to clean_conc_part_20261109.csv; run until 'DONE'. Columns asserted: Rd(NGC2841)=3.64.
import numpy as np, pandas as pd, os, io, contextlib, json, time, sys, warnings; warnings.filterwarnings('ignore')
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
assert abs(T.loc['NGC2841'].Rd-3.64)<1e-6 and abs(T.loc['DDO154'].Rd-0.37)<1e-6, 'MRT column error'
M=[json.loads(l) for l in open('marg_ns_20261030.jsonl')]; M=[r for r in M if r['arm']=='A' and r['sig_int']==0.03]
MG=pd.DataFrame([dict(gal=r['gal'],lrc_marg=r['post']['lrc'][1],w=r['post']['lrc'][2]-r['post']['lrc'][0]) for r in M]).set_index('gal')
gals=sorted(g for g in MG[MG.w<0.5].index if g in T.index and T.loc[g].Rd>0 and (P.gal==g).sum()>=5)
PART='clean_conc_part_20261109.csv'; done=set(pd.read_csv(PART).gal) if os.path.exists(PART) else set()
def refit(g,drop):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vo=d.vobs.values; w=1/(0.1*vo)**2
    m=((d.ev/d.vobs<0.1)&(d.vobs>0)&(d.vbar2>0)).values.copy(); m[:drop]=False
    if m.sum()<3: return np.nan
    def chi(l):
        o=solve_gal(g,l); return 1e12 if o is None else float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m]))
    return minimize_scalar(chi,bounds=(4.5,10.5),method='bounded',options={'xatol':0.02}).x
t0=time.time()
for g in [g for g in gals if g not in done]:
    if time.time()-t0>140: break
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); Rd=T.loc[g].Rd
    def idx(k):
        rr,vv=r[k:],vb2[k:]; Menc=vv*rr; half=0.5*Menc[-1]
        j=np.argmax(Menc>=half); rh=np.interp(half,[Menc[j-1],Menc[j]],[rr[j-1],rr[j]]) if j>0 else rr[0]
        return rr[np.argmax(vv)]/Rd, rh/Rd
    pk0,hm0=idx(0); pk2,hm2=idx(2)
    row=dict(gal=g,Q=int(T.loc[g].Q),npts=len(r),lrc_marg=MG.loc[g].lrc_marg,lrc_fix0=refit(g,0),lrc_fix2=refit(g,2),pk0=pk0,hm0=hm0,pk2=pk2,hm2=hm2,
             lMb=np.log10(0.5*T.loc[g].L*1e9+1.33*T.loc[g].MHI*1e9))
    pd.DataFrame([row]).to_csv(PART,mode='a',header=not os.path.exists(PART),index=False)
R=pd.read_csv(PART).drop_duplicates('gal'); print('progress %d / %d'%(len(R),len(gals)))
if len(R)>=len(gals):
    def resid(y,x): X=np.column_stack([np.ones(len(y)),x]); return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
    def line(lab,D,y,x):
        D=D[[y,x,'lMb']].dropna(); n=len(D); rho,p=spearmanr(D[x],D[y]); lrho,lp=spearmanr(np.log10(D[x]),D[y])
        rp,pp=spearmanr(resid(D[y].values,D.lMb.values),resid(np.log10(D[x].values),D.lMb.values))
        print('  %-52s N %3d  rho %+.2f (p %.3f) | at fixed M_b %+.2f (p %.3f)'%(lab,n,lrho,lp,rp,pp))
    print('DONE. Spearman of log(index) vs log rho_c | partial after log M_b')
    Q1=R[R.Q==1]
    print('A0 reproduce sec 7 (constrained, marg rho_c, all points):'); line('r_peak/Rd all pts, marg',R,'lrc_marg','pk0'); line('r_half/Rd all pts, marg',R,'lrc_marg','hm0')
    print('control: fixed-M/L refit, all points (method check):'); line('r_peak/Rd all pts, fix0',R,'lrc_fix0','pk0')
    print('A1 Q=1 only, marg rho_c, all points:'); line('r_peak/Rd all pts, marg',Q1,'lrc_marg','pk0'); line('r_half/Rd all pts, marg',Q1,'lrc_marg','hm0')
    print('A2 Q=1, inner two dropped from index AND fit (fixed-M/L refit):'); line('r_peak/Rd r[2:], fix2',Q1,'lrc_fix2','pk2'); line('r_peak/Rd r[2:], fix0 (drop index only)',Q1,'lrc_fix0','pk2')
    print('A3 Q=1, half-mass index, inner two dropped:'); line('r_half/Rd r[2:], fix2',Q1,'lrc_fix2','hm2'); line('r_half/Rd r[2:], marg',Q1,'lrc_marg','hm2')
    print('A2/A3 all Q (constrained), inner two dropped:'); line('r_peak/Rd r[2:], fix2',R,'lrc_fix2','pk2'); line('r_half/Rd r[2:], fix2',R,'lrc_fix2','hm2')
    print('refit vs marg: Spearman %.2f; median |fix0-marg| %.2f dex; median |fix2-fix0| %.2f dex'%(spearmanr(R.lrc_fix0,R.lrc_marg)[0],(R.lrc_fix0-R.lrc_marg).abs().median(),(R.lrc_fix2-R.lrc_fix0).abs().median()))
    print('bins by r_peak/Rd (r[2:]), Q=1, fix2 median rho_c:'); 
    for lo,hi in ((0,1.5),(1.5,3),(3,99)): s=Q1[(Q1.pk2>=lo)&(Q1.pk2<hi)]; print('   %4.1f-%4.1f N %3d  median %.2f'%(lo,hi,len(s),s.lrc_fix2.median()))
    print('Q counts in constrained set:',R.Q.value_counts().to_dict())
