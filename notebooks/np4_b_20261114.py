# np4_b_20261114.py [UNRATIFIED] Item 4, part B: (5) fixed-arm class medians on the 141 (table of record, V3); (6) SB_disk raw/partial, both arms (rhoc_corr
# method verbatim: Spearman; partial = Spearman of residuals after OLS on lMb, lSig, fgas); (7) marginalised-arm Spearmans of RESULT_RHOC_SIGMAD (P3, 1b);
# (8) matched pairs (diversity_20261105 method verbatim, marginalised MAP rho_c); (10) BS-AMEND-1 numbers. OLD = 147 / old table, must reproduce the record.
import json,itertools,numpy as np,pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
six=['CamB','D512-2','NGC4068','NGC6789','UGC05750','UGC07232']
f=lambda x:(np.percentile(x,84)-np.percentile(x,16))
T=pd.read_csv('rhoc_corr_20261105_table.csv').set_index('gal'); V=pd.read_csv('ep3_val_20261114_part.csv').set_index('key')
J=[json.loads(l) for l in open('marg_ns_20261030.jsonl')]; bul=pd.Series({r['gal']:r['bul'] for r in J if r['arm']=='A' and r['sig_int']==0.03})
tab=pd.read_csv('pass2_table2_part.csv').set_index('gal'); print('table of record lrc == ep3_val lrc_v3: max |diff| %.4f (N %d)'%((tab.lrc-V.lrc_v3.reindex(tab.index)).abs().max(),len(tab)))
def resid(y,X):
    X=np.column_stack([np.ones(len(y))]+[X[c].values for c in X.columns]); return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
for tag in ('OLD 147 (old fits)','V3 147','V3 141'):
    D=T.copy(); D['bul']=bul.reindex(D.index)
    if tag!='OLD 147 (old fits)': D['lrc']=tab.lrc.reindex(D.index)   # table of record (V3, md5 82dce260)
    if tag=='V3 141': D=D.drop(six)
    print('=== (5) fixed-arm class medians, %s: N, median, median - all, rob spread'%tag); m0=D.lrc.median()
    q=pd.qcut(D.lMb,4)
    cl=[('ALL',D.lrc==D.lrc),('bulge',D.bul==True),('no bulge',D.bul==False),('gas-rich f_gas>0.5',D.fgas>0.5),('gas-poor f_gas<=0.5',D.fgas<=0.5),('gas-poor f_gas<0.2',D.fgas<0.2),
        ('no V_f (slow risers)',~(D.Vf>0)),('V_f>0',D.Vf>0),('no V_f & lMb<9.5',(~(D.Vf>0))&(D.lMb<9.5))]+[('M_b quartile %d (%.1f-%.1f)'%(i+1,c.left,c.right),q==c) for i,c in enumerate(q.cat.categories)]
    for lab,m in cl:
        m=np.asarray(m); print('   %-30s N %3d  %.2f  %+.2f  rob %.2f'%(lab,m.sum(),D.lrc[m].median(),D.lrc[m].median()-m0,f(D.lrc[m])/2))
    sr=~(D.Vf>0); print('   slow risers vs V_f>0: MW p %.3f; vs no-bulge V_f>0: p %.3f'%(mannwhitneyu(D.lrc[sr],D.lrc[D.Vf>0]).pvalue,mannwhitneyu(D.lrc[sr],D.lrc[(D.Vf>0)&(D.bul==False)]).pvalue))
    print('=== (6) SB_disk, %s'%tag)
    for y,lab in (('lrc','fixed'),('lrc_marg','marginalised MAP')):
        d=D[[y,'lSBdisk','lMb','lSig','fgas']].dropna(); r,p=spearmanr(d[y],d.lSBdisk); rp,pp=spearmanr(resid(d[y].values,d[['lMb','lSig','fgas']]),resid(d.lSBdisk.values,d[['lMb','lSig','fgas']]))
        print('   %-18s N %d raw %+.3f (p %.4f) | partial %+.3f (p %.3f)'%(lab,len(d),r,p,rp,pp))
print('    record: fixed raw +0.27 [0.11,0.43], partial +0.08 (p 0.36); marginalised -0.02 (constrained 114, posterior medians)')
print('=== (7) marginalised arm, posterior medians (rhoc_sigmaD csvs)')
S=pd.read_csv('rhoc_sigmaD_20261031.csv').set_index('gal'); B=pd.read_csv('rhoc_sigmaD_1b_20261031.csv').set_index('gal')
for tag,drop in (('OLD',[]),('141',six)):
    ok=S.drop([g for g in drop if g in S.index]).dropna(subset=['lrcA','lMb','sb0']); ok=ok[ok.sb0>0]; ls=np.log10(ok.sb0)
    X=np.column_stack([np.ones(len(ok)),ok.lMb,ls]); b=np.linalg.lstsq(X,ok.lrcA.values,rcond=None)[0]
    o=B.drop([g for g in drop if g in B.index]).dropna(); d2=o.lrcA-o.lrc_nom; d1=o.lrcA-o.lrc_1b
    print('   %s N %d: rho_c vs lMb %+.2f (p %.3f); vs log Sigma_b0 %+.2f (p %.3f); vs 2-variable fit %+.2f (coef %+.3f, %+.3f; sd %.2f) | N %d: rho_c vs lud %+.2f (p %.3f), vs zD %+.2f (p %.4f) | RESIDUAL (A - rule_nom) vs lud %+.2f, vs zD %+.2f | (A - rule_1b) vs lud %+.2f, vs zD %+.2f'%(
        tag,len(ok),*spearmanr(ok.lrcA,ok.lMb),*spearmanr(ok.lrcA,ls),spearmanr(ok.lrcA,X@b)[0],b[1],b[2],ok.lrcA.std(),len(o),*spearmanr(o.lrcA,o.lud),*spearmanr(o.lrcA,o.zD),spearmanr(d2,o.lud)[0],spearmanr(d2,o.zD)[0],spearmanr(d1,o.lud)[0],spearmanr(d1,o.zD)[0]))
print('    record: -0.21 (M_b); 0.18 ("Spearman" after the 2-variable fit); -0.37 (lud), -0.29 (zD) quoted for "the residual (A - rule)"')
print('=== (8) matched pairs (marginalised MAP rho_c; M_b, Sigma_h, f_gas within 0.15)')
A=pd.read_csv('marg_fit_params_20261030.csv'); A=A[A.arm=='A'].set_index('gal')[['lrc']]; Z=pd.read_csv('marg_ns_20261030_pergal.csv').set_index('gal')[['logZ_A','logZ_B']]
TT=T[['lMb','lSig','fgas']].copy()
for tag,drop,src in (('OLD 147 MAP',[],'map'),('141 MAP',six,'map'),('141 FIXED V3',six,'fix')):
    D=TT.join(A,how='inner').join(Z,how='inner').dropna(subset=['lSig']); D['dZ']=D.logZ_A-D.logZ_B
    if src=='fix': D['lrc']=tab.lrc.reindex(D.index)
    D=D.drop([g for g in drop if g in D.index]); mt=[]; rn=[]
    for a,b in itertools.combinations(list(D.index),2):
        x,y=D.loc[a],D.loc[b]; d=abs(x.lrc-y.lrc)
        if abs(x.lMb-y.lMb)<=0.15 and abs(x.lSig-y.lSig)<=0.15 and abs(x.fgas-y.fgas)<=0.15: mt.append((d,max(x.dZ,y.dZ),min(x.dZ,y.dZ)))
        else: rn.append(d)
    mt=np.array(mt); k1=mt[mt[:,0]>=0.3]
    print('   %s: galaxies %d; matched %d, random %d; median |dlog rho_c| matched %.3f / random %.3f; pairs >=0.3: %d, >=0.5: %d; of >=0.3 with one member dZ>=+5: %d; clean fork (one dZ<=2, other >=10): %d'%(tag,len(D),len(mt),len(rn),np.median(mt[:,0]),np.median(rn),(mt[:,0]>=0.3).sum(),(mt[:,0]>=0.5).sum(),(k1[:,1]>=5).sum(),((k1[:,1]>=10)&(k1[:,2]<=2)).sum()))
print('    record: K3 0.49 / 0.51; 29 pairs')
print('=== (10) BS-AMEND-1: SPARC fixed arm, V3'); x=tab.lrc.drop(six).dropna()
print('   141: median %.2f; 16-84 half-width %.2f (full %.2f); std %.2f; fits at a scan edge: %d of 147 (%d of 141)'%(x.median(),f(x)/2,f(x),x.std(ddof=0),int(V.edge.sum()),int(V.edge.drop(six).sum())))
