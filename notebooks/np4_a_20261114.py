# np4_a_20261114.py [UNRATIFIED] Item 4 (numbers pass on V3 / 141), part A: B.1 class medians, quality splits, constrained-subset spread, bulge flag.
# JM 2026-11-14 17:34 "f1 accepted and start 4". Post-processing only; no fits. Inputs: marg_ns_20261030.jsonl (arm A, sig_int 0.03), marg_fit_params_20261030.csv,
# clean_conc_part_20261109.csv (old fix0), ep3_cc_20261114_part.csv (V3 fix0), rhoc_corr_20261105_table.csv (Vf, fD, MHI), ep3_val_20261114_part.csv (V3, 147).
# Each block prints OLD (must reproduce RESULT_CLEAN_CONC 5b/5e, RESULT_RHOC_CONSTRAINED secs 1, 6) then NEW (V3 fixed arm; EP3-4 cut: the six removed).
# Definitions recovered by matching the record: gas-rich = catalogue f_gas > 0.5; gas-poor no-bulge = f_gas < 0.2; MW test vs no-bulge V_f > 0; the record's Ns counted CamB (FIX NaN).
# Spread: full = p84 - p16; rob = full/2 (the note's "robust spread", EP-4).
import json,numpy as np,pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
six=['CamB','D512-2','NGC4068','NGC6789','UGC05750','UGC07232']
J=[json.loads(l) for l in open('marg_ns_20261030.jsonl')]; J=pd.DataFrame([dict(gal=r['gal'],med=r['post']['lrc'][1],w=r['post']['lrc'][2]-r['post']['lrc'][0],bul=r['bul'],npts=r['npts'],eDD=r['eD']/r['D0'],ei=r['ei']) for r in J if r['arm']=='A' and r['sig_int']==0.03]).set_index('gal')
P=pd.read_csv('marg_fit_params_20261030.csv'); P=P[P.arm=='A'].set_index('gal')
C=pd.read_csv('clean_conc_part_20261109.csv').drop_duplicates('gal').set_index('gal'); E=pd.read_csv('ep3_cc_20261114_part.csv').set_index('key')
T=pd.read_csv('rhoc_corr_20261105_table.csv').set_index('gal'); V=pd.read_csv('ep3_val_20261114_part.csv').set_index('key')
f=lambda x:(np.percentile(x,84)-np.percentile(x,16))
print('=== (1) RHOC_CONSTRAINED sec 1: marginalised arm, posterior width and spread of medians')
for lab,D in (('OLD 147',J),('NEW 141',J.drop(six))):
    print('  %s: width p50 %.2f; '%(lab,D.w.median())+'; '.join('w<%s N %d full %.2f rob %.2f std %.2f'%(c,(D.w<c).sum(),f(D.med[D.w<c]),f(D.med[D.w<c])/2,D.med[D.w<c].std(ddof=0)) for c in (0.3,0.5,9)))
print('    record: 114 of 147 < 0.5 (1.01, std 0.50); 77 < 0.3 (0.89); all 1.10, std 0.56; width p50 0.28')
print('=== (2) bulge flag, marginalised, constrained (w<0.5): residual after the M_b trend')
for lab,D in (('OLD',J),('NEW',J.drop(six))):
    D=D[D.w<0.5].join(P[['Mb']]); lM=np.log10(D.Mb); a=np.polyfit(lM,D.med,1); r=D.med-np.polyval(a,lM)
    X=np.column_stack([np.ones(len(D)),lM,D.bul.astype(float)]); b=np.linalg.lstsq(X,D.med.values,rcond=None)[0]; rs=D.med.values-X@b; cov=np.linalg.inv(X.T@X)*(rs@rs)/(len(D)-3)
    print('  %s N %d: bulge N %d median %.2f | no bulge N %d median %.2f | raw MW p %.3f | residual MEDIANS bulge %+.2f, none %+.2f -> %+.2f (record method) | residual MEANS %+.2f | MW on residuals p %.2f | OLS bulge term at fixed M_b %+.2f +- %.2f'%(lab,len(D),D.bul.sum(),D.med[D.bul].median(),(~D.bul).sum(),D.med[~D.bul].median(),mannwhitneyu(D.med[D.bul],D.med[~D.bul]).pvalue,r[D.bul].median(),r[~D.bul].median(),r[D.bul].median()-r[~D.bul].median(),r[D.bul].mean()-r[~D.bul].mean(),mannwhitneyu(r[D.bul],r[~D.bul]).pvalue,b[2],np.sqrt(cov[2,2])))
print('    record: bulge 29 7.25; none 85 7.61; p 0.03; -0.10 / +0.06 -> -0.16')
R=C.join(J[['bul','eDD','ei','w']]).join(T[['Vf','fD','MHI']]).join(P[['Mb']].rename(columns={'Mb':'MbP'}))
R['fg']=T.fgas.reindex(R.index)   # catalogue f_gas (reproduces the record's Ns)
def classes(R,col,fgcol,thr):
    fd=(~R.bul)&(R.Vf>0)&(R.lMb>=9.5)
    out=[('ALL',R.index==R.index),('bulge',R.bul),('no-bulge',~R.bul),('no-bulge flat discs lMb>=9.5',fd),('dwarfs lMb<9.5',R.lMb<9.5),('gas-rich dwarfs',(R.lMb<9.5)&(R[fgcol]>thr)),
         ('no V_f slow risers',~(R.Vf>0)),('no V_f & lMb<9.5',(~(R.Vf>0))&(R.lMb<9.5)),('gas-rich any mass',R[fgcol]>thr),('gas-poor any mass',R[fgcol]<=thr),('gas-poor no-bulge',(R[fgcol]<0.2)&(~R.bul)),
         ('r_peak/Rd < 1.5',R.pk0<1.5),('r_peak/Rd 1.5-3',(R.pk0>=1.5)&(R.pk0<3)),('r_peak/Rd > 3',R.pk0>=3),
         ('lMb 7-8.5',(R.lMb<8.5)),('lMb 8.5-9.5',(R.lMb>=8.5)&(R.lMb<9.5)),('lMb 9.5-10.5',(R.lMb>=9.5)&(R.lMb<10.5)),('lMb 10.5-12',R.lMb>=10.5)]
    return out,fd
if __name__=='__main__':
    import sys
    fgcol='fg'; thr=0.5
    for tag in ('OLD','NEW'):
        D=R.copy()
        if tag=='NEW': D['lrc_fix0']=E.fix0_v3.reindex(D.index); D=D.drop([g for g in six if g in D.index])
        D=D[D.lrc_fix0.notna()]
        cl,fd=classes(D,'lrc_fix0',fgcol,thr); fd2=(~D.bul)&(D.Vf>0)   # record's 'flat discs' for the MW test
        print('=== (3) class medians %s (N %d; fgas col %s > %s): MARG | FIX | FIX - 7.42'%(tag,len(D),fgcol,thr))
        for lab,m in cl:
            m=np.asarray(m); print('   %-30s N %3d  %.2f | %.2f | %+.2f   rob spread FIX %.2f'%(lab,m.sum(),D.lrc_marg[m].median(),D.lrc_fix0[m].median(),D.lrc_fix0[m].median()-7.42,f(D.lrc_fix0[m])/2 if m.sum()>5 else np.nan))
        sr=~(D.Vf>0); print('   slow risers vs flat discs MW p: FIX %.2f MARG %.2f | Spearman FIX vs lMb %+.2f (p %.2f), vs fgas %+.2f (p %.2f)'%(mannwhitneyu(D.lrc_fix0[sr],D.lrc_fix0[fd2]).pvalue,mannwhitneyu(D.lrc_marg[sr],D.lrc_marg[fd2]).pvalue,*spearmanr(D.lMb,D.lrc_fix0),*spearmanr(D[fgcol],D.lrc_fix0)))
        print('=== (4) quality splits %s: FIX full/rob/std | MARG full/rob/std'%tag)
        Q=pd.read_csv('clean_conc_part_20261109.csv').drop_duplicates('gal').set_index('gal').Q.reindex(D.index)
        for lab,m in (('all',Q>0),('Q=1',Q==1),('Q=2',Q==2),('eD/D<0.12',D.eDD<0.12),('eD/D>=0.25',D.eDD>=0.25),('direct D',D.fD!=1),('Hubble-flow D',D.fD==1),('Q1 & eD/D<0.12 & ei<=5',(Q==1)&(D.eDD<0.12)&(D.ei<=5)),('npts>=15',D.npts>=15)):
            m=np.asarray(m); print('   %-24s N %3d  %.2f/%.2f/%.2f | %.2f/%.2f/%.2f'%(lab,m.sum(),f(D.lrc_fix0[m]),f(D.lrc_fix0[m])/2,D.lrc_fix0[m].std(ddof=0),f(D.lrc_marg[m]),f(D.lrc_marg[m])/2,D.lrc_marg[m].std(ddof=0)))
    print('    record 5b: ALL 7.55|7.43; bulge(29) 7.25|7.48; flat discs(43) 7.57|7.42; dwarfs(38) 7.73|7.46; gas-rich dwarfs(32) 7.66|7.40; slow(13) 7.25|7.15 p .21/.16; noVf&low(9) 7.24|7.45; gas-rich(45) 7.69|7.33; gas-poor no-bulge(14) 7.39|7.48; Sp lMb -0.14 (.14), fgas -0.11 (.26)')
    print('    record 5e FIX|MARG: Q1(78) 0.94/0.52|0.86/0.50; Q2(35) 1.11/0.73|0.97/0.51; eD<.12(41) 1.13/0.55|1.07/0.50; eD>=.25(35) 1.15/0.60|0.92/0.46; direct(50) 0.95/0.63|0.89/0.44; HF(63) 1.12/0.56|1.05/0.54; clean19 0.86/0.55|0.97/0.54; npts>=15(71) 0.97/0.62|0.88/0.48')
