# ep3_n141_20261114.py -- EP3-4 RULED (JM 2026-11-14 02:18 "accepted"): the fits with < 3 fitted points (pass2 mask) are cut from every
# statistic of the fitted rho_c. Recomputes those statistics on the 141. Read-only on files of record; no solver call. [UNRATIFIED numbers]
# Inputs: pass2_table2_part.csv (V3, md5 82dce260), marg_fit_params_20261030.csv (arm A MAP), rhoc_corr_20261105_table.csv (baryon columns only),
# errprop2/errprop3 part csvs (tau), EP-3 scratch outputs under ~/ep3 (V3 runs of cdr_test, axisym_r1, rhoc_env), test1 + Burkert files (sec 3.3).
# EXPECTATIONS (written before the first run): fixed-arm median ~7.42, robust ~0.47 (EP-3 validation run); tau moves < 0.03; env partial slope
# stays negative; sec 3.3 numbers unchanged by the cut (none of the six is in TEST1); axisym unchanged (its solver already refuses < 3 points).
import os, io, contextlib, ast, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from scipy import odr
H=os.path.expanduser('~/ep3/'); SIX=['CamB','D512-2','NGC4068','NGC6789','UGC05750','UGC07232']
hw=lambda x:(np.percentile(np.asarray(x),84)-np.percentile(np.asarray(x),16))/2
V=pd.read_csv('pass2_table2_part.csv').set_index('gal')
ns0={}
with contextlib.redirect_stdout(io.StringIO()): exec(open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0],ns0)
P0=ns0['P']; nfit={g:int(((P0.gal==g)&(P0.ev/P0.vobs<0.1)&(P0.vobs>0)&(P0.vbar2>0)).sum()) for g in V.index}
cut=sorted(g for g,n in nfit.items() if n<3); assert cut==sorted(SIX), cut
print('0. cut (n_fit < 3): '+'; '.join('%s %d pt, lrc %.2f, chi2 %.2g'%(g,nfit[g],V.loc[g,'lrc'],V.loc[g,'chi_law']) for g in cut))
G=pd.read_csv('rhoc_corr_20261105_table.csv').set_index('gal'); assert abs(G.loc['NGC2841','Rd']-3.64)<1e-6
D=G[['lMb','lSig','fgas','fD','Q','btfr_res','lSBdisk']].copy(); D['fix']=V.lrc.reindex(D.index)
mg=pd.read_csv('marg_fit_params_20261030.csv'); D['marg']=mg[mg.arm=='A'].set_index('gal').lrc.reindex(D.index)
assert len(D)==147 and D.fix.notna().all() and D.marg.notna().all()
def ols(d,y,cols,B=2000,seed=0):
    X=np.column_stack([np.ones(len(d))]+[d[c].values for c in cols]); Y=d[y].values
    b=np.linalg.lstsq(X,Y,rcond=None)[0]; r=Y-X@b; rg=np.random.default_rng(seed); bs=[]
    for _ in range(B):
        i=rg.integers(0,len(d),len(d)); bs.append(np.linalg.lstsq(X[i],Y[i],rcond=None)[0])
    lo,hi=np.percentile(bs,[2.5,97.5],axis=0); return b,lo,hi,r
print('\n1. spread and slopes, both arms (std ddof 0)')
for arm in ('fix','marg'):
    for lab,d in (('147',D),('141',D.drop(SIX))):
        y=d[arm].values; print('  %-4s %s: median %.3f robust %.3f std %.3f max %.2f'%(arm,lab,np.median(y),hw(y),y.std(),y.max()))
        b,lo,hi,_=ols(d,arm,['lMb']); print('       marginal mass slope %+.3f [%+.2f,%+.2f] over %.1f dex'%(b[1],lo[1],hi[1],d.lMb.max()-d.lMb.min()))
        b,lo,hi,r=ols(d,arm,['lMb','lSig']); print('       2-var: lMb %+.3f lSig %+.3f [%+.2f,%+.2f] | resid robust %.3f std %.3f'%(b[1],b[2],lo[2],hi[2],hw(r),r.std()))
        b,lo,hi,r=ols(d,arm,['lMb','lSig','fgas']); print('       3-var: const %.3f lMb %+.3f [%+.2f,%+.2f] lSig %+.3f [%+.2f,%+.2f] fgas %+.3f [%+.2f,%+.2f] | resid robust %.3f std %.3f'%(b[0],b[1],lo[1],hi[1],b[2],lo[2],hi[2],b[3],lo[3],hi[3],hw(r),r.std()))
print('\n2. ERRPROP tau (BD-3 table), 147 vs 141')
ns={}
with contextlib.redirect_stdout(io.StringIO()): exec(open('errprop3_analysis_20261113.py').read(),ns)
Bt=ns['B']; print(Bt.loc[SIX,['lrc0','sD','sI','sU','stot','flag','fD']].round(3).to_string())
for lab,R in (('147',Bt),('141',Bt.drop(SIX))): ns['report'](lab,R)
print('\n3. rho_c vs Chae e_N,env (V3 table): rhoc_env Q1-Q3 code verbatim, rng re-seeded per run')
E_=pd.read_csv(H+'v3/notebooks/rhoc_env_20261105.csv',index_col=0); assert (E_.lrc-V.lrc.reindex(E_.index)).abs().max()<1e-9
def envrun(D,seed=20261105):
    rng=np.random.default_rng(seed)
    def boot(fn,n=1000):
        v=[fn(D.iloc[rng.integers(0,len(D),len(D))]) for _ in range(n)]; return np.percentile(v,[2.5,97.5])
    def ols1(d,cols):
        X=np.column_stack([np.ones(len(d))]+[d[c].values for c in cols]); return np.linalg.lstsq(X,d.lrc.values,rcond=None)[0][1]
    for col in ('log_eN_env_maxclust','log_eN_env_noclust'):
        r,p=spearmanr(D.lrc,D[col]); ci=boot(lambda d: spearmanr(d.lrc,d[col])[0])
        print('  N %d %-20s Q1 Spearman %+.3f (p %.3f) [%+.3f,%+.3f]'%(len(D),col,r,p,*ci))
        for cols,lab in (([col],'alone'),([col,'lMb','lSig','fgas'],'+ M_b, Sigma_h, f_gas')):
            b=ols1(D,cols); ci=boot(lambda d: ols1(d,cols)); print('       Q2 slope (%s) %+.3f [%+.3f,%+.3f]'%(lab,b,*ci))
        q=D[col].quantile([1/3,2/3]); lo=D[D[col]<=q.iloc[0]].lrc; hi=D[D[col]>=q.iloc[1]].lrc
        ci=boot(lambda d: d[d[col]>=d[col].quantile(2/3)].lrc.median()-d[d[col]<=d[col].quantile(1/3)].lrc.median())
        print('       Q3 tercile contrast %+.3f [%+.3f,%+.3f]'%(hi.median()-lo.median(),*ci))
envrun(E_); envrun(E_.drop([g for g in SIX if g in E_.index]))
print('\n4. central density relation (cdr_test eq10 + fit verbatim), V3 run points')
def eq10(p,ls): lS0,lSc,beta=p; x=10**ls/10**lSc; return lS0+(1-beta)*np.log10(1+x)+beta*np.log10(x)
def fit(ls,ld):
    d=odr.RealData(ls,ld,sx=0.14,sy=0.14); o=odr.ODR(d,odr.Model(eq10),beta0=[np.log10(1468),np.log10(1571),0.62]).run()
    res=ld-eq10(o.beta,ls); return o.beta,o.sd_beta,np.sqrt(np.mean(res**2))
Pc=pd.read_csv(H+'single/notebooks/cdr_points_20261030.csv')
for ex in ([],SIX):
    for arm in ('OBS','A0','U0','R0','B0'):
        a=Pc[(Pc.arm==arm)&np.isfinite(Pc.Sdyn)&(Pc.Sdyn>0)&~Pc.gal.isin(ex)]; b,sb,rms=fit(np.log10(a.Sstar),np.log10(a.Sdyn))
        print('  %-3s N %d beta %.3f +- %.3f resid rms %.3f'%(arm,len(a),b[2],sb[2],rms))
print('\n5. disc geometry (axisym_r1 V3 scratch files): fits of the six')
for mode in ('C0','full'):
    A=pd.read_csv(H+'single/notebooks/axisym_r1_%s_gals_part.csv'%mode).drop_duplicates('gal',keep='last').set_index('gal'); q=A[(A.Q==1)].lrc.dropna()
    print('  %-4s six: %s | Q1 N %d median %.3f robust %.3f std %.3f'%(mode,A.lrc.reindex(SIX).round(2).to_dict(),len(q),q.median(),hw(q),q.std(ddof=0)))
print('\n6. Table 2 per-point rms, with and without the six')
for col,lab in (('law','eq. (1), rho_c free'),('rfree','RAR fn, M/L free'),('rfix','RAR fn, fixed'),('r3','RAR fn, M/L D i free')):
    for ex,l2 in (([],'147'),(SIX,'141')):
        r=np.concatenate([np.array(ast.literal_eval(V.loc[g,col]),float) for g in V.index if g not in ex])
        print('  %-21s %s: points %d rms %.4f'%(lab,l2,len(r),np.sqrt(np.mean(r**2))))
print('\n7. sec 3.3 / Fig 4: fig4 statistics verbatim; v3 = TEST1 galaxies in the table of record, V3 values; 2-var fit = fig4_mu0_regress method')
src4=open('fig4_mu0_mond_20261027.py').read().split("fig,ax=plt.subplots")[0]
key="F=pd.read_csv('test1_capped_baryon_20261018.csv'); T=load_table()"; assert src4.count(key)==1
for MODE in ('t1','v3'):
    patch=key+"\nV3=pd.read_csv('pass2_table2_part.csv').set_index('gal').lrc\nif MODE=='v3':\n    F=F[F.gal.isin(V3.index)].reset_index(drop=True); F['l_tanh2']=V3.loc[F.gal].values\n"
    kB="B=pd.read_csv('test1c_burkert_data_20261018.csv',index_col=0); B=B[~B.r0_gt_rmax.astype(bool)]"; assert src4.count(kB)==1
    ns4={'MODE':MODE,'__name__':'x'}; print('  ==',MODE)
    exec(src4.replace(key,patch).replace(kB,kB+("; B=B[B.index.isin(F.gal)]" if MODE=='v3' else '')),ns4)  # v3: Burkert set restricted to the plotted galaxies
    F4=ns4['F']; Dd=pd.DataFrame(dict(gal=F4.gal,lMb=np.log10(F4.Mb.values),lSh=np.log10(F4.Sh.values),mu_t=ns4['mu_law'],mu_m=ns4['mu_mond'])).set_index('gal')
    B4=pd.read_csv('test1c_burkert_data_20261018.csv',index_col=0); B4=B4[~B4.r0_gt_rmax.astype(bool)]; B4=B4[B4.index.isin(Dd.index)]
    rg=np.random.default_rng(1)
    def fit2(x1,x2,y):
        X=np.c_[np.ones_like(x1),x1,x2]; c=np.linalg.lstsq(X,y,rcond=None)[0]; bs=[]
        for _ in range(2000):
            i=rg.integers(0,len(y),len(y)); bs.append(np.linalg.lstsq(X[i],y[i],rcond=None)[0])
        return c,np.std(bs,axis=0)
    for k,y in (('data (Burkert)',B4.mu_data.values),('this work',Dd.loc[B4.index,'mu_t'].values),('MOND',Dd.loc[B4.index,'mu_m'].values)):
        c,e=fit2(Dd.loc[B4.index,'lMb'].values,Dd.loc[B4.index,'lSh'].values,y); print('     2-var %-15s N=%3d  M_b^(%+.2f+-%.2f) Sigma_h^(%+.2f+-%.2f)'%(k,len(y),c[1],e[1],c[2],e[2]))
print('\n8. clean_conc without CamB (the only one of the six in its 112), V3 refits')
Cc=pd.read_csv('clean_conc_part_20261109.csv').drop_duplicates('gal').set_index('gal'); Ec=pd.read_csv('ep3_cc_20261114_part.csv').set_index('key')
Cc['lrc_fix0']=Ec.fix0_v3; Cc['lrc_fix2']=Ec.fix2_v3; Cc=Cc.drop('CamB')
def resid(y,x): X=np.column_stack([np.ones(len(y)),x]); return y-X@np.linalg.lstsq(X,y,rcond=None)[0]
for lab,d,y,x in (('all Q r_peak/Rd r[2:], fix2',Cc,'lrc_fix2','pk2'),('all Q r_half/Rd r[2:], fix2',Cc,'lrc_fix2','hm2'),('Q1 r_peak/Rd r[2:], fix2',Cc[Cc.Q==1],'lrc_fix2','pk2')):
    dd=d[[y,x,'lMb']].dropna(); r,p=spearmanr(np.log10(dd[x]),dd[y]); rp,pp=spearmanr(resid(dd[y].values,dd.lMb.values),resid(np.log10(dd[x].values),dd.lMb.values))
    print('  %-30s N %3d rho %+.2f (p %.3f) | fixed M_b %+.2f (p %.3f)'%(lab,len(dd),r,p,rp,pp))
print('\n9. rhoc_corr: the two Holm-significant variables, fixed arm V3')
for lab,d in (('147',D),('141',D.drop(SIX))):
    for v in ('btfr_res','lSBdisk'):
        dd=d[['fix',v]].dropna(); r,p=spearmanr(dd.fix,dd[v]); print('  %s %-9s N %3d Spearman %+.3f (p %.1e)'%(lab,v,len(dd),r,p))
