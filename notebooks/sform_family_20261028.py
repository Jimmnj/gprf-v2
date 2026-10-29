# sform_family_20261028.py -- run under PREREG_SFORM_20261028.md (STAMPED JM 2026-10-28 14:44; FROZEN). [UNRATIFIED result]
# Harness pass2_sparc_full_20261022.py unchanged except the global S. Selection as pass3 (147 gal, 2696 pts). Error injection as (15b): 50 realisations, seed 7.
import sys, io, contextlib, warnings, numpy as np, pandas as pd; warnings.filterwarnings('ignore')
from scipy.special import erf
CASE=sys.argv[1]   # 'univ' or 'rel3'
src=open('pass2_sparc_full_20261022.py').read().split("res=[o for o in")[0]
with contextlib.redirect_stdout(io.StringIO()): exec(src)
FORMS={'tanh':lambda x:2*np.tanh(x/2),'hard':lambda x:np.minimum(x,2.0),'erf':lambda x:2*erf(np.sqrt(np.pi)*x/4),
       'exp':lambda x:2*(1-np.exp(-x/2)),'rat':lambda x:x/(1+x/2)}
beta=np.array([11.063,-0.365,0.268,-0.929])
Lsel=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)&(P.vbar2>0); gals=sorted(P[Lsel].gal.unique())
def lrc(g):
    if CASE=='univ': return 7.45
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9
    if not t.Rd>0: return 7.45
    Sh=Mb/(2*np.pi*(1.678*t.Rd*1e3)**2); fg=1.33*t.MHI*1e9/Mb
    return beta[0]+beta[1]*np.log10(Mb)+beta[2]*np.log10(Sh)+beta[3]*fg
la0=np.log10(a0); edges=np.arange(-2.25,1.01,0.25)+la0
win=[(lo,hi) for lo,hi in zip(edges[:-1],edges[1:]) if lo-la0>=-1.751 and hi-la0<=-0.249]
rows=[]; summ={}
for name,Sf in FORMS.items():
    S=Sf; globals()['S']=Sf
    out=[]
    for g in gals:
        o=solve_gal(g,lrc(g))
        if o is None: continue
        out.append(o[(o.ev/o.vobs<0.1)&(o.vobs>0)&(o.vbar2>0)].copy())
    Lk0=pd.concat(out); Lk0['lgb']=np.log10(Lk0.vbar2/Lk0.r*1e6/kpc); Lk0['lgo']=np.log10(Lk0.vobs**2/Lk0.r*1e6/kpc); Lk0['lgm']=np.log10(Lk0.vmod2/Lk0.r*1e6/kpc)
    Lk0['lrar']=np.log10(10**Lk0.lgb/(1-np.exp(-np.sqrt(10**Lk0.lgb/a0)))); Lk0['res']=Lk0.lgo-Lk0.lgm; Lk0['resR']=Lk0.lgo-Lk0.lrar
    rng=np.random.default_rng(7); sc=[]
    for k in range(50):
        Lk=Lk0.copy(); gp={g:(rng.normal(0,1),rng.normal(0,1),rng.normal(0,0.1)) for g in Lk.gal.unique()}; eps=np.array([gp[g] for g in Lk.gal])
        dD=np.log10(1+eps[:,0]*Lk.eD/Lk.D).clip(-0.3,0.3); inc=np.radians(Lk.inc); inc2=np.radians((Lk.inc+eps[:,1]*Lk.einc).clip(20,90))
        di=2*np.log10(np.sin(inc)/np.sin(inc2)); dv=2*rng.normal(0,1,len(Lk))*Lk.ev/Lk.vobs/np.log(10)
        fst=(0.5*Lk.vd**2+0.7*Lk.vb**2)/Lk.vbar2; dML=np.log10(1+fst*(10**eps[:,2]-1))
        x=Lk.lgb+dML; y=Lk.lgm-dD+di+dv
        sc.append([((y-Lk.lrar)[(x>=lo)&(x<hi)]).std() for lo,hi in zip(edges[:-1],edges[1:])])
    sc=np.nanmean(np.array(sc),axis=0)
    for (lo,hi),s_m in zip(zip(edges[:-1],edges[1:]),sc):
        m=(Lk0.lgb>=lo)&(Lk0.lgb<hi)
        if m.sum()<10: continue
        rows.append(dict(form=name,lo=round(lo-la0,2),hi=round(hi-la0,2),N=int(m.sum()),med_res=Lk0.res[m].median(),own_sc=(Lk0.lgm-Lk0.lrar)[m].std(),
                         sc_err=s_m,sc_data=Lk0.resR[m].std(),inwin=bool((lo,hi) in win)))
    summ[name]=dict(rms=np.sqrt((Lk0.res**2).mean()),N=len(Lk0),ngal=Lk0.gal.nunique())
B=pd.DataFrame(rows); B.to_csv('sform_family_bins_%s_20261028.csv'%CASE,index=False)
print('CASE %s; forms %s'%(CASE,list(FORMS)))
print('per-bin median log(g_obs/g_model):'); print(B.pivot(index=['lo','hi','N'],columns='form',values='med_res')[list(FORMS)].round(3).to_string())
print('\nwindow: model own scatter | model+errors | data'); W=B[B.inwin]
print(W.pivot(index=['lo','hi'],columns='form',values='own_sc')[list(FORMS)].round(3).to_string())
print(W.pivot(index=['lo','hi'],columns='form',values='sc_err')[list(FORMS)].round(3).to_string())
print('data sd:',W[W.form=='tanh'].sc_data.round(3).tolist())
print('\nform  N    rms    M1(max|dmed|)  M3(max|d own sc|)  M4 bins with (err-data)<=-0.01 [of %d]   D1'%len(win))
base=B[B.form=='tanh'].set_index(['lo','hi'])
for f in FORMS:
    b=B[B.form==f].set_index(['lo','hi']); M1=(b.med_res-base.med_res).abs().max()
    w=b[b.inwin]; M3=(w.own_sc-base[base.inwin].own_sc).abs().max(); M4=((w.sc_err-w.sc_data)<=-0.01).sum()
    print('%-5s %d %.4f   %.3f          %.3f              %d                                %s'%(f,summ[f]['N'],summ[f]['rms'],M1,M3,M4,'ok' if (M1<=0.03 and M3<=0.02) else 'BREAKS'))
