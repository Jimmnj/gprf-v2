# pass2_sparc_full_20261022.py [UNRATIFIED] 2026-10-22 17:10 (JM 'yes' to Pass 2; simulated Lelli review items 13-15).
# (13) Lelli sample: Q<=2, inc>=30, point cut dV/V<0.1 (McGaugh+2016/Lelli+2017 RAR selection, reported 153 galaxies / 2693 points) vs ours (Q<=2, inc>30, V_f>0).
#      Solve EVERY galaxy at universal rho_c = 10^7.45 (sigma from M_b, no V_f needed). Binned RAR vs data for (all) and (galaxies without V_f).
# (14) V(R_last): model/data at the last point vs log Sigma_h and log M_b; also V_model(R_last)/V_f vs V_obs(R_last)/V_f.
# (15) residual correlations (per point) with log r, log Sigma_h, f_gas, log M_b; error-injected model scatter vs data scatter per bin.
import numpy as np, pandas as pd, glob, os, sys, warnings; warnings.filterwarnings('ignore')
from scipy.integrate import solve_ivp
from scipy.stats import spearmanr
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2)
MRT='../sparc_galaxy_data/SPARC_Lelli2016c.mrt'; SD='../sparc_galaxy_data/sparc_database'
T=[]
for line in open(MRT):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: T.append(dict(gal=p[0],D=float(p[2]),eD=float(p[3]),inc=float(p[5]),einc=float(p[6]),L=float(p[7]),Rd=float(p[11]),MHI=float(p[13]),Vf=float(p[15]),Q=int(p[17])))
        except ValueError: pass
T=pd.DataFrame(T).set_index('gal')
rows=[]
for f in glob.glob(os.path.join(SD,'*_rotmod.dat')):
    g=os.path.basename(f).replace('_rotmod.dat','')
    if g not in T.index: continue
    d=np.atleast_2d(np.loadtxt(f,comments='#'))
    rows.append(pd.DataFrame(dict(gal=g,r=d[:,0],vobs=d[:,1],ev=d[:,2],vg=d[:,3],vd=d[:,4],vb=d[:,5])))
P=pd.concat(rows).join(T,on='gal')
P['vbar2']=np.sign(P.vg)*P.vg**2+0.5*P.vd**2+0.7*P.vb**2
selL=(P.Q<=2)&(P.inc>=30)&(P.ev/P.vobs<0.1)&(P.vobs>0)
selO=(P.Q<=2)&(P.inc>30)&(P.vbar2>0)&(P.vobs>0)
for lab,m in (('Lelli-style (Q<=2, inc>=30, dV/V<0.1)',selL),('ours before V_f cut (Q<=2, inc>30, vbar2>0)',selO),('Lelli-style & vbar2>0',selL&(P.vbar2>0))):
    print('%-48s galaxies %3d points %4d'%(lab,P[m].gal.nunique(),m.sum()))
print('galaxies with Q<=2, inc>=30:',T[(T.Q<=2)&(T.inc>=30)].shape[0],' of which V_f>0:',T[(T.Q<=2)&(T.inc>=30)&(T.Vf>0)].shape[0])
gals=sorted(T[(T.Q<=2)&(T.inc>=30)].index)
def solve_gal(g,lrc=7.45):
    d=P[P.gal==g].sort_values('r'); r=d.r.values; vb2=np.clip(d.vbar2.values,0,None); rl=r[-1]
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2])
    vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    t=T.loc[g]; Mb=0.5*t.L*1e9+1.33*t.MHI*1e9; s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**lrc
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200)
    if not s.success or s.y.shape[1]!=len(r): return None
    out=d.copy(); out['vmod2']=vb2+G*s.y[2]/r; out['Mb']=Mb
    out['Sh']=Mb/(2*np.pi*(1.678*t.Rd*1e3)**2) if t.Rd>0 else np.nan; out['fgas']=1.33*t.MHI*1e9/Mb
    return out
res=[o for o in (solve_gal(g) for g in gals) if o is not None]
A=pd.concat(res); A.to_csv('pass2_sparc_full_points_20261022.csv',index=False)
A=A[(A.vbar2>0)]
A['lgb']=np.log10(A.vbar2/A.r*1e6/kpc); A['lgo']=np.log10(A.vobs**2/A.r*1e6/kpc); A['lgm']=np.log10(A.vmod2/A.r*1e6/kpc)
A['lrar']=np.log10(10**A.lgb/(1-np.exp(-np.sqrt(10**A.lgb/a0))))
L=A[(A.ev/A.vobs<0.1)]
print('\nsolved galaxies %d; Lelli-style points used %d in %d galaxies'%(A.gal.nunique(),len(L),L.gal.nunique()))
la0=np.log10(a0); edges=np.arange(-3,1.01,0.25)+la0
def table(D,lab):
    print('\n'+lab+'\n  lg(gb/a0)    N  data-model(med)  data-RARfn(med)  scat(data-model)  scat(data-RARfn)')
    for lo,hi in zip(edges[:-1],edges[1:]):
        m=(D.lgb>=lo)&(D.lgb<hi)
        if m.sum()<10: continue
        x=(D.lgo-D.lgm)[m]; y=(D.lgo-D.lrar)[m]
        print('  %+5.2f..%+5.2f %4d   %+.3f            %+.3f            %.3f             %.3f'%(lo-la0,hi-la0,m.sum(),x.median(),y.median(),x.std(),y.std()))
    x=D.lgo-D.lgm; y=D.lgo-D.lrar; print('  ALL: rms data-model %.3f, data-RARfn %.3f; median %+.3f / %+.3f'%(np.sqrt((x**2).mean()),np.sqrt((y**2).mean()),x.median(),y.median()))
table(L,'(13a) ALL galaxies, Lelli-style points, universal rho_c 10^7.45')
nof=L[L.Vf<=0]; print('\ngalaxies without V_f in sample: %d (%d points)'%(nof.gal.nunique(),len(nof)))
table(nof,'(13b) galaxies WITHOUT V_f (dropped before)')
# (14)
G2=[]
for g,d in L.groupby('gal'):
    d=d.sort_values('r'); last=d.iloc[-1]
    G2.append(dict(gal=g,Mb=last.Mb,Sh=last.Sh,fgas=last.fgas,Vf=last.Vf,ratio=np.sqrt(last.vmod2)/last.vobs,
                   vm_vf=np.sqrt(last.vmod2)/last.Vf if last.Vf>0 else np.nan, vo_vf=last.vobs/last.Vf if last.Vf>0 else np.nan))
G2=pd.DataFrame(G2); G2.to_csv('pass2_vlast_20261022.csv',index=False)
print('\n(14) V_model(R_last)/V_obs(R_last), universal rho_c: median %.3f, 16-84%% %.3f-%.3f, N=%d'%(G2.ratio.median(),*G2.ratio.quantile([.16,.84]),len(G2)))
for col in ('Sh','Mb','fgas'):
    x=np.log10(G2[col]) if col!='fgas' else G2[col]; ok=np.isfinite(x)&np.isfinite(G2.ratio)
    rho,pv=spearmanr(x[ok],np.log10(G2.ratio[ok])); print('    vs %-5s Spearman %+.2f (p %.1e)'%(col,rho,pv))
qs=pd.qcut(np.log10(G2.Sh),4); print('    by Sigma_h quartile (median ratio):',G2.groupby(qs).ratio.median().round(3).to_dict())
q2=pd.qcut(np.log10(G2.Mb),4); print('    by M_b quartile (median ratio):',G2.groupby(q2).ratio.median().round(3).to_dict())
h=G2.dropna(subset=['vm_vf']); print('    galaxies with V_f: V_model(R_last)/V_f median %.3f vs V_obs(R_last)/V_f median %.3f (N=%d)'%(h.vm_vf.median(),h.vo_vf.median(),len(h)))
# (15) residual correlations per point
L=L.copy(); L['res']=L.lgo-L.lgm; L['resR']=L.lgo-L.lrar
print('\n(15a) per-point residual correlations (Spearman), law (univ rho_c) vs RAR function:')
for lab,x in (('log r',np.log10(L.r)),('log Sigma_h',np.log10(L.Sh)),('f_gas',L.fgas),('log M_b',np.log10(L.Mb))):
    ok=np.isfinite(x); print('    %-12s law %+.2f   RARfn %+.2f'%(lab,spearmanr(x[ok],L.res[ok])[0],spearmanr(x[ok],L.resR[ok])[0]))
# (15b) error injection
rng=np.random.default_rng(7); N=50; sc=[]
for k in range(N):
    Lk=L.copy()
    gpert={g:(rng.normal(0,1),rng.normal(0,1),rng.normal(0,0.1)) for g in Lk.gal.unique()}
    eps=np.array([gpert[g] for g in Lk.gal])
    dD=np.log10(1+eps[:,0]*Lk.eD/Lk.D).clip(-0.3,0.3)          # g_obs ∝ 1/D at fixed angles
    inc=np.radians(Lk.inc); inc2=np.radians((Lk.inc+eps[:,1]*Lk.einc).clip(20,90))
    di=2*np.log10(np.sin(inc)/np.sin(inc2))                     # g_obs ∝ 1/sin^2 i
    dv=2*rng.normal(0,1,len(Lk))*Lk.ev/Lk.vobs/np.log(10)
    fst=(0.5*Lk.vd**2+0.7*Lk.vb**2)/Lk.vbar2; dML=np.log10(1+fst*(10**eps[:,2]-1))  # M/L error moves g_bar
    Lk['x']=Lk.lgb+dML; Lk['y']=Lk.lgm-dD+di+dv
    sc.append([((Lk.y-Lk.lrar)[(Lk.x>=lo)&(Lk.x<hi)]).std() for lo,hi in zip(edges[:-1],edges[1:])])
sc=np.nanmean(np.array(sc),axis=0)
print('\n(15b) scatter of log(g/g_RARfn) per bin: data vs law+injected errors (D, i, 0.1 dex M/L, velocity errors; 50 realisations)')
for (lo,hi),s_m in zip(zip(edges[:-1],edges[1:]),sc):
    m=(L.lgb>=lo)&(L.lgb<hi)
    if m.sum()<10: continue
    print('  %+5.2f..%+5.2f  N %4d   data %.3f   law+errors %.3f   law alone %.3f'%(lo-la0,hi-la0,m.sum(),(L.resR[m]).std(),s_m,(L.lgm-L.lrar)[m].std()))
