# dms_fb_20261105.py [RATIFIED 2026-11-05 21:31; T1 scored as the velocity ratio V_b/V_c, the DMS definition; T2 void] -- PREREG_DMS_20261105: eq. (1) at universal rho_c vs DiskMass F_b(2.2 h_R) and BTFR.
import re, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from scipy.integrate import solve_ivp
from scipy.special import i0,k0,i1,k1
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2); LRC=7.45
B='../data/diskmass/'
def rows(path):
    out=[]
    for l in open(path,encoding='utf-8',errors='ignore'):
        if l.lstrip().startswith('%') or '&' not in l: continue
        c=[re.sub(r'\\(it|bf|phantom\{[^}]*\})','',x).replace('{','').replace('}','').replace('$','').strip() for x in l.replace('\\\\','').split('&')]
        if re.match(r'^\d{3,5}$',c[0]): out.append(c)
    return out
def num(s):  # first number in a cell like "10.25_{-0.51}^{+0.23}" or "3.9 \pm 0.2" or "-23.93"
    s=s.replace('−','-').replace('$-$','-'); m=re.search(r'[-+]?\d+\.?\d*',s); return float(m.group(0)) if m else np.nan
M={r[0]:dict(lMd=num(r[1]),lMbul=num(r[2]),lMat=num(r[3]),lMmol=num(r[4]),lMb=num(r[5])) for r in rows(B+'arxiv_1308.0336_src/tables/Chap5_tabMasses.tex')}
H={}
for r in rows(B+'arxiv_1308.0336_src/tables/tab_hRhz.tex'):
    H[r[0]]=dict(hR=num(r[1]),hz=num(r[2]))
    if len(r)>=6 and re.match(r'^\d{3,5}$',r[3]): H[r[3]]=dict(hR=num(r[4]),hz=num(r[5]))
F={}
for r in rows(B+'arxiv_1308.0336_src/tables/Chap5_tabF.tex'):
    F[r[0]]=dict(Fb=num(r[1]),eFb=num(r[1].split('pm')[-1]) if 'pm' in r[1] else np.nan,Fbmax=num(r[2]))
    if len(r)>=6 and re.match(r'^\d{3,5}$',r[3]): F[r[3]]=dict(Fb=num(r[4]),eFb=num(r[4].split('pm')[-1]) if 'pm' in r[4] else np.nan,Fbmax=num(r[5]))
DM={r[0]:dict(UpsK=num(r[1]),Rc=num(r[3]),lrho0=num(r[4])) for r in rows(B+'arxiv_1308.0336_src/tables/Chap5_tabDM.tex')}
V={r[0]:dict(VflatTF=num(r[5]),Vcsini=num(r[6]),iTF=num(r[7]),Vstar=num(r[8]),Voiii=num(r[10]) if len(r)>10 else np.nan) for r in rows(B+'arxiv_1307.8130_src/tables/Chap2_tab5.tex')}
D=pd.DataFrame(M).T.join(pd.DataFrame(H).T).join(pd.DataFrame(F).T).join(pd.DataFrame(DM).T).join(pd.DataFrame(V).T).astype(float)
print('galaxies with masses %d, hR %d, Fb %d, V %d'%(len(M),len(H),len(F),len(V))); print(D[['lMb','hR','Fb','Voiii','iTF','UpsK']].describe().loc[['count','50%']].round(2).to_string())
def vexp2(r,Mass,l):
    y=np.clip(r/(2*l),1e-6,None); return 4*np.pi*G*(Mass/(2*np.pi*l**2))*l*y**2*(i0(y)*k0(y)-i1(y)*k1(y))
out=[]
for g,t in D.iterrows():
    if not np.isfinite(t.lMb) or not np.isfinite(t.hR): continue
    Md,Mbul,Mat,Mmol=10**t.lMd,10**t.lMbul if np.isfinite(t.lMbul) else 0.0,10**t.lMat,10**t.lMmol if np.isfinite(t.lMmol) else 0.0
    Mb=Md+Mbul+Mat+Mmol; rl=max(8*t.hR,30.0); r=np.linspace(0.05,rl,400)
    vb2=vexp2(r,Md,t.hR)+G*Mbul/r+vexp2(r,Mat,2*t.hR)+vexp2(r,Mmol,t.hR)  # DMS M_atom, M_mol already x1.4 for He+metals (DMS VII l.499, 523)
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2]); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vb2[-1]*rl/x
    xs=np.linspace(1e-4,rl,20001); phib0=np.trapezoid(np.interp(xs,rr,vv)/xs,xs)+vb2[-1]
    s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rc=10**LRC
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    s=solve_ivp(f,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/400)
    vt2=vb2+G*s.y[2]/r; i22=np.argmin(abs(r-2.2*t.hR))
    out.append(dict(gal=g,lMb=np.log10(Mb),hR=t.hR,Fb_dms=t.Fb,eFb=t.eFb,Fb_pred=vb2[i22]/vt2[i22],Vpred=(G*a0k*Mb)**0.25,Voiii=t.Voiii,Vstar=t.Vstar,Vmod_last=np.sqrt(vt2[-1]),V22=np.sqrt(vt2[i22]),UpsK=t.UpsK,Rc=t.Rc,lrho0=t.lrho0))
R=pd.DataFrame(out).set_index('gal'); R['Fb_pred_v']=np.sqrt(R.Fb_pred); R['dF']=R.Fb_pred_v-R.Fb_dms  # DMS F_b is V_b/V_c (velocity ratio); Fb_pred column is V_b^2/V_c^2
R=R; R['dlogV']=np.log10(R.Voiii/R.Vpred)
R['v_inf2_piso']=4*np.pi*G*10**R.lrho0*R.Rc**2; R['two_sigma2']=2*np.sqrt(a0k*G*10**R.lMb)/2
pd.set_option('display.width',220); print(R.round(3).to_string())
print('\nT1 F_b (V_b/V_c): N=%d median dF %+.3f, rms %.3f, 16-84 [%+.3f,%+.3f]'%(R.dF.notna().sum(),R.dF.median(),np.sqrt((R.dF**2).mean()),*R.dF.quantile([.16,.84])))
m=R.dlogV.notna(); print('T2 BTFR: N=%d median dlogV %+.3f, rms %.3f'%(m.sum(),R.dlogV[m].median(),np.sqrt((R.dlogV[m]**2).mean())))
m3=R.v_inf2_piso.notna(); print('T3: median log(v_inf^2 pISO / 2 sigma^2) = %+.3f (N=%d)'%(np.log10(R.v_inf2_piso[m3]/R.two_sigma2[m3]).median(),m3.sum()))
R.to_csv('dms_fb_20261105.csv')
