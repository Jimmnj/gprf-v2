# df1_trail_20261114.py [UNRATIFIED] DF-1 (queue item 3; JM 2026-11-14 14:43 "yes"). NGC1052-DF2/DF4/DF9 at the June 2026 distances under ONE-FIELD SCOPED (A').
# RULE 3a: the dwarf protocol's machinery is used VERBATIM: n9_dwarfs_20261003.py up to 'def run(' (G, a0k, rho_bg, host_rho, solve, sigma_los), exec'd exactly as in
#   dwarf_prereg_20261020.py / onefield_dwarfs_20261102.py / rhoc_sigma_b_oos_20261019.py; rule (3) regression verbatim from dwarf_prereg; trend-only (M_b, Sigma_h)
#   regression verbatim from rhoc_sigma_b_oos. Stage 'rep' reproduces the record (rhoc_sigma_b_oos_df_rerun_20261020.txt) before any new number is used.
# A': the body is solved in its own phi_b with its own sigma (eq. 2, M_b = M*, gas-free); the host enters only as ambient density.
#   FIELD = rho_bg 35.5 (DF2 is 1.8-3.8 Mpc in front of NGC 1052 at 17.6 Mpc); HOST bracket = NGC 1052 250 km/s at 80 kpc + rho_bg (the old protocol).
# RADIUS CONVENTION (flag): n9 sets the Plummer scale a = rh/1.305, i.e. reads rh as the 3-D half-light radius, and averages sigma_los over projected R <= rh.
#   C = consistent with the note's stated protocol ('projected dispersion within the half-light radius'): a = R_e (projected), aperture R <= R_e.
#   P = as coded for sec 4.3: rh fed as before (DF2 r_1/2 = 2.7 kpc 3-D [Danieli+2019 = 4/3 x R_e,c 2.0]; DF4 R_e 1.6 [Shen+2023]); DF9 R_e 1.1 as Appendix B rows.
# INPUTS (LIT_TRITONSTATION sec 8t; Danieli+2019 checked 2026-11-14 15:0x): M* and R_e scale as D^2 and D.
#   DF2: M* 2.0e8, R_e,c 2.0 kpc at 20 Mpc (vD+2018; Cohen+2018 via Danieli+2019); sigma 8.5 +2.3/-3.1 (Danieli+2019 ApJL 874 L12; arXiv v1 8.4 +- 2.1); D 17.6 (JWST TRGB) | alt 21.7 (HST TRGB)
#   DF4: M* 1.5e8 (2 x 7.5e7 inside r_1/2, Shen+2023), R_e 1.6 at 20 Mpc; sigma 8.0 +2.3/-1.9 (Shen+2023); D 18.8 (SBF, Tang+2026) | alt 20.0 (TRGB, Danieli+2020)
#   DF9: M* 1.42e8 (2 x 0.71e8 inside r_e), R_e 1.1 at 20.6 Mpc; sigma 6.4 +4.0/-4.3 (Keim+2026); D 18.7 (SBF, Tang+2026) | alt 20.6 (Keim's assumed)
import sys, time, warnings, numpy as np, pandas as pd; warnings.filterwarnings('ignore')
from multiprocessing import Pool
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
STAGE=sys.argv[1]; OUT='df1_trail_20261114_report.txt'
def say(s):
    print(s); open(OUT,'a').write(s+'\n')
LMED=7.42                      # V3 first-minimum, fixed arm, 141 galaxies (RESULT_EP3_20261114; EP3-4 cut)
FIELD=rho_bg; HOST=host_rho(250.,80.)+rho_bg
def sig_one(args):
    Mb,a,ap,amb,lrc=args
    s2=0.5*np.sqrt(G*a0k*Mb); rmax=max(8*1.305*a,3*np.sqrt(s2/(2*np.pi*G*amb)))   # = protocol rmax max(8 rh, 3 r_t) with rh = 1.305 a
    s=solve(Mb,a,s2,amb,lrc,rmax); return sigma_los(s,a,ap,rmax)
def stars_only(Mb,a,ap):
    s2=0.5*np.sqrt(G*a0k*Mb); s=solve(Mb,a,s2,rho_bg,-99,60*a); return sigma_los(s,a,ap,60*a)
def rules():
    from gate2c_robustness_20260925 import load_table
    T=load_table(); F=pd.read_csv('test1_capped_baryon_20261018.csv').set_index('gal').join(T[['Rd','L','MHI']]); F=F[F.Rd>0]
    lM=np.log10(F.Mb.values); lS=np.log10(F.Mb.values/(2*np.pi*(1.678*F.Rd.values*1e3)**2)); fg=1.33*F.MHI.values*1e9/F.Mb.values; y=F.l_tanh2.values
    X=np.column_stack([np.ones(len(y)),lM,lS,fg]); beta=np.linalg.lstsq(X,y,rcond=None)[0]
    X2=np.column_stack([np.ones(len(y)),lM,lS]); beta2=np.linalg.lstsq(X2,y,rcond=None)[0]
    return beta,beta2,len(y)
LR=np.round(np.arange(4.5,8.501,0.1),2)
OBJ=[('DF2',2.0e8,2.0,2.7,20.0,[(8.5,2.3,3.1,'ApJL'),(8.4,2.1,2.1,'v1')],17.6,21.7),
     ('DF4',1.5e8,1.6,1.6,20.0,[(8.0,2.3,1.9,'Shen23')],18.8,20.0),
     ('DF9',1.42e8,1.1,1.1,20.6,[(6.4,4.0,4.3,'Keim26')],18.7,20.6)]
def req(sig,v):
    if v<=sig[0]: return -np.inf
    if v>=sig[-1]: return np.inf
    i=np.argmax(sig>=v); return LR[i-1]+(v-sig[i-1])/(sig[i]-sig[i-1])*(LR[i]-LR[i-1])
fmt=lambda x: ('<%.1f'%LR[0]) if x==-np.inf else (('>%.1f'%LR[-1]) if x==np.inf else '%.2f'%x)
if STAGE=='rep':
    t=time.time(); say('=== df1_trail_20261114 [UNRATIFIED] stage rep %s'%time.strftime('%Y-%m-%d %H:%M'))
    beta,beta2,n=rules(); say('rule (3) refit (dwarf_prereg verbatim, N=%d): log rho_c = %.3f %+.3f log Mb %+.3f log Sigma_h %+.3f f_gas  [note eq: 11.06 -0.365 +0.268 -0.93]'%(n,*beta))
    say('trend-only (rhoc_sigma_b_oos verbatim): %.3f %+.3f log Mb %+.3f log Sigma_h'%tuple(beta2))
    LR25=np.arange(5.0,8.01,0.25)
    for nm,Mb,rh,so,ep,em,rec in (('DF2',2.0e8,2.7,8.4,2.1,2.1,'22.4 | req 5.82 [nan,6.22] | stars 6.6'),('DF4',1.5e8,1.6,8.0,2.3,1.9,'17.3 | req 5.44 [nan,6.52] | stars 7.8')):
        a=rh/1.305; sc=np.array([sig_one((Mb,a,rh,HOST,l)) for l in LR25])
        f=lambda v: np.interp(v,sc,LR25) if sc.min()<v<sc.max() else np.nan
        s2=0.5*np.sqrt(G*a0k*Mb); sN_host=sigma_los(solve(Mb,a,s2,HOST,-99,8*rh),a,rh,8*rh)
        say('REP %s: at 7.45 %.2f | req %.2f [%.2f,%.2f] | stars-only (host amb, onefield form) %.2f, (field, pure) %.2f || RECORD %s'%(nm,sig_one((Mb,a,rh,HOST,7.45)),f(so),f(so-em),f(so+ep),sN_host,stars_only(Mb,a,rh),rec))
    say('rep time %.0fs'%(time.time()-t))
if STAGE in ('main','alt'):
    t=time.time(); beta,beta2,n=rules(); jobs=[]; keys=[]
    for nm,M0,Re0,rhP0,D0,obs,Dp,Da in OBJ:
        runs=[(Dp,'C','FIELD'),(Dp,'C','HOST')] if STAGE=='main' else [(Da,'C','FIELD'),(Dp,'P','FIELD')]
        for D,cv,am in runs:
            k=D/D0; Mb=M0*k*k; Re=Re0*k; rhP=rhP0*k
            a,ap=(Re,Re) if cv=='C' else (rhP/1.305,rhP)
            amb=FIELD if am=='FIELD' else HOST
            keys.append((nm,D,cv,am,Mb,Re,a,ap,amb,obs)); jobs+= [(Mb,a,ap,amb,l) for l in LR]
    with Pool(2) as p: S=np.array(p.map(sig_one,jobs)).reshape(len(keys),len(LR))
    say('=== stage %s %s (scan log rho_c %.1f..%.1f step 0.1; %d solves; %.0fs)'%(STAGE,time.strftime('%H:%M'),LR[0],LR[-1],len(jobs),time.time()-t))
    rows=[]
    for (nm,D,cv,am,Mb,Re,a,ap,amb,obs),sig in zip(keys,S):
        mono=bool(np.all(np.diff(sig)>0)); s2=0.5*np.sqrt(G*a0k*Mb)
        lSh=np.log10(Mb/(2*np.pi*(Re*1e3)**2)); l3=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh; lt=beta2[0]+beta2[1]*np.log10(Mb)+beta2[2]*lSh
        sAt=lambda l: np.interp(l,LR,sig) if LR[0]<=l<=LR[-1] else np.nan
        sN=stars_only(Mb,a,ap); sMi=(4*G*a0k*Mb/81)**0.25; sW=np.sqrt(G*(Mb/2)/(3*4/3*Re))
        say('%s D %.1f %s %-5s | Mb %.2e R_e %.2f a %.2f ap %.2f | eq2 sigma %.1f | stars-only %.1f (Wolf %.1f) | MOND iso %.1f | at median 7.42: %.1f | trend %.2f -> %.1f | rule(3) f_gas=0 %.2f -> %s | log Sigma_h %.2f | mono %s'%(
            nm,D,cv,am,Mb,Re,a,ap,np.sqrt(s2),sN,sW,sMi,sAt(LMED),lt,sAt(lt),l3,('%.1f'%sAt(l3)) if np.isfinite(sAt(l3)) else '>scan',lSh,mono))
        for so,ep,em,lab in obs:
            r0,rl,rhh=req(sig,so),req(sig,so-em),req(sig,so+ep)
            say('      obs %s %.1f +%.1f/-%.1f -> REQUIRED log rho_c %s [%s, %s] | median - required = %s dex'%(lab,so,ep,em,fmt(r0),fmt(rl),fmt(rhh),('%.2f'%(LMED-r0)) if np.isfinite(r0) else 'n/a'))
            rows.append(dict(stage=STAGE,obj=nm,D=D,conv=cv,amb=am,Mb=Mb,Re=Re,a=a,ap=ap,sig_eq2=np.sqrt(s2),s_stars=sN,s_wolf=sW,s_mond=sMi,s_med=sAt(LMED),l_trend=lt,s_trend=sAt(lt),l_rule3=l3,
                             obs=lab,s_obs=so,ep=ep,em=em,req=r0,req_lo=rl,req_hi=rhh))
        say('      scan: '+' '.join('%.1f:%.1f'%(l,s) for l,s in zip(LR[::5],sig[::5])))
    pd.DataFrame(rows).to_csv('df1_trail_20261114_%s.csv'%STAGE,index=False)
if STAGE=='appB':
    t=time.time(); beta,beta2,n=rules()
    ns={}; exec('dw=['+open('dwarf_prereg_20261020.py').read().split('\ndw=[')[1].split('\nrows=[]')[0],globals(),ns); dw=ns['dw']
    rec=pd.read_csv('dwarf_prereg_20261020.csv').set_index('name')
    jobs=[]
    for name,LV,rh,so,ep,em,vc,Dh,cls in dw:
        Mb=1.6*LV; lSh=np.log10(Mb/(2*np.pi*(rh*1e3)**2)); lrc=beta[0]+beta[1]*np.log10(Mb)+beta[2]*lSh+beta[3]*0.0; amb=host_rho(vc,Dh)+rho_bg
        jobs+=[(Mb,rh/1.305,rh,amb,lrc),(Mb,rh,rh,amb,lrc)]
    with Pool(2) as p: S=np.array(p.map(sig_one,jobs)).reshape(len(dw),2)
    say('=== stage appB %s: Appendix B rows (dwarf_prereg protocol, rule (3) f_gas=0, host ambient) under P (as coded, a = rh/1.305) and C (a = rh = projected R_e) (%.0fs)'%(time.strftime('%H:%M'),time.time()-t))
    out=[]
    for (name,LV,rh,so,ep,em,vc,Dh,cls),(sP,sC) in zip(dw,S):
        bar=np.sqrt(0.10**2+(0.5*(ep+em)/so/np.log(10))**2); dP=np.log10(sP/so); dC=np.log10(sC/so)
        out.append(dict(name=name,cls=cls,s_obs=so,s_P=sP,s_rec=rec.loc[name,'s_GPRF'],s_C=sC,dP=dP,dC=dC,bar=bar,pP=abs(dP)<=bar,pC=abs(dC)<=bar))
        say('  %-10s %-3s obs %5.1f | P %5.2f (record %5.2f) %s | C %5.2f %s | C/P %.3f'%(name,cls,so,sP,rec.loc[name,'s_GPRF'],'PASS' if abs(dP)<=bar else 'fail',sC,'PASS' if abs(dC)<=bar else 'fail',sC/sP))
    O=pd.DataFrame(out); O.to_csv('df1_trail_20261114_appB.csv',index=False)
    for lab,d,pp in (('P',O.dP,O.pP),('C',O.dC,O.pC)):
        say('  %s: pass %d/%d | rms %.3f dex | mean %+.3f | embedded %+.3f (N=%d) isolated %+.3f (N=%d)'%(lab,pp.sum(),len(O),np.sqrt((d**2).mean()),d.mean(),d[O.cls!='ISO'].mean(),(O.cls!='ISO').sum(),d[O.cls=='ISO'].mean(),(O.cls=='ISO').sum()))
    say('  max |P - record| = %.3f km/s'%np.abs(O.s_P-O.s_rec).max())
