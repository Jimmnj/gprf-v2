# pb_solve_canon_20260925.py  [UNRATIFIED]
# Canon-form spherical Poisson-Boltzmann solve (SNAP_20260922_001 sec 0c, 6b):
#   rho_nb = rho_amb [exp(phi_loc/sigma^2) - 1],  phi_loc = phi_loc(0) - psi(r)
#        = rho_c exp(-psi/sigma^2) - rho_amb,     rho_c = rho_amb exp(phi_loc(0)/sigma^2)
# Differs from pb_solve_20260925.py ONLY by the uniform -rho_amb in the source.
# Free parameter per galaxy: rho_c (log). rho_amb scanned: 0 (=9/25 solve), 1e3 (Mistele bound), 3e4 (0c demo).
# Everything else identical: sigma modes, 10% weights, cuts, bootstrap (nb=200), seed.
# Usage (from notebooks\): python pb_solve_canon_20260925.py [ideal|btfr_meas|vflat] [rho_amb,rho_amb,...]
# See RESULT_PB_SOLVE_CANON_20260925.md
import sys, numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import gate2c_robustness_20260925 as g2
from gate2c_robustness_20260925 import load, boot, G, a0k, kpc, MRT
def table_col(i):
    out={}
    for line in open(MRT):
        p=line.split()
        if len(p)>=18 and p[1].isdigit():
            try: out[p[0]]=float(p[i])
            except ValueError: pass
    return out
def solve_gal(d, s2, rho_amb, wmode='frac'):
    r=d.r.values; vb2=np.clip(d.vbar2.values,0,None)
    rr=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2])
    vbar=lambda x: np.interp(x,rr,vv) if x<=r[-1] else vb2[-1]*r[-1]/x
    def model(l):
        rc=10**l
        f=lambda x,y:[(vbar(max(x,1e-4))+G*y[1]/max(x,1e-4))/max(x,1e-4),
                      4*np.pi*x*x*(rc*np.exp(-y[0]/s2)-rho_amb)]
        s=solve_ivp(f,[1e-4,r[-1]],[0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=r[-1]/200)
        return vb2+G*s.y[1]/r
    w=1/(0.1*d.vobs.values)**2 if wmode=='frac' else 1/np.maximum(d.ev.values,2)**2
    chi=lambda l: np.sum(w*(np.sqrt(np.clip(model(l),1,None))-d.vobs.values)**2)
    l=minimize_scalar(chi,bounds=(4,11),method='bounded',options={'xatol':0.02}).x
    r0=np.sqrt(9*s2/(4*np.pi*G*10**l))
    return model(l), l, r0
def run(P, mode, rho_amb, smult=1.0, wmode='frac', btfr=None):
    pts=[]; cores=[]
    for g,d in P.groupby('gal'):
        Mb=float(0.5*d.L.iloc[0]*1e9+1.33*d.MHI.iloc[0]*1e9)
        if mode=='ideal': s2=smult*np.sqrt(a0k*G*Mb)/2
        elif mode=='btfr_meas': s2=(10**btfr[1]*Mb**btfr[0])**2/2
        else: s2=d.Vf.iloc[0]**2/2
        v2,l,r0=solve_gal(d,s2,rho_amb,wmode)
        dd=d.copy(); dd['lgPB']=np.log10(np.clip(v2,1e-3,None)/d.r.values*1e6/kpc); pts.append(dd)
        cores.append(dict(gal=g,lrhoc=l,r0=r0,Rd=d.Rd.iloc[0],s2=s2))
    Pp=pd.concat(pts).reset_index(drop=True); Pp['resPB']=Pp.lgo-Pp.lgPB
    return Pp, pd.DataFrame(cores)
if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'ideal'
    ambs=[float(x) for x in sys.argv[2].split(',')] if len(sys.argv)>2 else [0,1e3,3e4]
    P=load(); P['Vf']=P.gal.map(table_col(15))
    if mode!='ideal': P=P[P.Vf>0].reset_index(drop=True)
    g1=P.groupby('gal').first(); Mbg=0.5*g1.L*1e9+1.33*g1.MHI*1e9
    btfr=np.polyfit(np.log10(Mbg),np.log10(g1.Vf),1) if (g1.Vf>0).all() else None
    print('mode',mode,'N',len(P),'gal',P.gal.nunique())
    prev=None
    for amb in ambs:
        g2.rng=np.random.default_rng(1)          # same bootstrap draws for every rho_amb
        Pp,C=run(P,mode,amb,btfr=btfr)
        clean=(Pp.Q==1)&(Pp.inc<=80)&(Pp.rdr>=4)
        print(f'rho_amb={amb:.0e}: scatter %.3f mean %+.3f'%(Pp.resPB.std(),Pp.resPB.mean()))
        for k,m in {'ALL':Pp.Q>0,'Q1 clean':clean}.items():
            dd=Pp[m].reset_index(drop=True)
            e,lo,hi,_,_=boot(dd,['lrRd'],y='resPB',nb=200); print(f'  {k}: r/Rd slope {e[0]:+.3f} [{lo[0]:+.3f},{hi[0]:+.3f}]')
        print('  median r0/Rd %.2f ; galaxies at lower bound %d ; median log rho_c %.2f'%(np.median(C.r0/C.Rd),(C.lrhoc<4.05).sum(),C.lrhoc.median()))
        if prev is not None:
            dl=np.abs(Pp.lgPB.values-prev.lgPB.values)
            print('  max |dlog g_PB| vs previous rho_amb: %.4f dex ; median %.5f'%(dl.max(),np.median(dl)))
        prev=Pp
