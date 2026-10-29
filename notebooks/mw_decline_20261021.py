# mw_decline_20261021.py [UNRATIFIED] PRE-REGISTERED 2026-10-21 16:20 (JM 'run', after Triton Station comment: 'The RC of the MW is declining and it
# should be declining in MOND because that is what the mass distribution says it should do').
# TEST: frozen law (tanh2 cap, BTFR sigma, rho_amb 35.5) on the Ou et al. 2024 (arXiv:2303.12838) baryon model (their Table 2, from de Salas+2019 /
# Misiriotis+2006): Hernquist bulge 1.55e10 a=0.70; exp disc 3.65e10 Rd 2.35; HI 8.2e9 Rd 18.24; H2 1.3e9 Rd 2.57. Disc terms by the Freeman thin-disc
# formula; nb response spherical (as in the SPARC pipeline). Compare v(R) 6-27 kpc with their Table 1 (Gaia DR3 Jeans curve) at UNIVERSAL rho_c=10^7.45
# and at FITTED rho_c; against Newton (baryons only) and spherical MOND (RAR nu, same baryons). DECLARED outputs: rms residual (km/s) 8-27 kpc; slope
# dv/dR 10-25 kpc (km/s/kpc) for data/law/MOND; asymptotic v_flat. Error model for fits: quadrature of the table errors and a 5 km/s floor (inner bins
# have 0.5 km/s errors that no spherical model should be judged by). PASS criterion (declared): law at universal rho_c within 2x the MOND rms, and both
# slopes negative within 50% of the data slope 10-25 kpc.
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.special import i0,i1,k0,k1
from scipy.optimize import minimize_scalar
G=4.30091e-6; a0=1.2e-10; kpc=3.0857e19; a0k=a0*kpc/1e6; rho_amb=35.5
S=lambda x:2*np.tanh(x/2)
rows="""6.27 231.07 1.28 1.00
6.78 230.93 0.97 0.97
7.28 232.87 1.13 0.79
7.86 234.14 0.63 0.62
8.19 232.83 0.56 0.55
8.71 231.27 0.54 0.66
9.23 230.47 0.44 0.46
9.72 229.16 0.54 0.43
10.24 229.37 0.25 0.50
10.74 227.95 0.53 0.34
11.23 227.09 0.44 0.46
11.73 227.15 0.35 0.45
12.23 226.90 0.39 0.41
12.73 225.61 0.64 0.53
13.22 224.95 0.66 0.69
13.72 222.79 0.68 0.53
14.22 222.13 1.04 0.64
14.73 220.08 0.65 0.89
15.24 218.25 1.14 0.79
15.72 221.16 1.06 1.24
16.24 218.30 2.17 1.88
16.77 217.07 1.39 1.35
17.21 219.56 1.80 1.69
17.77 215.49 1.98 2.08
18.23 214.62 2.24 1.59
18.73 210.89 1.42 1.32
19.22 208.48 2.20 1.65
19.71 205.97 1.20 1.55
20.22 202.97 1.55 2.25
20.72 195.16 2.60 1.97
21.22 200.20 2.84 1.45
21.72 201.11 2.72 3.81
22.27 196.79 4.82 6.25
22.71 218.65 14.93 17.54
23.40 192.49 4.25 4.8
25.02 191.48 6.41 9.61
27.31 172.98 15.82 17.07"""
D=pd.DataFrame([list(map(float,l.split())) for l in rows.splitlines()],columns=['R','v','ep','em']); D['e']=np.sqrt(((D.ep+D.em)/2)**2+5.0**2)
def vdisc2(R,M,Rd):
    y=R/(2*Rd); S0=M/(2*np.pi*Rd**2); return 4*np.pi*G*S0*Rd*y*y*(i0(y)*k0(y)-i1(y)*k1(y))
def vhern2(R,M,a): return G*M*R/(R+a)**2
comps=dict(bulge=(1.55e10,0.70,'h'),disc=(3.65e10,2.35,'d'),HI=(8.2e9,18.24,'d'),H2=(1.3e9,2.57,'d'))
Mb=sum(c[0] for c in comps.values())
def vb2(R): return sum(vhern2(R,M,s) if k=='h' else vdisc2(R,M,s) for M,s,k in comps.values())
rl=150.0; rr=np.linspace(1e-3,rl,30001); vv=vb2(rr)
vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vv[-1]*rl/x
phib0=np.trapezoid(vv/rr,rr)+vv[-1]; s2b=np.sqrt(a0k*G*Mb)/2
print('M_b=%.3e Msun  sigma=%.1f km/s  v_flat(BTFR)=%.1f km/s  phi_b0/sigma^2=%.2f  r_t=%.2f Mpc'%(Mb,np.sqrt(s2b),np.sqrt(2*s2b),phib0/s2b,np.sqrt(s2b)/np.sqrt(2*np.pi*G*rho_amb)/1e3))
Rg=np.concatenate([D.R.values,[30,35,40,50,60,80,100]])
def law(l,R=Rg):
    rc=10**l; xb0=phib0/s2b; Sb0=S(xb0)
    f=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*(rc*np.exp(S(xb0-y[0]/s2b)-Sb0-y[1]/s2b)-rho_amb)]
    s=solve_ivp(f,[1e-4,R[-1]],[0,0,0],t_eval=R,rtol=1e-7,atol=1e-4,max_step=0.2)
    return np.sqrt(vb2(R)+G*s.y[2]/R)
def mond(R=Rg):
    gb=vb2(R)/R; nu=1/(1-np.exp(-np.sqrt(gb/a0k))); return np.sqrt(gb*nu*R)
def chi(v): return float(np.sum(((v[:len(D)]-D.v.values)/D.e.values)**2))
vN=np.sqrt(vb2(Rg)); vM=mond(); vU=law(7.45)
fit=minimize_scalar(lambda l:chi(law(l)),bounds=(5.5,9.5),method='bounded'); vF=law(fit.x)
n=len(D); m=(D.R>=8)
def rms(v): return float(np.sqrt(np.mean((v[:n][m]-D.v.values[m])**2)))
def slope(R,v):
    k=(R>=10)&(R<=25); return float(np.polyfit(R[k],v[k],1)[0])
print('fitted log rho_c = %.2f  (chi2 %.1f, N=%d)'%(fit.x,fit.fun,n))
print('%-28s %8s %8s %8s %8s'%('model','rms8-27','slope10-25','v(8.2)','v(27.3)'))
for name,v,c in (('Newton (baryons only)',vN,chi(vN)),('MOND spherical (RAR nu)',vM,chi(vM)),('law, rho_c=10^7.45',vU,chi(vU)),('law, rho_c fitted %.2f'%fit.x,vF,chi(vF))):
    print('%-28s %8.1f %+8.2f %8.1f %8.1f   chi2=%.1f'%(name,rms(v),slope(Rg[:n],v[:n]),v[4],v[n-1],c))
print('%-28s %8s %+8.2f %8.1f %8.1f'%('DATA (Ou+2024)','-',slope(D.R.values,D.v.values),D.v[4],D.v[n-1]))
print('\nR(kpc)  data   Newton  MOND   law7.45  lawfit')
for i,R in enumerate(Rg):
    d='%6.1f'%D.v[i] if i<n else '   -  '
    if i<n and i%3!=0 and R<22: continue
    print('%5.1f  %s %6.1f %6.1f %6.1f %6.1f'%(R,d,vN[i],vM[i],vU[i],vF[i]))
pd.DataFrame(dict(R=Rg,newton=vN,mond=vM,law745=vU,lawfit=vF)).to_csv('mw_decline_20261021.csv',index=False)
