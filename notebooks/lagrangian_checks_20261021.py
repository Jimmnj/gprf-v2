# lagrangian_checks_20261021.py [UNRATIFIED] 2026-10-21 19:40 (JM 'ok' to LAGRANGIAN? recs (a)+(b)).
# (A) MOMENTUM: the law places each body's nb mass on that body; the force of B on A's nb mass has no recipient. For two isolated bodies (point baryons,
#     spherical truncated isothermal tails, shell theorem, d < r_t) net force on the pair F_net = G/d^2 [M_A M_hB(<d) - M_B M_hA(<d)] (toward B).
#     Report F_net / (mutual force on the lighter body's baryons) and the implied centre-of-mass acceleration (inertia = baryons + nb inside d of each)
#     and Delta v over 10 Gyr. Pairs: MW-M31 (6.15e10, 1.0e11, 780 kpc), MW-LMC (6.15e10, 3e9, 50 kpc), equal pair (control: must be 0).
# (B) REACTION in an action: if L = -|grad Phi|^2/8piG - rho_b Phi + P(Phi, W[phi_b]) with P = sigma^2 rho_c exp[-(Phi-W)/sigma^2],
#     W = phi_b + sigma^2 S(-phi_b/sigma^2), baryons feel rho_b + rho_nb S'(x_b) (S' = sech^2(x_b/2)), not rho_b + rho_nb. Compute the MW curve the
#     baryons would trace under that action (Ou baryons, universal and fitted rho_c) against Ou's data; compare rms with the constitutive law.
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.special import i0,i1,k0,k1
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6; ramb=35.5; Gyr=3.156e16; kms_per_kpc2=1/3.0857e16  # (km/s)^2/kpc -> km/s^2
def Mh(M,d):
    s2=np.sqrt(G*a0k*M)/2; rt=np.sqrt(s2)/np.sqrt(2*np.pi*G*ramb); r=min(d,rt); return 2*s2*r/G-4/3*np.pi*ramb*r**3
print('(A) PAIR MOMENTUM (point baryons, truncated isothermal tails, rho_amb 35.5)')
for name,MA,MB,d in (('MW-M31',6.15e10,1.0e11,780.),('MW-LMC',6.15e10,3e9,50.),('equal 6e10 @500',6e10,6e10,500.)):
    hA,hB=Mh(MA,d),Mh(MB,d); Fnet=G/d**2*(MA*hB-MB*hA); Fmut=G*MB*(MA+hA)/d**2  # force on B's baryons from all of A
    acm=abs(Fnet)/(MA+MB+hA+hB)*kms_per_kpc2  # km/s^2
    dv=acm*10*Gyr
    print('  %-16s M_hA(<d) %.2e M_hB(<d) %.2e | F_net/F_on_B %.3f | a_cm %.2e km/s^2 (%.1e a0) | dv(10 Gyr) %.0f km/s'%(name,hA,hB,Fnet/Fmut,acm,acm*1e3/1.2e-10,dv))
# (B)
rows=open('mw_decline_20261021.py').read().split('rows="""')[1].split('"""')[0]
D=pd.DataFrame([list(map(float,l.split())) for l in rows.splitlines()],columns=['R','v','ep','em']); D['e']=np.sqrt(((D.ep+D.em)/2)**2+25)
def vdisc2(R,M,Rd):
    y=R/(2*Rd); S0=M/(2*np.pi*Rd**2); return 4*np.pi*G*S0*Rd*y*y*(i0(y)*k0(y)-i1(y)*k1(y))
comps=[(1.55e10,0.70,'h'),(3.65e10,2.35,'d'),(8.2e9,18.24,'d'),(1.3e9,2.57,'d')]; Mb=sum(c[0] for c in comps)
vb2=lambda R: sum(G*M*R/(R+s)**2 if k=='h' else vdisc2(R,M,s) for M,s,k in comps)
rl=150.; rr=np.linspace(1e-3,rl,30001); vv=vb2(rr); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vv[-1]*rl/x
phib0=np.trapezoid(vv/rr,rr)+vv[-1]; s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; S=lambda x:2*np.tanh(x/2); Sp=lambda x:1/np.cosh(x/2)**2
def run(l):
    rc=10**l
    def f(x,y):
        x_=max(x,1e-4); xb=xb0-y[0]/s2; rho=rc*np.exp(S(xb)-S(xb0)-y[1]/s2)-ramb
        return [vbar(x_)/x_, G*y[2]/x_**2, 4*np.pi*x*x*rho, 4*np.pi*x*x*rho*Sp(xb)]
    s=solve_ivp(f,[1e-4,D.R.max()],[0,0,0,0],t_eval=D.R.values,rtol=1e-7,atol=1e-4,max_step=0.2)
    R=D.R.values; vlaw=np.sqrt(vb2(R)+G*s.y[2]/R); vact=np.sqrt(np.maximum(vb2(R)+G*s.y[3]/R,0)); xb=xb0-s.y[0]/s2
    return R,vlaw,vact,s.y[2],s.y[3],xb
print('\n(B) ACTION REACTION, Milky Way (Ou+2024 baryons/data; x_b0 = %.1f, S\'(x_b0) = %.1e)'%(xb0,Sp(xb0)))
m=D.R>=8
for l in (7.45,7.88):
    R,vl,va,Mn,Me,xb=run(l)
    rms=lambda v: np.sqrt(np.mean((v[m]-D.v[m])**2))
    print('  log rho_c %.2f: constitutive law rms %.1f km/s | action (baryons feel rho_nb S\') rms %.1f km/s'%(l,rms(vl),rms(va)))
    for i in (4,12,24,36): print('     R %5.1f  x_b %5.2f  S\' %.2f  M_nb %.2e  M_nb,felt %.2e (%.0f%%)  v_law %5.1f  v_action %5.1f  data %5.1f'%(R[i],xb[i],Sp(xb[i]),Mn[i],Me[i],100*Me[i]/Mn[i],vl[i],va[i],D.v[i]))
