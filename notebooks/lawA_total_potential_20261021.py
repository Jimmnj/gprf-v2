# lawA_total_potential_20261021.py [UNRATIFIED] 2026-10-21 20:50 (JM 'ok' to the action attempt).
# DERIVED (DERIVATION_ACTION_20261021.md): with canonical auxiliary fields and local couplings (QUMOND class), an action exists only if
# rho_nb = F(psi + kappa chi); the BTFR isothermal tail needs kappa = 1, so rho_nb = F(total potential): the single-field perfect-fluid action
# L = -|grad Phi|^2/8piG - rho_b Phi + P(Phi), dP/dPhi = -rho_nb (matter couples to Phi only; momentum conserved).
# LAW-A: rho_nb = rho_c exp[S(x) - S(x0)] - rho_amb, x = (total bump potential)/sigma^2 (baryons + nb), x0 central. Identical to the frozen law
# where S is linear; differs where saturated (the nb self-potential is saturated too). Absolute x needs the nb central potential: iterate.
# TEST (declared): Milky Way (Ou+2024) rms 8-27 kpc at universal and fitted rho_c vs the frozen law (19.7 / 8.6). Pass = within 1.5x.
import numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.special import i0,i1,k0,k1
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6; ramb=35.5; S=lambda x:2*np.tanh(x/2)
rows=open('mw_decline_20261021.py').read().split('rows="""')[1].split('"""')[0]
D=pd.DataFrame([list(map(float,l.split())) for l in rows.splitlines()],columns=['R','v','ep','em']); D['e']=np.sqrt(((D.ep+D.em)/2)**2+25)
def vdisc2(R,M,Rd):
    y=R/(2*Rd); S0=M/(2*np.pi*Rd**2); return 4*np.pi*G*S0*Rd*y*y*(i0(y)*k0(y)-i1(y)*k1(y))
comps=[(1.55e10,0.70,'h'),(3.65e10,2.35,'d'),(8.2e9,18.24,'d'),(1.3e9,2.57,'d')]; Mb=sum(c[0] for c in comps)
vb2=lambda R: sum(G*M*R/(R+s)**2 if k=='h' else vdisc2(R,M,s) for M,s,k in comps)
rl=150.; rr=np.linspace(1e-3,rl,30001); vv=vb2(rr); vbar=lambda x: np.interp(x,rr,vv) if x<=rl else vv[-1]*rl/x
phib0=np.trapezoid(vv/rr,rr)+vv[-1]; s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2
def solve(l,law,xn0=0.0,rmax=6000.):
    rc=10**l; x0=xb0+xn0
    def f(x,y):
        x_=max(x,1e-4); xb=xb0-y[0]/s2
        if law=='frozen': e=S(xb)-S(xb0)-y[1]/s2
        else: e=S(xb+xn0-y[1]/s2)-S(x0)
        return [vbar(x_)/x_, G*y[2]/x_**2, 4*np.pi*x*x*max(rc*np.exp(e)-ramb,-ramb)]
    R=np.concatenate([D.R.values,np.geomspace(30,rmax,200)])
    s=solve_ivp(f,[1e-4,rmax],[0,0,0],t_eval=R,rtol=1e-7,atol=1e-4,max_step=1.0)
    # nb central bump potential = int_0^inf G M_nb/r^2 dr (Mnb const beyond rmax)
    Rr=s.t; M=np.maximum(s.y[2],0); xn=(np.trapezoid(G*M/Rr**2,Rr)+G*M[-1]/Rr[-1])/s2
    return s, xn
def curve(l,law):
    xn=0.0
    if law=='A':
        for it in range(25):
            s,xn_new=solve(l,'A',xn)
            if abs(xn_new-xn)<1e-3: break
            xn=0.5*xn+0.5*xn_new
    s,_=solve(l,law,xn); n=len(D); R=D.R.values
    return np.sqrt(vb2(R)+G*s.y[2][:n]/R), xn
m=(D.R>=8).values
rms=lambda v: float(np.sqrt(np.mean((v[m]-D.v.values[m])**2)))
chi=lambda v: float(np.sum(((v-D.v.values)/D.e.values)**2))
print('Milky Way (Ou+2024). x_b0 = %.2f'%xb0)
for law in ('frozen','A'):
    vu,xn=curve(7.45,law)
    fit=minimize_scalar(lambda l: chi(curve(l,law)[0]),bounds=(6.0,9.5),method='bounded',options={'xatol':0.02})
    vf,xnf=curve(fit.x,law)
    print('  %-7s universal 7.45: rms %5.1f (x_nb0 %.2f) | fitted %.2f: rms %5.1f (x_nb0 %.2f) | v(8.2) %5.1f/%5.1f  v(27.3) %5.1f/%5.1f'%(law,rms(vu),xn,fit.x,rms(vf),xnf,vu[4],vf[4],vu[-1],vf[-1]))
