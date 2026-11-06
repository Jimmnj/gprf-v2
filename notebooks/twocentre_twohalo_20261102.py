# Comparison column (ii): two-halo reading (pre-ONE-FIELD). LMC own-sigma isothermal solve, rho_amb = MW nb at 50 kpc + 35.5,
# rho_c 10^7.45, Plummer baryons; tail cut where its density meets the MW's (density match); plus nesting factor 1.36 on the MW nb inside r_t.
import numpy as np, json
from scipy.integrate import solve_ivp
G=4.30091e-6; a0k=1.2e-10*3.0857e19/1e6; RHO_BG=35.5; S=lambda x:2*np.tanh(x/2)
j=json.load(open('twocentre_primary_20261102.json')); rho_mw=j['rho_mw_at_D']; amb=rho_mw+RHO_BG
for (Mb,a,lab) in ((3.0e9,1.5,'Plummer 3e9'),):
  s2=0.5*np.sqrt(G*a0k*Mb); rc=10**7.45; xb0=G*Mb/a/s2
  phib=lambda r:-G*Mb/np.sqrt(r*r+a*a); Mbar=lambda r: Mb*r**3/(r*r+a*a)**1.5
  def f(r,y):
    xb=-phib(r)/s2; rho=max(rc*np.exp(np.clip(S(xb)-S(xb0)-y[0]/s2,-700,50))-amb,-amb); return [G*y[1]/r**2,4*np.pi*r*r*rho]
  sol=solve_ivp(f,[1e-4,200.],[0.,0.],rtol=1e-9,atol=1e-3,dense_output=True)
  rho=lambda r: rc*np.exp(S(-phib(r)/s2)-S(xb0)-sol.sol(r)[0]/s2)-amb
  rr=np.logspace(-2,2.3,4000); rt=rr[np.argmax(rho(rr)<=0)]   # density match: own excess reaches zero against the MW density
  print(f'{lab}: sigma_own {np.sqrt(s2):.1f}  rho_amb {amb:.3e}  r_t(density match) {rt:.1f} kpc  M_own,nb(<r_t) {sol.sol(rt)[1]:.3e}')
  out={}
  for R in (10.,30.,50.):
    Rc=min(R,rt); Mnb=sol.sol(Rc)[1]; Mmw=j['M_mw'][str(R)]*(min(R,rt)/R)**3*0+0  # placeholder
    out[str(R)]=dict(M_own_nb=Mnb,M_b=Mbar(R))
  json.dump(dict(rt=rt,sigma=np.sqrt(s2),rho_amb=amb,cols=out),open('twocentre_twohalo_20261102.json','w'),indent=1)
  for R in (10.,30.,50.): print(f'  R {R:3.0f}: M_b {Mbar(R):.3e}  M_own,nb(<min(R,r_t)) {out[str(R)]["M_own_nb"]:.3e}')
