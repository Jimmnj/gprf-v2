# V1 control: 1D spherical shooting of the same law on the same spherical-equivalent MW baryons (no LMC); compare with twocentre MW-alone 2D.
import numpy as np, json, sys
from scipy.integrate import solve_ivp
exec(open('twocentre_20261102.py').read().split("if __name__")[0].replace("args=dict(a.split('=') for a in sys.argv[1:])","args={}"))
rc=10**LRC
def Mb(r): return sum(m*r*r/(r+a)**2 for m,a in MW)
phib0=phi_mw(1e-6); xb0=-phib0/S2
def f(r,y):  # y = [psi_nb, M_nb]
    psi=y[0]; Mn=y[1]; xb=-phi_mw(r)/S2
    rho=max(rc*np.exp(np.clip(S(xb)-S(xb0)-psi/S2,-700,50))-RHO_BG,-RHO_BG)
    return [G*Mn/r**2, 4*np.pi*r*r*rho]
sol=solve_ivp(f,[1e-4,3000.],[0.0,0.0],rtol=1e-9,atol=1e-3,dense_output=True)
Rv=[8.,12.,16.,20.,27.]
vc=[np.sqrt(G*(Mb(R)+sol.sol(R)[1])/R) for R in Rv]
tag=sys.argv[1] if len(sys.argv)>1 else 'primary'
j=json.load(open(f'twocentre_{tag}_20261102.json'))
print('R    v_1D    v_2D(MWalone)  diff')
for R,v,v2 in zip(Rv,vc,j['vc_mw']): print(f'{R:4.0f}  {v:7.2f}  {v2:7.2f}  {v2-v:+.3f}')
# M_nb within 50 kpc of MW centre (1D) vs 2D monopole mass
d=np.load(f'twocentre_{tag}_20261102.npz'); rr=d['rr']; mu=d['mu']; w=np.polynomial.legendre.leggauss(len(mu))[1]
rho0=d['rho_mw']@w/2; M2=np.trapz(4*np.pi*rr**2*rho0*(rr<=50),rr)
print(f'M_nb(<50 kpc, MW centre): 1D {sol.sol(50.)[1]:.4e}   2D {M2:.4e}   ratio {M2/sol.sol(50.)[1]:.5f}')
rho1D=lambda r: rc*np.exp(S(-phi_mw(r)/S2)-S(xb0)-sol.sol(r)[0]/S2)-RHO_BG
print(f'rho_nb at 50 kpc: 1D {rho1D(50.):.5e}  2D {j["rho_mw_at_D"]:.5e}')
