# df1_sat_20261114.py [UNRATIFIED] DF-1 diagnostic: the dwarf solver (n9) implements eq. (1) with S(x) = x; the frozen R1 has S = 2 tanh(x/2) (pass2 solve_gal).
# Variant: same ODE layout as n9 solve (y = [psi_nb, M_tot], M_tot incl. Plummer stars) but exponent S(x_b) - S(x_b0) - psi_nb/s2, x_b = G M_b/(s2 sqrt(r^2+a^2)).
# Check 1: with S = identity it must reproduce n9 solve. Then DF2/DF4/DF9 (convention C, FIELD, primary distances): sigma at 7.42 and required rho_c.
import warnings, numpy as np; warnings.filterwarnings('ignore')
from scipy.integrate import solve_ivp
src=open('n9_dwarfs_20261003.py').read().split('def run(')[0]; exec(src)
Sid=lambda x:x; Sth=lambda x:2*np.tanh(x/2)
def solve_S(Mb,a,s2,rho_amb,lrc,rmax,S):
    rc=10**lrc; xb0=G*Mb/(a*s2); Sb0=S(xb0)
    rhob=lambda r: 3*Mb/(4*np.pi*a**3)*(1+r*r/(a*a))**-2.5
    Mbr=lambda r: Mb*r**3/(r*r+a*a)**1.5
    def f(r,y):
        pn,Mt=y; r=max(r,1e-5); xb=G*Mb/(s2*np.sqrt(r*r+a*a))
        return [G*(Mt-Mbr(r))/r**2, 4*np.pi*r*r*(rhob(r)+max(rc*np.exp(S(xb)-Sb0-pn/s2)-rho_amb,-rho_amb))]
    return solve_ivp(f,[1e-5,rmax],[0,0],dense_output=True,rtol=1e-8,atol=1e-3,max_step=rmax/2000)
LR=np.round(np.arange(4.5,8.501,0.1),2)
def scan(Mb,a,amb,S):
    s2=0.5*np.sqrt(G*a0k*Mb); rmax=max(8*1.305*a,3*np.sqrt(s2/(2*np.pi*G*amb)))
    return np.array([sigma_los(solve_S(Mb,a,s2,amb,l,rmax,S),a,a,rmax) for l in LR]),G*Mb/(a*s2)
def req(sig,v):
    if v<=sig[0]: return -np.inf
    if v>=sig[-1]: return np.inf
    i=np.argmax(sig>=v); return LR[i-1]+(v-sig[i-1])/(sig[i]-sig[i-1])*(LR[i]-LR[i-1])
out=open('df1_trail_20261114_report.txt','a')
def say(s): print(s); out.write(s+'\n')
# check 1 (identity vs n9 solve, DF2 C FIELD 17.6 at 7.42)
Mb=2.0e8*0.88**2; a=2.0*0.88; s2=0.5*np.sqrt(G*a0k*Mb); amb=rho_bg; rmax=max(8*1.305*a,3*np.sqrt(s2/(2*np.pi*G*amb)))
s_n9=sigma_los(solve(Mb,a,s2,amb,7.42,rmax),a,a,rmax); s_id=sigma_los(solve_S(Mb,a,s2,amb,7.42,rmax,Sid),a,a,rmax)
say('=== saturation check (df1_sat_20261114.py): identity-S variant vs n9 solve, DF2 at 7.42: %.3f vs %.3f km/s'%(s_id,s_n9))
for nm,Mb,a,so,ep,em in (('DF2',2.0e8*0.88**2,2.0*0.88,8.5,2.3,3.1),('DF4',1.5e8*0.94**2,1.6*0.94,8.0,2.3,1.9),('DF9',1.42e8*(18.7/20.6)**2,1.1*18.7/20.6,6.4,4.0,4.3)):
    sl,xb0=scan(Mb,a,rho_bg,Sid); st,_=scan(Mb,a,rho_bg,Sth)
    f=lambda sig: (np.interp(7.42,LR,sig),req(sig,so),req(sig,so+ep))
    L,T=f(sl),f(st)
    say('  %s x_b0 %.2f | S=x: at 7.42 %.2f, req %.2f, +1sd bound %.2f | S=2tanh(x/2): at 7.42 %.2f, req %.2f, +1sd bound %.2f | d(req) %+.3f dex'%(nm,xb0,L[0],L[1],L[2],T[0],T[1],T[2],(T[1]-L[1]) if np.isfinite(T[1]-L[1]) else np.nan))
