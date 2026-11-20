# feature_axisym_20261109.py -- PREREG_FEATURE_AXISYM_20261109 (STAMPED JM 2026-11-09 13:20). NGC 1560 in disc geometry with the frozen law.
# Step 1: NNLS ring reconstruction of Sigma_gas, Sigma_disk from Gentile's V_gas, V_disk (fidelity gate sec 2).
# Step 2: control A0 (sphere on reconstructed V_bar) vs RESULT_NGC1560 R1.  Step 3: A1/A2/A3 full-mode solves; criteria sec 4. [UNRATIFIED result]
import numpy as np, io, contextlib, warnings, sys; warnings.filterwarnings('ignore')
from scipy.optimize import nnls, minimize_scalar
from scipy.special import eval_legendre
from scipy.integrate import cumulative_trapezoid as ctz, solve_ivp
import axisym_r1_20261109 as ax
from gate2b_axisym_20261004 import ring_sources, ring_phi, NMU, MU_ALL, NR, NRC
G=ax.G; a0k=ax.a0k; ramb=ax.ramb; S=ax.S
d=np.loadtxt('../data/ngc1560/ngc1560_gentile2010_rc_extracted_20261108.csv',delimiter=',')
r,vobs,lo,hi,vg,vd=d[:,0],d[:,1],d[:,2],d[:,3],d[:,4],d[:,5]; ev=(hi-lo)/2
vb2=np.sign(vg)*vg**2+np.sign(vd)*vd**2; vb2c=np.clip(vb2,0,None)
MHI=2.36e5*3.45**2*294.6*1.33; MSTAR={'primary':1.3e9,'lo':1.0e9,'hi':1.6e9}; RD=1.3; HZ=0.196*RD**0.633; HG=0.2
# ---- step 1: ring reconstruction
Rr=np.geomspace(0.1,25,120)   # AMENDMENT FA-4 (13:30): 120 rings to 25 kpc + soft total-gas-mass constraint
def ring_v2_matrix(Rr,h):
    """in-plane v^2 at data radii r from unit-Sigma rings (each ring = ring_sources with Sigma=1 on that ring only)"""
    A=np.zeros((len(r),len(Rr)))
    for j in range(len(Rr)):
        sig=np.zeros(len(Rr)); sig[j]=1.0
        Rs,zs,ms,_=ring_sources(Rr,sig,h)
        rq=np.concatenate([r*0.999,r*1.001])
        phi=ring_phi(rq,np.zeros_like(rq),Rs,zs,ms,h/2)
        A[:,j]=r*(phi[len(r):]-phi[:len(r)])/(r*0.002)
    return A
edges=np.concatenate([[0],np.sqrt(Rr[1:]*Rr[:-1]),[Rr[-1]*np.sqrt(Rr[-1]/Rr[-2])]]); ring_area=np.pi*np.diff(edges**2)
def recon(v,h,label,Mtarget=None):
    A=ring_v2_matrix(Rr,h); y=np.sign(v)*v**2; rows=[A]; ys=[y]
    if Mtarget: rows.append((ring_area*1e-7)[None,:]*5); ys.append(np.array([Mtarget*1e-7*5]))
    sig,res=nnls(np.vstack(rows),np.concatenate(ys))
    vfit=np.sqrt(np.clip(A@sig,0,None)); m=(r>=1)&(r<=9.5)
    rms=np.sqrt(np.mean((vfit[m]-np.abs(v[m]))**2))
    return sig,vfit,rms
sg,vgf,rg=recon(vg,HG,'gas',Mtarget=MHI); sd,vdf,rd=recon(vd,HZ,'disk')
pk_in=r[(r>1)&(r<6)][np.argmax(vg[(r>1)&(r<6)])]; pk_fit=r[(r>1)&(r<6)][np.argmax(vgf[(r>1)&(r<6)])]

Mg_rec=np.sum(sg*ring_area)
print('STEP 1 ring reconstruction: gas rms %.2f km/s (peak in %.2f -> fit %.2f kpc) ; disk rms %.2f km/s ; M_gas(rings) %.2e vs 1.33 M_HI %.2e'%(rg,pk_in,pk_fit,rd,Mg_rec,MHI))
GATE = rg<1.0 and rd<1.0 and abs(pk_fit-pk_in)<0.3
print('FIDELITY GATE:', 'PASS' if GATE else 'FAIL')
np.savetxt('feature_axisym_sigma_20261109.csv',np.column_stack([Rr,sg,sd]),delimiter=',',header='R_kpc,Sigma_gas_Msun_pc2(x1e-6 of Msun/kpc2 units: values are Msun/kpc2),Sigma_disk_Msun_kpc2',comments='')
if not GATE: sys.exit('baryon model failed -- nothing scored')
# ---- geometry: baryon potential drop on the (r, mu) grid, phi_b(0)
R_out=2.0*r[-1]; rr=np.concatenate([[0],np.geomspace(0.005,R_out,NR-1)])
rc=np.concatenate([[0],np.geomspace(0.01,R_out,NRC)]); Rf=(rc[:,None]*np.sqrt(1-MU_ALL**2)).ravel(); zf=(rc[:,None]*MU_ALL).ravel()
Rg_,zg_,mg_,_=ring_sources(Rr,sg,HG)
def disc_pot(mstar):
    scale=mstar/np.sum(sd*ring_area)
    Rd_,zd_,md_,_=ring_sources(Rr,sd*scale,HZ)
    phi=(ring_phi(Rf,zf,Rd_,zd_,md_,HZ/2)+ring_phi(Rf,zf,Rg_,zg_,mg_,0.1)).reshape(len(rc),len(MU_ALL))
    psic=phi-phi[0,0]; psib=np.column_stack([np.interp(rr,rc,psic[:,j]) for j in range(len(MU_ALL))])
    return psib,-phi[0,0],scale
# in-plane baryon v^2 for the model: Gentile's V_bar^2 with the stellar part rescaled to M_* (PREREG sec 1: M_* bracket moves sigma AND v_disk)
def vbar2_of(mstar): return np.clip(np.sign(vg)*vg**2+np.sign(vd)*vd**2*(mstar/1.3e9),0,None)
def solve(mode,lrc,mstar):
    Mb=MHI+mstar; s2=np.sqrt(a0k*G*Mb)/2; vb2m=vbar2_of(mstar)
    if mode=='C0':   # sphere on the RING-MODEL in-plane v_bar (control A0 uses Gentile's vb2c directly below)
        psib,phib0,_=disc_pot(mstar); lmax=0
    else:
        psib,phib0,_=disc_pot(mstar); lmax=12
    M,p0,Phil,dPhil,Sb,psi,it=ax.solve_halo(rr,psib,phib0,s2,lrc,ramb,lmax)
    vnb2=G*np.interp(r,rr,M)/r
    ls=np.arange(2,lmax+1,2)
    if len(ls): vnb2=vnb2+r*np.interp(r,rr,(dPhil*eval_legendre(ls,0.0)[:,None]).sum(0))
    q2=np.nan
    if lmax:
        ps=-(Sb-psi/s2); R0=2*RD; tgt=np.interp(R0,rr,ps[:,NMU]); axv=np.maximum.accumulate(ps[:,NMU+1]); q2=np.interp(tgt,axv,rr)/R0
    return np.sqrt(np.clip(vb2m+vnb2,1,None)),it,q2
# ---- step 2: control. A0 = spherical solve on the reconstructed in-plane v_bar (same code path as ngc1560_solve relation())
def sphere(lrc,Mb,vb2in):
    rl=r[-1]; rr2=np.concatenate([[0],r]); vv=np.concatenate([[0],vb2in])
    vbar=lambda x: np.interp(x,rr2,vv) if x<=rl else vb2in[-1]*rl/x
    xs=np.linspace(1e-4,rl,4001); phib0=np.trapezoid(np.interp(xs,rr2,vv)/xs,xs)+vb2in[-1]
    s2=np.sqrt(a0k*G*Mb)/2; xb0=phib0/s2; Sb0=S(xb0); rcc=10**lrc
    fn=lambda x,y:[vbar(max(x,1e-4))/max(x,1e-4), G*y[2]/max(x,1e-4)**2, 4*np.pi*x*x*max(rcc*np.exp(S(xb0-y[0]/s2)-Sb0-y[1]/s2)-ramb,-ramb)]
    so=solve_ivp(fn,[1e-4,rl],[0,0,0],t_eval=r,rtol=1e-6,atol=1e-3,max_step=rl/200); return np.sqrt(np.clip(vb2in+G*so.y[2]/r,1,None))
ref=np.loadtxt('ngc1560_curves_20261108.csv',delimiter=',',skiprows=1); vR1=ref[:,4]; vR3=ref[:,6]
vb2_rings=np.clip(vgf**2+vdf**2,0,None)
vA0=sphere(7.45,MHI+1.3e9,vb2_rings); vR1_re=sphere(7.45,MHI+1.3e9,vb2c)
print('STEP 2 control: R1 re-solved vs stored max|dv| %.3f km/s ; A0 (sphere on ring-model v_bar) vs R1 max|dv| %.3f km/s (gate < 0.5)'%(np.abs(vR1_re-vR1).max(),np.abs(vA0-vR1).max()))
# ---- step 3
w=(r>=2)&(r<=8); fw=(r[w]>=3)&(r[w]<=6.5)
def detrend(v,w): c=np.polyfit(np.log10(r[w]),np.log10(v[w]),2); return np.log10(v[w])-np.polyval(c,np.log10(r[w]))
def feature(v): f=detrend(v,w); return f[fw].max()-f[fw].min(), r[w][fw][f[fw].argmax()]
def rms(v): return np.sqrt(np.mean(np.log10(vobs/v)**2))
def chi2(v): return np.sum(((v-vobs)/ev)**2)
vA1,it1,q2=solve('full',7.45,1.3e9)
fit=minimize_scalar(lambda l: chi2(solve('full',l,1.3e9)[0]),bounds=(6,10),method='bounded',options={'xatol':0.01}).x
vA2,_,_=solve('full',fit,1.3e9)
vA3lo,_,_=solve('full',7.45,1.0e9); vA3hi,_,_=solve('full',7.45,1.6e9)
print('\nSTEP 3 (full, l<=12; iters A1 %d; nb isodensity q(2Rd) %.2f)'%(it1,q2))
print('  %-16s rms(dex)  chi2    A(feature)  r_peak'%'curve')
for nm,v in [('V_obs',vobs),('R3 RARfn',vR3),('R1 sphere',vR1),('A0 sphere/rings',vA0),('A1 disc 7.45',vA1),('A2 disc fit %.2f'%fit,vA2),('A3 disc M*1.0e9',vA3lo),('A3 disc M*1.6e9',vA3hi)]:
    A,rp=feature(v); print('  %-16s %.4f    %7.1f  %.4f      %.2f'%(nm,rms(v) if nm!='V_obs' else 0,chi2(v) if nm!='V_obs' else 0,A,rp))
fo=detrend(vobs,w); print('  Pearson detrended: A1 vs obs %.2f ; R1 vs obs %.2f ; R3 vs obs %.2f'%(np.corrcoef(detrend(vA1,w),fo)[0,1],np.corrcoef(detrend(vR1,w),fo)[0,1],np.corrcoef(detrend(vR3,w),fo)[0,1]))
AA1=feature(vA1)[0]; AR1=feature(vR1)[0]; AR3=feature(vR3)[0]
verdict='(i) GEOMETRY RESTORES' if AA1>=0.8*AR3 else ('(ii) HELPS, LAW STILL DAMPS' if AA1>=1.5*AR1 else '(iii) GEOMETRY IRRELEVANT')
print('  A(A1)/A(R3) = %.2f ; A(A1)/A(R1 sph) = %.2f  -> %s'%(AA1/AR3,AA1/AR1,verdict))
np.savetxt('feature_axisym_curves_20261109.csv',np.column_stack([r,vobs,ev,vR1,vR3,vA0,vA1,vA2,vA3lo,vA3hi]),delimiter=',',header='r,vobs,ev,R1_sph,R3_RAR,A0_sph_rings,A1_disc745,A2_disc_fit%.2f,A3_lo,A3_hi'%fit,comments='')
