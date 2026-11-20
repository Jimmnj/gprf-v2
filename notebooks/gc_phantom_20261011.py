#!/usr/bin/env python3
"""gc_phantom_20261011.py -- QUMOND phantom mass around Plummer GCs / dwarfs in the MW external field. [UNRATIFIED]
2026-10-11 Claude (Fable 5.1). rho_ph = (nu-1) rho_b - nu'(y) (g.T.g)/(4 pi G a0 |g|), y = |g_N|/a0, RAR nu.
g_N = Plummer field + uniform g_e along -z. Units pc, Msun, (km/s)^2/pc. Inputs RECALLED (VERIFY), as in
GC_vs_dwarfs_regimes_20261011.xlsx. Axisymmetric shells: 400-pt Gauss-Legendre in cos(theta), log radial grid."""
import numpy as np
G=4.30091e-3; A0=1.2e-10*3.0857e16/1e6
nu=lambda y:1/(1-np.exp(-np.sqrt(y)))
def nup(y):
    u=np.sqrt(y); e=np.exp(-u); return -e/((1-e)**2*2*u)
def rho_ph(x,z,M,a,ge):
    d2=x*x+z*z+a*a; d3=d2**1.5; d5=d2**2.5
    gx=-G*M*x/d3; gz=-G*M*z/d3-ge
    Txx=-G*M*(1/d3-3*x*x/d5); Tzz=-G*M*(1/d3-3*z*z/d5); Txz=G*M*3*x*z/d5
    g=np.sqrt(gx*gx+gz*gz); gTg=gx*gx*Txx+2*gx*gz*Txz+gz*gz*Tzz
    y=g/A0; rb=3*M/(4*np.pi*a**3)*(1+(x*x+z*z)/a**2)**-2.5
    return (nu(y)-1)*rb - nup(y)*gTg/(4*np.pi*G*A0*g), rb
objs=[("Typical inner GC",2e5,3,8,220),("NGC 2419",9e5,20,90,190),("Pal 14",1.3e4,25,70,190),
      ("Pal 4",3e4,23,110,185),("Segue 1 (UFD)",6e2,29,23,220),("Draco (dSph)",3e5,220,76,190)]
ct,wt=np.polynomial.legendre.leggauss(400); st=np.sqrt(1-ct*ct)
rows=[]
for n,M,rh,R,v in objs:
    a=rh/1.305; ge=v*v/(R*1e3); r0=np.sqrt(G*M/ge); cush=v*v/(4*np.pi*G*(R*1e3)**2)
    rmax=max(30*rh,30*r0); rs=np.geomspace(1e-3*rh,rmax,1500); dr=np.gradient(rs)
    cn=cp=cng=0; out={}
    Mb=lambda r: M*r**3/(r*r+a*a)**1.5
    mins=np.inf; minpos=None
    for r,d in zip(rs,dr):
        rp,rb=rho_ph(r*st,r*ct,M,a,ge); sh=2*np.pi*r*r*wt*rp
        cn+=sh.sum()*d; cp+=sh[sh>0].sum()*d; cng+=sh[sh<0].sum()*d
        if rp.min()<mins: mins=rp.min(); minpos=r
        for key,rr in (("rh",rh),("3rh",3*rh),("r0",r0)):
            if key not in out and r>=rr: out[key]=(cn,Mb(r))
    out["max"]=(cn,M)
    boost=lambda k:(out[k][1]+out[k][0])/out[k][1]
    rows.append((n,ge/A0,r0,out["rh"][0]/M,boost("rh"),boost("3rh"),out["max"][0]/M,cng/M,mins,minpos/rh,cush,abs(mins)/cush,nu(ge/A0)))
hd=("object","g_e/a0","r0 pc","Mph(<rh)/M","boost<rh","boost<3rh","Mph net/M","gross neg/M","min rho","at r/rh","cushion","|min|/cush","nu_e")
print(("{:18s}"+"{:>11s}"*12).format(*hd))
for r in rows: print(("{:18s}"+"{:11.3g}"*12).format(*r))
np.save("gc_phantom_rows.npy",np.array([r[1:] for r in rows]))
