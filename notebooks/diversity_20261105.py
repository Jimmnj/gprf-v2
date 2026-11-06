# diversity_20261105.py [UNRATIFIED] -- PREREG_DIVERSITY_20261105
import numpy as np, pandas as pd, itertools
T=[]
for line in open('../sparc_galaxy_data/SPARC_Lelli2016c.mrt'):
    p=line.split()
    if len(p)>=18 and p[1].isdigit():
        try: T.append(dict(gal=p[0],L=float(p[7]),Rd=float(p[11]),MHI=float(p[13]),Q=int(p[17]),inc=float(p[5])))
        except ValueError: pass
T=pd.DataFrame(T).set_index('gal'); T['Mb']=0.5*T.L*1e9+1.33*T.MHI*1e9; T['lMb']=np.log10(T.Mb); T['fgas']=1.33*T.MHI*1e9/T.Mb; T['lSig']=np.log10(T.Mb/(2*np.pi*(1.678*T.Rd*1e3)**2))
A=pd.read_csv('marg_fit_params_20261030.csv'); A=A[A.arm=='A'].set_index('gal')[['lrc','lud','D','D0','i0','chi2']]
Z=pd.read_csv('marg_ns_20261030_pergal.csv').set_index('gal')[['logZ_A','logZ_B','logZ_U']]
D=T.join(A,how='inner').join(Z,how='inner'); D['dZ']=D.logZ_A-D.logZ_B; print('galaxies',len(D))
g=list(D.index); matched=[]; rnd=[]
for a,b in itertools.combinations(g,2):
    x,y=D.loc[a],D.loc[b]; d=abs(x.lrc-y.lrc)
    if abs(x.lMb-y.lMb)<=0.15 and abs(x.lSig-y.lSig)<=0.15 and abs(x.fgas-y.fgas)<=0.15: matched.append(dict(a=a,b=b,dlrc=d,lrc_a=x.lrc,lrc_b=y.lrc,dZa=x.dZ,dZb=y.dZ,dZmax=max(x.dZ,y.dZ),lMb=x.lMb,lSig=x.lSig,fgas=x.fgas))
    else: rnd.append(d)
M=pd.DataFrame(matched).sort_values('dlrc',ascending=False); rnd=np.array(rnd)
print('matched pairs %d; random pairs %d'%(len(M),len(rnd)))
print('K1: pairs with dlog rho_c >= 0.3: %d; >= 0.5: %d'%((M.dlrc>=0.3).sum(),(M.dlrc>=0.5).sum()))
k1=M[M.dlrc>=0.3]; print('K2: of those, with one member dZ >= +5: %d of %d'%((k1.dZmax>=5).sum(),len(k1)))
print('K3: median dlog rho_c matched %.3f / random %.3f = %.2f'%(M.dlrc.median(),np.median(rnd),M.dlrc.median()/np.median(rnd)))
pd.set_option('display.width',220); print(M.head(15).round(2).to_string(index=False)); M.to_csv('diversity_pairs_20261105.csv',index=False)
print('per-galaxy dZ (A-B) over all: median %+.1f; fraction dZ>=5: %.2f; fraction dZ<=-5: %.2f'%(D.dZ.median(),(D.dZ>=5).mean(),(D.dZ<=-5).mean()))
