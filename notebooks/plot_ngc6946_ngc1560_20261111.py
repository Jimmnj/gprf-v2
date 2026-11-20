# plot_ngc6946_ngc1560_20261111.py -- figure for JM 2026-11-11 14:57 ("run both and gen a plot"). [UNRATIFIED; display]
# Inputs: ngc6946_curves_20261111.csv; ngc1560_curves_supp_D2.94_ML0.98_20261111.csv; ngc1560_curves_supp_D3.45_ML1.43_20261111.csv.
# Also prints the PREREG_NGC1560 sec 13 checks K3, K4 (post-hoc shape metric of CHECK_NP17 H5 on the same points).
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def load(f):
    h=open(f).readline().strip().split(','); d=np.loadtxt(f,delimiter=',',skiprows=1); return h,{k:d[:,i] for i,k in enumerate(h)}
INK='#0b0b0b'; MUT='#898781'; GRID='#e1e0d9'; BLUE='#2a78d6'; ORANGE='#eb6834'; AQUA='#1baf7a'
panels=[('ngc6946_curves_20261111.csv','NGC 6946\nSPARC; M/L 0.5 disc, 0.7 bulge; catalogue D, i',(0,21),(0,215),(-30,50)),
        ('ngc1560_curves_supp_D2.94_ML0.98_20261111.csv','NGC 1560\n2.94 Mpc, M/L 0.98 (Gentile et al. MOND-fit values)',(0,8.6),(0,92),(-12,18)),
        ('ngc1560_curves_supp_D3.45_ML1.43_20261111.csv','NGC 1560\n3.45 Mpc, M/L 1.43 (pre-registered)',(0,10.1),(0,92),(-12,18))]
plt.rcParams.update({'font.size':10,'axes.edgecolor':'#c3c2b7','axes.labelcolor':INK,'xtick.color':MUT,'ytick.color':MUT,'axes.titlesize':10.5,'font.family':'DejaVu Sans'})
fig,ax=plt.subplots(2,3,figsize=(15.5,8.2),gridspec_kw={'height_ratios':[2.1,1],'hspace':0.08,'wspace':0.2},sharex='col')
for j,(f,title,xl,yl,rl) in enumerate(panels):
    h,d=load(f); r=d['r']; fitk=[k for k in h if k.startswith('R2_fit')][0]; lrc=fitk.split('lrc')[1]
    sel=d['sel'].astype(bool) if 'sel' in d else np.ones(len(r),bool)
    a=ax[0,j]; b=ax[1,j]
    a.errorbar(r[sel],d['vobs'][sel],yerr=d['ev'][sel],fmt='o',ms=4.2,mfc='white',mec=INK,mew=0.9,ecolor=MUT,elinewidth=0.8,capsize=0,zorder=3,label='Data')
    if (~sel).any(): a.errorbar(r[~sel],d['vobs'][~sel],yerr=d['ev'][~sel],fmt='o',ms=4.2,mfc='white',mec=MUT,mew=0.9,ecolor=GRID,elinewidth=0.8,zorder=2,label='Data, cut (error > 10%)')
    a.plot(r,d['vbar'],color=MUT,lw=1.4,label='Newtonian baryons',zorder=1)
    a.plot(r,d['R3_RAR'],color=BLUE,lw=2,solid_capstyle='round',label='Acceleration relation',zorder=4)
    a.plot(r,d['R1_745'],color=ORANGE,lw=2,ls=(0,(5,2.5)),label=r'Eq. (1), universal $\rho_c=10^{7.45}$',zorder=5)
    a.plot(r,d[fitk],color=AQUA,lw=2,ls=(0,(1.2,1.8)),label=r'Eq. (1), fitted $\rho_c=10^{%s}$'%lrc,zorder=6)
    for k,c,ls in [('R3_RAR',BLUE,'-'),('R1_745',ORANGE,(0,(5,2.5))),(fitk,AQUA,(0,(1.2,1.8)))]:
        b.plot(r,d[k]-d['vobs'],color=c,lw=2,ls=ls)
    b.axhline(0,color='#c3c2b7',lw=1.2,zorder=0)
    a.set_title(title,loc='left',color=INK); a.set_xlim(*xl); a.set_ylim(*yl); b.set_ylim(*rl)
    for q in (a,b):
        q.grid(axis='y',color=GRID,lw=0.8); q.set_axisbelow(True)
        for s in ('top','right'): q.spines[s].set_visible(False)
    b.set_xlabel('Radius (kpc)')
    if j==0: a.set_ylabel('Speed (km/s)'); b.set_ylabel('Model minus data (km/s)')
    a.legend(loc='lower right',frameon=False,fontsize=9,labelcolor='#52514e',handlelength=2.6)
fig.savefig('fig_ngc6946_ngc1560_curves_20261111.png',dpi=130,bbox_inches='tight',facecolor='white')
# --- PREREG_NGC1560 sec 13 checks
h0,d0=load('ngc1560_curves_supp_D3.45_ML1.43_20261111.csv'); idx=(d0['r']>=4.0)&(d0['r']<=6.03); idx2=(d0['r']>=6.0)&(d0['r']<=6.91)
for f,tag in [('ngc1560_curves_supp_D3.45_ML1.43_20261111.csv','3.45,1.43'),('ngc1560_curves_supp_D2.94_ML0.98_20261111.csv','2.94,0.98')]:
    h,d=load(f); r=d['r']; fitk=[k for k in h if k.startswith('R2_fit')][0]
    sl=lambda v,m: np.polyfit(r[m],v[m],1)[0]
    p=np.polyfit(r[idx],d['vobs'][idx],1,w=1/d['ev'][idx],cov=True)
    print('K3 (%s) shelf slope km/s/kpc, same %d points: obs %.2f (weighted %.2f +- %.2f) | R3 %.2f | R1 %.2f | R2 %.2f || step: obs %.2f R3 %.2f R1 %.2f R2 %.2f'%(tag,idx.sum(),sl(d['vobs'],idx),p[0][0],np.sqrt(p[1][0,0]),sl(d['R3_RAR'],idx),sl(d['R1_745'],idx),sl(d[fitk],idx),sl(d['vobs'],idx2),sl(d['R3_RAR'],idx2),sl(d['R1_745'],idx2),sl(d[fitk],idx2)))
    o=d0['r']>2.0
    for k in ('R3_RAR','R1_745',fitk):
        df=(d[k]-d['vobs'])[o]; print('   K4 (%s) %-16s model-obs beyond 2 kpc(at 3.45): median %+.1f  range %+.1f..%+.1f km/s'%(tag,k,np.median(df),df.min(),df.max()))
