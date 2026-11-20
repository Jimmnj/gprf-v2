# errprop_grid_20261112.py -- PREREG_ERRPROP_20261112 sec 7 grid check (declared 21:55 before the run). Solver and chi^2 verbatim via errprop_20261112.fit internals. [UNRATIFIED]
import sys, numpy as np, pandas as pd
sys.argv=['x','--validate','0']
exec(open('errprop_20261112.py').read().split("args=sys.argv[1:]")[0])
def chi_curve(g,grid):
    d0=P0[P0.gal==g].sort_values('r'); t0=T0.loc[g]; m=((d0.ev/d0.vobs<0.1)&(d0.vobs>0)&(d0.vbar2>0)).values
    ns['P']=d0; ns['T']=pd.DataFrame([t0],index=[g]); vo=d0.vobs.values; w=1/(0.1*vo)**2; out=[]
    for l in grid:
        o=solve_gal(g,l); out.append(1e12 if o is None else float(np.sum((w*(np.sqrt(o.vmod2.values)-vo)**2)[m])))
    return np.array(out)
grid=np.round(np.arange(4.5,10.5001,0.05),2)
for g in sys.stdin.read().split():
    c=chi_curve(g,grid); i=np.argmin(c); ref=pd.read_csv('pass2_table2_part.csv').set_index('gal').lrc[g]
    j=np.argmin(abs(grid-ref)); lm=[grid[k] for k in range(1,len(c)-1) if c[k]<c[k-1] and c[k]<=c[k+1]]
    print('%-9s grid min %.2f chi2 %.1f | at pass2 value %.2f chi2 %.1f | chi2(7.45) %.1f | local minima %s'%(g,grid[i],c[i],ref,c[j],c[np.argmin(abs(grid-7.45))],lm[:6]),flush=True)
