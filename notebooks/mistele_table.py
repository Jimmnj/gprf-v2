# mistele_table.py -- byte-position parser for Mistele+2024 Table 1 (MRT). Replaces whitespace split(), which shifted values into the wrong
# mass bins on rows with blank columns (1836.51 and 2496.47 kpc). Fixed 2026-10-22 after the simulated Mistele review.
import pandas as pd
COLS=[(33,40,47),(54,61,68),(75,82,89),(96,103,110)]   # start of Vc, stat, syst per bin (0-based); each field 6 chars
def load(sample='All',path='../data/mistele2024_table1_mrt.txt'):
    rows=[]
    for line in open(path):
        if line[:24].strip()!=sample: continue
        try: R=float(line[25:32])
        except ValueError: continue
        rec={'R':R}
        for b,(c,s,y) in enumerate(COLS,1):
            f=line[c:c+6].strip()
            if f: rec['v%d'%b]=float(f); rec['s%d'%b]=float(line[s:s+6]); rec['y%d'%b]=float(line[y:y+6])
        rows.append(rec)
    return pd.DataFrame(rows)
if __name__=='__main__':
    T=load(); print(T.tail(3).to_string())
