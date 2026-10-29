#!/usr/bin/env python3
"""Fetch the SPARC database (Lelli, McGaugh & Schombert 2016, AJ 152, 157) into ../sparc_galaxy_data/,
the layout the scripts in notebooks/ expect. SPARC is not redistributed here; cite Lelli et al. (2016a) when using it.
Sources tried in order: the Zenodo deposit (10.5281/zenodo.16284118, CC BY 4.0), then the SPARC site.
MD5 sums are those of the copies used for v2.0.0; a mismatch is reported, not fatal."""
import os, sys, hashlib, zipfile, urllib.request
root=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','sparc_galaxy_data')
FILES={'SPARC_Lelli2016c.mrt':'6181df386bfc05868a3700c196e800da',
       'sparc_database.zip':'471e2b38a151c3b63b05021bf985b45f',
       'Rotmod_LTG.zip':'e4c8b92766026770ed35e5889064e12b',
       'MassModels_Lelli2016c.mrt':'fe6188538c3f5504f70f486ff6b4d29c'}
BASES=['https://zenodo.org/records/16284118/files/%s?download=1','http://astroweb.cwru.edu/SPARC/%s']
os.makedirs(root,exist_ok=True)
for name,md5 in FILES.items():
    dst=os.path.join(root,name)
    if not os.path.exists(dst):
        for b in BASES:
            try: urllib.request.urlretrieve(b%name,dst); print('got',name,'from',b.split('/%s')[0]); break
            except Exception as e: print('  failed',b%name,'--',e)
        else: sys.exit('could not fetch '+name)
    h=hashlib.md5(open(dst,'rb').read()).hexdigest()
    print('  %-28s md5 %s %s'%(name,h,'OK' if h==md5 else 'DIFFERS from the v2.0.0 copy'))
z=os.path.join(root,'sparc_database.zip'); out=os.path.join(root,'sparc_database')
if not os.path.isdir(out) or not os.listdir(out):
    os.makedirs(out,exist_ok=True)
    with zipfile.ZipFile(z) as f:
        for m in f.namelist():
            if m.endswith('_rotmod.dat'): open(os.path.join(out,os.path.basename(m)),'wb').write(f.read(m))
print(len([x for x in os.listdir(out) if x.endswith('_rotmod.dat')]),'rotation-curve files in sparc_galaxy_data/sparc_database (expect 175)')
