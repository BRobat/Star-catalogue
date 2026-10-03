"""Build stars.bin + meta.json for the viewer from the AT-HYG HYGLike catalogue.

Usage: python tools/prep.py path/to/hyglike_from_athyg_v32.csv.gz   (run from the repo root)
Writes stars.bin (13 float32 columns, column-major), stars.b64.txt and meta.json.
"""
import sys, base64
import pandas as pd, numpy as np, json
src = sys.argv[1] if len(sys.argv) > 1 else 'hyglike_from_athyg_v32.csv.gz'
df = pd.read_csv(src, low_memory=False)
print(len(df)); print(df.dist_src.value_counts().head(8))
df = df[(df.dist < 99999) & df.dist.notna()]
df = df[df.absmag.notna() & df.mag.notna()]
# equatorial xyz (pc) -> galactic
R = np.array([[-0.0548755604,-0.8734370902,-0.4838350155],
              [ 0.4941094279,-0.4448296300, 0.7469822445],
              [-0.8676661490,-0.1980763734, 0.4559837762]])
eq = df[['x','y','z']].to_numpy()
g = eq @ R.T   # gx toward GC, gy toward rotation (l=90), gz NGP
ci = df.ci.fillna(0.65).clip(-0.4,2.2).to_numpy()
gaia = df.dist_src.astype(str).str.startswith('G').to_numpy().astype(np.float32)
hip = df.hip.fillna(0).to_numpy()
spects = df.spect.fillna('').astype(str)
uniq = sorted(spects.unique()); idx = {s:i for i,s in enumerate(uniq)}
sp = spects.map(idx).to_numpy()
v = df[['vx','vy','vz']].fillna(0).to_numpy() @ R.T * 1e6   # pc/Myr, galactic
hasrv = df.rv.notna().to_numpy().astype(np.float32)
cols = [g[:,0],g[:,1],g[:,2],df.mag.to_numpy(),df.absmag.to_numpy(),ci,hip,gaia,sp,v[:,0],v[:,1],v[:,2],hasrv]
# closest approaches (linear)
tmin = -(g*v).sum(1)/np.maximum((v*v).sum(1),1e-12); dmin = np.linalg.norm(g+v*tmin[:,None],axis=1)
ok=(tmin>-5)&(tmin<5)&(hasrv>0)
for k in np.argsort(np.where(ok,dmin,1e9))[:8]: print(df.iloc[k][['proper','hip','gl']].tolist(), round(tmin[k],4),'Myr', round(dmin[k],3),'pc')
arr = np.stack([c.astype(np.float32) for c in cols])  # column-major
arr.tofile('stars.bin')
names = {}
for i,(p,b,f,c,hd,gl) in enumerate(zip(df.proper, df.bayer, df.flam, df.con, df.hd, df.gl)):
    n = p if isinstance(p,str) and p else None
    alt = None
    if isinstance(b,str) and b and isinstance(c,str): alt = f"{b} {c}"
    elif isinstance(f,str) and f and isinstance(c,str): alt = f"{f} {c}"
    elif pd.notna(f) and isinstance(c,str) and str(f).strip(): alt = f"{int(float(f))} {c}"
    if not alt and isinstance(gl,str) and gl.strip(): alt = gl.strip()
    if n or alt: names[i] = [n or '', alt or '']
RENAME = {'Rigil Kentaurus': ['Alpha Centauri A', 'Rigil Kentaurus'], 'Toliman': ['Alpha Centauri B', 'Toliman']}
for k, v in names.items():
    if v[0] in RENAME: names[k] = RENAME[v[0]]
meta = {"n": len(df), "cols": ["x","y","z","mag","absmag","ci","hip","gaia","spect","vx","vy","vz","hasrv"], "spect": uniq, "names": names}
json.dump(meta, open('meta.json','w'), separators=(',',':'))
open('stars.b64.txt','w').write(base64.b64encode(open('stars.bin','rb').read()).decode())
print(len(df), len(names), len(uniq), gaia.mean(), np.percentile(df.dist,[50,90,99]))
