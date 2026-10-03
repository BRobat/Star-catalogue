"""Export a star-lane graph (JSON) from stars.bin + meta.json.

Usage: python tools/lanes.py [radius_pc] [max_jump_pc] [links_per_star]   (run from the repo root)
"""
import sys
import numpy as np, json
from scipy.spatial import cKDTree
m=json.load(open('meta.json'));n=m['n'];a=np.fromfile('stars.bin',dtype=np.float32).reshape(-1,n)
x,y,z,mag,absm,ci,hip,gaia,sp,vx,vy,vz,hasrv=a
R=float(sys.argv[1]) if len(sys.argv)>1 else 25
J=float(sys.argv[2]) if len(sys.argv)>2 else 3.5
K=int(sys.argv[3]) if len(sys.argv)>3 else 3
d=np.sqrt(x*x+y*y+z*z); idx=np.where(d<=R)[0]
P=np.stack([x[idx],y[idx],z[idx]],1); tree=cKDTree(P)
dd,nn=tree.query(P,k=K+1,distance_upper_bound=J)
edges={}
for i in range(len(idx)):
    for dist,j in zip(dd[i][1:],nn[i][1:]):
        if np.isfinite(dist): edges[(min(i,j),max(i,j))]=round(float(dist),3)
# components
par=list(range(len(idx)))
def f(u):
    while par[u]!=u: par[u]=par[par[u]]; u=par[u]
    return u
for (u,v) in edges: par[f(u)]=f(v)
roots=[f(i) for i in range(len(idx))]; rid={r:k for k,r in enumerate(dict.fromkeys(roots))}
def T(bv): return round(4600*(1/(0.92*bv+1.7)+1/(0.92*bv+0.62))/10)*10
stars=[]
for k,i in enumerate(idx):
    nm=m['names'].get(str(i),['',''])
    name=nm[0] or nm[1] or (f"HIP {int(hip[i])}" if hip[i]>0 else f"AT-HYG #{i}")
    o={"id":k,"name":name}
    if nm[0] and nm[1]: o["alias"]=nm[1]
    if hip[i]>0: o["hip"]=int(hip[i])
    o.update(x=round(float(x[i]),4),y=round(float(y[i]),4),z=round(float(z[i]),4),
      vx=round(float(vx[i]/1.02271),2),vy=round(float(vy[i]/1.02271),2),vz=round(float(vz[i]/1.02271),2),
      absmag=round(float(absm[i]),2),tempK=T(float(ci[i])),cluster=rid[roots[k]])
    s=m['spect'][int(sp[i])]
    if s: o["spect"]=s
    stars.append(o)
out={"meta":{"source":"AT-HYG v3.2 (HYGLike subset), Gaia DR3 distances and motions where available","units":{"position":"parsec","velocity":"km/s"},
 "frame":"heliocentric galactic: +x toward galactic centre, +y toward galactic rotation, +z toward north galactic pole","epoch":"J2016",
 "regionRadiusPc":R,"maxJumpPc":J,"linksPerStar":K,"note":"Each star linked to its K nearest neighbours within maxJumpPc. Catalogue misses most faint red dwarfs beyond ~20 pc."},
 "stars":stars,"lanes":[[int(u),int(v),dd_] for (u,v),dd_ in edges.items()]}
json.dump(out,open(f'exports/star_lanes_{R:g}pc.json','w'),separators=(',',':'))
from collections import Counter
c=Counter(roots); print(len(stars),len(edges),len(c),max(c.values())/len(stars))
