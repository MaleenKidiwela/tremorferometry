#!/usr/bin/env python
"""COMBINED explorer: slow-slip model products (Gualandi 2025 via cascadia-sse-strain, in our 30-d trailing frame) AND
the interface delta-beta/beta movie (2-4 s daily joint inversion; production operator and the corrected S-pinned
differential-kernel operator) + near-station site terms, on ONE Leaflet map with ONE weekly time slider.
Layers: model variable (grid rectangles, 24 products) / interface db/b cells (production or G_diff) / site dv/v
(stations); any combination. Click a cell or station -> Plotly series of its dv/v (daily + 90-d mean) with the selected
model variable at the nearest model node on a twin axis; click a model node -> model series with the nearest cell's db/b.
Era lines (whitelisted instrument changes) on station series. -> figures/sse_dvv_combined.html"""
import sys, json, gzip, base64, numpy as np, pandas as pd, xarray as xr
from scipy.spatial import cKDTree
REPO = "/home/jovyan/cascadia-sse-strain"; sys.path.insert(0, f"{REPO}/src")
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
# ---------------- model grid products (as in build_sse_model_map.py) ----------------
ds = xr.open_dataset(f"{SD}/megathrust_strain_trailing30.nc"); lat, lon = ds.lat.values, ds.lon.values; times = pd.to_datetime(ds.time.values)
LA, LO = np.meshgrid(lat, lon, indexing="ij"); valid = np.isfinite(ds["dilatation"].isel(depth=0).values).any(0); ii, jj = np.where(valid); npts = len(ii)
UNITS = {"strain": ("1e-8", 1e-8), "stress": ("kPa", 1e3), "dilatation": ("1e-8", 1e-8), "eps_eq": ("1e-8", 1e-8), "max_shear_strain": ("1e-8", 1e-8), "mean_stress": ("kPa", 1e3), "pressure": ("kPa", 1e3),
         "von_mises_stress": ("kPa", 1e3), "max_shear_stress": ("kPa", 1e3), "coulomb": ("kPa", 1e3), "shear_stress": ("kPa", 1e3), "normal_stress": ("kPa", 1e3)}
def unit_of(v):
    for k, u in UNITS.items():
        if v == k or v.startswith(k + "_"): return u
    return ("", 1.0)
order = ["slip_thrust_inc30", "slip_thrust_cum", "coulomb", "shear_stress", "normal_stress", "dilatation", "eps_eq", "max_shear_strain", "pressure", "mean_stress", "von_mises_stress", "max_shear_stress"] + \
        [f"strain_{c}" for c in ("exx", "eyy", "ezz", "exy", "exz", "eyz")] + [f"stress_s{c}" for c in ("xx", "yy", "zz", "xy", "xz", "yz")]
from sse_strain import load_fault_ascii, load_solution
from sse_strain.geodesy import llh2localxy
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False); Ls, Ld = sol.component_slip_patterns()
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); A_trail = zc["A_trail"].astype(float); A_raw = zc["A_raw"].astype(float); adays = pd.to_datetime(zc["days"])
cll = fault.centroid_llh; px, py = llh2localxy(cll[:, 1], cll[:, 0], fault.origin); gx, gy = llh2localxy(LA[valid], LO[valid], fault.origin)
_, ip = cKDTree(np.column_stack([px, py])).query(np.column_stack([gx, gy]))
thr = -(Ld @ np.nan_to_num(A_trail).T)[ip]; thr_raw = -(Ld @ np.nan_to_num(A_raw).T)[ip]
tidx = adays.get_indexer(times, method="nearest"); cum = thr[:, tidx]; inc = np.full_like(cum, np.nan); t30 = adays.get_indexer(times - pd.Timedelta(days=30), method="nearest"); ok = t30 >= 0; inc[:, ok] = thr_raw[:, tidx[ok]] - thr_raw[:, t30[ok]]
extra = {"slip_thrust_cum": (cum.T * 1000, "mm"), "slip_thrust_inc30": (inc.T * 1000, "mm / 30 d")}
LABEL = {"slip_thrust_inc30": "thrust slip within the trailing 30-d window", "slip_thrust_cum": "cumulative thrust slip (rel. 2019-24)"}
def pack_i16(arr):
    vmax = float(np.nanpercentile(np.abs(arr), 99.5)) or 1.0; q = np.round(np.clip(arr / vmax, -1, 1) * 32000).astype(np.int16); q[np.isnan(arr)] = -32768
    return base64.b64encode(gzip.compress(q.tobytes(), 6)).decode(), vmax
VARS, PAY = {}, {}
for v in order:
    if v in extra: arr, unit = extra[v]; arr = arr.astype(float)
    else: unit, sc = unit_of(v); arr = ds[v].isel(depth=0).values[:, ii, jj].astype(float) / sc
    PAY[v], vmax = pack_i16(arr); VARS[v] = dict(unit=unit, vmax=vmax, label=LABEL.get(v, v.replace("_", " ")), sym=bool(v != "slip_thrust_inc30"))
print(f"model: {npts} nodes x {len(times)} weekly frames, {len(order)} variables", flush=True)
# ---------------- interface db/b products on the 872 cells: weekly frames for the map + daily series for clicks ----------------
zp = np.load(f"{RD}/inversion_4d_daily.npz", allow_pickle=True); zq = np.load("fault_tomography/inversion/splits/Gdiff24/inversion_4d_daily.npz", allow_pickle=True)
mdays = pd.to_datetime(zp["days"]); cells = zp["cell"]; clat, clon, cdep = zp["clat"], zp["clon"], zp["depth"]; M = len(cells)
dgrid = pd.date_range(mdays.min(), mdays.max(), freq="D"); DPROD = {}
# per-cell transmission correction (raw / T_cell), the SAME construction as interface_dvv_explorer.html: T_cell = pattern-
# weighted mean of the 20 patch transmissions (time-median of transmission_T.npz T[t,p], +1 % Gaussian-taper blob through
# day_solve). Production: the table recovered from that page (data/interface_cell_T_from_explorer.csv; reproduced to 1.2 %
# by the construction). G_diff re-inversion: splits/Gdiff24/cell_T_patchmapped.csv (same construction on its own
# transmission_T.npz, 2026-09-08). Cells with T_cell <= 0.2 are left uncorrected. NOTE: single-cell self-recovery
# (cell_recovery.csv, cell_recovery_table.py) is a different, much smaller quantity and is NOT what the explorer applies.
def tcell(path):
    Tcsv = pd.read_csv(path).set_index("cell").T_cell; Tc = np.array([Tcsv.get(str(c), np.nan) for c in cells]); return np.where(np.isfinite(Tc) & (Tc > 0.2), Tc, 1.0)
TCELL = {"production": tcell("data/interface_cell_T_from_explorer.csv"), "gdiff": tcell("fault_tomography/inversion/splits/Gdiff24/cell_T_patchmapped.csv")}
for name, z, Tc in (("production, transmission-corrected (as interface explorer)", zp, TCELL["production"]), ("corrected kernel (G_diff), transmission-corrected", zq, TCELL["gdiff"])):
    D = pd.DataFrame(np.asarray(z["MODEL"], float), index=pd.to_datetime(z["days"])).reindex(dgrid).values / Tc[None, :]
    wk = pd.DataFrame(D, index=dgrid).reindex(times, method="nearest", tolerance=pd.Timedelta(days=4)).values
    DPROD[name] = dict(daily=pack_i16(D), weekly=pack_i16(wk))
S = pd.DataFrame(np.asarray(zp["SITES"], float), index=mdays, columns=list(zp["sta_tags"])).reindex(dgrid)
sll = pd.read_csv(f"{RD}/pairs.csv").groupby("tag")[["sta_lat", "sta_lon"]].first(); stas = [t for t in S.columns if t in sll.index]; S = S[stas]
Swk = S.reindex(times, method="nearest", tolerance=pd.Timedelta(days=4)).values
SITE = dict(daily=pack_i16(S.values), weekly=pack_i16(Swk), tags=stas, lat=[round(float(sll.sta_lat[t]), 4) for t in stas], lon=[round(float(sll.sta_lon[t]), 4) for t in stas])
et = pd.read_csv(f"{RD}/era_table.csv"); wl = {x.strip().upper() for x in open("data/era_split_whitelist.txt") if x.strip()}
ERAS = {r.tag: [x for x in str(r.boundaries).split(";") if x and x != "nan"] for _, r in et.iterrows() if str(r.tag).upper().replace("P90F40", "") in wl and isinstance(r.boundaries, str)}
# nearest model node per cell and nearest cell per model node (for paired series)
mxy = np.column_stack([(LO[valid] + 124) * 111 * np.cos(np.radians(47)), (LA[valid] - 47) * 111]); cxy = np.column_stack([(clon + 124) * 111 * np.cos(np.radians(47)), (clat - 47) * 111])
c2m = cKDTree(mxy).query(cxy)[1]; m2c = cKDTree(cxy).query(mxy)[1]; s2m = cKDTree(mxy).query(np.column_stack([(np.array(SITE["lon"]) + 124) * 111 * np.cos(np.radians(47)), (np.array(SITE["lat"]) - 47) * 111]))[1]
meta = dict(lat=[round(float(x), 3) for x in LA[valid]], lon=[round(float(x), 3) for x in LO[valid]], times=[t.strftime("%Y-%m-%d") for t in times], iface=[round(float(x), 1) for x in ds.interface_depth.values[valid]],
            vars=VARS, order=order, n=npts, dlat=float(np.diff(lat).mean()), dlon=float(np.diff(lon).mean()),
            cells=dict(id=[str(c) for c in cells], lat=[round(float(x), 3) for x in clat], lon=[round(float(x), 3) for x in clon], dep=[round(float(x), 1) for x in cdep], n=M, c2m=[int(x) for x in c2m], m2c=[int(x) for x in m2c],
                       days=[d.strftime("%Y-%m-%d") for d in dgrid], nd=len(dgrid), prods={k: dict(vmax_daily=v["daily"][1], vmax_weekly=v["weekly"][1]) for k, v in DPROD.items()}),
            sites=dict(tags=SITE["tags"], lat=SITE["lat"], lon=SITE["lon"], n=len(stas), vmax_daily=SITE["daily"][1], vmax_weekly=SITE["weekly"][1], s2m=[int(x) for x in s2m], eras=ERAS))
payload = dict(model=PAY, cells={k: dict(daily=v["daily"][0], weekly=v["weekly"][0]) for k, v in DPROD.items()}, sites=dict(daily=SITE["daily"][0], weekly=SITE["weekly"][0]))
HTML = r'''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia: slow-slip model + interface δβ/β</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/><script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script><script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
<style>:root{--bg:#f1f5f9;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb}
*{box-sizing:border-box}html,body{margin:0;height:100%;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,sans-serif;color:var(--ink);background:var(--bg)}
#wrap{display:flex;height:100vh}#map{flex:1;height:100%;min-width:200px}#drag{width:7px;flex:none;cursor:col-resize;background:var(--line)}#drag:hover,#drag.on{background:var(--accent)}
#side{width:600px;flex:none;min-width:380px;display:flex;flex-direction:column;border-left:1px solid var(--line);background:var(--bg);overflow:auto}
#hdr{padding:10px 16px;background:var(--card);border-bottom:1px solid var(--line)}.title{font-size:15px;font-weight:700}.sub{font-size:11.5px;color:var(--muted);margin-top:3px;line-height:1.45}
.panel{margin:10px 12px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}.lbl{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
select,button{font:inherit;font-size:12.5px;padding:4px 8px;border:1px solid var(--line);border-radius:6px;background:#fff}button.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.row{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-top:6px}#plot{flex:1;min-height:320px;margin:10px 12px;background:#fff;border:1px solid var(--line);border-radius:10px}
#info{font-size:12px;color:var(--muted);padding:6px 16px;line-height:1.5}.cbar{height:12px;border-radius:4px;margin-top:4px;border:1px solid var(--line)}.cbl{display:flex;justify-content:space-between;font-size:10.5px;color:var(--muted)}
.leaflet-container{background:#e5e9f0}</style></head><body><div id="wrap"><div id="map"></div><div id="drag" title="drag to widen the panel"></div><aside id="side">
<div id="hdr"><div class="title">Cascadia: slow-slip model (Gualandi 2025) + LFE-coda interface δβ/β</div>
<div class="sub">Both in the dv/v measurement frame: 30-day trailing mean, reference 2019&ndash;24. Model on its own grid (10 km above the interface, 21-km low-pass); δβ/β on the 872 LFE cells (2&ndash;4 s daily joint inversion); site dv/v at stations. One weekly slider; click any node, cell or station for its series paired with the other dataset.</div></div>
<div class="panel"><div class="row"><span class="lbl">layers</span><button id="lyM" class="on">model</button><button id="lyC" class="on">interface δβ/β</button><button id="lyS">site dv/v</button>
<span class="lbl" style="margin-left:8px">δβ/β product</span><select id="cprod"></select></div>
<div class="row"><span class="lbl">model variable</span><select id="var"></select><span class="lbl">date</span><b id="date"></b></div>
<div class="row"><button id="play">&#9654;</button><input type="range" id="t" min="0" max="0" value="0" step="1" style="flex:1;accent-color:var(--accent)"><button id="prev">&lsaquo;</button><button id="next">&rsaquo;</button></div>
<div class="row" style="gap:16px"><div style="flex:1"><div class="lbl">model scale</div><div id="cbarM" class="cbar"></div><div class="cbl"><span id="mmin"></span><span id="mmax"></span></div></div>
<div style="flex:1"><div class="lbl">δβ/β scale <select id="cscale"><option value="0.3">±0.3%</option><option value="0.5" selected>±0.5%</option><option value="1">±1%</option><option value="2">±2%</option></select></div><div id="cbarC" class="cbar"></div><div class="cbl"><span id="cmin"></span><span id="cmax"></span></div></div></div></div>
<div id="plot"></div><div id="info">Click a model node (large square), a δβ/β cell (0.1° square) or a station (triangle).</div></aside></div><script>
const META=__META__, PAY=__PAY__, cache={};
function dec(key,b64,n){if(cache[key])return cache[key];const b=atob(b64);const u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);const raw=pako.inflate(u);cache[key]=new Int16Array(raw.buffer,raw.byteOffset,raw.byteLength/2);return cache[key];}
const N=META.n,NT=META.times.length,C=META.cells,S=META.sites;let cur={v:META.order[0],t:NT-1,sel:null,timer:null,cprod:Object.keys(C.prods)[0],ly:{M:true,C:true,S:false}};
function mval(v,t,i){const q=dec('m:'+v,PAY.model[v])[t*N+i];return q===-32768?null:q/32000*META.vars[v].vmax;}
function cweek(t,i){const p=C.prods[cur.cprod];const q=dec('cw:'+cur.cprod,PAY.cells[cur.cprod].weekly)[t*C.n+i];return q===-32768?null:q/32000*p.vmax_weekly;}
function cdaily(i){const p=C.prods[cur.cprod];const a=dec('cd:'+cur.cprod,PAY.cells[cur.cprod].daily);const out=new Array(C.nd);for(let d=0;d<C.nd;d++){const q=a[d*C.n+i];out[d]=q===-32768?null:q/32000*p.vmax_daily;}return out;}
function sweek(t,j){const q=dec('sw',PAY.sites.weekly)[t*S.n+j];return q===-32768?null:q/32000*S.vmax_weekly;}
function sdaily(j){const a=dec('sd',PAY.sites.daily);const out=new Array(C.nd);for(let d=0;d<C.nd;d++){const q=a[d*S.n+j];out[d]=q===-32768?null:q/32000*S.vmax_daily;}return out;}
function mseries(v,i){const a=dec('m:'+v,PAY.model[v]);const out=new Array(NT);for(let t=0;t<NT;t++){const q=a[t*N+i];out[t]=q===-32768?null:q/32000*META.vars[v].vmax;}return out;}
const vsel=document.getElementById('var');META.order.forEach(v=>{const o=document.createElement('option');o.value=v;o.textContent=META.vars[v].label+' ('+META.vars[v].unit+')';vsel.appendChild(o);});
const psel=document.getElementById('cprod');Object.keys(C.prods).forEach(k=>{const o=document.createElement('option');o.value=k;o.textContent=k;psel.appendChild(o);});
const sl=document.getElementById('t');sl.max=NT-1;sl.value=NT-1;
const map=L.map('map').setView([46.8,-123.5],6);L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',{attribution:'&copy; OpenStreetMap',maxZoom:19}).addTo(map);
function lerp(a,b,x){return a.map((c,i)=>Math.round(c+(b[i]-c)*x));}const RDBU=[[5,48,97],[67,147,195],[247,247,247],[214,96,77],[103,0,31]];const MAG=[[0,0,4],[81,18,124],[183,55,121],[252,137,97],[252,253,191]];
function col(st,x){x=Math.max(0,Math.min(1,x));const k=(st.length-1)*x,i=Math.min(st.length-2,Math.floor(k));return 'rgb('+lerp(st[i],st[i+1],k-i).join(',')+')';}
function colM(x){const V=META.vars[cur.v];if(x==null)return '#cbd5e1';return V.sym?col(RDBU,(x/V.vmax+1)/2):col(MAG,x/V.vmax);}
function colC(x){if(x==null)return '#cbd5e1';const s=+document.getElementById('cscale').value;return col(RDBU,(x/s+1)/2);}
const grpM=L.layerGroup().addTo(map),grpC=L.layerGroup().addTo(map),grpS=L.layerGroup();const rects=[],circs=[],tris=[];const h=META.dlat/2,w=META.dlon/2;let selmk=null;
function mark(la,lo){if(selmk)map.removeLayer(selmk);selmk=L.circleMarker([la,lo],{radius:9,color:'#000',weight:2.5,fill:false}).addTo(map);}
for(let i=0;i<N;i++){const la=META.lat[i],lo=META.lon[i];const r=L.rectangle([[la-h,lo-w],[la+h,lo+w]],{weight:0,fillOpacity:0.75,fillColor:'#cbd5e1'}).addTo(grpM);
 r.on('click',()=>{cur.sel={type:'m',i:i};mark(la,lo);series();});r.on('mouseover',()=>{const x=mval(cur.v,cur.t,i);r.bindTooltip('model node: '+(x==null?'n/a':x.toFixed(3))+' '+META.vars[cur.v].unit+' &middot; interface '+META.iface[i]+' km',{sticky:true}).openTooltip();});rects.push(r);}
const ch=0.05,cw=0.05;for(let i=0;i<C.n;i++){const la=C.lat[i],lo=C.lon[i];const c=L.rectangle([[la-ch,lo-cw],[la+ch,lo+cw]],{weight:0.4,color:'#334155',opacity:0.6,fillOpacity:0.95,fillColor:'#cbd5e1'}).addTo(grpC);
 c.on('click',()=>{cur.sel={type:'c',i:i};mark(la,lo);series();});c.on('mouseover',()=>{const x=cweek(cur.t,i);c.bindTooltip('cell '+C.id[i]+' ('+C.dep[i]+' km): '+(x==null?'n/a':x.toFixed(3)+' %'),{sticky:true}).openTooltip();});circs.push(c);}
for(let j=0;j<S.n;j++){const la=S.lat[j],lo=S.lon[j];const t=L.marker([la,lo],{icon:L.divIcon({className:'',html:'<div style="width:0;height:0;border-left:6px solid transparent;border-right:6px solid transparent;border-bottom:11px solid #94a3b8"></div>',iconSize:[12,11],iconAnchor:[6,6]})}).addTo(grpS);
 t.on('click',()=>{cur.sel={type:'s',i:j};mark(la,lo);series();});t.bindTooltip(S.tags[j]);tris.push(t);}
function paint(){const t=cur.t,V=META.vars[cur.v];document.getElementById('date').textContent=META.times[t];
 if(cur.ly.M)for(let i=0;i<N;i++)rects[i].setStyle({fillColor:colM(mval(cur.v,t,i))});
 if(cur.ly.C)for(let i=0;i<C.n;i++)circs[i].setStyle({fillColor:colC(cweek(t,i))});
 if(cur.ly.S)for(let j=0;j<S.n;j++){const x=sweek(t,j);tris[j].setIcon(L.divIcon({className:'',html:'<div style="width:0;height:0;border-left:6px solid transparent;border-right:6px solid transparent;border-bottom:11px solid '+colC(x)+'"></div>',iconSize:[12,11],iconAnchor:[6,6]}));}
 const st=V.sym?RDBU:MAG;document.getElementById('cbarM').style.background='linear-gradient(90deg,'+st.map((c,k)=>'rgb('+c.join(',')+') '+(100*k/(st.length-1))+'%').join(',')+')';document.getElementById('mmin').textContent=(V.sym?'-'+V.vmax.toPrecision(3):'0');document.getElementById('mmax').textContent='+'+V.vmax.toPrecision(3)+' '+V.unit;
 const s=+document.getElementById('cscale').value;document.getElementById('cbarC').style.background='linear-gradient(90deg,'+RDBU.map((c,k)=>'rgb('+c.join(',')+') '+(100*k/4)+'%').join(',')+')';document.getElementById('cmin').textContent='-'+s+'%';document.getElementById('cmax').textContent='+'+s+'% δβ/β';}
function roll90(days,y){const out=new Array(y.length).fill(null);let acc=[];for(let i=0;i<y.length;i++){acc.push(y[i]);if(acc.length>90)acc.shift();const v=acc.filter(x=>x!=null);out[i]=v.length>=20?v.reduce((a,b)=>a+b,0)/v.length:null;}return out;}
function series(){const s=cur.sel;if(!s)return;const V=META.vars[cur.v];let tr=[],title='',mi;
 if(s.type==='m'){mi=s.i;const ci=C.m2c[mi];const y=cdaily(ci);tr.push({x:C.days,y:y,type:'scattergl',mode:'lines',line:{color:'#cbd5e1',width:0.5},name:'δβ/β daily (cell '+C.id[ci]+')',yaxis:'y2'});tr.push({x:C.days,y:roll90(C.days,y),type:'scattergl',mode:'lines',line:{color:'#1d4ed8',width:1.5},name:'δβ/β 90-d mean (nearest cell)',yaxis:'y2'});title='model node '+META.lat[mi].toFixed(2)+'N '+META.lon[mi].toFixed(2)+'E ('+META.iface[mi]+' km) + nearest δβ/β cell '+C.id[ci];}
 else if(s.type==='c'){mi=C.c2m[s.i];const y=cdaily(s.i);tr.push({x:C.days,y:y,type:'scattergl',mode:'lines',line:{color:'#cbd5e1',width:0.5},name:'δβ/β daily',yaxis:'y2'});tr.push({x:C.days,y:roll90(C.days,y),type:'scattergl',mode:'lines',line:{color:'#1d4ed8',width:1.5},name:'δβ/β 90-d mean',yaxis:'y2'});title='cell '+C.id[s.i]+' ('+C.dep[s.i]+' km, '+cur.cprod+') + nearest model node';}
 else{mi=S.s2m[s.i];const y=sdaily(s.i);tr.push({x:C.days,y:y,type:'scattergl',mode:'lines',line:{color:'#fde68a',width:0.5},name:'site dv/v daily',yaxis:'y2'});tr.push({x:C.days,y:roll90(C.days,y),type:'scattergl',mode:'lines',line:{color:'#b45309',width:1.5},name:'site dv/v 90-d mean',yaxis:'y2'});title='station '+S.tags[s.i]+' site dv/v + nearest model node';}
 tr.unshift({x:META.times,y:mseries(cur.v,mi),type:'scatter',mode:'lines',line:{color:'#dc2626',width:1.3},name:V.label+' ('+V.unit+')'});
 tr.unshift({x:META.times,y:mseries('slip_thrust_inc30',mi),type:'scatter',mode:'lines',fill:'tozeroy',line:{color:'rgba(245,158,11,0.9)',width:0.6},fillcolor:'rgba(245,158,11,0.2)',name:'thrust slip / 30 d (mm)',yaxis:'y3'});
 const shapes=[{type:'line',xref:'x',x0:META.times[cur.t],x1:META.times[cur.t],yref:'paper',y0:0,y1:1,line:{color:'#000',width:1,dash:'dot'}}];
 if(s.type==='s'&&S.eras[S.tags[s.i]])S.eras[S.tags[s.i]].forEach(d=>shapes.push({type:'line',xref:'x',x0:d,x1:d,yref:'paper',y0:0,y1:1,line:{color:'#dc2626',width:1,dash:'dot'}}));
 Plotly.react('plot',tr,{margin:{l:52,r:70,t:34,b:30},title:{text:title,font:{size:11}},showlegend:true,legend:{orientation:'h',y:1.14,font:{size:9}},
  yaxis:{title:V.unit,zeroline:true,titlefont:{color:'#dc2626'}},yaxis2:{title:'dv/v (%)',overlaying:'y',side:'right',zeroline:true,titlefont:{color:'#1d4ed8'}},yaxis3:{overlaying:'y',side:'right',position:1.0,showticklabels:false,showgrid:false,rangemode:'tozero'},
  shapes:shapes,paper_bgcolor:'#fff',plot_bgcolor:'#fff'},{responsive:true});
 document.getElementById('info').innerHTML='<b>'+title+'</b>. Red = selected model variable; orange = modeled thrust slip in the trailing 30-d window; blue/brown = coda dv/v (right axis). Black dotted = slider date; red dotted = instrument era boundary.';}
vsel.onchange=()=>{cur.v=vsel.value;paint();series();};psel.onchange=()=>{cur.cprod=psel.value;paint();series();};document.getElementById('cscale').onchange=paint;
sl.oninput=()=>{cur.t=+sl.value;paint();if(cur.sel)series();};document.getElementById('prev').onclick=()=>{sl.value=Math.max(0,cur.t-1);sl.oninput();};document.getElementById('next').onclick=()=>{sl.value=Math.min(NT-1,cur.t+1);sl.oninput();};
document.getElementById('play').onclick=function(){if(cur.timer){clearInterval(cur.timer);cur.timer=null;this.innerHTML='&#9654;';return;}this.innerHTML='&#10074;&#10074;';cur.timer=setInterval(()=>{sl.value=(cur.t+1)%NT;sl.oninput();},150);};
function toggle(k,btn,grp){cur.ly[k]=!cur.ly[k];btn.classList.toggle('on',cur.ly[k]);if(cur.ly[k]){grp.addTo(map);paint();}else map.removeLayer(grp);}
document.getElementById('lyM').onclick=function(){toggle('M',this,grpM);};document.getElementById('lyC').onclick=function(){toggle('C',this,grpC);};document.getElementById('lyS').onclick=function(){toggle('S',this,grpS);};
(function(){const dr=document.getElementById('drag'),sd=document.getElementById('side');let on=false;dr.addEventListener('mousedown',e=>{on=true;dr.classList.add('on');document.body.style.cursor='col-resize';e.preventDefault();});
 window.addEventListener('mousemove',e=>{if(!on)return;const w=Math.min(window.innerWidth-220,Math.max(380,window.innerWidth-e.clientX));sd.style.width=w+'px';});
 window.addEventListener('mouseup',()=>{if(!on)return;on=false;dr.classList.remove('on');document.body.style.cursor='';map.invalidateSize();Plotly.Plots.resize('plot');});})();
window.addEventListener('resize',()=>{map.invalidateSize();Plotly.Plots.resize('plot');});paint();
</script></body></html>'''
out = HTML.replace("__META__", json.dumps(meta, separators=(",", ":"))).replace("__PAY__", json.dumps(payload, separators=(",", ":")))
open("figures/sse_dvv_combined.html", "w").write(out); print(f"wrote figures/sse_dvv_combined.html ({len(out) / 1e6:.1f} MB)")
