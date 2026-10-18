#!/usr/bin/env python
"""Interactive map + time slider of the slow-slip model products in OUR measurement frame (30-d trailing mean on
the vbICA amplitudes, gauge 2019-24): reads res_catalog_sse/megathrust_strain_trailing30.nc (repo grid 0.1 deg,
weekly) and adds thrust slip (cumulative + 30-d increment) at the grid nodes from the nearest mesh patch.
Variable selector (all 22 fields + 2 slip), play/slider, per-variable fixed robust colour scale, hover values,
click a node -> time series (Plotly) with the slider date marked. -> figures/sse_model_map.html
Data: per variable int16-quantised, gzip+base64, decoded on demand in the browser (pako)."""
import sys, json, gzip, base64, numpy as np, pandas as pd, xarray as xr
REPO = "/home/jovyan/cascadia-sse-strain"; sys.path.insert(0, f"{REPO}/src")
SD = "fault_tomography/inversion/res_catalog_sse"
ds = xr.open_dataset(f"{SD}/megathrust_strain_trailing30.nc")
lat, lon = ds.lat.values, ds.lon.values; times = pd.to_datetime(ds.time.values)
LA, LO = np.meshgrid(lat, lon, indexing="ij")
valid = np.isfinite(ds["dilatation"].isel(depth=0).values).any(0)
ii, jj = np.where(valid); npts = len(ii)
print(f"grid nodes {npts} | epochs {len(times)} ({times[0].date()}..{times[-1].date()})")
UNITS = {"strain": ("1e-8", 1e-8), "stress": ("kPa", 1e3), "dilatation": ("1e-8", 1e-8), "eps_eq": ("1e-8", 1e-8),
         "max_shear_strain": ("1e-8", 1e-8), "mean_stress": ("kPa", 1e3), "pressure": ("kPa", 1e3), "von_mises_stress": ("kPa", 1e3),
         "max_shear_stress": ("kPa", 1e3), "coulomb": ("kPa", 1e3), "shear_stress": ("kPa", 1e3), "normal_stress": ("kPa", 1e3)}
def unit_of(v):
    for k, u in UNITS.items():
        if v == k or v.startswith(k + "_"): return u
    return ("", 1.0)
VARS = {}
order = ["slip_thrust_inc30", "slip_thrust_cum", "coulomb", "shear_stress", "normal_stress", "dilatation", "eps_eq", "max_shear_strain",
         "pressure", "mean_stress", "von_mises_stress", "max_shear_stress"] + [f"strain_{c}" for c in ("exx", "eyy", "ezz", "exy", "exz", "eyz")] + \
        [f"stress_s{c}" for c in ("xx", "yy", "zz", "xy", "xz", "yz")]
# slip at grid nodes from the mesh (nearest patch), trailing frame
from sse_strain import load_fault_ascii, load_solution
from sse_strain.geodesy import llh2localxy
from scipy.spatial import cKDTree
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True)
A_trail = zc["A_trail"].astype(float); A_raw = zc["A_raw"].astype(float); adays = pd.to_datetime(zc["days"])
cll = fault.centroid_llh; px, py = llh2localxy(cll[:, 1], cll[:, 0], fault.origin)
gx, gy = llh2localxy(LA[valid], LO[valid], fault.origin)
_, ip = cKDTree(np.column_stack([px, py])).query(np.column_stack([gx, gy]))
thr = -(Ld @ np.nan_to_num(A_trail).T)[ip]                                       # (npts, days) thrust-positive, m
thr_raw = -(Ld @ np.nan_to_num(A_raw).T)[ip]                                     # un-averaged cumulative
tidx = adays.get_indexer(times, method="nearest")
cum = thr[:, tidx]; inc = np.full_like(cum, np.nan)
t30 = adays.get_indexer(times - pd.Timedelta(days=30), method="nearest"); ok = t30 >= 0
inc[:, ok] = thr_raw[:, tidx[ok]] - thr_raw[:, t30[ok]]                            # slip within (t-30, t]
extra = {"slip_thrust_cum": (cum.T * 1000, "mm"), "slip_thrust_inc30": (inc.T * 1000, "mm / 30 d")}
LABEL = {"slip_thrust_inc30": "thrust slip within the trailing 30-d window", "slip_thrust_cum": "cumulative thrust slip (rel. 2019-24)"}
payload = {}
for v in order:
    if v in extra:
        arr, unit = extra[v]; arr = arr.astype(float)
    else:
        unit, sc = unit_of(v); arr = ds[v].isel(depth=0).values[:, ii, jj].astype(float) / sc
    vmax = float(np.nanpercentile(np.abs(arr), 99.5)) or 1.0
    q = np.round(np.clip(arr / vmax, -1, 1) * 32000).astype(np.int16); q[np.isnan(arr)] = -32768
    payload[v] = base64.b64encode(gzip.compress(q.tobytes(), 6)).decode()
    VARS[v] = dict(unit=unit, vmax=vmax, label=LABEL.get(v, v.replace("_", " ")), sym=bool(v != "slip_thrust_inc30"))
    print(f"  {v}: scale ±{vmax:.3g} {unit}, {len(payload[v]) / 1e6:.1f} MB")
meta = dict(lat=[round(float(x), 3) for x in LA[valid]], lon=[round(float(x), 3) for x in LO[valid]], times=[t.strftime("%Y-%m-%d") for t in times],
            iface=[round(float(x), 1) for x in ds.interface_depth.values[valid]], vars=VARS, order=order, n=npts, dlat=float(np.diff(lat).mean()), dlon=float(np.diff(lon).mean()))
HTML = r'''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia slow-slip model — dv/v frame</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/><script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script><script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
<style>:root{--bg:#f1f5f9;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb}
*{box-sizing:border-box}html,body{margin:0;height:100%;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,sans-serif;color:var(--ink);background:var(--bg)}
#wrap{display:flex;height:100vh}#map{flex:1;height:100%;min-width:200px}#drag{width:7px;flex:none;cursor:col-resize;background:var(--line)}#drag:hover,#drag.on{background:var(--accent)}#side{width:560px;flex:none;min-width:360px;display:flex;flex-direction:column;border-left:1px solid var(--line);background:var(--bg);overflow:auto}
#hdr{padding:12px 16px;background:var(--card);border-bottom:1px solid var(--line)}.title{font-size:15px;font-weight:700}.sub{font-size:11.5px;color:var(--muted);margin-top:3px;line-height:1.45}
.panel{margin:10px 12px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.lbl{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
select,button{font:inherit;font-size:12.5px;padding:4px 8px;border:1px solid var(--line);border-radius:6px;background:#fff}
#plot{flex:1;min-height:320px;margin:10px 12px;background:#fff;border:1px solid var(--line);border-radius:10px}
#info{font-size:12px;color:var(--muted);padding:6px 16px;line-height:1.5}#cbar{height:14px;border-radius:4px;margin-top:6px;border:1px solid var(--line)}
.cbl{display:flex;justify-content:space-between;font-size:10.5px;color:var(--muted)}
.leaflet-container{background:#e5e9f0}</style></head><body><div id="wrap"><div id="map"></div><div id="drag" title="drag to widen the panel"></div><aside id="side">
<div id="hdr"><div class="title">Cascadia slow-slip model (Gualandi 2025) in the dv/v measurement frame</div>
<div class="sub">30-day trailing mean (t&minus;30, t] on the vbICA time functions &middot; reference 2019&ndash;2024 &middot; evaluated 10 km above the Slab2 interface, 21-km low-pass &middot; weekly frames 2007&ndash;2026. Click a node for its time series.</div></div>
<div class="panel"><span class="lbl">variable</span> <select id="var"></select> <span class="lbl" style="margin-left:10px">date</span> <b id="date"></b>
<div style="display:flex;align-items:center;gap:8px;margin-top:8px"><button id="play">&#9654;</button><input type="range" id="t" min="0" max="0" value="0" step="1" style="flex:1;accent-color:var(--accent)"><button id="prev">&lsaquo;</button><button id="next">&rsaquo;</button></div>
<div id="cbar"></div><div class="cbl"><span id="cmin"></span><span id="cmid"></span><span id="cmax"></span></div></div>
<div id="plot"></div><div id="info">Click a grid node.</div></aside></div><script>
const META=__META__, PAY=__PAY__, cache={};
function decode(v){if(cache[v])return cache[v];const b=atob(PAY[v]);const u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);const raw=pako.inflate(u);cache[v]=new Int16Array(raw.buffer,raw.byteOffset,raw.byteLength/2);return cache[v];}
const N=META.n,NT=META.times.length;let cur={v:META.order[0],t:NT-1,pt:null,timer:null};
function val(v,t,i){const q=decode(v)[t*N+i];return q===-32768?null:q/32000*META.vars[v].vmax;}
const sel=document.getElementById('var');META.order.forEach(v=>{const o=document.createElement('option');o.value=v;o.textContent=META.vars[v].label+' ('+META.vars[v].unit+')';sel.appendChild(o);});
const sl=document.getElementById('t');sl.max=NT-1;sl.value=NT-1;
const map=L.map('map').setView([46.5,-123.5],6);L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',{attribution:'&copy; OpenStreetMap',maxZoom:19}).addTo(map);
// colour scales: diverging (RdBu) for signed fields, sequential (magma-like) for the slip increment
function lerp(a,b,x){return a.map((c,i)=>Math.round(c+(b[i]-c)*x));}
const RDBU=[[5,48,97],[67,147,195],[247,247,247],[214,96,77],[103,0,31]];const MAG=[[0,0,4],[81,18,124],[183,55,121],[252,137,97],[252,253,191]];
function col(stops,x){x=Math.max(0,Math.min(1,x));const k=(stops.length-1)*x,i=Math.min(stops.length-2,Math.floor(k));const c=lerp(stops[i],stops[i+1],k-i);return 'rgb('+c.join(',')+')';}
function colorFor(v,x){const V=META.vars[v];if(x==null)return '#cbd5e1';return V.sym?col(RDBU,(x/V.vmax+1)/2):col(MAG,x/V.vmax);}
const rects=[];let selmk=null;const h=META.dlat/2,w=META.dlon/2;
for(let i=0;i<N;i++){const la=META.lat[i],lo=META.lon[i];const r=L.rectangle([[la-h,lo-w],[la+h,lo+w]],{weight:0,fillOpacity:0.85,fillColor:'#cbd5e1'}).addTo(map);
 r.on('click',()=>{if(cur.pt!=null)rects[cur.pt].setStyle({weight:0});cur.pt=i;r.setStyle({weight:3,color:'#000',opacity:1});r.bringToFront();if(selmk)map.removeLayer(selmk);selmk=L.circleMarker([la,lo],{radius:7,color:'#000',weight:2,fillColor:'#fff',fillOpacity:0.9}).addTo(map);series();});r.on('mouseover',()=>{const x=val(cur.v,cur.t,i);r.bindTooltip((x==null?'n/a':x.toFixed(3))+' '+META.vars[cur.v].unit+' &middot; '+la.toFixed(2)+'N '+lo.toFixed(2)+'E &middot; interface '+META.iface[i]+' km',{sticky:true}).openTooltip();});rects.push(r);}
function paint(){const v=cur.v,t=cur.t,V=META.vars[v];for(let i=0;i<N;i++)rects[i].setStyle({fillColor:colorFor(v,val(v,t,i))});document.getElementById('date').textContent=META.times[t];
 const cb=document.getElementById('cbar');const st=V.sym?RDBU:MAG;cb.style.background='linear-gradient(90deg,'+st.map((c,k)=>'rgb('+c.join(',')+') '+(100*k/(st.length-1))+'%').join(',')+')';
 document.getElementById('cmin').textContent=(V.sym?'-':'0')+(V.sym?V.vmax.toPrecision(3):'');document.getElementById('cmid').textContent=V.sym?'0':(V.vmax/2).toPrecision(3);document.getElementById('cmax').textContent='+'+V.vmax.toPrecision(3)+' '+V.unit;}
function series(){const i=cur.pt;if(i==null)return;const v=cur.v,V=META.vars[v];const y=[];for(let t=0;t<NT;t++)y.push(val(v,t,i));const y2=[];for(let t=0;t<NT;t++)y2.push(val('slip_thrust_inc30',t,i));
 Plotly.react('plot',[{x:META.times,y:y2,type:'scatter',mode:'lines',fill:'tozeroy',line:{color:'rgba(245,158,11,0.9)',width:0.8},fillcolor:'rgba(245,158,11,0.25)',name:'thrust slip / 30 d (mm)',yaxis:'y2'},
  {x:META.times,y:y,type:'scatter',mode:'lines',line:{color:'#1d4ed8',width:1.4},name:V.label}],
  {margin:{l:50,r:50,t:30,b:30},title:{text:V.label+' @ '+META.lat[i].toFixed(2)+'N '+META.lon[i].toFixed(2)+'E (interface '+META.iface[i]+' km)',font:{size:12}},showlegend:true,legend:{orientation:'h',y:1.12,font:{size:10}},
   yaxis:{title:V.unit,zeroline:true},yaxis2:{title:'mm/30d',overlaying:'y',side:'right',showgrid:false,rangemode:'tozero',titlefont:{color:'#b45309',size:10},tickfont:{color:'#b45309',size:9}},
   shapes:[{type:'line',xref:'x',x0:META.times[cur.t],x1:META.times[cur.t],yref:'paper',y0:0,y1:1,line:{color:'#dc2626',width:1.2,dash:'dot'}}],paper_bgcolor:'#fff',plot_bgcolor:'#fff'},{responsive:true});
 document.getElementById('info').innerHTML='<b>'+V.label+'</b> at node '+i+' &mdash; red dotted = slider date. Orange = modeled thrust slip in the trailing 30-day window (the same window the dv/v uses).';}
sel.onchange=()=>{cur.v=sel.value;paint();series();};sl.oninput=()=>{cur.t=+sl.value;paint();if(cur.pt!=null)series();};
document.getElementById('prev').onclick=()=>{sl.value=Math.max(0,cur.t-1);sl.oninput();};document.getElementById('next').onclick=()=>{sl.value=Math.min(NT-1,cur.t+1);sl.oninput();};
document.getElementById('play').onclick=function(){if(cur.timer){clearInterval(cur.timer);cur.timer=null;this.innerHTML='&#9654;';return;}this.innerHTML='&#10074;&#10074;';cur.timer=setInterval(()=>{sl.value=(cur.t+1)%NT;sl.oninput();},120);};
paint();
(function(){const dr=document.getElementById('drag'),sd=document.getElementById('side');let on=false;
 dr.addEventListener('mousedown',e=>{on=true;dr.classList.add('on');document.body.style.cursor='col-resize';e.preventDefault();});
 window.addEventListener('mousemove',e=>{if(!on)return;const w=Math.min(window.innerWidth-220,Math.max(360,window.innerWidth-e.clientX));sd.style.width=w+'px';});
 window.addEventListener('mouseup',()=>{if(!on)return;on=false;dr.classList.remove('on');document.body.style.cursor='';map.invalidateSize();Plotly.Plots.resize('plot');});})();
window.addEventListener('resize',()=>{map.invalidateSize();Plotly.Plots.resize('plot');});
</script></body></html>'''
out = HTML.replace("__META__", json.dumps(meta, separators=(",", ":"))).replace("__PAY__", json.dumps(payload, separators=(",", ":")))
open("figures/sse_model_map.html", "w").write(out)
print(f"wrote figures/sse_model_map.html ({len(out) / 1e6:.1f} MB)")
