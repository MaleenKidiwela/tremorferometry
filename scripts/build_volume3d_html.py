#!/usr/bin/env python
"""Standalone interactive HTML for the coarse 3D volume dv/v: pick a DEPTH ZONE (shallow/mid/deep),
slide TIME (months). Leaflet basemap + colored 3D cells, recoloured on depth/time change. Reads
fault_tomography/inversion/fault_volume3d.npz -> fault_tomography/cascadia_volume3d.html (self-contained).
RESOLUTION IS COARSE (~150-200 km horizontal, 3 depth zones) -- stated in the page.
"""
import json
import numpy as np, pandas as pd

z = np.load('fault_tomography/inversion/fault_volume3d.npz', allow_pickle=True)
lon = z['lon'].astype(float); lat = z['lat'].astype(float); dep = z['depth'].astype(float)
MV = z['MV']; ok = z['ok'].astype(bool); months = z['months'].astype(str); DEPTHS = z['DEPTHS'].astype(float)
DLAB = ['shallow ~%.0f km' % DEPTHS[0], 'mid ~%.0f km' % DEPTHS[1], 'deep ~%.0f km' % DEPTHS[2]]
mok = months[ok]
# 3-month smoothed version for cleaner scrubbing
df = pd.DataFrame(MV.T, index=pd.PeriodIndex(months, freq='M'))
MVs = df.rolling(3, center=True, min_periods=2).mean().values.T

layers = []
for di, d0 in enumerate(DEPTHS):
    m = np.abs(dep - d0) < 1
    cells = [{'lat': round(float(a), 3), 'lon': round(float(o), 3)} for a, o in zip(lat[m], lon[m])]
    def pack(A):
        return [[None if not np.isfinite(v) else round(float(v), 3) for v in A[m, j]] for j in np.where(ok)[0]]
    layers.append({'label': DLAB[di], 'depth': float(d0), 'cells': cells, 'vals': pack(MV), 'svals': pack(MVs)})
DATA = {'months': list(mok), 'layers': layers, 'grid': 0.4}
js = json.dumps(DATA, separators=(',', ':'))
print('depths:', DEPTHS.tolist(), '| months:', len(mok), '| cells/depth:', [len(l['cells']) for l in layers],
      '| data %.2f MB' % (len(js)/1e6))

HTML = '''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia 3D dv/v volume</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<style>
:root{--bg:#f1f5f9;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb;--accent-soft:#eff6ff;--shadow:0 1px 3px rgba(15,23,42,.08)}
*{box-sizing:border-box}html,body{margin:0;height:100%;font-family:-apple-system,system-ui,Segoe UI,Roboto,Arial,sans-serif;color:var(--ink);background:var(--bg)}
#wrap{display:flex;height:100vh}#map{flex:1;height:100%}
#side{width:340px;flex:none;border-left:1px solid var(--line);background:var(--card);padding:16px 18px;overflow-y:auto}
h1{font-size:16px;margin:0 0 2px}.sub{font-size:12px;color:var(--muted);line-height:1.45;margin-bottom:14px}
.lbl{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:.04em;color:var(--muted);margin:14px 0 7px}
.seg{display:flex;background:var(--bg);border:1px solid var(--line);border-radius:10px;padding:3px;gap:2px}
.seg button{flex:1;font-size:12px;font-weight:600;color:var(--muted);background:none;border:none;padding:8px 6px;border-radius:8px;cursor:pointer}
.seg button.on{background:var(--card);color:var(--accent);box-shadow:var(--shadow)}
#tlabel{font-size:15px;font-weight:700;margin:4px 0 6px}
input[type=range]{width:100%}
.row{display:flex;align-items:center;gap:8px;font-size:12.5px;color:#475569;margin-top:10px}
.cbar{height:14px;border-radius:4px;margin-top:8px;background:linear-gradient(90deg,#1e5fbf,#9ecae1,#fff,#fcae91,#cb1b1b)}
.cbw{display:flex;justify-content:space-between;font-size:11px;color:var(--muted)}
.note{font-size:11px;color:var(--muted);line-height:1.5;margin-top:16px;border-top:1px solid var(--line);padding-top:12px}
.leaflet-container{font-family:inherit}
</style></head>
<body><div id="wrap"><div id="map"></div>
<div id="side">
<h1>Cascadia 3D dv/v volume</h1>
<div class="sub">Coarse 3D inversion (multi-window 1&ndash;3 / 2&ndash;4 / 3&ndash;5 s). Pick a depth zone, slide time. Red = velocity drop, blue = increase.</div>
<div class="lbl">depth zone</div><div class="seg" id="depseg"></div>
<div class="lbl">month</div><div id="tlabel"></div>
<input type="range" id="tslider" min="0" value="0">
<label class="row"><input type="checkbox" id="smooth" checked> 3-month smoothed</label>
<div class="lbl">&delta;&beta;/&beta; (%)</div><div class="cbar"></div>
<div class="cbw"><span>&minus;0.2 (slower)</span><span>0</span><span>+0.2 (faster)</span></div>
<div class="note"><b>Resolution is coarse</b> &mdash; ~150&ndash;200 km horizontal, 3 depth zones (validated by a 3D checkerboard: coarse corr 0.86, shallow vs deep separable). Fine detail is NOT resolved; read large patches, not single cells. Stations = black triangles.</div>
</div></div><script>
const DATA=__DATA__;
const map=L.map('map').setView([45.5,-123.3],6);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',{attribution:'&copy; OSM, &copy; CARTO',subdomains:'abcd',maxZoom:12}).addTo(map);
const STA=[[49.22,-124.81],[49.23,-123.99],[48.83,-125.13],[48.65,-123.45],[48.2,-124.43],[47.81,-122.91],[47.65,-123.05],[47.56,-122.83],[47.51,-123.81],[46.99,-122.22],[46.98,-123.02],[46.38,-123.84],[46.3,-122.26],[46.14,-122.17],[46.11,-123.08],[45.95,-123.93],[45.31,-123.82],[45.17,-122.44],[44.59,-123.3],[44.49,-122.96],[43.97,-122.77],[43.67,-123.39],[43.29,-123.12],[42.51,-123.38],[41.83,-122.42],[41.47,-122.48],[40.48,-123.57]];
STA.forEach(s=>L.marker(s,{icon:L.divIcon({className:'',html:'<div style=\"color:#111;font-size:13px\">&#9650;</div>',iconSize:[12,12]})}).addTo(map));
function cmap(v){if(v==null||isNaN(v))return 'rgba(180,180,180,0.25)';var t=Math.max(-1,Math.min(1,v/0.2));var r,g,b;
 if(t>=0){r=255;g=Math.round(255-180*t);b=Math.round(255-200*t);}else{r=Math.round(255+200*t);g=Math.round(255+180*t);b=255;}return 'rgb('+r+','+g+','+b+')';}
let curD=2, curT=0;
const g=DATA.grid/2; let rects=[];
function smooth(){return document.getElementById('smooth').checked;}
function drawDepth(di){rects.forEach(r=>map.removeLayer(r));rects=[];const L1=DATA.layers[di];
 L1.cells.forEach(c=>{const r=L.rectangle([[c.lat-g,c.lon-g],[c.lat+g,c.lon+g]],{stroke:false,fillOpacity:0.78}).addTo(map);rects.push(r);});}
function recolor(){const L1=DATA.layers[curD];const vals=(smooth()?L1.svals:L1.vals)[curT];
 rects.forEach((r,i)=>r.setStyle({fillColor:cmap(vals[i])}));
 document.getElementById('tlabel').textContent=DATA.months[curT];}
// depth buttons
const ds=document.getElementById('depseg');
DATA.layers.forEach((L1,di)=>{const b=document.createElement('button');b.textContent=L1.label;b.onclick=()=>{curD=di;[...ds.children].forEach((x,k)=>x.classList.toggle('on',k===di));drawDepth(di);recolor();};ds.appendChild(b);});
ds.children[curD].classList.add('on');
// time slider
const ts=document.getElementById('tslider');ts.max=DATA.months.length-1;ts.value=DATA.months.length-1;curT=DATA.months.length-1;
ts.addEventListener('input',e=>{curT=+e.target.value;recolor();});
document.getElementById('smooth').addEventListener('change',recolor);
drawDepth(curD);recolor();
</script></body></html>'''
with open('fault_tomography/cascadia_volume3d.html', 'w') as f:
    f.write(HTML.replace('__DATA__', js))
print('wrote fault_tomography/cascadia_volume3d.html')
