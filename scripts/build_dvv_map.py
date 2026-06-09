#!/usr/bin/env python
"""Standalone interactive Cascadia dv/v map.
TEMPORAL AXIS = daily 30-day causal rolling (NO monthly binning, NO median over time).
The only median is the cross-family median PER DAY (a spatial median across patches), shown as a
reference line.  Click a station -> all families' daily lines + cross-family median; click a family ->
that family's daily line + the median.  Daily data is gzip-compressed per station and inflated on click
(pako) so the file stays openable.  Window selector (2-4/1-4/1-3) + deseason switch + local tremor overlay.
Reads data/daily_dvv_{S}_{2to4,1to4,1to3}_{roll,des}.csv -> fault_tomography/cascadia_dvv_map.html
"""
import json, gzip, base64
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

STA = {
 'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732),
}
CCMIN = 0.7; MINDAYS = 60
WINS = [('24','2to4_cal'),('24d','2to4_cal_des'),('13','1to3_cal'),('13d','1to3_cal_des'),('35','3to5_cal'),('35d','3to5_cal_des')]

def famll(pid):
    a = pid.split('__')[0].split('_'); return float(a[0]), float(a[1])

def daily_block(csv):
    """Daily 30-day-rolling per family (NO temporal median) + per-day cross-family median."""
    try:
        d = pd.read_csv(csv)
    except FileNotFoundError:
        return None
    d = d[d.cc_max > CCMIN].copy()
    if not len(d):
        return None
    d['ds'] = pd.to_datetime(d['date']).dt.strftime('%Y-%m-%d')
    piv = d.pivot_table(index='ds', columns='patch', values='dvv', aggfunc='first') * 100.0  # one value/day, no agg
    days = sorted(piv.index); piv = piv.reindex(days)
    byfam = {}
    for fid in piv.columns:
        col = piv[fid]
        if col.notna().sum() < MINDAYS:
            continue
        byfam[fid] = [None if pd.isna(v) else round(float(v), 3) for v in col.values]
    if not byfam:
        return None
    med = piv[list(byfam)].median(axis=1)            # cross-family median PER DAY (spatial median, not over time)
    med = [None if pd.isna(v) else round(float(v), 3) for v in med.values]
    return {'dates': days, 'byfam': byfam, 'med': med, 'cc': round(float(d.cc_max.mean()), 2)}

# --- tremor context (monthly counts near families): total=footprint, byfam=local within R km ---
TREMOR_R_KM = 20.0; LAT0, LON0 = 45.5, -123.0
def _km(lat, lon):
    return ((np.asarray(lon, float) - LON0) * 111.0 * np.cos(np.radians(LAT0)), (np.asarray(lat, float) - LAT0) * 111.0)
_tr = pd.read_csv('catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
_tr['t'] = pd.to_datetime(_tr['time'], errors='coerce'); _tr = _tr.dropna(subset=['t', 'lat', 'lon'])
_tr_ym = _tr['t'].dt.to_period('M').astype(str).values
_trx, _try = _km(_tr['lat'].values, _tr['lon'].values); _tr_xy = np.column_stack([_trx, _try]); _tr_tree = cKDTree(_tr_xy)
print(f'tremor catalog {len(_tr):,} detections (R={TREMOR_R_KM:.0f} km)')
def tremor_for(fam):
    if not fam: return None
    fx, fy = _km([f['lat'] for f in fam], [f['lon'] for f in fam]); fxy = np.column_stack([fx, fy])
    dmin, _ = cKDTree(fxy).query(_tr_xy, k=1); near = dmin < TREMOR_R_KM
    if int(near.sum()) == 0: return None
    months = sorted(set(_tr_ym[near])); mi = {m: i for i, m in enumerate(months)}
    total = np.zeros(len(months), int)
    for ym in _tr_ym[near]: total[mi[ym]] += 1
    byfam = {}
    for k, f in enumerate(fam):
        cnt = np.zeros(len(months), int); hit = False
        for j in _tr_tree.query_ball_point(fxy[k], TREMOR_R_KM):
            m = _tr_ym[j]
            if m in mi: cnt[mi[m]] += 1; hit = True
        if hit: byfam[f['id']] = [int(v) for v in cnt]
    return {'dates': months, 'total': [int(v) for v in total], 'byfam': byfam}

DATA, DAILY = {}, {}
for s, (la, lo) in STA.items():
    blocks, ccs = {}, {}
    for k, suf in WINS:
        b = daily_block(f'data/daily_dvv_{s}_{suf}.csv')
        if b:
            ccs[k] = b.pop('cc'); blocks['s' + k] = b
    if 's24' not in blocks and 's13' not in blocks:
        print('  skip', s); continue
    fam_ids = sorted(set().union(*[set(blocks[w]['byfam']) for w in blocks]))
    fam = [{'id': fid, 'lat': round(famll(fid)[0], 3), 'lon': round(famll(fid)[1], 3)} for fid in fam_ids]
    DATA[s] = {'lat': la, 'lon': lo, 'fam': fam, 'tremor': tremor_for(fam)}
    for k, _ in WINS:
        DATA[s]['cc' + k] = ccs.get(k)
    raw = json.dumps(blocks, separators=(',', ':')).encode()
    DAILY[s] = base64.b64encode(gzip.compress(raw, 6)).decode()
    print(f'  {s}: {len(fam)} fam, daily {sum(len(b["dates"]) for b in blocks.values()):,} days -> {len(DAILY[s])/1e6:.2f} MB gz')

djs = json.dumps(DATA, separators=(',', ':'))
djson = json.dumps(DAILY, separators=(',', ':'))
print(f'DATA {len(djs)/1e6:.1f} MB | DAILY (compressed) {len(djson)/1e6:.1f} MB')

HTML = '''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia LFE-coda dv/v map</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
<style>
:root{--bg:#f1f5f9;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb;--accent-soft:#eff6ff;--nb:#ea580c;--nb-soft:#fff7ed;--shadow:0 1px 3px rgba(15,23,42,.08),0 1px 2px rgba(15,23,42,.04)}
*{box-sizing:border-box}
html,body{margin:0;height:100%;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,Helvetica,Arial,sans-serif;color:var(--ink);background:var(--bg)}
#wrap{display:flex;height:100vh}
#map{flex:1;height:100%;min-width:240px}
#drag{width:6px;flex:none;cursor:col-resize;background:var(--line);transition:background .12s}
#drag:hover,#drag.on{background:var(--accent)}
#side{width:520px;flex:none;min-width:380px;display:flex;flex-direction:column;background:var(--bg)}
#hdr{padding:14px 18px 12px;background:var(--card);border-bottom:1px solid var(--line)}
#hdr .title{font-size:16px;font-weight:700;letter-spacing:-.01em}
#hdr .sub{font-size:12px;color:var(--muted);margin-top:3px;line-height:1.45}
.panel{margin:12px 14px 0;background:var(--card);border:1px solid var(--line);border-radius:12px;box-shadow:var(--shadow)}
.lbl{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:.04em;color:var(--muted);margin-bottom:7px}
#ctrl{padding:12px 14px}
.seg{display:inline-flex;background:var(--bg);border:1px solid var(--line);border-radius:10px;padding:3px;gap:2px}
.seg input{display:none}
.seg label{font-size:12.5px;font-weight:600;color:var(--muted);padding:6px 14px;border-radius:8px;cursor:pointer;transition:.12s;user-select:none}
.seg label:hover{color:var(--ink)}
.seg input:checked+label{background:var(--card);color:var(--accent);box-shadow:var(--shadow)}
#dswrap{display:flex;align-items:center;gap:9px;margin-top:12px;cursor:pointer;font-size:12.5px;color:#475569;user-select:none}
#dswrap input{display:none}
.switch{width:34px;height:19px;border-radius:999px;background:#cbd5e1;position:relative;transition:.15s;flex:none}
.switch::after{content:'';position:absolute;top:2px;left:2px;width:15px;height:15px;border-radius:50%;background:#fff;box-shadow:var(--shadow);transition:.15s}
#dswrap input:checked+.switch{background:var(--accent)}
#dswrap input:checked+.switch::after{transform:translateX(15px)}
.dslbl{font-weight:500}
#plot{height:42vh;min-height:300px;padding:8px 8px 2px;overflow:hidden}
#info{font-size:12.5px;color:#475569;line-height:1.55;padding:10px 16px 12px}
#info b.k{color:var(--ink)}
#stawrap{flex:1;min-height:0;display:flex;flex-direction:column;padding:12px 14px 14px}
.nbkey{color:var(--nb);font-weight:600}
#stafilter{width:100%;padding:7px 11px;border:1px solid var(--line);border-radius:8px;font-size:12.5px;margin-bottom:8px;outline:none;font-family:inherit}
#stafilter:focus{border-color:var(--accent);box-shadow:0 0 0 3px var(--accent-soft)}
#stalist{flex:1;overflow-y:auto}
.chip{display:inline-block;padding:3px 9px;margin:2px 3px;font-size:12px;font-weight:600;border:1px solid var(--line);border-radius:999px;cursor:pointer;background:var(--card);color:#334155;transition:.12s}
.chip:hover{border-color:var(--accent);color:var(--accent);background:var(--accent-soft)}
.chip.nb{border-color:#fed7aa;color:var(--nb);background:var(--nb-soft)}
.chip.active{background:var(--accent);color:#fff;border-color:var(--accent)}
.chip.nb.active{background:var(--nb);border-color:var(--nb);color:#fff}
.nblabel{background:rgba(255,247,237,.95);border:1px solid var(--nb);font-weight:700;font-size:10px;padding:0 4px;color:var(--nb);border-radius:4px}
.leaflet-container{font-family:inherit}
</style></head>
<body><div id="wrap"><div id="map"></div><div id="drag" title="drag to resize"></div>
<aside id="side">
<div id="hdr"><div class="title">Cascadia LFE-coda dv/v</div><div class="sub">daily, true 30-calendar-day trailing stack &middot; click a <b>station</b> (all patches + median) or a <b>patch</b> (that patch + median)</div></div>
<div id="ctrl" class="panel"><div class="lbl">coda window</div>
<div class="seg" id="winseg">
<input type="radio" name="win" id="w13" value="13"><label for="w13">1&ndash;3 s</label>
<input type="radio" name="win" id="w24" value="24" checked><label for="w24">2&ndash;4 s</label>
<input type="radio" name="win" id="w35" value="35"><label for="w35">3&ndash;5 s</label>
</div>
<label id="dswrap"><input type="checkbox" id="deseason"><span class="switch"></span><span class="dslbl">deseason &mdash; remove annual cycle</span></label>
</div>
<div id="plot" class="panel"></div>
<div id="info">Each line is the <b class=k>daily</b> dv/v from a 30-day causal rolling stack (each day = the 30 days ending that date). Orange = tremor on the right axis.</div>
<div id="stawrap"><div class="lbl">jump to station &nbsp;<span class="nbkey">&#9679; non-borehole</span></div>
<input id="stafilter" placeholder="filter stations&hellip;" autocomplete="off">
<div id="stalist"></div></div>
</aside></div><script>
const DATA=__DATA__;
const DAILY=__DAILY__;
const _dc={};
function dailyOf(s){if(!(s in _dc)){try{const b=atob(DAILY[s]);const u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);_dc[s]=JSON.parse(pako.inflate(u,{to:'string'}));}catch(e){_dc[s]=null;}}return _dc[s];}
const PAL=['#e6194B','#3cb44b','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#fabed4','#469990','#dcbeff','#9A6324','#800000','#aaffc3','#808000','#000075','#a9a9a9','#e6194B','#3cb44b','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#469990','#9A6324','#800000','#000075'];
const names=Object.keys(DATA).sort((a,b)=>DATA[b].lat-DATA[a].lat);
const col={};names.forEach((s,i)=>col[s]=PAL[i%PAL.length]);
const NB=new Set(['GNW','HDW','NLLB','COR','PGC','COLT']);
const map=L.map('map').setView([46.0,-123.0],6);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',{attribution:'&copy; OSM, &copy; CARTO',subdomains:'abcd',maxZoom:19}).addTo(map);
let cur={s:null,fam:null};
function win(){return document.querySelector('input[name=win]:checked').value;}
function dson(){return document.getElementById('deseason').checked;}
function wkey(){return win()+(dson()?'d':'');}
function wlab(){var w=win();var l=w==='13'?'1&ndash;3 s':w==='35'?'3&ndash;5 s':'2&ndash;4 s';return l+(dson()?' deseason':'');}
function cc(s){return DATA[s]['cc'+wkey()];}
function blk(s){var DD=dailyOf(s);return DD?DD['s'+wkey()]:null;}
function noData(s){Plotly.purge('plot');document.getElementById('info').innerHTML='<b>'+s+'</b>: no data in the '+wlab()+' set.';}
const Y2={title:'tremor / mo',overlaying:'y',side:'right',showgrid:false,rangemode:'tozero',titlefont:{color:'#ea580c',size:11},tickfont:{color:'#ea580c',size:10}};
function showStation(s){cur={s:s,fam:null};hiStation(s);setActiveChip(s);const D=blk(s);if(!D)return noData(s);
 const tr=[];const T=DATA[s].tremor;
 if(T)tr.push({x:T.dates,y:T.total,type:'bar',yaxis:'y2',marker:{color:'rgba(234,88,12,0.6)'},name:'tremor/mo',hoverinfo:'x+y'});
 const ids=Object.keys(D.byfam);let X=[],Y=[];
 for(const id of ids){const a=D.byfam[id];for(let i=0;i<a.length;i++){X.push(D.dates[i]);Y.push(a[i]);}X.push(null);Y.push(null);}
 tr.push({x:X,y:Y,type:'scattergl',mode:'lines',line:{color:'rgba(71,85,105,0.33)',width:0.7},hoverinfo:'skip'});
 tr.push({x:D.dates,y:D.med,type:'scattergl',mode:'lines',line:{color:col[s],width:2.4},name:'cross-family median'});
 Plotly.react('plot',tr,{margin:{l:54,r:48,t:36,b:34},showlegend:false,bargap:0.1,title:{text:s+'  '+wlab()+'  &mdash; '+ids.length+' patches (daily) + cross-family median  &middot;  cc '+(cc(s)||'?'),font:{size:13}},yaxis:{title:'dv/v (%)',zeroline:true,zerolinecolor:'#e33',range:[-0.6,0.6]},yaxis2:Y2,xaxis:{title:'date'}},{responsive:true});
 document.getElementById('info').innerHTML='<b class=k>'+s+'</b> '+DATA[s].lat.toFixed(2)+'&deg;N &mdash; all '+ids.length+' patches (daily 30-day rolling) + the per-day cross-family median (bold, '+wlab()+'). Click one patch to isolate it.';}
function showFam(s,id){cur={s:s,fam:id};hiPatch(s,id);setActiveChip(s);const D=blk(s);if(!D)return noData(s);
 const y=D.byfam[id];const f=DATA[s].fam.find(x=>x.id===id);
 if(!y){Plotly.purge('plot');document.getElementById('info').innerHTML='<b>'+s+'</b> patch '+id+': not in the '+wlab()+' set.';return;}
 const ftr=[];const T=DATA[s].tremor;
 if(T&&T.byfam[id])ftr.push({x:T.dates,y:T.byfam[id],type:'bar',yaxis:'y2',marker:{color:'rgba(234,88,12,0.65)'},name:'tremor/mo',hoverinfo:'x+y'});
 ftr.push({x:D.dates,y:D.med,type:'scattergl',mode:'lines',line:{color:'#94a3b8',width:1.4,dash:'dot'},name:'cross-family median'});
 ftr.push({x:D.dates,y:y,type:'scattergl',mode:'lines',line:{color:col[s],width:2},name:'this patch'});
 Plotly.react('plot',ftr,{margin:{l:54,r:48,t:36,b:34},showlegend:false,bargap:0.1,title:{text:s+' &mdash; patch '+id+'  ('+wlab()+')',font:{size:13}},yaxis:{title:'dv/v (%)',zeroline:true,zerolinecolor:'#e33',autorange:true},yaxis2:Y2,xaxis:{title:'date'}},{responsive:true});
 document.getElementById('info').innerHTML='<b class=k>'+s+'</b> &mdash; patch @ '+f.lat+'&deg;N '+f.lon+'&deg; (daily 30-day rolling, bold) vs the cross-family median (dotted). Orange = tremor within 20 km of this patch. Click the station marker for all patches.';}
function refresh(){if(cur.s){cur.fam?showFam(cur.s,cur.fam):showStation(cur.s);}}
const layers={},chips={};
names.forEach(s=>{const d=DATA[s];const c=col[s];const grp=L.layerGroup();const fams={};
 d.fam.forEach(f=>{
  const line=L.polyline([[d.lat,d.lon],[f.lat,f.lon]],{color:c,weight:1,opacity:0.18}).on('click',e=>{L.DomEvent.stop(e);showFam(s,f.id);}).addTo(grp);
  const dot=L.circleMarker([f.lat,f.lon],{radius:2.5,color:c,fillColor:c,weight:0,fillOpacity:0.55}).on('click',e=>{L.DomEvent.stop(e);showFam(s,f.id);}).bindTooltip(s+' patch '+f.id,{direction:'top'}).addTo(grp);
  fams[f.id]={line:line,dot:dot};});
 grp.addTo(map);
 const isnb=NB.has(s);
 const mstyle={radius:isnb?8:6,color:isnb?'#b35900':'#222',weight:isnb?2.8:1.4,fillColor:c,fillOpacity:0.95};
 const marker=L.circleMarker([d.lat,d.lon],mstyle).on('click',e=>{L.DomEvent.stop(e);showStation(s);}).bindTooltip(s+(isnb?' (broadband)':'')+' - '+d.lat.toFixed(1)+'N',{direction:'top',permanent:isnb,className:isnb?'nblabel':''}).addTo(map);
 layers[s]={marker:marker,mstyle:mstyle,c:c,fams:fams};});
function dimAll(){for(const k in layers){const o=layers[k];
 o.marker.setStyle(Object.assign({},o.mstyle,{fillOpacity:0.4}));
 for(const id in o.fams){o.fams[id].line.setStyle({color:o.c,opacity:0.05,weight:1});o.fams[id].dot.setStyle({fillColor:o.c,color:o.c,weight:0,radius:2,fillOpacity:0.12});}}}
function emph(o){o.marker.setStyle(Object.assign({},o.mstyle,{radius:o.mstyle.radius+2,fillOpacity:1}));}
function hiStation(s){if(!layers[s])return;dimAll();const o=layers[s];emph(o);
 for(const id in o.fams){o.fams[id].line.setStyle({color:o.c,opacity:0.6,weight:2});o.fams[id].dot.setStyle({fillColor:o.c,color:o.c,weight:0,radius:4,fillOpacity:0.9});}}
function hiPatch(s,id){if(!layers[s])return;dimAll();const o=layers[s];emph(o);
 if(o.fams[id]){o.fams[id].line.setStyle({color:o.c,opacity:0.95,weight:4});o.fams[id].dot.setStyle({fillColor:o.c,color:'#000',weight:1.8,radius:6,fillOpacity:1});}}
const sl=document.getElementById('stalist');
names.forEach(s=>{const b=document.createElement('span');b.textContent=s;b.className='chip'+(NB.has(s)?' nb':'');b.title=DATA[s].lat.toFixed(2)+'N';b.onclick=()=>showStation(s);sl.appendChild(b);chips[s]=b;});
function setActiveChip(s){for(const k in chips)chips[k].classList.toggle('active',k===s);}
document.getElementById('stafilter').addEventListener('input',function(e){const q=e.target.value.toUpperCase();for(const k in chips)chips[k].style.display=k.toUpperCase().indexOf(q)>=0?'':'none';});
document.querySelectorAll('input[name=win]').forEach(r=>r.addEventListener('change',refresh));
document.getElementById('deseason').addEventListener('change',refresh);
// draggable divider: widen the side panel by dragging left
(function(){const dr=document.getElementById('drag'),sd=document.getElementById('side');let on=false;
 dr.addEventListener('mousedown',e=>{on=true;dr.classList.add('on');document.body.style.cursor='col-resize';e.preventDefault();});
 window.addEventListener('mousemove',e=>{if(!on)return;const w=Math.min(window.innerWidth-260,Math.max(380,window.innerWidth-e.clientX));sd.style.width=w+'px';});
 window.addEventListener('mouseup',()=>{if(!on)return;on=false;dr.classList.remove('on');document.body.style.cursor='';map.invalidateSize();Plotly.Plots.resize('plot');});})();
showStation('GNW');
</script></body></html>'''

out = HTML.replace('__DATA__', djs).replace('__DAILY__', djson)
with open('fault_tomography/cascadia_dvv_map.html', 'w') as f:
    f.write(out)
print(f'wrote fault_tomography/cascadia_dvv_map.html ({len(out)/1e6:.2f} MB), {len(DATA)} stations')
