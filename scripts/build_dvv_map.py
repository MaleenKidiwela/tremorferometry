#!/usr/bin/env python
"""Rebuild the standalone interactive Cascadia dv/v map with the corrected rolling-stack data.
Window toggle = 1-4 s vs 2-4 s (SAME 30-day causal rolling + SVD-Wiener method, only the coda
window differs) so clicking a family shows the S-dilution directly.
Station click -> all families (bold) + cross-family median.  Family click -> that family, bold + isolated.
Reads data/daily_dvv_{S}_{1to4,2to4}_roll.csv -> fault_tomography/cascadia_dvv_map.html
"""
import json
import numpy as np, pandas as pd

STA = {  # station: (lat, lon)
 'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732),
}
CCMIN = 0.7; MINMON = 6

def famll(pid):
    a = pid.split('__')[0].split('_'); return float(a[0]), float(a[1])

def window_block(csv):
    try:
        d = pd.read_csv(csv)
    except FileNotFoundError:
        return None
    d = d[d.cc_max > CCMIN].copy()
    if not len(d): return None
    d['date'] = pd.to_datetime(d['date']); d['ym'] = d['date'].dt.to_period('M').astype(str)
    piv = d.groupby(['ym', 'patch'])['dvv'].median().unstack('patch') * 100.0
    months = sorted(piv.index); piv = piv.reindex(months)
    byfam = {}
    for fid in piv.columns:
        col = piv[fid]
        if col.notna().sum() < MINMON: continue
        byfam[fid] = [None if pd.isna(v) else round(float(v), 3) for v in col.values]
    if not byfam: return None
    med = piv[list(byfam)].median(axis=1)
    med = [None if pd.isna(v) else round(float(v), 3) for v in med.values]
    return {'dates': months, 'byfam': byfam, 'med': med, 'cc': round(float(d.cc_max.mean()), 2)}

DATA = {}
for s, (la, lo) in STA.items():
    s14 = window_block(f'data/daily_dvv_{s}_1to4_roll.csv')
    s24 = window_block(f'data/daily_dvv_{s}_2to4_roll.csv')
    s24d = window_block(f'data/daily_dvv_{s}_2to4_des.csv')   # deseasoned 2-4 s
    if s24 is None and s14 is None and s24d is None: print('  skip', s); continue
    fam_ids = sorted(set(list((s24 or {}).get('byfam', {})) + list((s14 or {}).get('byfam', {}))
                        + list((s24d or {}).get('byfam', {}))))
    fam = [{'id': fid, 'lat': round(famll(fid)[0], 3), 'lon': round(famll(fid)[1], 3)} for fid in fam_ids]
    DATA[s] = {'lat': la, 'lon': lo, 'fam': fam,
               'cc14': (s14 or {}).get('cc'), 'cc24': (s24 or {}).get('cc'), 'cc24d': (s24d or {}).get('cc'),
               's14': {k: s14[k] for k in ('dates', 'byfam', 'med')} if s14 else None,
               's24': {k: s24[k] for k in ('dates', 'byfam', 'med')} if s24 else None,
               's24d': {k: s24d[k] for k in ('dates', 'byfam', 'med')} if s24d else None}
    print(f'  {s}: {len(fam)} families  cc24 {DATA[s]["cc24"]}  cc24d {DATA[s]["cc24d"]}')

js = json.dumps(DATA, separators=(',', ':'))
with open('fault_tomography/dvv_map_data.json', 'w') as f:
    f.write(js)

HTML = '''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia LFE-coda dv/v map</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>html,body{margin:0;height:100%;font-family:system-ui,Arial,sans-serif}
#wrap{display:flex;height:100vh}#map{flex:1.2;height:100%}
#side{flex:1;display:flex;flex-direction:column;border-left:1px solid #ccc;min-width:420px}
#bar{padding:8px 12px;background:#f4f4f4;border-bottom:1px solid #ddd}#bar h3{margin:0 0 6px;font-size:14px}
.win{margin-right:14px;cursor:pointer}
#plot{height:48vh;width:100%;border-bottom:1px solid #eee}
#info{padding:8px 12px;font-size:12px;color:#444;line-height:1.5}
.hint{color:#888;font-style:italic}b.k{color:#222}
#stalist{margin-top:6px;line-height:1.95}
.chip{display:inline-block;padding:1px 6px;margin:1px 2px;font-size:11px;border:1px solid #bbb;border-radius:9px;cursor:pointer;background:#fff}
.chip:hover{background:#eef}
.chip.nb{border-color:#e07000;color:#b35900;font-weight:bold;background:#fff6ec}
.nblabel{background:rgba(255,246,236,0.92);border:1px solid #e07000;font-weight:bold;font-size:10px;padding:0 3px;color:#b35900}</style></head>
<body><div id="wrap"><div id="map"></div><div id="side">
<div id="bar"><h3>Cascadia LFE-coda dv/v &mdash; click a <b>station</b> (all families) or a <b>family path/dot</b> (that family). 30-day causal rolling stack + SVD-Wiener.</h3>
<label class="win"><input type="radio" name="win" value="24" checked> 2&ndash;4 s</label>
<label class="win"><input type="radio" name="win" value="24d"> 2&ndash;4 s deseason</label>
<label class="win"><input type="radio" name="win" value="14"> 1&ndash;4 s (diluted)</label>
<div id="stalist"></div></div>
<div id="plot"></div>
<div id="info" class="hint">Toggle 2&ndash;4 s vs 1&ndash;4 s on the same family to see the direct-S flatten the dv/v ~50&times;.</div>
</div></div><script>
const DATA=__DATA__;
const PAL=['#e6194B','#3cb44b','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#fabed4','#469990','#dcbeff','#9A6324','#800000','#aaffc3','#808000','#000075','#a9a9a9','#e6194B','#3cb44b','#4363d8','#f58231','#911eb4','#42d4f4','#f032e6','#bfef45','#469990','#9A6324','#800000','#000075'];
const names=Object.keys(DATA).sort((a,b)=>DATA[b].lat-DATA[a].lat);
const col={};names.forEach((s,i)=>col[s]=PAL[i%PAL.length]);
const NB=new Set(['GNW','HDW','NLLB','COR','PGC','COLT']);
const map=L.map('map').setView([46.0,-123.0],6);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',{attribution:'&copy; OSM, &copy; CARTO',subdomains:'abcd',maxZoom:19}).addTo(map);
let cur={s:null,fam:null};
function win(){return document.querySelector('input[name=win]:checked').value;}
function wlab(){var w=win();return w==='14'?'1&ndash;4 s':w==='24d'?'2&ndash;4 s deseason':'2&ndash;4 s';}
function S(s){var w=win();return w==='14'?DATA[s].s14:w==='24d'?DATA[s].s24d:DATA[s].s24;}
function cc(s){var w=win();return w==='14'?DATA[s].cc14:w==='24d'?DATA[s].cc24d:DATA[s].cc24;}
function noData(s){Plotly.purge('plot');document.getElementById('info').innerHTML='<b>'+s+'</b>: no data in the '+wlab()+' set.';}
function showStation(s){cur={s:s,fam:null};const D=S(s);if(!D)return noData(s);
 const tr=[];Object.values(D.byfam).forEach(c=>tr.push({x:D.dates,y:c,mode:'lines',line:{color:'rgba(70,70,70,0.55)',width:1.1},hoverinfo:'skip',showlegend:false,connectgaps:false}));
 tr.push({x:D.dates,y:D.med,mode:'lines',line:{color:col[s],width:3.0},connectgaps:false});
 Plotly.react('plot',tr,{margin:{l:55,r:15,t:36,b:36},showlegend:false,title:{text:s+'  '+wlab()+'  &mdash; '+Object.keys(D.byfam).length+' families + median (bold)  cc '+(cc(s)||'?'),font:{size:13}},yaxis:{title:'dv/v (%)',zeroline:true,zerolinecolor:'#e33',range:[-0.5,0.5]},xaxis:{title:'date'}},{responsive:true});
 document.getElementById('info').innerHTML='<b class=k>'+s+'</b> '+DATA[s].lat.toFixed(2)+'&deg;N &mdash; all '+Object.keys(D.byfam).length+' families (bold grey) + cross-family median ('+wlab()+'). Click one family to isolate it.';}
function showFam(s,id){cur={s:s,fam:id};const D=S(s);if(!D)return noData(s);
 const y=D.byfam[id];const f=DATA[s].fam.find(x=>x.id===id);
 if(!y){Plotly.purge('plot');document.getElementById('info').innerHTML='<b>'+s+'</b> family '+id+': not present in the '+wlab()+' set.';return;}
 Plotly.react('plot',[{x:D.dates,y:y,mode:'lines',line:{color:col[s],width:3.2},connectgaps:false}],
  {margin:{l:55,r:15,t:36,b:36},showlegend:false,title:{text:s+' &mdash; family '+id+'  ('+wlab()+')',font:{size:13}},yaxis:{title:'dv/v (%)',zeroline:true,zerolinecolor:'#e33',autorange:true},xaxis:{title:'date'}},{responsive:true});
 document.getElementById('info').innerHTML='<b class=k>'+s+'</b> &mdash; single family @ '+f.lat+'&deg;N '+f.lon+'&deg; (id '+id+'), bold. Toggle the window to see the S-dilution. Click the station marker for all families.';}
function refresh(){if(cur.s){cur.fam?showFam(cur.s,cur.fam):showStation(cur.s);}}
names.forEach(s=>{const d=DATA[s];const c=col[s];const grp=L.layerGroup();
 d.fam.forEach(f=>{L.polyline([[d.lat,d.lon],[f.lat,f.lon]],{color:c,weight:1,opacity:0.18}).on('click',e=>{L.DomEvent.stop(e);showFam(s,f.id);}).addTo(grp);
  L.circleMarker([f.lat,f.lon],{radius:2.5,color:c,weight:0,fillOpacity:0.55}).on('click',e=>{L.DomEvent.stop(e);showFam(s,f.id);}).bindTooltip(s+' fam '+f.id,{direction:'top'}).addTo(grp);});
 grp.addTo(map);
 const isnb=NB.has(s);
 L.circleMarker([d.lat,d.lon],{radius:isnb?8:6,color:isnb?'#b35900':'#222',weight:isnb?2.8:1.4,fillColor:c,fillOpacity:0.95}).on('click',e=>{L.DomEvent.stop(e);showStation(s);}).bindTooltip(s+(isnb?' (broadband)':'')+' - '+d.lat.toFixed(1)+'N',{direction:'top',permanent:isnb,className:isnb?'nblabel':''}).addTo(map);});
const sl=document.getElementById('stalist');
sl.innerHTML='<b style="font-size:11px">jump to station &mdash; <span style="color:#b35900">orange = non-borehole</span>:</b> ';
names.forEach(s=>{const b=document.createElement('span');b.textContent=s;b.className='chip'+(NB.has(s)?' nb':'');b.title=DATA[s].lat.toFixed(2)+'N';b.onclick=()=>showStation(s);sl.appendChild(b);});
document.querySelectorAll('input[name=win]').forEach(r=>r.addEventListener('change',refresh));
showStation('GNW');
</script></body></html>'''

out = HTML.replace('__DATA__', js)
with open('fault_tomography/cascadia_dvv_map.html', 'w') as f:
    f.write(out)
print(f'wrote fault_tomography/cascadia_dvv_map.html ({len(out)/1e6:.2f} MB), {len(DATA)} stations')
