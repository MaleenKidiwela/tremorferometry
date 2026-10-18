#!/usr/bin/env python
"""Interface dv/v explorer for a COARSE inversion product (0.5 x 0.5 deg boxes): same idea as interface_dvv_explorer.html —
interface mode (boxes) / site mode (stations), weekly slider + play, click -> daily series (30-d trailing frames) with
era lines; per-box recovery correction (divide by the box's mean self-transmission T_diag from the +1 % injection),
boxes with T_diag < 0.7 greyed (poorly resolved); reference 2019-24; 2009-10 excluded. Leaflet + Plotly.
Usage: python scripts/build_coarse_dvv_explorer.py fault_tomography/inversion/coarse05x05 figures/interface_dvv_explorer_05x05.html"""
import sys, json, gzip, base64, numpy as np, pandas as pd
RESDIR, OUT = sys.argv[1], sys.argv[2]; tag = RESDIR.rstrip("/").split("/")[-1]; RD = "fault_tomography/inversion/res_catalog"
z = np.load(f"{RESDIR}/inversion_4d_daily.npz", allow_pickle=True); days = pd.to_datetime(z["days"]); cells = [str(c) for c in z["cell"]]; M = len(cells)
T = pd.read_csv(f"{RESDIR}/{tag}_cells_summary.csv").set_index("cell"); Tc = np.array([T.T_diag.get(c, np.nan) for c in cells]); weak = (Tc < 0.7) | ~np.isfinite(Tc)
corr = np.where(weak, 1.0, np.clip(Tc, 0.3, None))
dgrid = pd.date_range("2011-01-01", days.max(), freq="D"); D = pd.DataFrame(np.asarray(z["MODEL"], float), index=days).reindex(dgrid).values / corr[None, :]
wk = pd.date_range("2011-01-06", dgrid.max(), freq="7D"); W = pd.DataFrame(D, index=dgrid).reindex(wk).values
S = pd.DataFrame(np.asarray(z["SITES"], float), index=days, columns=list(z["sta_tags"])).reindex(dgrid)
sll = pd.read_csv(f"{RD}/pairs.csv").groupby("tag")[["sta_lat", "sta_lon"]].first(); stas = [t for t in S.columns if t in sll.index]; S = S[stas]; SW = S.reindex(wk).values
et = pd.read_csv(f"{RD}/era_table.csv"); wl = {x.strip().upper() for x in open("data/era_split_whitelist.txt") if x.strip()}
ERAS = {r.tag: [x for x in str(r.boundaries).split(";") if x and x != "nan"] for _, r in et.iterrows() if str(r.tag).upper().replace("P90F40", "") in wl and isinstance(r.boundaries, str)}
def pack(a, vmax):
    q = np.round(np.clip(a / vmax, -1, 1) * 32000).astype(np.int16); q[np.isnan(a)] = -32768; return base64.b64encode(gzip.compress(q.tobytes(), 6)).decode()
VM = 3.0
boxes = [dict(id=c, la=float(c.split("_")[0]), lo=float(c.split("_")[1]), T=round(float(Tc[j]), 2) if np.isfinite(Tc[j]) else None, weak=bool(weak[j]), dep=round(float(z["depth"][j]), 1)) for j, c in enumerate(cells)]
meta = dict(tag=tag, boxes=boxes, n=M, days=[d.strftime("%Y-%m-%d") for d in dgrid], nd=len(dgrid), wk=[d.strftime("%Y-%m-%d") for d in wk], nw=len(wk), vmax=VM,
            sites=dict(tags=stas, lat=[round(float(sll.sta_lat[t]), 4) for t in stas], lon=[round(float(sll.sta_lon[t]), 4) for t in stas], n=len(stas), eras=ERAS))
pay = dict(cw=pack(W, VM), cd=pack(D, VM), sw=pack(SW, VM), sd=pack(S.values, VM))
HTML = r'''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Interface dv/v explorer — 0.5° boxes</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/><script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script><script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
<style>:root{--bg:#f1f5f9;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#0e7490}*{box-sizing:border-box}html,body{margin:0;height:100%;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,sans-serif;color:var(--ink);background:var(--bg)}
#wrap{display:flex;height:100vh}#map{flex:1;height:100%;min-width:200px}#drag{width:7px;flex:none;cursor:col-resize;background:var(--line)}#drag:hover{background:var(--accent)}#side{width:560px;flex:none;min-width:360px;display:flex;flex-direction:column;border-left:1px solid var(--line);background:var(--bg);overflow:auto}
#hdr{padding:10px 16px;background:var(--card);border-bottom:1px solid var(--line)}.title{font-size:15px;font-weight:700}.sub{font-size:11.5px;color:var(--muted);margin-top:3px;line-height:1.45}.panel{margin:10px 12px 0;background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 12px}
.lbl{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}select,button{font:inherit;font-size:12.5px;padding:4px 8px;border:1px solid var(--line);border-radius:6px;background:#fff}button.on{background:var(--accent);color:#fff;border-color:var(--accent)}
.row{display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-top:6px}#plot{flex:1;min-height:320px;margin:10px 12px;background:#fff;border:1px solid var(--line);border-radius:10px}#info{font-size:12px;color:var(--muted);padding:6px 16px;line-height:1.5}
.cbar{height:12px;border-radius:4px;margin-top:4px;border:1px solid var(--line)}.cbl{display:flex;justify-content:space-between;font-size:10.5px;color:var(--muted)}.leaflet-container{background:#e5e9f0}</style></head><body><div id="wrap"><div id="map"></div><div id="drag"></div><aside id="side">
<div id="hdr"><div class="title">Interface dv/v explorer — 0.5° × 0.5° boxes (__TAG__)</div><div class="sub"><b>Interface mode:</b> δβ/β per 0.5° box from the coarse daily joint inversion (corrected S-pinned kernel), divided by each box's self-recovery of a +1 % injection (T ≈ 0.98 typical); grey boxes = poorly resolved (T &lt; 0.7); reference 2019&ndash;24; 2011&rarr;. <b>Site mode:</b> per-station near-surface dv/v (site term, true sign). Weekly frames of the 30-d trailing stack; click a box or station for its daily series. Single-frame values are noise-dominated: read multi-year shapes; amplitudes provisional.</div></div>
<div class="panel"><div class="row"><button id="mInt" class="on">interface δβ/β</button><button id="mSit">site dv/v</button><span class="lbl" style="margin-left:8px">scale</span><select id="cscale"><option value="0.3">±0.3%</option><option value="0.5">±0.5%</option><option value="1" selected>±1%</option><option value="2">±2%</option></select><span class="lbl">date</span><b id="date"></b></div>
<div class="row"><button id="play">&#9654;</button><input type="range" id="t" min="0" max="0" value="0" step="1" style="flex:1;accent-color:var(--accent)"><button id="prev">&lsaquo;</button><button id="next">&rsaquo;</button></div>
<div id="cbar" class="cbar"></div><div class="cbl"><span id="cmin"></span><span>0</span><span id="cmax"></span></div></div><div id="plot"></div><div id="info">Click a box or a station.</div></aside></div><script>
const META=__META__,PAY=__PAY__,cache={};function dec(k,b){if(cache[k])return cache[k];const s=atob(b);const u=new Uint8Array(s.length);for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);const r=pako.inflate(u);cache[k]=new Int16Array(r.buffer,r.byteOffset,r.byteLength/2);return cache[k];}
const NB=META.n,NS=META.sites.n,NW=META.nw,ND=META.nd,VM=META.vmax;let cur={mode:'int',t:NW-1,sel:null,timer:null};
function cw(t,i){const q=dec('cw',PAY.cw)[t*NB+i];return q===-32768?null:q/32000*VM;}function sw(t,j){const q=dec('sw',PAY.sw)[t*NS+j];return q===-32768?null:q/32000*VM;}
function cd(i){const a=dec('cd',PAY.cd);const o=new Array(ND);for(let d=0;d<ND;d++){const q=a[d*NB+i];o[d]=q===-32768?null:q/32000*VM;}return o;}function sd(j){const a=dec('sd',PAY.sd);const o=new Array(ND);for(let d=0;d<ND;d++){const q=a[d*NS+j];o[d]=q===-32768?null:q/32000*VM;}return o;}
const RDBU=[[5,48,97],[67,147,195],[247,247,247],[214,96,77],[103,0,31]];function lerp(a,b,x){return a.map((c,i)=>Math.round(c+(b[i]-c)*x));}function col(x){x=Math.max(0,Math.min(1,x));const k=4*x,i=Math.min(3,Math.floor(k));return 'rgb('+lerp(RDBU[i],RDBU[i+1],k-i).join(',')+')';}
function colv(v){if(v==null)return '#e5e7eb';const s=+document.getElementById('cscale').value;return col((v/s+1)/2);}
const map=L.map('map').setView([46.5,-123.5],6);L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',{attribution:'&copy; OpenStreetMap',maxZoom:19}).addTo(map);
const gB=L.layerGroup().addTo(map),gS=L.layerGroup();const rects=[],tris=[];let selmk=null;function mark(la,lo){if(selmk)map.removeLayer(selmk);selmk=L.circleMarker([la,lo],{radius:9,color:'#000',weight:2.5,fill:false}).addTo(map);}
META.boxes.forEach((b,i)=>{const r=L.rectangle([[b.la,b.lo],[b.la+0.5,b.lo+0.5]],{weight:0.7,color:'#334155',fillOpacity:0.85,fillColor:'#e5e7eb'}).addTo(gB);r.on('click',()=>{cur.sel={type:'b',i:i};mark(b.la+0.25,b.lo+0.25);series();});
 r.on('mouseover',()=>{const x=cw(cur.t,i);r.bindTooltip('box '+b.id+' (mean depth '+b.dep+' km, T='+b.T+(b.weak?', poorly resolved':'')+'): '+(x==null?'n/a':x.toFixed(3)+' %'),{sticky:true}).openTooltip();});rects.push(r);});
META.sites.tags.forEach((s,j)=>{const la=META.sites.lat[j],lo=META.sites.lon[j];const m=L.circleMarker([la,lo],{radius:5,weight:0.6,color:'#334155',fillOpacity:0.95,fillColor:'#e5e7eb'}).addTo(gS);m.on('click',()=>{cur.sel={type:'s',i:j};mark(la,lo);series();});m.bindTooltip(s);tris.push(m);});
function paint(){const t=cur.t;document.getElementById('date').textContent=META.wk[t];if(cur.mode==='int'){META.boxes.forEach((b,i)=>rects[i].setStyle({fillColor:b.weak?'#d1d5db':colv(cw(t,i)),fillOpacity:b.weak?0.5:0.85}));}else{for(let j=0;j<NS;j++)tris[j].setStyle({fillColor:colv(sw(t,j))});}
 const s=+document.getElementById('cscale').value;document.getElementById('cbar').style.background='linear-gradient(90deg,'+RDBU.map((c,k)=>'rgb('+c.join(',')+') '+(25*k)+'%').join(',')+')';document.getElementById('cmin').textContent='-'+s+'%';document.getElementById('cmax').textContent='+'+s+'%';}
function roll(y,n){const o=new Array(y.length).fill(null);let acc=[];for(let i=0;i<y.length;i++){acc.push(y[i]);if(acc.length>n)acc.shift();const v=acc.filter(x=>x!=null);o[i]=v.length>=Math.max(10,n/4)?v.reduce((a,b)=>a+b,0)/v.length:null;}return o;}
function series(){const s=cur.sel;if(!s)return;let y,title,shapes=[{type:'line',xref:'x',x0:META.wk[cur.t],x1:META.wk[cur.t],yref:'paper',y0:0,y1:1,line:{color:'#000',width:1,dash:'dot'}}];
 if(s.type==='b'){const b=META.boxes[s.i];y=cd(s.i);title='box '+b.id+' — δβ/β (mean depth '+b.dep+' km, T='+b.T+(b.weak?', POORLY RESOLVED':'')+')';}
 else{const tg=META.sites.tags[s.i];y=sd(s.i);title='station '+tg+' — site dv/v';(META.sites.eras[tg]||[]).forEach(d=>shapes.push({type:'line',xref:'x',x0:d,x1:d,yref:'paper',y0:0,y1:1,line:{color:'#dc2626',width:1,dash:'dot'}}));}
 Plotly.react('plot',[{x:META.days,y:y,type:'scattergl',mode:'lines',line:{color:'#cbd5e1',width:0.5},name:'daily (30-d trailing frames)'},{x:META.days,y:roll(y,90),type:'scattergl',mode:'lines',line:{color:'#0e7490',width:1.6},name:'90-d mean'},{x:META.days,y:roll(y,365),type:'scattergl',mode:'lines',line:{color:'#7c2d12',width:1.4},name:'365-d mean'}],
  {margin:{l:50,r:20,t:34,b:30},title:{text:title,font:{size:11}},showlegend:true,legend:{orientation:'h',y:1.12,font:{size:9}},yaxis:{title:'dv/v (%)',zeroline:true},shapes:shapes,paper_bgcolor:'#fff',plot_bgcolor:'#fff'},{responsive:true});
 document.getElementById('info').innerHTML='<b>'+title+'</b>. Grey daily, teal 90-d mean, brown 365-d mean; black dotted = slider date; red dotted = instrument era boundary.';}
const sl=document.getElementById('t');sl.max=NW-1;sl.value=NW-1;sl.oninput=()=>{cur.t=+sl.value;paint();if(cur.sel)series();};document.getElementById('prev').onclick=()=>{sl.value=Math.max(0,cur.t-1);sl.oninput();};document.getElementById('next').onclick=()=>{sl.value=Math.min(NW-1,cur.t+1);sl.oninput();};
document.getElementById('play').onclick=function(){if(cur.timer){clearInterval(cur.timer);cur.timer=null;this.innerHTML='&#9654;';return;}this.innerHTML='&#10074;&#10074;';cur.timer=setInterval(()=>{sl.value=(cur.t+1)%NW;sl.oninput();},120);};document.getElementById('cscale').onchange=paint;
function setMode(m){cur.mode=m;document.getElementById('mInt').classList.toggle('on',m==='int');document.getElementById('mSit').classList.toggle('on',m==='sit');if(m==='int'){map.removeLayer(gS);gB.addTo(map);}else{map.removeLayer(gB);gS.addTo(map);}paint();}
document.getElementById('mInt').onclick=()=>setMode('int');document.getElementById('mSit').onclick=()=>setMode('sit');
(function(){const dr=document.getElementById('drag'),sd_=document.getElementById('side');let on=false;dr.addEventListener('mousedown',e=>{on=true;document.body.style.cursor='col-resize';e.preventDefault();});window.addEventListener('mousemove',e=>{if(!on)return;sd_.style.width=Math.min(window.innerWidth-220,Math.max(360,window.innerWidth-e.clientX))+'px';});window.addEventListener('mouseup',()=>{if(!on)return;on=false;document.body.style.cursor='';map.invalidateSize();Plotly.Plots.resize('plot');});})();
window.addEventListener('resize',()=>{map.invalidateSize();Plotly.Plots.resize('plot');});paint();
</script></body></html>'''
out = HTML.replace("__TAG__", tag).replace("__META__", json.dumps(meta, separators=(",", ":"))).replace("__PAY__", json.dumps(pay, separators=(",", ":")))
open(OUT, "w").write(out); print(f"wrote {OUT} ({len(out) / 1e6:.1f} MB); boxes {M}, weak {int(weak.sum())}, stations {len(stas)}")
