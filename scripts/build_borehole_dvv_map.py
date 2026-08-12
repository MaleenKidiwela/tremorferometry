#!/usr/bin/env python
"""Interactive borehole LFE-coda dv/v map (current causality workflow).
Stations = boreholes + broadband/fleet with a finalized summary json; families = CERTIFIED only (causality
cert, else fwd-vs-rev ratio>1.5). Reads data/daily_dvv_<TAG>_Z_2to4.csv. Click a station -> all certified
families' daily dv/v + cross-family median; a family SLIDER above the plot steps through each family (median
= black dotted, family = station color). y-axis autoscales+LOCKS at station-select (deseason & the slider
never rescale, so raw vs deseasoned compare like-for-like); x-axis clips to the actual dv/v span (no empty
2010-2026 padding). Deseason toggle (annual+semiannual removed server-side) + local tremor overlay.
-> figures/borehole_dvv_map.html"""
import os, json, gzip, base64, glob
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import f as fdist

COORD = {r.sta: (r.lat, r.lon) for _, r in pd.read_csv("data/candidate_stations_post2020.csv").iterrows()}
# broadband/land stations absent from the post-2020 candidate file (permanent networks).
COORD.update({"PGC": (48.6498, -123.4521), "SHB": (49.5990, -123.8770), "CLRS": (48.8203, -124.1309)})
CCMIN, MINDAYS = 0.6, 60


def famll(pid):
    a = pid.split('__')[0].split('_'); return float(a[0]), float(a[1])


def certfams(s):
    if os.path.exists(f"data/{s}_causality_cert.csv"):
        c = pd.read_csv(f"data/{s}_causality_cert.csv"); return set(c[c.reliable].fam)
    if os.path.exists(f"data/{s}_fwd_vs_rev_coda.csv"):
        return set(pd.read_csv(f"data/{s}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    return None


def deseason_series(dates, vals):
    v = np.array([np.nan if x is None else x for x in vals], float)
    ok = ~np.isnan(v)
    if ok.sum() < 90: return vals
    t = (pd.to_datetime(dates) - pd.to_datetime(dates[0])).days.values.astype(float)
    yr = 365.25
    X = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr),
                         np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    beta, *_ = np.linalg.lstsq(X[ok], v[ok], rcond=None)
    # CONDITIONAL deseason (per family -- seasonality is location-specific): only remove the cycle where
    # it is statistically present. F-test the 4 harmonic terms vs the trend-only (intercept+t) model; if
    # not significant (p>=0.05) the family has NO real seasonality -> leave it UNCHANGED (don't fit noise).
    n = int(ok.sum())
    rss_full = float(((v[ok] - X[ok] @ beta) ** 2).sum())
    bt, *_ = np.linalg.lstsq(X[ok, :2], v[ok], rcond=None)         # trend-only reduced model
    rss_red = float(((v[ok] - X[ok, :2] @ bt) ** 2).sum())
    F = ((rss_red - rss_full) / 4) / (rss_full / (n - 6)) if rss_full > 0 and n > 6 else 0.0
    if 1 - fdist.cdf(F, 4, n - 6) >= 0.05:
        return vals                                                # no significant seasonality -> unchanged
    # subtract ONLY the annual+semiannual harmonics (cols 2:); KEEP the mean (col 0) and the linear
    # trend (col 1) -- the trend term is in the fit so harmonics don't absorb it, but deseason must not
    # remove the long-term signal or the level (else raw vs deseasoned won't compare like-for-like).
    r = v - X[:, 2:] @ beta[2:]
    return [None if np.isnan(x) else round(float(x), 3) for x in r]


def blocks_for(tag, fams):
    # dv/v preference: UNIFORM eps_max=0.05 re-run (2026-08-12, certified families only, rail verified cleared)
    # > the Jul-28 _eps05 pass (inconsistent per-station ranges 0.023-0.049) > the original +/-2% clipped file.
    ucsv = f"data/daily_dvv_{tag}_Z_2to4_unif05.csv"
    ecsv = f"data/daily_dvv_{tag}_Z_2to4_eps05.csv"
    csv = ucsv if os.path.exists(ucsv) else (ecsv if os.path.exists(ecsv) else f"data/daily_dvv_{tag}_Z_2to4.csv")
    if not os.path.exists(csv): return None
    d = pd.read_csv(csv)
    if csv in (ucsv, ecsv):   # cc-gate the extended large-stretch TAIL: drop |dv/v|>2% with cc<0.8 (cycle-skip runaways)
        d = d[~((d.dvv.abs() > 0.02) & (d.cc_max < 0.80))].copy()
    return _blocks(d, fams, mcsv=f"data/daily_dvv_{tag}_MIRROR_2to4.csv")


def _blocks(d, fams, mcsv=None):
    d = d[(d.cc_max > CCMIN) & (d.patch.isin(fams))].copy()
    if not len(d): return None
    d['ds'] = pd.to_datetime(d['date']).dt.strftime('%Y-%m-%d')
    piv = d.pivot_table(index='ds', columns='patch', values='dvv', aggfunc='first') * 100.0
    days = sorted(piv.index); piv = piv.reindex(days)
    keep = [fid for fid in piv.columns if piv[fid].notna().sum() >= MINDAYS]
    if not keep: return None
    piv = piv[keep]

    def blk(P):
        byfam = {fid: [None if pd.isna(v) else round(float(v), 3) for v in P[fid].values] for fid in P.columns}
        md = P.median(axis=1)
        return {'dates': days, 'byfam': byfam, 'med': [None if pd.isna(v) else round(float(v), 3) for v in md.values]}

    def deseason_blk(P):
        Pd = pd.DataFrame({fid: deseason_series(days, [None if pd.isna(v) else float(v) for v in P[fid].values])
                           for fid in P.columns}, index=days)
        return blk(Pd)

    raw = blk(piv); raw['cc'] = round(float(d.cc_max.mean()), 2)
    des = deseason_blk(piv)
    # mirror-cleaned: coda - beta*mirror  (beta from the station's deseasoned common-mode regression)
    if mcsv and os.path.exists(mcsv):
        m = pd.read_csv(mcsv); m = m[(m.cc_max > CCMIN) & (m.patch.isin(keep))]
        m['ds'] = pd.to_datetime(m['date']).dt.strftime('%Y-%m-%d')
        mpiv = m.pivot_table(index='ds', columns='patch', values='dvv', aggfunc='first') * 100.0
        mpiv = mpiv.reindex(index=days, columns=keep)
        cc_c = pd.Series(deseason_series(days, list(piv.median(axis=1).values)), index=days)
        cc_m = pd.Series(deseason_series(days, list(mpiv.median(axis=1).values)), index=days)
        jj = pd.concat([cc_c, cc_m], axis=1).dropna()
        beta = float(np.polyfit(jj.iloc[:, 1], jj.iloc[:, 0], 1)[0]) if len(jj) > 30 else 0.0
        cpiv = piv - beta * mpiv
        clean = blk(cpiv); cleand = deseason_blk(cpiv)
    else:
        clean, cleand = None, None   # no mirror dv/v yet -> toggle falls back to raw in JS
    return raw, des, clean, cleand


# tremor context
LAT0, LON0, TREMOR_R_KM = 47.5, -124.0, 20.0
_tr = pd.read_csv('catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
_tr['t'] = pd.to_datetime(_tr['time'], errors='coerce'); _tr = _tr.dropna(subset=['t', 'lat', 'lon'])
_tr_ym = _tr['t'].dt.to_period('M').astype(str).values
def _km(lat, lon): return ((np.asarray(lon, float)-LON0)*111*np.cos(np.radians(LAT0)), (np.asarray(lat, float)-LAT0)*111)
_trx, _try = _km(_tr['lat'].values, _tr['lon'].values); _tr_xy = np.column_stack([_trx, _try])
def tremor_for(fam):
    if not fam: return None
    fx, fy = _km([f['lat'] for f in fam], [f['lon'] for f in fam]); fxy = np.column_stack([fx, fy])
    dmin, _ = cKDTree(fxy).query(_tr_xy, k=1); near = dmin < TREMOR_R_KM
    if int(near.sum()) == 0: return None
    months = sorted(set(_tr_ym[near])); mi = {m: i for i, m in enumerate(months)}
    total = np.zeros(len(months), int)
    for ym in _tr_ym[near]: total[mi[ym]] += 1
    return {'dates': months, 'total': [int(v) for v in total]}


DATA, DAILY = {}, {}
# NEW-pipeline stations only: identified by the p90f40 product naming (old margin workflow used _cal).
_stations = sorted({os.path.basename(f).split("daily_dvv_")[1].split("p90f40")[0].upper()
                    for f in glob.glob("data/daily_dvv_*p90f40_Z_2to4.csv")})
print(f"new-pipeline stations: {_stations}")
for s in _stations:
    if s not in COORD: continue
    fams = certfams(s.lower())
    if not fams: continue
    bl = blocks_for(f"{s}p90f40", fams)   # dv/v CSVs are named with UPPERCASE station tag (B001p90f40)
    if not bl: continue
    raw, des, clean, cleand = bl
    la, lo = COORD[s]
    fam_ids = sorted(raw['byfam'])
    fam = [{'id': fid, 'lat': round(famll(fid)[0], 3), 'lon': round(famll(fid)[1], 3)} for fid in fam_ids]
    cc = raw.pop('cc')
    DATA[s] = {'lat': la, 'lon': lo, 'fam': fam, 'cc': cc, 'nfam': len(fam), 'tremor': tremor_for(fam)}
    payload = json.dumps({'s24': raw, 's24d': des},
                         separators=(',', ':')).encode()
    DAILY[s] = base64.b64encode(gzip.compress(payload, 6)).decode()
    print(f"  {s}: {len(fam)} certified fam -> {len(DAILY[s])/1e6:.2f} MB gz")

# --- broadband/land entries (borehole loop above discovers p90f40 products only). PGC = ONE marker, both
#     instrument eras (BHZ 2010-2017 + HHZ 2018-2026) in ONE figure. Each era keeps its OWN per-era reference
#     (unjoined: never averaged across the Aug-2017 sensor swap); a synthetic all-NaN day at the swap breaks
#     the line so the two segments read as two instruments, not one continuous series. Fleet broadbands add here.
def bb_entry(parts, brk):
    """parts=[(dvv_tag, cert_stem), ...] concatenated into ONE figure. brk=True inserts a synthetic all-NaN seam
    between different-INSTRUMENT eras (PGC/SHB: per-era refs, unjoined); False merges same-instrument cohorts
    (CLRS original+new families, one continuous series)."""
    dfs, allf = [], set()
    for i, (tag, stem) in enumerate(parts):
        f = certfams(stem); csv = f"data/daily_dvv_{tag}_Z_2to4.csv"
        if not f or not os.path.exists(csv): continue
        d = pd.read_csv(csv); d = d[d.patch.isin(f)]
        allf |= f; dfs.append(d)
        if brk and i < len(parts) - 1 and len(d):
            seam = (pd.to_datetime(d.date).max() + pd.Timedelta(days=30)).strftime('%Y-%m-%d')
            dfs.append(pd.DataFrame({'patch': [sorted(f)[0]], 'date': [seam], 'dvv': [np.nan], 'cc_max': [1.0], 'n_stack': [1]}))
    if not dfs: return None, None
    return _blocks(pd.concat(dfs, ignore_index=True), allf, mcsv=None), allf

BB = {                                                                   # key -> ([(dvv_tag, cert_stem)...], break_eras)
    "PGC":  ([("PGCbhz", "pgcbhz"), ("PGChhz", "pgchhz")], True),        # BHZ 2010-2017 + HHZ 2018-2026
    "SHB":  ([("SHB", "shb"), ("SHBhhz", "shbhhz")], True),             # BHZ 2010-2019 + HHZ 2020-2026
    "CLRS": ([("CLRS", "clrs"), ("CLRSnew", "clrsnew")], False),        # 2 cohorts, same HHZ -> merged
}
for key, (parts, brk) in BB.items():
    if key not in COORD: continue
    res, _ = bb_entry(parts, brk)
    if not res: continue
    raw, des, clean, cleand = res
    la, lo = COORD[key]
    fam_ids = sorted(raw['byfam'])
    fam = [{'id': fid, 'lat': round(famll(fid)[0], 3), 'lon': round(famll(fid)[1], 3)} for fid in fam_ids]
    cc = raw.pop('cc')
    DATA[key] = {'lat': la, 'lon': lo, 'fam': fam, 'cc': cc, 'nfam': len(fam), 'tremor': tremor_for(fam)}
    payload = json.dumps({'s24': raw, 's24d': des},
                         separators=(',', ':')).encode()
    DAILY[key] = base64.b64encode(gzip.compress(payload, 6)).decode()
    print(f"  {key}: {len(fam)} certified fam -> {len(DAILY[key])/1e6:.2f} MB gz")

# --- fleet broadband stations (single-band, auto-discovered from their dv/v products). INCLUDE bar: >=20
#     certified reliable families. Coords from the fleet order. Picks up every completed fleet station.
_fleet_coord = {r.sta: (r.lat, r.lon) for _, r in pd.read_csv("data/broadband_fleet_order.csv").iterrows()}
_anchor_tags = {"PGCbhz", "PGChhz", "SHB", "SHBhhz", "CLRS", "CLRSnew"}
for csv in sorted(glob.glob("data/daily_dvv_*_Z_2to4.csv")):
    tag = os.path.basename(csv).split("daily_dvv_")[1].split("_Z_2to4")[0]
    if "p90f40" in tag or "MIRROR" in tag or tag in _anchor_tags or tag in DATA or tag not in _fleet_coord:
        continue
    fams = certfams(tag.lower())
    if not fams or len(fams) < 20:          # INCLUDE bar
        continue
    bl = blocks_for(tag, fams)
    if not bl:
        continue
    raw, des, clean, cleand = bl
    la, lo = _fleet_coord[tag]
    fam_ids = sorted(raw['byfam'])
    fam = [{'id': fid, 'lat': round(famll(fid)[0], 3), 'lon': round(famll(fid)[1], 3)} for fid in fam_ids]
    cc = raw.pop('cc')
    DATA[tag] = {'lat': la, 'lon': lo, 'fam': fam, 'cc': cc, 'nfam': len(fam), 'tremor': tremor_for(fam)}
    payload = json.dumps({'s24': raw, 's24d': des},
                         separators=(',', ':')).encode()
    DAILY[tag] = base64.b64encode(gzip.compress(payload, 6)).decode()
    print(f"  {tag}: {len(fam)} certified fam (fleet broadband) -> {len(DAILY[tag])/1e6:.2f} MB gz")

djs = json.dumps(DATA, separators=(',', ':')); djson = json.dumps(DAILY, separators=(',', ':'))
print(f"{len(DATA)} stations | DATA {len(djs)/1e6:.2f} MB | DAILY {len(djson)/1e6:.2f} MB")

HTML = r'''<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia borehole LFE-coda dv/v</title>
<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>
<script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/pako@2.1.0/dist/pako.min.js"></script>
<style>
:root{--bg:#f1f5f9;--card:#fff;--ink:#0f172a;--muted:#64748b;--line:#e2e8f0;--accent:#2563eb;--soft:#eff6ff;--shadow:0 1px 3px rgba(15,23,42,.08)}
*{box-sizing:border-box}html,body{margin:0;height:100%;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",system-ui,Roboto,sans-serif;color:var(--ink);background:var(--bg)}
#wrap{display:flex;height:100vh}#map{flex:1;height:100%;min-width:200px}
#drag{width:7px;flex:none;cursor:col-resize;background:var(--line);transition:background .12s}#drag:hover,#drag.on{background:var(--accent)}
#side{width:540px;flex:none;min-width:340px;display:flex;flex-direction:column;background:var(--bg);border-left:1px solid var(--line)}
#hdr{padding:14px 18px;background:var(--card);border-bottom:1px solid var(--line)}#hdr .title{font-size:16px;font-weight:700}#hdr .sub{font-size:12px;color:var(--muted);margin-top:3px;line-height:1.45}
.panel{margin:12px 14px 0;background:var(--card);border:1px solid var(--line);border-radius:12px}
.lbl{font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:.04em;color:var(--muted)}
#ctrl{padding:12px 14px}.dswrap{display:flex;align-items:center;gap:9px;cursor:pointer;font-size:12.5px;color:var(--ink);user-select:none}
.dswrap input{display:none}.switch{width:34px;height:19px;border-radius:999px;background:#cbd5e1;position:relative;transition:.15s;flex:none}
.switch::after{content:'';position:absolute;top:2px;left:2px;width:15px;height:15px;border-radius:50%;background:#fff;transition:.15s}
.dswrap input:checked+.switch{background:var(--accent)}.dswrap input:checked+.switch::after{transform:translateX(15px)}
#plot{height:44vh;min-height:300px;padding:8px}#info{font-size:12.5px;color:var(--muted);line-height:1.55;padding:10px 16px}
#stawrap{flex:1;min-height:0;display:flex;flex-direction:column;padding:12px 14px}
.chip{display:inline-block;padding:4px 10px;margin:3px;font-size:12px;font-weight:600;border:1px solid var(--line);border-radius:999px;cursor:pointer;background:var(--card);color:#cbd5e1}
.chip:hover,.chip.active{border-color:var(--accent);color:var(--accent)}.chip.active{background:var(--soft)}
.leaflet-container{font-family:inherit;background:#e5e9f0}
.chip{box-shadow:var(--shadow)}
</style></head><body><div id="wrap"><div id="map"></div><div id="drag" title="drag to widen the panel"></div>
<aside id="side"><div id="hdr"><div class="title">Cascadia LFE-coda dv/v (boreholes + broadband)</div>
<div class="sub">causality-certified families &middot; daily 30-day trailing stack, 2&ndash;4 s coda &middot; click a <b>station</b> or a <b>family</b></div></div>
<div id="ctrl" class="panel"><label class="dswrap"><input type="checkbox" id="deseason"><span class="switch"></span><span>deseason &mdash; remove annual + semiannual cycle</span></label></div>
<div id="sliderwrap" class="panel" style="display:none;padding:11px 14px"><div style="display:flex;align-items:center;gap:12px"><span class="lbl">family</span><input type="range" id="famslider" min="0" max="0" value="0" step="1" style="flex:1;accent-color:var(--accent)"><span id="slabel" class="lbl" style="min-width:132px;text-align:right;color:var(--ink)"></span></div></div>
<div id="plot" class="panel"></div><div id="info">Click a station marker.</div>
<div id="stawrap"><div class="lbl">stations (certified family count)</div><div id="stalist"></div></div></aside></div><script>
const DATA=__DATA__,DAILY=__DAILY__,_dc={};
function dailyOf(s){if(!(s in _dc)){try{const b=atob(DAILY[s]);const u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);_dc[s]=JSON.parse(pako.inflate(u,{to:'string'}));}catch(e){_dc[s]=null;}}return _dc[s];}
const names=Object.keys(DATA).sort((a,b)=>DATA[b].lat-DATA[a].lat);
const PAL=['#2563eb','#dc2626','#059669','#d97706','#7c3aed','#0891b2','#db2777','#65a30d','#9333ea','#c2410c'];
const col={};names.forEach((s,i)=>col[s]=PAL[i%PAL.length]);
const map=L.map('map').setView([48.4,-124.2],7);
L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png',{attribution:'&copy; OSM, &copy; CARTO',subdomains:'abcd',maxZoom:19}).addTo(map);
let cur={s:null,fi:0,yr:[-0.7,0.7]};   // fi: 0=all-families overview, 1..nfam=that family; yr locked at station-select
function dson(){return document.getElementById('deseason').checked;}
function wkey(){return 's24'+(dson()?'d':'');}
function tagTxt(){return dson()?' (deseasoned)':'';}
function blk(s){const DD=dailyOf(s);return DD?DD[wkey()]:null;}
function rawblk(s){const DD=dailyOf(s);return DD?DD['s24']:null;}
const Y2={title:'tremor/mo',overlaying:'y',side:'right',showgrid:false,rangemode:'tozero',titlefont:{color:'#f59e0b',size:10},tickfont:{color:'#f59e0b',size:9}};
const LAY={paper_bgcolor:'#fff',plot_bgcolor:'#fff',font:{color:'#334155'}};
// symmetric robust y-range from the RAW family cloud + median; locked at station-select so deseason & the
// family slider never rescale (compare like-for-like). 1/99 pct trims wild single-family excursions.
function computeYR(D){const v=[];for(const id in D.byfam){const a=D.byfam[id];for(let i=0;i<a.length;i++)if(a[i]!=null)v.push(a[i]);}for(let i=0;i<D.med.length;i++)if(D.med[i]!=null)v.push(D.med[i]);if(!v.length)return[-0.7,0.7];v.sort((a,b)=>a-b);const q=p=>v[Math.min(v.length-1,Math.max(0,Math.round(p*(v.length-1))))];let m=Math.max(Math.abs(q(0.01)),Math.abs(q(0.99)));m=Math.max(m*1.1,0.05);return[-m,m];}
// x-range = the actual span where the shown series exists (0.5% stray-trim). No empty 2010-2026 padding.
function xrFromVals(dates,vals){const idx=[];for(let i=0;i<vals.length;i++)if(vals[i]!=null)idx.push(i);if(!idx.length)return[dates[0],dates[dates.length-1]];const lo=idx[Math.floor(0.005*(idx.length-1))],hi=idx[Math.ceil(0.995*(idx.length-1))];return[dates[lo],dates[hi]];}
function selectStation(s){const R=rawblk(s);cur={s:s,fi:0,yr:R?computeYR(R):[-0.7,0.7]};const sl2=document.getElementById('famslider');document.getElementById('sliderwrap').style.display='';sl2.max=DATA[s].nfam;sl2.value=0;render();}
function render(){const s=cur.s;if(!s)return;setChip(s);const D=blk(s);const slab=document.getElementById('slabel');if(!D){Plotly.purge('plot');return;}
 if(cur.fi===0){hiStation(s);
  const tr=[];const T=DATA[s].tremor;if(T)tr.push({x:T.dates,y:T.total,type:'bar',yaxis:'y2',marker:{color:'rgba(245,158,11,0.4)'},hoverinfo:'x+y',name:'tremor'});
  let X=[],Y=[];for(const id in D.byfam){const a=D.byfam[id];for(let i=0;i<a.length;i++){X.push(D.dates[i]);Y.push(a[i]);}X.push(null);Y.push(null);}
  tr.push({x:X,y:Y,type:'scattergl',mode:'lines',line:{color:'rgba(148,163,184,0.28)',width:0.6},hoverinfo:'skip'});
  tr.push({x:D.dates,y:D.med,type:'scattergl',mode:'lines',line:{color:col[s],width:2.4},name:'median'});
  Plotly.react('plot',tr,Object.assign({margin:{l:52,r:46,t:34,b:34},showlegend:false,title:{text:s+' &mdash; '+DATA[s].nfam+' certified families + median'+tagTxt(),font:{size:13}},yaxis:{title:'dv/v (%)',zeroline:true,zerolinecolor:'#e33',range:cur.yr.slice()},yaxis2:Y2,xaxis:{title:'date',range:xrFromVals(D.dates,D.med)}},LAY),{responsive:true});
  if(slab)slab.textContent='all families ('+DATA[s].nfam+')';
  document.getElementById('info').innerHTML='<b>'+s+'</b> '+DATA[s].lat.toFixed(2)+'&deg;N &mdash; '+DATA[s].nfam+' certified families (thin) + per-day cross-family median (bold). Orange = tremor/mo within 20 km. Slide above to step through families.';
 }else{const f=DATA[s].fam[cur.fi-1];if(!f)return;const id=f.id,yv=D.byfam[id]||[];hiPatch(s,id);
  const ft=[];
  ft.push({x:D.dates,y:D.med,type:'scattergl',mode:'lines',line:{color:'#000',width:1.3,dash:'dot'},name:'median'});
  ft.push({x:D.dates,y:yv,type:'scattergl',mode:'lines',line:{color:col[s],width:2},name:'family'});
  Plotly.react('plot',ft,Object.assign({margin:{l:52,r:46,t:34,b:34},showlegend:false,title:{text:s+' &mdash; family '+cur.fi+'/'+DATA[s].nfam+tagTxt(),font:{size:12}},yaxis:{title:'dv/v (%)',zeroline:true,zerolinecolor:'#e33',range:cur.yr.slice()},xaxis:{title:'date',range:xrFromVals(D.dates,yv)}},LAY),{responsive:true});
  if(slab)slab.textContent='family '+cur.fi+' / '+DATA[s].nfam;
  document.getElementById('info').innerHTML='<b>'+s+'</b> &mdash; family '+cur.fi+'/'+DATA[s].nfam+' @ '+f.lat+'&deg;N '+f.lon+'&deg; (bold, station color) vs station median (black dotted).';}}
function showStation(s){selectStation(s);}
function showFam(s,id){if(cur.s!==s)selectStation(s);const i=DATA[s].fam.findIndex(x=>x.id===id);if(i<0)return;cur.fi=i+1;document.getElementById('famslider').value=cur.fi;render();}
const layers={},chips={};
names.forEach(s=>{const d=DATA[s],c=col[s],grp=L.layerGroup(),fams={};
 d.fam.forEach(f=>{const line=L.polyline([[d.lat,d.lon],[f.lat,f.lon]],{color:c,weight:1,opacity:0.2}).on('click',e=>{L.DomEvent.stop(e);showFam(s,f.id);}).addTo(grp);
  const dot=L.circleMarker([f.lat,f.lon],{radius:2.5,color:c,fillColor:c,weight:0,fillOpacity:0.6}).on('click',e=>{L.DomEvent.stop(e);showFam(s,f.id);}).bindTooltip(s+' '+f.id,{direction:'top'}).addTo(grp);
  fams[f.id]={line,dot};});grp.addTo(map);
 const marker=L.circleMarker([d.lat,d.lon],{radius:7,color:'#1e293b',weight:1.6,fillColor:c,fillOpacity:0.95}).on('click',e=>{L.DomEvent.stop(e);showStation(s);}).bindTooltip(s+' ('+d.nfam+' fam)',{direction:'top',permanent:false,className:'',sticky:true}).addTo(map);
 layers[s]={marker,c,fams};});
function dimAll(){for(const k in layers){const o=layers[k];o.marker.setStyle({fillOpacity:0.45});for(const id in o.fams){o.fams[id].line.setStyle({opacity:0.06});o.fams[id].dot.setStyle({radius:2,fillOpacity:0.15});}}}
function hiStation(s){dimAll();const o=layers[s];if(!o)return;o.marker.setStyle({fillOpacity:1,radius:9});for(const id in o.fams){o.fams[id].line.setStyle({opacity:0.6});o.fams[id].dot.setStyle({radius:4,fillOpacity:0.9});}}
function hiPatch(s,id){dimAll();const o=layers[s];if(!o)return;o.marker.setStyle({fillOpacity:1,radius:9});if(o.fams[id]){o.fams[id].line.setStyle({opacity:0.95,weight:3});o.fams[id].dot.setStyle({radius:6,color:'#fff',weight:1.5,fillOpacity:1});}}
const sl=document.getElementById('stalist');
names.forEach(s=>{const b=document.createElement('span');b.textContent=s+' ('+DATA[s].nfam+')';b.className='chip';b.onclick=()=>showStation(s);sl.appendChild(b);chips[s]=b;});
function setChip(s){for(const k in chips)chips[k].classList.toggle('active',k===s);}
document.getElementById('deseason').addEventListener('change',render);
document.getElementById('famslider').addEventListener('input',function(){cur.fi=+this.value;render();});
(function(){const dr=document.getElementById('drag'),sd=document.getElementById('side');let on=false;
 dr.addEventListener('mousedown',e=>{on=true;dr.classList.add('on');document.body.style.cursor='col-resize';e.preventDefault();});
 window.addEventListener('mousemove',e=>{if(!on)return;const w=Math.min(window.innerWidth-220,Math.max(340,window.innerWidth-e.clientX));sd.style.width=w+'px';});
 window.addEventListener('mouseup',()=>{if(!on)return;on=false;dr.classList.remove('on');document.body.style.cursor='';map.invalidateSize();Plotly.Plots.resize('plot');});})();
window.addEventListener('resize',()=>{map.invalidateSize();Plotly.Plots.resize('plot');});
if(names.length)showStation(names[0]);
</script></body></html>'''
out = HTML.replace('__DATA__', djs).replace('__DAILY__', djson)
os.makedirs('figures', exist_ok=True)
with open('figures/borehole_dvv_map.html', 'w') as f:
    f.write(out)
print(f"wrote figures/borehole_dvv_map.html ({len(out)/1e6:.2f} MB), {len(DATA)} stations")
