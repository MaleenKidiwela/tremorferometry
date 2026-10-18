#!/usr/bin/env python
"""Web (claude.ai artifact) version of figures/borehole_dvv_map.html (2026-10-16).
The local map is a single 420 MB file (all station data inline) that cannot be shared or opened remotely. This
converter re-packs it, without re-running the map build:
- per-station data become separate text files d/<i>.txt, loaded only when a station is clicked (synchronous XHR, so
  the page logic is unchanged). Artifacts serve no raw binary, so each file is base64 of gzip'd JSON; to fit the size
  cap the per-family daily series (byfam, x1000, and nst) are delta-coded integers, decoded in the page (lossless:
  the values are stored at 3 decimals);
- the 3-5 s coda window (s35, s35d) is dropped to stay under the artifact's 256 MB per-version limit (2-4 s kept,
  raw + deseasoned, plus stack images); the toggle then falls back to 2-4 s, as it already does for stations
  without a 3-5 s product;
- the artifact sandbox blocks external images and stylesheets, so map tiles are replaced by an embedded Natural
  Earth 10m coastline + US-Canada border, Leaflet's CSS is inlined, and Plotly loads from jsDelivr.
Input  : figures/borehole_dvv_map.html, a local copy of leaflet@1.9.4 leaflet.css (argv[1]).
Output : figures/borehole_dvv_web/{index.html, d/*.txt, files.json}"""
import os, sys, json, gzip, base64
import shapefile                                   # pyshp
from shapely.geometry import shape, box, mapping

SRC, CSS = 'figures/borehole_dvv_map.html', sys.argv[1]
OUT = 'figures/borehole_dvv_web'
NE = os.path.expanduser('~/.local/share/cartopy/shapefiles/natural_earth')
BBOX = box(-132.0, 37.5, -116.0, 52.5)

s = open(SRC).read()
i = s.find('const DATA='); dec = json.JSONDecoder()
DATA, e = dec.raw_decode(s, i + len('const DATA='))
assert s[e:e + 7] == ',DAILY='
DAILY, e2 = dec.raw_decode(s, e + 7)
head, tail = s[:i], s[e2:]

def delta(arr, scale):
    out, prev = [], 0
    for v in arr:
        if v is None: out.append(None); continue
        q = round(v * scale); out.append(q - prev); prev = q
    return out


os.makedirs(f'{OUT}/d', exist_ok=True)
FILES, tot = {}, 0
for n, (k, v) in enumerate(sorted(DAILY.items())):
    p = json.loads(gzip.decompress(base64.b64decode(v)))
    p.pop('s35', None); p.pop('s35d', None)
    for w in ('s24', 's24d'):
        if w not in p: continue
        p[w]['byfam'] = {f: delta(a, 1000) for f, a in p[w]['byfam'].items()}
        if 'nst' in p[w]: p[w]['nst'] = {f: delta(a, 1) for f, a in p[w]['nst'].items()}
    p['_delta'] = 1
    b = base64.b64encode(gzip.compress(json.dumps(p, separators=(',', ':')).encode(), 9))
    open(f'{OUT}/d/{n}.txt', 'wb').write(b)
    FILES[k] = f'd/{n}.txt'; tot += len(b)
print(f'{len(FILES)} station files, {tot / 1e6:.1f} MB, max {max(os.path.getsize(f"{OUT}/{f}") for f in FILES.values()) / 1e6:.1f} MB')


def lines(path, tol):
    feats = []
    for sr in shapefile.Reader(path).shapeRecords():
        g = shape(sr.shape.__geo_interface__)
        if not g.intersects(BBOX): continue
        g = g.intersection(BBOX).simplify(tol, preserve_topology=False)
        if g.is_empty: continue
        m = mapping(g)
        m = json.loads(json.dumps(m), parse_float=lambda x: round(float(x), 3))
        feats.append({'type': 'Feature', 'properties': {}, 'geometry': m})
    return feats


coast = lines(f'{NE}/physical/ne_10m_coastline.shp', 0.004)
border = lines(f'{NE}/cultural/ne_10m_admin_0_boundary_lines_land.shp', 0.004)
GEO = json.dumps({'coast': {'type': 'FeatureCollection', 'features': coast},
                  'border': {'type': 'FeatureCollection', 'features': border}}, separators=(',', ':'))
print(f'coastline {len(coast)} features, border {len(border)}, {len(GEO) / 1e6:.2f} MB')

css = open(CSS).read()
head = head.replace('<!DOCTYPE html><html><head><meta charset="utf-8"><title>Cascadia borehole LFE-coda dv/v</title>',
                    '<title>Borehole dv/v Map</title>')
head = head.replace('<link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css"/>', f'<style>{css}</style>')
head = head.replace('https://cdn.plot.ly/plotly-2.27.0.min.js', 'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.27.0/plotly.min.js')
head = head.replace('</style></head><body>',
                    '@media (max-width:820px){#wrap{flex-direction:column;height:auto}#map{height:55vh;flex:none}'
                    '#drag{display:none}#side{width:auto!important;min-width:0;border-left:0}#plot{height:auto}}\n</style>')
head = head.replace('toggle for 3&ndash;5 s (same S-anchored estimator; stations without a 3&ndash;5 s product fall back to 2&ndash;4 s)',
                    '3&ndash;5 s is in the full local map only (left out here to fit the web size limit)')
assert '__' not in head[-200:]
tail = tail.replace(',_dc={};', ',_dc={},GEO=' + GEO + ';', 1)
tail = tail.replace(
    "function dailyOf(s){if(!(s in _dc)){try{const b=atob(DAILY[s]);const u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);",
    "function undelta(a,sc){let p=0;return a.map(v=>{if(v===null)return null;p+=v;return sc===1?p:p/sc;});}"
    "function _unpack(P){if(P&&P._delta){for(const w of ['s24','s24d']){const W=P[w];if(!W)continue;"
    "for(const f in W.byfam)W.byfam[f]=undelta(W.byfam[f],1000);if(W.nst)for(const f in W.nst)W.nst[f]=undelta(W.nst[f],1);}}return P;}"
    "function dailyOf(s){if(!(s in _dc)){try{const x=new XMLHttpRequest();x.open('GET',DAILY[s],false);x.send();"
    "if(x.status!==200&&x.status!==0)throw new Error(x.status);"
    "const b=atob(x.responseText.trim());const u=new Uint8Array(b.length);for(let i=0;i<b.length;i++)u[i]=b.charCodeAt(i);")
tail = tail.replace("_dc[s]=JSON.parse(pako.inflate(u,{to:'string'}));", "_dc[s]=_unpack(JSON.parse(pako.inflate(u,{to:'string'})));", 1)
assert 'XMLHttpRequest' in tail and '_unpack(JSON.parse' in tail
TILE = ("L.tileLayer('https://{s}.basemaps.cartocdn.com/rastertiles/voyager/{z}/{x}/{y}{r}.png',"
        "{attribution:'&copy; OpenStreetMap contributors &copy; CARTO',subdomains:'abcd',maxZoom:19}).addTo(map);")
assert TILE in tail
tail = tail.replace(TILE, "L.geoJSON(GEO.coast,{style:{color:'#475569',weight:1.1,fill:false},interactive:false}).addTo(map);"
                          "L.geoJSON(GEO.border,{style:{color:'#94a3b8',weight:1,dashArray:'4 3'},interactive:false}).addTo(map);"
                          "map.attributionControl.addAttribution('Coastline: Natural Earth');")
tail = tail.replace('</body></html>', '')
open(f'{OUT}/index.html', 'w').write(head + 'const DATA=' + json.dumps(DATA, separators=(',', ':')) + ',DAILY='
                                     + json.dumps(FILES) + tail)
json.dump(list(FILES.values()), open(f'{OUT}/files.json', 'w'))
print(f'index.html {os.path.getsize(f"{OUT}/index.html") / 1e6:.2f} MB')
