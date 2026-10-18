"""Digitize Royer & Bostock (2014, EPSL 402:247-256) LFE template locations
from the vector graphics of the paper and its supplement.

  Fig. 4  (main p.6)  Hyp2000 epicentres, panels A (overview) B (NVI) C (SVI) D (NW+SW)
  Fig. S1 (supp p.3)  hypoDD  epicentres, same layout
  Fig. S2 (supp p.4)  Hyp2000 depth sections H,I,J,A,B (left col) C,D,E,F,G (right col)
  Fig. 5  (main p.7)  hypoDD  depth sections, same layout

Outputs (this folder):
  royer2014_templates_digitized.csv
  royer2014_georef_ticks.csv       per-tick residuals of every map panel fit
  royer2014_profiles.csv           profile endpoints (lon/lat) + fitted distance offset
  royer2014_section_points.csv     every section diamond (profile, d_km, z_km, matched row)
  qc_*.png                         overlays on rasterised figures
Run:  python digitize_royer2014.py
"""
import os
import sys
import numpy as np
import pandas as pd
import pymupdf
from pyproj import Geod
from scipy.optimize import linear_sum_assignment

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import pdfvec  # noqa: E402
import georef  # noqa: E402
import panels_cfg as pc  # noqa: E402

NOTES = os.path.abspath(os.path.join(HERE, "../../../../notes"))
MAIN = os.path.join(NOTES, "1-s2.0-S0012821X13004780-main.pdf")
SUPP = os.path.join(NOTES, "royer2014_supp.pdf")
GEOD = Geod(ellps="WGS84")

MAPS = {  # fig -> (pdf, page index, locator, colour table)
    "Fig4": (MAIN, 5, "Hyp2000", {(0.93, 0.91, 0.25): "SVI", (0.8, 0.21, 0.19): "NW",
                                   (0.51, 0.73, 0.31): "SW", (0.89, 0.63, 0.22): "NVI"}),
    "FigS1": (SUPP, 2, "hypoDD", {(0.93, 0.91, 0.13): "SVI", (0.93, 0.19, 0.18): "NW",
                                  (0.38, 0.74, 0.27): "SW", (0.98, 0.64, 0.17): "NVI"}),
}
SECTIONS = {  # locator -> (pdf, page, fig name, colour table, map fig providing bars)
    "Hyp2000": (SUPP, 3, "FigS2", MAPS["FigS1"][3], "Fig4"),
    "hypoDD": (MAIN, 6, "Fig5", MAPS["Fig4"][3], "FigS1"),
}
PRIMARY_PANEL = {"NVI": "B", "SVI": "C", "NW": "D", "SW": "D"}
EXPECTED = {"Hyp2000": {"NW": 122, "SW": 54, "NVI": 47},
            "hypoDD": {"NW": round(0.77 * 122), "SW": round(0.72 * 54), "NVI": round(0.91 * 47)}}
GEOREF_KIND = "aea"          # best tick fit (free standard parallels); see README
ALT_KINDS = [("tmerc", False), ("lcc", False), ("eqdc", True)]  # model-spread check
PROFILE_ORDER_LEFT = ["H", "I", "J", "A", "B"]
PROFILE_ORDER_RIGHT = ["C", "D", "E", "F", "G"]
# map bars in panel A, identified from the text labels X / X' of Fig. S1 (sorted
# by the SW end's page y, then x); see identify_bars()
SWATH_KM = 25.0


def sorted_panels(page):
    return sorted(pdfvec.panels(page), key=lambda b: (round(b[1] / 50), b[0]))


def fit_panel(page, box, name, kind=GEOREF_KIND, free_std=True):
    cfg = pc.PANELS[name]
    tk = pdfvec.ticks(page, box)
    T = pc.build_ticks(cfg, box, tk)
    g = georef.PanelGeoref(kind, cfg["lonlim"], cfg["latlim"], T, free_std=free_std).fit(box)
    return g, T


def identify_bars(page, boxA):
    """White profile bars in panel A -> centre-line end points (page coords).
    Returns list of (xSW, ySW, xNE, yNE) sorted as in PROFILE_ORDER via label logic."""
    bars = []
    for g in page.get_drawings():
        r = g["rect"]
        if not (r.x0 >= boxA[0] - 1 and r.x1 <= boxA[2] + 1 and r.y0 >= boxA[1] - 1 and r.y1 <= boxA[3] + 1):
            continue
        c = pdfvec.rgb(g)
        if c is None or min(c) < 0.97 or max(r.width, r.height) < 20 or len(g["items"]) > 6:
            continue
        sp = pdfvec.subpaths(g)
        a = np.array(sp[0][:4])
        # short sides: (0,3)->? find the two shortest edges
        e = [(i, (i + 1) % 4, np.hypot(*(a[(i + 1) % 4] - a[i]))) for i in range(4)]
        e.sort(key=lambda t: t[2])
        m1 = (a[e[0][0]] + a[e[0][1]]) / 2
        m2 = (a[e[1][0]] + a[e[1][1]]) / 2
        sw, ne = (m1, m2) if m1[0] < m2[0] else (m2, m1)
        if abs(ne[0] - sw[0]) > 0.9 * (boxA[2] - boxA[0]):
            continue  # full-width white path = not a profile bar
        bars.append((sw[0], sw[1], ne[0], ne[1]))
    return bars


def label_bars(bars):
    """Assign profile letters. Geometry (both figures): H,I,J are the three
    parallel NVI bars (SW end x<200 in S1 / <180 in Fig4 frame units relative to
    panel), ordered N->S; A (short, steep), B (short, less steep); C,D,E parallel
    NW bars ordered N->S; F ~horizontal; G dips to the SE (NE end lower)."""
    bars = [np.array(b) for b in bars]
    out = {}
    assert len(bars) == 10, len(bars)
    xs = np.array([b[0] for b in bars])
    nvi = [b for b in bars if b[0] < xs.min() + 10]
    rest = [b for b in bars if b[0] >= xs.min() + 10]
    assert len(nvi) == 3
    nvi.sort(key=lambda b: b[1])
    for k, b in zip("HIJ", nvi):
        out[k] = b
    slope = lambda b: (b[3] - b[1]) / (b[2] - b[0])  # page y down => negative = rises to NE
    length = lambda b: np.hypot(b[2] - b[0], b[3] - b[1])
    short = sorted(rest, key=length)[:2]
    short.sort(key=slope)  # steeper (more negative) first -> A
    out["A"], out["B"] = short
    long_ = [b for b in rest if not any(b is s for s in short)]
    g = [b for b in long_ if slope(b) > 0.05]
    f = [b for b in long_ if abs(slope(b)) <= 0.05]
    cde = sorted([b for b in long_ if slope(b) < -0.05], key=lambda b: b[1])
    out["G"], out["F"] = g[0], f[0]
    for k, b in zip("CDE", cde):
        out[k] = b
    return out


def extract_maps():
    rows, tickrows, georefs = [], [], {}
    for fig, (pdf, pg, loc, cols) in MAPS.items():
        page = pymupdf.open(pdf)[pg]
        boxes = dict(zip("ABCD", sorted_panels(page)))
        for name, box in boxes.items():
            g, T = fit_panel(page, box, name)
            alts = [fit_panel(page, box, name, k, fs)[0] for k, fs in ALT_KINDS]
            georefs[(fig, name)] = (g, alts, box)
            for (typ, val, edge, tc), r in zip(T, g.res_pt):
                tickrows.append(dict(fig=fig, panel=name, tick_type=typ, value=val, page_edge=edge,
                                     page_tick=tc, resid_pt=r, resid_km=r * g.km_per_pt))
        q = pdfvec.filled_quads(page, set(cols))
        for d in q:
            inside = [n for n, b in boxes.items()
                      if b[0] - 3 <= d["cx"] <= b[2] + 3 and b[1] - 3 <= d["cy"] <= b[3] + 3]
            if not inside:
                continue
            pan = inside[0]
            g, alts, box = georefs[(fig, pan)]
            lon, lat = g.inv(g.p, d["cx"], d["cy"])
            spread = 0.0
            for a in alts:
                lo2, la2 = a.inv(a.p, d["cx"], d["cy"])
                spread = max(spread, GEOD.inv(lon, lat, lo2, la2)[2] / 1000)
            rms_km = float(np.sqrt(np.mean(g.res_pt ** 2)) * g.km_per_pt)
            rows.append(dict(fig=fig, locator=loc, panel=pan, subarray=cols[d["colour"]],
                             page_x=d["cx"], page_y=d["cy"], sym_w=d["w"], sym_h=d["h"],
                             lon=float(lon), lat=float(lat), tick_rms_km=rms_km,
                             model_spread_km=spread, km_per_pt=g.km_per_pt,
                             symbol_style=f"diamond fill={d['colour']} {d['w']:.1f}x{d['h']:.1f}pt"))
    return pd.DataFrame(rows), pd.DataFrame(tickrows), georefs


def crosscheck_overview(df):
    """Panel-A symbols vs zoom-panel symbols of the same fig/subarray (1-1 Hungarian)."""
    out = []
    for (fig, sub), grp in df.groupby(["fig", "subarray"]):
        a = grp[grp.panel == "A"]
        z = grp[grp.panel == PRIMARY_PANEL[sub]]
        if len(a) == 0 or len(z) == 0:
            continue
        D = np.array([[GEOD.inv(r1.lon, r1.lat, r2.lon, r2.lat)[2] / 1000 for r2 in z.itertuples()]
                      for r1 in a.itertuples()])
        i, j = linear_sum_assignment(D)
        d = D[i, j]
        out.append(dict(fig=fig, subarray=sub, nA=len(a), nZoom=len(z), med_km=np.median(d),
                        p90_km=np.percentile(d, 90), max_km=d.max(),
                        A_kmpt=a.km_per_pt.iloc[0], zoom_kmpt=z.km_per_pt.iloc[0]))
    return pd.DataFrame(out)


def extract_sections(loc):
    pdf, pg, fig, cols, _ = SECTIONS[loc]
    page = pymupdf.open(pdf)[pg]
    dr = page.get_drawings()
    boxes = sorted({tuple(round(v, 2) for v in g["rect"]) for g in dr
                    if pdfvec.rgb(g) == (1.0, 1.0, 1.0) and g["rect"].width > 150},
                   key=lambda b: (b[0] > page.rect.width / 2, b[1]))
    assert len(boxes) == 10, boxes
    names = PROFILE_ORDER_LEFT + PROFILE_ORDER_RIGHT
    pts = []
    for g in dr:
        c = pdfvec.rgb(g)
        if c not in cols:
            continue
        r = g["rect"]
        if not (2.3 < r.width < 4.5 and 3.0 < r.height < 5.0):
            continue
        cx, cy = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        for nm, b in zip(names, boxes):
            if b[0] - 2 <= cx <= b[2] + 2 and b[1] - 2 <= cy <= b[3] + 2:
                d = -100 + 300 * (cx - b[0]) / (b[2] - b[0])
                z = 80 * (cy - b[1]) / (b[3] - b[1])
                pts.append(dict(locator=loc, fig=fig, profile=nm, subarray=cols[c], d_km=d, z_km=z,
                                page_x=cx, page_y=cy, kmpt=300 / (b[2] - b[0])))
    return pd.DataFrame(pts), dict(zip(names, boxes)), page


def along_cross(lon, lat, p0, p1):
    """Along-track / cross-track km of points relative to great circle p0->p1."""
    az01, _, L = GEOD.inv(p0[0], p0[1], p1[0], p1[1])
    az0p, _, dist = GEOD.inv(np.full(len(lon), p0[0]), np.full(len(lon), p0[1]), lon, lat)
    th = np.radians(az0p - az01)
    return dist * np.cos(th) / 1000, dist * np.sin(th) / 1000, L / 1000


def join_depths(maps, loc, sec):
    """For each profile: project map epicentres of the matching locator onto
    the bar, fit the along-profile offset c (d_sec = d_map + c) by 1-D Hungarian
    matching, and attach section depths to map symbols."""
    barfig = SECTIONS[loc][4]
    prim = maps[(maps.locator == loc) & (maps.panel == maps.subarray.map(PRIMARY_PANEL))].copy()
    prof_rows, match_rows = [], []
    for nm, (x0, y0, x1, y1, lo0, la0, lo1, la1) in BARS[barfig].items():
        s = sec[sec.profile == nm]
        for sub, ss in s.groupby("subarray"):
            m = prim[prim.subarray == sub]
            if len(m) == 0:
                continue
            al, cr, L = along_cross(m.lon.values, m.lat.values, (lo0, la0), (lo1, la1))
            best = None
            for c in np.arange(-200, 200.01, 0.25):
                D = np.abs(ss.d_km.values[:, None] - (al[None, :] + c))
                D = np.where(np.abs(cr)[None, :] <= SWATH_KM + 5, D, 1e3)
                i, j = linear_sum_assignment(D)
                cost = np.sum(np.minimum(D[i, j], 3.0))
                if best is None or cost < best[0]:
                    best = (cost, c, i, j, D[i, j])
            cost, c, i, j, dd = best
            n_in = int(np.sum(np.abs(cr) <= SWATH_KM))
            prof_rows.append(dict(locator=loc, profile=nm, subarray=sub, n_section=len(ss),
                                  n_map_in_swath25=n_in, offset_c_km=c,
                                  median_abs_dd_km=float(np.median(dd)), n_dd_gt1km=int(np.sum(dd > 1)),
                                  bar_len_km=L))
            for ii, jj, e in zip(i, j, dd):
                srow = ss.iloc[ii]
                # ambiguity: 2nd-nearest candidate in along distance
                dall = np.abs(srow.d_km - (al + c))
                dall = np.where(np.abs(cr) <= SWATH_KM + 5, dall, np.inf)
                srt = np.sort(dall)
                match_rows.append(dict(locator=loc, profile=nm, subarray=sub, map_index=m.index[jj],
                                       d_sec=srow.d_km, z_km=srow.z_km, d_map=al[jj] + c,
                                       cross_km=cr[jj], dd_km=e,
                                       gap2_km=float(srt[1]) if len(srt) > 1 else np.nan))
    return pd.DataFrame(prof_rows), pd.DataFrame(match_rows)


BARS = {}


def thesis_crosscheck(prim):
    """Thesis Fig 4.1 families vs EPSL Hyp2000 / hypoDD symbols, per subarray:
    nearest-neighbour distance from each thesis family to the EPSL set, and the
    Hungarian one-to-one matched distances."""
    out = []
    th = prim[prim.source_fig == "Thesis_Fig4.1"]
    for loc in ["Hyp2000", "hypoDD"]:
        for sub in ["SVI", "NW", "SW"]:
            a = th[th.subarray == sub]
            b = prim[(prim.locator == loc) & (prim.subarray == sub)]
            D = GEOD.inv(np.repeat(a.lon.values, len(b)), np.repeat(a.lat.values, len(b)),
                         np.tile(b.lon.values, len(a)), np.tile(b.lat.values, len(a)))[2].reshape(len(a), len(b)) / 1000
            nn = D.min(1)
            i, j = linear_sum_assignment(D)
            out.append(dict(epsl_set=loc, subarray=sub, n_thesis=len(a), n_epsl=len(b),
                            nn_med_km=np.median(nn), nn_p90_km=np.percentile(nn, 90),
                            frac_nn_lt0p3km=np.mean(nn < 0.3), frac_nn_lt1km=np.mean(nn < 1.0),
                            hung_med_km=np.median(D[i, j])))
    return pd.DataFrame(out)


def thesis_duplicate_figs():
    """Thesis Fig 3.4 (p.57) and A.1 (p.125) reproduce EPSL Fig 4 / S1: compare
    coloured-diamond counts per colour."""
    doc = pymupdf.open(os.path.join(NOTES, "royer_sources/royer2014_phd_thesis.pdf"))
    res = {}
    for nm, pg in (("Thesis_Fig3.4", 57), ("Thesis_FigA.1", 125)):
        import collections
        c = collections.Counter()
        for g in doc[pg].get_drawings():
            r = g["rect"]
            col = pdfvec.rgb(g)
            if col and len(set(col)) > 1 and len(g["items"]) == 4 and 2.8 < min(r.width, r.height) < 4.4 \
                    and max(r.width, r.height) < 5.8:
                c[col] += 1
        res[nm] = dict(c)
    return res


def main():
    maps, ticks, georefs = extract_maps()
    ticks.to_csv(os.path.join(HERE, "royer2014_georef_ticks.csv"), index=False)
    print(ticks.groupby(["fig", "panel"]).agg(n=("resid_pt", "size"),
                                              rms_pt=("resid_pt", lambda r: np.sqrt(np.mean(r ** 2))),
                                              max_pt=("resid_pt", lambda r: np.abs(r).max()),
                                              rms_km=("resid_km", lambda r: np.sqrt(np.mean(r ** 2))),
                                              max_km=("resid_km", lambda r: np.abs(r).max())).round(3))
    print(maps.groupby(["fig", "panel", "subarray"]).size().unstack(fill_value=0))
    cc = crosscheck_overview(maps)
    print("overview-vs-zoom cross-check (km):\n", cc.round(2))
    cc.to_csv(os.path.join(HERE, "royer2014_overview_crosscheck.csv"), index=False)

    # profile bars from panel A of each map figure
    prof_out = []
    for fig in MAPS:
        page = pymupdf.open(MAPS[fig][0])[MAPS[fig][1]]
        g, _, boxA = georefs[(fig, "A")]
        lab = label_bars(identify_bars(page, boxA))
        BARS[fig] = {}
        for k, b in lab.items():
            lo0, la0 = g.inv(g.p, b[0], b[1])
            lo1, la1 = g.inv(g.p, b[2], b[3])
            BARS[fig][k] = (*b, float(lo0), float(la0), float(lo1), float(la1))
            prof_out.append(dict(fig=fig, profile=k, page_x0=b[0], page_y0=b[1], page_x1=b[2], page_y1=b[3],
                                 lon0=float(lo0), lat0=float(la0), lon1=float(lo1), lat1=float(la1)))
    prof_df = pd.DataFrame(prof_out)

    all_prof, all_match, all_sec = [], [], []
    secboxes = {}
    for loc in SECTIONS:
        sec, sboxes, spage = extract_sections(loc)
        secboxes[loc] = (sboxes, spage)
        pr, mt = join_depths(maps, loc, sec)
        all_prof.append(pr)
        all_match.append(mt)
        all_sec.append(sec)
    pr = pd.concat(all_prof)
    mt = pd.concat(all_match)
    sec = pd.concat(all_sec)
    print(pr.round(2).to_string())
    prof_df.to_csv(os.path.join(HERE, "royer2014_profiles.csv"), index=False)
    pr.to_csv(os.path.join(HERE, "royer2014_profile_fits.csv"), index=False)
    mt.to_csv(os.path.join(HERE, "royer2014_section_matches.csv"), index=False)
    sec.to_csv(os.path.join(HERE, "royer2014_section_points.csv"), index=False)

    # depth per map symbol: median over profiles where matched with |dd|<=1 km
    good = mt[mt.dd_km <= 1.0]
    dep = good.groupby("map_index").agg(depth_km=("z_km", "median"), depth_n=("z_km", "size"),
                                        depth_spread_km=("z_km", lambda z: z.max() - z.min()),
                                        depth_profiles=("profile", lambda p: "".join(sorted(p))),
                                        min_gap2_km=("gap2_km", "min"))
    out = maps.join(dep)
    out["depth_n"] = out["depth_n"].fillna(0).astype(int)
    out["template_id"] = np.nan
    out["source_fig"] = out.fig + out.panel
    out["georef_residual_km"] = np.sqrt(out.tick_rms_km ** 2 + out.model_spread_km ** 2)
    prim = out[out.panel == out.subarray.map(PRIMARY_PANEL)].copy()
    cols = ["subarray", "source_fig", "locator", "lon", "lat", "depth_km", "template_id", "symbol_style",
            "georef_residual_km", "tick_rms_km", "model_spread_km", "km_per_pt", "depth_n",
            "depth_spread_km", "depth_profiles", "min_gap2_km", "page_x", "page_y"]
    prim = prim.sort_values(["locator", "subarray", "lat", "lon"])[cols]
    # thesis Fig 4.1 families (with IDs) appended as extra rows
    import digitize_thesis_fig41 as th
    tdf, tinfo = th.extract()
    print("thesis Fig4.1:", tinfo)
    tdf["depth_n"] = 0
    for c in cols:
        if c not in tdf:
            tdf[c] = np.nan
    prim = pd.concat([prim, tdf.sort_values("template_id")[cols]], ignore_index=True)
    xc = thesis_crosscheck(prim)
    print(xc.round(2).to_string())
    xc.to_csv(os.path.join(HERE, "royer2014_thesis_vs_epsl.csv"), index=False)
    print(thesis_duplicate_figs())
    # join thesis family IDs onto EPSL Hyp2000 templates (one-to-one, <= 0.3 km)
    prim["thesis_match_km"] = np.nan
    th_ = prim[prim.source_fig == "Thesis_Fig4.1"]
    for sub in ["NW", "SW"]:  # SVI thesis families use Bostock 2012 locations: median NN 2.2 km, no join
        a = th_[th_.subarray == sub]
        b = prim[(prim.locator == "Hyp2000") & (prim.subarray == sub)]
        D = GEOD.inv(np.repeat(a.lon.values, len(b)), np.repeat(a.lat.values, len(b)),
                     np.tile(b.lon.values, len(a)), np.tile(b.lat.values, len(a)))[2].reshape(len(a), len(b)) / 1000
        i, j = linear_sum_assignment(D)
        for ii, jj in zip(i, j):
            if D[ii, jj] <= 0.3:
                prim.loc[b.index[jj], "template_id"] = a.template_id.iloc[ii]
                prim.loc[b.index[jj], "thesis_match_km"] = D[ii, jj]
                prim.loc[a.index[ii], "thesis_match_km"] = D[ii, jj]
        print(sub, "Hyp2000 templates given thesis family ID:",
              int(prim.loc[b.index, "template_id"].notna().sum()), "of", len(b))
    prim.to_csv(os.path.join(HERE, "royer2014_templates_digitized.csv"), index=False, float_format="%.5f")
    # also keep overview-panel duplicates for reference
    out[cols + ["panel"]].to_csv(os.path.join(HERE, "royer2014_all_panel_symbols.csv"), index=False,
                                 float_format="%.5f")
    summ = prim.groupby(["locator", "subarray"]).agg(
        n=("lon", "size"), with_depth=("depth_km", lambda z: int(z.notna().sum())),
        with_id=("template_id", lambda z: int(z.notna().sum())),
        zmin=("depth_km", "min"), zmed=("depth_km", "median"), zmax=("depth_km", "max"),
        georef_km_med=("georef_residual_km", "median"), georef_km_max=("georef_residual_km", "max"))
    print(summ.round(2))
    summ.to_csv(os.path.join(HERE, "royer2014_summary.csv"))
    make_qc(maps, georefs, sec, secboxes, prim)


def make_qc(maps, georefs, sec, secboxes, prim):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image
    import io
    COLS = {"NW": "red", "SW": "lime", "NVI": "orange", "SVI": "yellow"}
    dpi = 400
    for fig, (pdf, pg, loc, _) in MAPS.items():
        page = pymupdf.open(pdf)[pg]
        f, axs = plt.subplots(2, 2, figsize=(16, 13))
        for ax, name in zip(axs.flat, "ABCD"):
            g, _, box = georefs[(fig, name)]
            clip = pymupdf.Rect(box[0] - 14, box[1] - 4, box[2] + 4, box[3] + 12)
            pix = page.get_pixmap(dpi=dpi, clip=clip)
            img = Image.open(io.BytesIO(pix.tobytes("png")))
            ext = [clip.x0, clip.x1, clip.y1, clip.y0]
            ax.imshow(img, extent=ext)
            # fitted graticule
            cfg = pc.PANELS[name]
            lons = sorted(set(cfg["top"] + cfg["bottom"]))
            lats = sorted(set(cfg["left"]))
            la = np.linspace(cfg["latlim"][0] - 1, cfg["latlim"][1] + 1, 200)
            for L in lons:
                x, y = g.fwd(g.p, np.full_like(la, L), la)
                ax.plot(x, y, "c-", lw=0.4)
            lo = np.linspace(cfg["lonlim"][0] - 1, cfg["lonlim"][1] + 1, 200)
            for P in lats:
                x, y = g.fwd(g.p, lo, np.full_like(lo, P))
                ax.plot(x, y, "c-", lw=0.4)
            m = maps[(maps.fig == fig) & (maps.panel == name)]
            for sub, mm in m.groupby("subarray"):
                x, y = g.fwd(g.p, mm.lon.values, mm.lat.values)  # round trip through lon/lat
                ax.plot(x, y, "+", color="k", ms=4, mew=0.6)
                ax.plot(x, y, "o", mfc="none", mec=COLS[sub], ms=6, mew=0.6)
            ax.set_xlim(clip.x0, clip.x1)
            ax.set_ylim(clip.y1, clip.y0)
            ax.set_title(f"{fig} panel {name}: cyan = fitted graticule, +/o = digitized (lon,lat) re-projected")
            ax.set_xticks([])
            ax.set_yticks([])
        f.tight_layout()
        f.savefig(os.path.join(HERE, f"qc_map_{fig}.png"), dpi=130)
        plt.close(f)
    # sections: overlay extracted points + the map-derived projection
    for loc, (boxes, page) in secboxes.items():
        f, axs = plt.subplots(5, 2, figsize=(16, 16))
        names = PROFILE_ORDER_LEFT + PROFILE_ORDER_RIGHT
        for ax, nm in zip(axs.T.flat, names):
            b = boxes[nm]
            clip = pymupdf.Rect(b[0] - 2, b[1] - 2, b[2] + 2, b[3] + 2)
            pix = page.get_pixmap(dpi=dpi, clip=clip)
            img = Image.open(io.BytesIO(pix.tobytes("png")))
            # show in data coordinates
            kx = 300 / (b[2] - b[0])
            kz = 80 / (b[3] - b[1])
            ax.imshow(img, extent=[-100 - 2 * kx, 200 + 2 * kx, 80 + 2 * kz, -2 * kz], aspect="auto")
            s = sec[(sec.locator == loc) & (sec.profile == nm)]
            ax.plot(s.d_km, s.z_km, "k+", ms=5, mew=0.6)
            ax.set_title(f"{SECTIONS[loc][2]} {nm}-{nm}'  ({loc}): + = digitized diamond centres", fontsize=9)
            ax.set_xlim(-100, 200)
            ax.set_ylim(80, 0)
        f.tight_layout()
        f.savefig(os.path.join(HERE, f"qc_sections_{SECTIONS[loc][2]}.png"), dpi=110)
        plt.close(f)
    # final truth table map
    f, axs = plt.subplots(1, 3, figsize=(22, 8))
    for ax, loc in zip(axs, ["Hyp2000", "hypoDD", "thesis_family"]):
        p = prim[prim.locator == loc]
        if loc == "thesis_family":
            sc = ax.scatter(p.lon, p.lat, c=p.template_id, cmap="gist_rainbow_r", s=25,
                            edgecolors=[COLS[s] for s in p.subarray], linewidths=0.8)
            for r in p.itertuples():
                ax.text(r.lon + 0.01, r.lat, str(int(r.template_id)), fontsize=4)
            ax.set_xlim(-125.1, -121.9)
            ax.set_ylim(46.5, 49.2)
            ax.set_aspect(1 / np.cos(np.radians(48)))
            ax.set_title("Thesis Fig 4.1 families (label = decoded family ID)")
            plt.colorbar(sc, ax=ax, label="family ID")
            continue
        sc = ax.scatter(p.lon, p.lat, c=p.depth_km, cmap="viridis_r", vmin=28, vmax=50, s=25,
                        edgecolors=[COLS[s] for s in p.subarray], linewidths=0.8)
        nd = p[p.depth_km.isna()]
        ax.plot(nd.lon, nd.lat, "kx", ms=5, label="no depth")
        for fig_, bars in BARS.items():
            if MAPS[fig_][2] != loc:
                continue
            for k, b in bars.items():
                ax.plot([b[4], b[6]], [b[5], b[7]], "-", color="0.6", lw=0.8)
                ax.text(b[4], b[5], k, fontsize=8)
        ax.set_xlim(-128.2, -121.8)
        ax.set_ylim(46.4, 50.7)
        ax.set_aspect(1 / np.cos(np.radians(48.5)))
        ax.set_title(f"{loc}: digitized templates (edge colour = subarray, fill = depth)")
        ax.legend(loc="lower left")
        plt.colorbar(sc, ax=ax, label="depth km")
    f.tight_layout()
    f.savefig(os.path.join(HERE, "qc_truth_table.png"), dpi=110)
    plt.close(f)


if __name__ == "__main__":
    main()
