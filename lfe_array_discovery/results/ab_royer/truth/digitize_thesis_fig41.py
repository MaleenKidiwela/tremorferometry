"""Royer PhD thesis (UBC 2014) Fig. 4.1 (pdf page index 73): 242 LFE *family*
locations (93 SVI + 100 NW + 49 SW regrouped 'geographically independent'
families of the tidal chapter), GMT Mercator map, one filled diamond per family,
colour = family ID on a continuous colour bar (embedded 2244x1 image).

Family ID recovery: each diamond fill colour is located on the colour-bar image
(nearest RGB), converted to a continuous ID via the colour-bar tick labels
(10..240), then ranked -> integer ID 1..242. Validated because (a) the caption
says IDs are ordered by latitude N->S within each catalogue and the decoded IDs
have ZERO latitude inversions inside 1-93 / 94-193 / 194-242; (b) the three
families labelled 49, 123, 242 in the figure (leader lines) decode to those IDs.
"""
import os
import numpy as np
import pandas as pd
import pymupdf
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
THESIS = os.path.abspath(os.path.join(HERE, "../../../../notes/royer_sources/royer2014_phd_thesis.pdf"))
PAGE = 73
CB_XREF, CB_X0, CB_X1 = 1371, 104.82, 481.81
# frame ticks (stroked GMT frame; read from page.get_drawings(), see README)
LON_TICKS = {-125.0: 139.54, -124.5: 198.01, -124.0: 256.64, -123.5: 315.10,
             -123.0: 373.74, -122.5: 432.20, -122.0: 490.83}
LAT_TICKS = {46.5: 563.93, 47.0: 478.76, 47.5: 392.74, 48.0: 306.05, 48.5: 218.36, 49.0: 129.82}


def merc(lat):
    return np.log(np.tan(np.pi / 4 + np.radians(lat) / 2))


def georef():
    lo = np.array(list(LON_TICKS)); xs = np.array(list(LON_TICKS.values()))
    la = np.array(list(LAT_TICKS)); ys = np.array(list(LAT_TICKS.values()))
    A = np.vstack([np.ones_like(lo), lo]).T
    cx, *_ = np.linalg.lstsq(A, xs, rcond=None)
    B = np.vstack([np.ones_like(la), merc(la)]).T
    cy, *_ = np.linalg.lstsq(B, ys, rcond=None)
    rx = A @ cx - xs
    ry = B @ cy - ys
    # equirectangular (linear lat) alternative for comparison
    C = np.vstack([np.ones_like(la), la]).T
    cl, *_ = np.linalg.lstsq(C, ys, rcond=None)
    ry_lin = C @ cl - ys
    return cx, cy, rx, ry, ry_lin


def extract():
    doc = pymupdf.open(THESIS)
    page = doc[PAGE]
    pix = pymupdf.Pixmap(doc, CB_XREF)
    cb = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.width, pix.n)[:, :3].astype(float)
    cx, cy, rx, ry, ry_lin = georef()
    km_x = 111.19 * np.cos(np.radians(47.75)) / cx[1]          # km per pt E-W
    km_y = 111.19 * 0.5 / (LAT_TICKS[47.5] - LAT_TICKS[48.0])   # km per pt N-S (approx)
    rms_km = float(np.sqrt(np.mean(np.concatenate([(rx * km_x) ** 2, (ry * km_y) ** 2]))))
    rows = []
    for g in page.get_drawings():
        if g["type"] != "f" or len(g["items"]) != 3 or abs(g["rect"].width - 6.05) > 0.3:
            continue
        r = g["rect"]
        pts = [(it[1].x, it[1].y) for it in g["items"]] + [(g["items"][-1][2].x, g["items"][-1][2].y)]
        px, py = (r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2
        c = np.array(g["fill"]) * 255
        d = np.linalg.norm(cb - c, axis=1)
        idx = np.where(d <= d.min() + 0.75)[0]
        xcb = CB_X0 + (idx.mean() + 0.5) / pix.width * (CB_X1 - CB_X0)
        id_cont = 10 + (xcb - 118.75) * 230 / (477.1 - 118.75)
        lon = (px - cx[0]) / cx[1]
        m = (py - cy[0]) / cy[1]
        lat = np.degrees(2 * np.arctan(np.exp(m)) - np.pi / 2)
        rows.append(dict(page_x=px, page_y=py, lon=lon, lat=lat, id_cont=id_cont, cb_rgb_dist=d.min(),
                         fill=tuple(round(v, 3) for v in g["fill"])))
    df = pd.DataFrame(rows)
    df["template_id"] = df.id_cont.rank(method="first").astype(int)
    df["subarray"] = np.where(df.template_id <= 93, "SVI", np.where(df.template_id <= 193, "NW", "SW"))
    # checks
    inv = {}
    for sub, g in df.groupby("subarray"):
        yy = g.sort_values("template_id").page_y.values
        inv[sub] = int(np.sum(np.diff(yy) < -0.01))
    # labelled families: leader lines end at the highlighted diamonds
    lab = {}
    for t in page.get_text("dict")["blocks"]:
        for l in t.get("lines", []):
            for s in l["spans"]:
                if s["text"].strip() in ("49", "123", "242") and s["bbox"][1] < 560:
                    lab[int(s["text"])] = s["bbox"]
    leader = {}
    for g in page.get_drawings():
        if g["type"] == "s" and len(g["items"]) == 1 and g["rect"].width > 30:
            p0, p1 = g["items"][0][1], g["items"][0][2]
            for k, bb in lab.items():
                for a, b in ((p0, p1), (p1, p0)):
                    if abs(a.x - (bb[0] + bb[2]) / 2) < 20 and abs(a.y - (bb[1] + bb[3]) / 2) < 12:
                        dd = np.hypot(df.page_x - b.x, df.page_y - b.y)
                        leader[k] = int(df.template_id.iloc[int(dd.argmin())])
    df["source_fig"] = "Thesis_Fig4.1"
    df["locator"] = "thesis_family"
    df["depth_km"] = np.nan
    df["tick_rms_km"] = rms_km
    df["model_spread_km"] = float(np.abs(ry - ry_lin).max() * km_y)  # Mercator vs linear-lat
    df["georef_residual_km"] = rms_km
    df["km_per_pt"] = km_y
    df["symbol_style"] = ["diamond fill=%s 6.0x6.0pt" % (f,) for f in df.fill]
    info = dict(n=len(df), counts=df.subarray.value_counts().to_dict(), lat_inversions=inv,
                labelled_check=leader, tick_rms_pt=float(np.sqrt(np.mean(np.concatenate([rx, ry]) ** 2))),
                tick_rms_km=rms_km, lin_lat_rms_pt=float(np.sqrt(np.mean(ry_lin ** 2))),
                max_cb_rgb_dist=float(df.cb_rgb_dist.max()),
                min_id_gap=float(np.diff(np.sort(df.id_cont)).min()))
    return df, info


if __name__ == "__main__":
    df, info = extract()
    print(info)
