#!/usr/bin/env python
"""Precompute the 'correlation-function' stack image for every certified family, for the interactive map.

The map builder only reads dv/v CSVs; these images need the daily-stack waveforms (283 GB of npz), so
they are built ONCE into a cache and the map just embeds it.

Per (station, certified family): rebuild exactly what dvv_roll30cal.py measures -- 30-calendar-day
trailing stacks (>=5 daily stacks), peak-normalised -- then
  * regrid rows onto a real daily calendar so OUTAGES stay visible (rendered light grey),
  * bin to <= ROWS rows, crop to lapse LO..HI s,
  * colour-scale on the 2-4 s CODA (the direct S is ~16x larger and would saturate everything),
  * encode RGB WebP q75 (~1-3 KB/family; PNG is ~10x larger for the same picture).
Writes data/stack_img/<TAG>.json.gz = {fam_id: {"img": <b64 webp>, "t0","t1", "lo","hi", "n","ngap"}}
"""
import os, sys, io, json, gzip, base64, glob, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from PIL import Image
import matplotlib; matplotlib.use("Agg"); import matplotlib.cm as cm

ROWS, LO, HI, Q = 200, 1.0, 5.0, 80
# Gaps must be UNMISTAKEABLE. RdBu_r is near-white (~247) at zero, so a light grey is invisible;
# use a saturated colour that cannot occur in the colormap. Row binning (~30 d/row) would also hide
# short gaps entirely, so each row additionally carries a COVERAGE gutter on the left.
GAP_RGB = (242, 242, 242)      # same flat grey the matplotlib stack figures use for 'no data'
# Embedding 12,436 families forces row-binning (daily resolution would be ~1.8 GB fleet-wide), and a
# binned row containing even ONE day of data would otherwise render as full-strength colour. So each
# row is blended toward GAP_RGB in proportion to how much of it is missing: complete rows keep full
# colour, half-empty rows wash out, fully-empty rows are flat grey exactly as in the daily figures.
OUT = "data/stack_img"; os.makedirs(OUT, exist_ok=True)


def cert_fams(stem):
    f1 = f"data/{stem}_causality_cert.csv"
    if os.path.exists(f1):
        return set(pd.read_csv(f1).query("reliable").fam)
    f2 = f"data/{stem}_fwd_vs_rev_coda.csv"
    if os.path.exists(f2):
        return set(pd.read_csv(f2).query("ratio>1.5").fam)
    return set()


def station(TAG, stem, force=False):
    dst = f"{OUT}/{TAG}.json.gz"
    if os.path.exists(dst) and not force:
        return None
    npz = f"data/long_window_daily_{TAG}_Z.npz"
    fams = cert_fams(stem)
    if not (os.path.exists(npz) and fams):
        return None
    d = np.load(npz, allow_pickle=True)
    t = d["t"]; pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]); St = d["stacks"]
    lap = (t >= LO) & (t < HI)
    cod = (t >= 2.0) & (t < 4.0)
    out = {}
    for fam in fams:
        m = pat == fam
        if m.sum() < 60:
            continue
        S = St[m].astype(np.float64); ds = dates[m]
        o = np.argsort(ds.values); S, ds = S[o], ds[o]
        Dn = ds.values.astype("datetime64[D]"); n = len(S); idx = np.arange(n)
        cs = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(S, axis=0)])
        los = np.searchsorted(Dn, Dn - np.timedelta64(30, "D"), side="right")
        c = idx + 1 - los; keep = c >= 5
        if keep.sum() < 30:
            continue
        R = (cs[idx + 1] - cs[los]) / c[:, None]
        R, Dr = R[keep], ds[keep]
        pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0
        Rn = R / pk
        v = float(np.percentile(np.abs(Rn[:, cod]), 99))     # scale on the CODA, not the direct S
        if not np.isfinite(v) or v <= 0:
            continue
        full = pd.date_range(Dr.min(), Dr.max(), freq="D")
        G = np.full((len(full), int(lap.sum())), np.nan)
        G[full.get_indexer(Dr)] = Rn[:, lap]
        k = max(1, len(full) // ROWS); nb = len(full) // k
        H = G[:nb * k].reshape(nb, k, -1)
        B = np.nanmean(H, axis=1)
        cov = 1.0 - np.isnan(H[:, :, 0]).mean(axis=1)        # fraction of this row's days that HAVE a stack
        Qq = np.clip(B / v, -1, 1)
        rgb = (cm.RdBu_r((Qq + 1) / 2)[..., :3] * 255).astype(float)
        rgb[np.isnan(Qq)] = GAP_RGB
        a = cov[:, None, None]                               # blend toward grey by the missing fraction
        rgb = (rgb * a + np.asarray(GAP_RGB, float) * (1 - a)).astype(np.uint8)
        # TRANSPOSE to (lapse, date) so the image plots with DATE on x -- same axis as the dv/v curve
        # above it -- and lapse on y. Row 0 must be the LARGEST lapse (Plotly anchors a layout image
        # at its top-left), hence the vertical flip.
        rgb = rgb.transpose(1, 0, 2)[::-1]
        buf = io.BytesIO(); Image.fromarray(rgb).save(buf, "WEBP", quality=Q)
        out[fam] = {"img": base64.b64encode(buf.getvalue()).decode(),
                    "t0": f"{Dr.min():%Y-%m-%d}", "t1": f"{Dr.max():%Y-%m-%d}",
                    "lo": LO, "hi": HI, "n": int(keep.sum()),
                    "ngap": int(np.isnan(B[:, 0]).sum()),
                    "cov": round(float(cov.mean()), 3)}
    del d, St
    if not out:
        return None
    blob = gzip.compress(json.dumps(out, separators=(",", ":")).encode(), 6)
    open(dst, "wb").write(blob)
    return len(out), len(blob)


if __name__ == "__main__":
    args = sys.argv[1:]
    force = "--force" in args
    only = [a for a in args if not a.startswith("--")]
    tags = []
    for f in sorted(glob.glob("data/long_window_daily_*_Z.npz")):
        tag = os.path.basename(f).replace("long_window_daily_", "").replace("_Z.npz", "")
        if "rev" in tag.lower() or "MIRROR" in tag:
            continue
        tags.append((tag, tag.lower().replace("p90f40", "")))
    if only:
        tags = [x for x in tags if x[0] in only]
    tot = 0
    for tag, stem in tags:
        try:
            r = station(tag, stem, force)
        except Exception as e:
            print(f"  {tag}: FAILED {type(e).__name__}: {e}", flush=True); continue
        if r:
            tot += r[1]
            print(f"  {tag}: {r[0]} families -> {r[1]/1024:.0f} KB   (cache total {tot/2**20:.1f} MB)", flush=True)
    print(f"DONE  cache {tot/2**20:.1f} MB in {OUT}/")
