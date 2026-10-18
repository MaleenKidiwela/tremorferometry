#!/usr/bin/env python
"""s11_gates -- S11 step 6, gates (c)-(g) on both data paths (rules: svi_science/results/s11_depth/PLAN.md D7), plus the
step-4 precision check. Headline path from tables/headline.json. eps = -dv/v; values in %.
  (c) scrambled control: D from control cells within 2 sigma of 0 in >= 90 % of measurable episodes (>= 5 control cells)
  (d) split network {PGC site} vs {B926, VGZ, SNB}; families odd vs even: Theil-Sen slope 0.7-1.3 and Pearson r > 0.5
  (e) leave-one-site-out D vs site-only D same sign >= 70 % (mean over sites); |Spearman| < 0.2 of cell dv/v (series-
      demeaned) vs pre-arrival RMS change and vs (1 - cc)
  (f) D(3.5-6 s) vs D(2-3.5 s), each with its own kernels: Theil-Sen slope 0.7-1.3 (kernel-predicted ratio 1.0)
  (g) annual amplitude of D < 2 sigma while S carries the cycle
    python s11_gates.py   -> tables/gates.csv, tables/gates.json, GATES_DONE
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
TAG = os.environ.get("S11_TAG", "")


def Dser(par):
    D = par[par.kind == "D"].set_index("episode")
    return D


def compare(Da, Db, nmin=8):
    j = pd.concat([Da.value.rename("a"), Db.value.rename("b")], axis=1).dropna()
    if len(j) < nmin:
        return dict(n=len(j), slope=np.nan, r=np.nan, sign=np.nan)
    sl = stats.theilslopes(j.b, j.a)[0]
    return dict(n=len(j), slope=float(sl), r=float(np.corrcoef(j.a, j.b)[0, 1]), sign=float(np.mean(np.sign(j.a) == np.sign(j.b))))


def run_path(path, log):
    out = []
    head = pd.read_csv(f"{C.TB}/inv_{path}_head{TAG}.csv")
    cells = pd.read_csv(f"{C.TB}/cells_{path}_head{TAG}.csv")
    # ---------------- (c) scrambled control
    cc_ = I.load_cells(path, C.HEAD, TAG, control=True)
    cc_ = I.prepare(path, C.HEAD, TAG, cells=cc_)
    if len(cc_) >= 10:
        pc, sc = I.invert(cc_, extras=False)
        D = pc[(pc.kind == "D") & (pc.n_cells >= 5)]
        frac = float(np.mean(np.abs(D.value / D.sig) < 2)) if len(D) else np.nan
        pc.to_csv(f"{C.TB}/inv_{path}_control{TAG}.csv", index=False)
        out.append(dict(gate="c", path=path, stat="frac_within_2sig", value=frac, n=len(D), n_cells=len(cc_),
                        passed=bool(frac >= 0.9) if len(D) else None,
                        note="measurable = episodes with >= 5 control cells"))
    else:
        out.append(dict(gate="c", path=path, stat="frac_within_2sig", value=np.nan, n=0, n_cells=len(cc_), passed=None,
                        note="fewer than 10 usable control cells: not measurable"))
    # ---------------- (d) split network and odd/even families
    def sub(mask):
        c = cells[mask].reset_index(drop=True)
        key = c.fam.astype(str) + "|" + c.netsta + "|" + c.epoch
        c = c[key.map(key.value_counts()) >= 3].reset_index(drop=True)
        return Dser(I.invert(c, extras=False)[0]) if len(c) > 20 else None
    for name, ma, mb in (("network", cells.site == "PGC", cells.site != "PGC"),
                         ("families", cells.fam % 2 == 1, cells.fam % 2 == 0)):
        Da, Db = sub(ma), sub(mb)
        r = compare(Da, Db) if Da is not None and Db is not None else dict(n=0, slope=np.nan, r=np.nan, sign=np.nan)
        out.append(dict(gate="d", path=path, stat=f"split_{name}", value=r["slope"], r=r["r"], n=r["n"],
                        passed=bool(0.7 <= r["slope"] <= 1.3 and r["r"] > 0.5) if np.isfinite(r["slope"]) else None))
    # ---------------- (e) per-site sign agreement + cell correlations
    signs = []
    for s in C.SITES:
        Ds, Do = sub(cells.site == s), sub(cells.site != s)
        if Ds is None or Do is None:
            continue
        r = compare(Ds, Do, nmin=5); signs.append(r["sign"])
        out.append(dict(gate="e", path=path, stat=f"site_{s}_vs_rest_sign", value=r["sign"], r=r["r"], n=r["n"], passed=None))
    ms = float(np.nanmean(signs)) if signs else np.nan
    out.append(dict(gate="e", path=path, stat="mean_site_sign_agreement", value=ms, n=len(signs), passed=bool(ms >= 0.7)))
    w = 1 / cells.sig ** 2
    key = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    dm = cells.d - (cells.d * w).groupby(key).transform("sum") / w.groupby(key).transform("sum")
    for nm, x in (("pre_rms_change", np.log(cells.pre_rel)), ("one_minus_cc", 1 - cells.cc)):
        ok = np.isfinite(x) & np.isfinite(dm)
        rho = float(stats.spearmanr(dm[ok], x[ok])[0])
        out.append(dict(gate="e", path=path, stat=f"rho_dvv_vs_{nm}", value=rho, n=int(ok.sum()), passed=bool(abs(rho) < 0.2)))
    # ---------------- (f) lapse
    try:
        a = Dser(pd.read_csv(f"{C.TB}/inv_{path}_w2035{TAG}.csv")); b = Dser(pd.read_csv(f"{C.TB}/inv_{path}_w3560{TAG}.csv"))
        r = compare(a, b)
        K = pd.read_parquet(f"{C.TB}/kernel_coef.parquet"); K = K[(K.variant == "ell40") & (K.comp == "3c")]
        ad = K.pivot_table(index="series", columns="win", values="a_d")
        ratio = float(np.nanmedian(ad["3560"] / ad["2035"]))
        out.append(dict(gate="f", path=path, stat="D3560_vs_D2035_slope", value=r["slope"], r=r["r"], n=r["n"],
                        passed=bool(0.7 <= r["slope"] <= 1.3) if np.isfinite(r["slope"]) else None,
                        note=f"raw a_d ratio 3.5-6/2-3.5 median {ratio:.2f}; kernel-corrected predicted ratio 1.0"))
    except FileNotFoundError:
        pass
    # ---------------- (g) seasonal
    D = head[head.kind == "D"]
    f = C.annual_fit(pd.to_datetime(D.mid), D.value.to_numpy(), D.sig.to_numpy())
    out.append(dict(gate="g", path=path, stat="D_annual_amp", value=f["amp"], sig=f["amp_sig"], peak_doy=f["peak_doy"],
                    n=f["n"], passed=bool(f["amp"] < 2 * f["amp_sig"])))
    for s in C.SITES:
        q = head[(head.kind == "S") & (head.site == s)]
        fs = C.annual_fit(pd.to_datetime(q.mid), q.value.to_numpy(), q.sig.to_numpy())
        out.append(dict(gate="g", path=path, stat=f"S_{s}_annual_amp", value=fs["amp"], sig=fs["amp_sig"],
                        peak_doy=fs["peak_doy"], n=fs["n"], passed=None))
    # pooled cell annual cycle (for reference: what the data carry, kernel-free)
    t = cells.episode.map(C.episodes().mid)
    fc = C.annual_fit(t, dm.to_numpy(), cells.sig.to_numpy())
    out.append(dict(gate="g", path=path, stat="cells_annual_amp_kernel_free", value=fc["amp"], sig=fc["amp_sig"],
                    peak_doy=fc["peak_doy"], n=fc["n"], passed=None))
    for r in out:
        log(f"{path} gate {r['gate']} {r['stat']}: {r['value']:.4f} " + " ".join(f"{k}={v}" for k, v in r.items()
                                                                              if k not in ('gate', 'path', 'stat', 'value')))
    return out


def main():
    log = C.Log(f"{C.LOGD}/gates{TAG}.log")
    rows = []
    # step-4 precision check (raw 3C 2-5 s bootstrap sigma at n >= 100, boreholes)
    E = pd.read_parquet(f"{C.TB}/episode_raw{TAG}.parquet")
    bh = E[(E.n >= 100) & E.netsta.str.startswith("PB.")]
    s = 100 * bh.groupby("netsta").sig3c_25.median()
    rows.append(dict(gate="step4", path="raw", stat="median_sig3c25_n100_boreholes_max", value=float(s.max()),
                     passed=bool(s.max() <= 0.25), note="; ".join(f"{k} {v:.3f} %" for k, v in s.items())))
    log(f"step-4: borehole per-cell sigma (n >= 100): {s.round(3).to_dict()}")
    for path in ("raw", "svd"):
        try:
            rows += run_path(path, log)
        except Exception as e:
            import traceback
            log(f"ERROR {path}: {e} {traceback.format_exc()[-800:]}")
    G = pd.DataFrame(rows); G.to_csv(f"{C.TB}/gates{TAG}.csv", index=False)
    head = json.load(open(f"{C.TB}/headline{TAG}.json"))["headline"]
    js = {f"{r['gate']}|{r['path']}|{r['stat']}": dict(value=r["value"], passed=r.get("passed")) for r in rows}
    json.dump(dict(headline=head, gates=js, written=C.pt_now()), open(f"{C.TB}/gates{TAG}.json", "w"), indent=1, default=str)
    open(f"{C.OUT}/GATES{TAG}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
