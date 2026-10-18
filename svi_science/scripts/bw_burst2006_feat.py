#!/usr/bin/env python
"""bw_burst2006_feat -- who fires and who carries the correlation on the early-2006 burst days (first-look follow-up).

Read-only on mf/feat and p4/envflag/win.  Single process, light.  No Bostock events read (only family positions from
data/catalogue/families_v1_3.csv, SVI:<k> = detector family k, pos_lat/pos_lon).
Out: svi_science/results/backward/first_look/burst2006/{day_summary.csv, station_share.csv, hourly.csv,
     fam_counts.csv, env_windows.csv, fig_*.png}
"""
import os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"
BW = f"{ROOT}/svi_science/results/backward"
FEAT = f"{BW}/mf/feat"
ENV = f"{BW}/p4/envflag/win"
OUT = f"{BW}/first_look/burst2006"
T = 9.0                                               # frozen B2 T_core; era shown at the same provisional 9.0
DAYS = ["2006-01-15", "2006-01-16", "2006-03-27", "2006-03-28", "2006-03-29", "2006-03-30", "2006-04-24",
        "2006-05-01", "2006-05-03", "2006-05-11", "2006-05-25", "2006-05-26", "2006-06-18"]
REF = [("p5f_2006", "2006-03-20", "quiet ref"), ("p5f_2007", "2007-01-28", "ETS ref"),
       ("p5e_2005", "2005-09-05", "ETS ref (era)"), ("p5e_2006", "2006-03-20", "quiet ref (era)")]
STRIKE = 310.0                                        # along-strike axis azimuth (NW positive); down-dip = 40 deg
ENV_C = {"CN4": 0.42, "all": 0.32}                    # frozen calibration (ENVFLAG.md)

F = pd.read_csv(f"{ROOT}/data/catalogue/families_v1_3.csv")
F = F[F.region == "S_Vancouver_Is"].copy()
F["fam"] = F.family_uid.str.split(":").str[1].astype(int)
F = F.set_index("fam")[["pos_lat", "pos_lon", "depth_km"]].assign(in_v13=True)
# 15 detector families (of 130) are not in v1.3 (dropped there); their position falls back to e5r_families_2012 (same bkey)
_E = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv").set_index("bkey")
_miss = _E.index.difference(F.index)
F = pd.concat([F, pd.DataFrame(dict(pos_lat=_E.lat[_miss], pos_lon=_E.lon[_miss], depth_km=_E.dep[_miss], in_v13=False),
                               index=_miss)]).sort_index()
LAT0, LON0 = F.pos_lat.mean(), F.pos_lon.mean()


def xy(lat, lon):
    x = (np.asarray(lon) - LON0) * 111.19 * np.cos(np.radians(LAT0)); y = (np.asarray(lat) - LAT0) * 111.19
    a = np.radians(STRIKE)
    s = x * np.sin(a) + y * np.cos(a)                 # along strike, km, NW +
    a2 = np.radians(STRIKE - 270)                     # 40 deg
    d = x * np.sin(a2) + y * np.cos(a2)               # down dip, km, NE +
    return s, d


F["s_km"], F["d_km"] = xy(F.pos_lat, F.pos_lon)


def load(run, d):
    p = f"{FEAT}/{run}/{d}.parquet"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None, None
    A = pd.read_parquet(p)
    J = json.load(open(f"{FEAT}/{run}/{d}.json"))
    return A, J


def summarise(run, d, label):
    A, J = load(run, d)
    if A is None:
        return None, [], [], []
    stas = J["stations"]
    rows, share, hourly, famc = {}, [], [], []
    out = dict(run=run, date=d, label=label, stations=",".join(stas), n_sta=len(stas))
    for s in ("data", "c1"):
        a = A[(A.series == s) & (A.net_mad >= T)].copy()
        out[f"n_{s}"] = len(a)
        if not len(a):
            continue
        out[f"n_fam_{s}"] = a.fam.nunique()
        out[f"dom_top_med_{s}"] = float(a.dom_top.median())
        out[f"frac_z3_med_{s}"] = float(a.frac_z3.median())
        out[f"n_z3_med_{s}"] = float(a.n_z3.median())
        out[f"hfr_rel_med_{s}"] = float(a.hfr_rel.median())
        out[f"net_cc_med_{s}"] = float(a.net_cc.median())
        zc = [c for c in a.columns if c.startswith("z_") and c[2:] in stas]
        Z = a[zc].to_numpy(float)
        Zn = np.where(np.isfinite(Z), Z, -np.inf)
        top = np.array([zc[i][2:] for i in Zn.argmax(1)])
        for c in zc:
            st = c[2:]
            share.append(dict(run=run, date=d, series=s, sta=st, z_med=float(np.nanmedian(a[c])),
                              cc_med=float(np.nanmedian(a[f"cc_{st}"])), frac_top=float((top == st).mean()),
                              frac_z3=float((a[c] >= 3).mean()), frac_nan=float(a[c].isna().mean())))
        # timing
        h = ((a.t.to_numpy() - pd.Timestamp(d).value / 1e9) // 3600).astype(int)
        for k, v in zip(*np.unique(h, return_counts=True)):
            hourly.append(dict(run=run, date=d, series=s, hour_utc=int(k), n=int(v)))
        # families / positions
        vc = a.fam.value_counts()
        for fam, n in vc.items():
            famc.append(dict(run=run, date=d, series=s, fam=int(fam), n=int(n)))
        pos = F.reindex(a.fam.to_numpy())
        sk, dk = pos.s_km.to_numpy(), pos.d_km.to_numpy()
        out[f"s_q10_{s}"], out[f"s_q50_{s}"], out[f"s_q90_{s}"] = np.nanpercentile(sk, [10, 50, 90])
        out[f"d_q10_{s}"], out[f"d_q50_{s}"], out[f"d_q90_{s}"] = np.nanpercentile(dk, [10, 50, 90])
        out[f"lat_med_{s}"], out[f"lon_med_{s}"] = float(np.nanmedian(pos.pos_lat)), float(np.nanmedian(pos.pos_lon))
        tt = pd.to_datetime(a.t, unit="s")
        out[f"t_q10_{s}"], out[f"t_q50_{s}"], out[f"t_q90_{s}"] = [str(x.strftime("%H:%M")) for x in tt.quantile([.1, .5, .9])]
        topst = pd.Series(top).value_counts(normalize=True)
        out[f"top_sta_{s}"] = f"{topst.index[0]}:{topst.iloc[0]:.2f}"
    return out, share, hourly, famc


def env_windows(d):
    y = d[:4]
    p = f"{ENV}/env_{y}.parquet"
    if not os.path.exists(p):
        return pd.DataFrame()
    W = pd.read_parquet(p, filters=[("date", "==", d)])
    W = W[W.series == "data"].copy()
    W["C"] = W.sset.map(ENV_C)
    W["tremor"] = (W.n_valid >= 3) & (W.coh >= W.C) & (W.s3 >= W.C) & (W.dfam_km <= 30)
    return W


def main():
    os.makedirs(OUT, exist_ok=True)
    S, SH, H, FC = [], [], [], []
    jobs = [(r, d, "burst") for d in DAYS for r in ("p5f_2006", "p5e_2006")] + REF
    for r, d, lab in jobs:
        o, sh, h, fc = summarise(r, d, lab)
        if o is None:
            continue
        S.append(o); SH += sh; H += h; FC += fc
    S = pd.DataFrame(S); SH = pd.DataFrame(SH); H = pd.DataFrame(H); FC = pd.DataFrame(FC)
    S.to_csv(f"{OUT}/day_summary.csv", index=False, float_format="%.3f")
    SH.to_csv(f"{OUT}/station_share.csv", index=False, float_format="%.3f")
    H.to_csv(f"{OUT}/hourly.csv", index=False)
    FC.to_csv(f"{OUT}/fam_counts.csv", index=False)
    # envelope tremor windows
    E = []
    for d in DAYS + ["2006-03-20"]:
        W = env_windows(d)
        if len(W):
            E.append(W)
    E = pd.concat(E, ignore_index=True)
    E.to_csv(f"{OUT}/env_windows.csv", index=False, float_format="%.4f")
    Et = E[E.tremor]
    print(Et.groupby(["date", "sset"]).agg(n=("win", "size"), lat=("lat", "median"), lon=("lon", "median"),
                                           t0=("win", "min"), t1=("win", "max"), imp=("imp", "median")).to_string())

    # ---- figure 1: hourly data vs control (frozen + era) and env tremor windows, per burst day
    bd = ["2006-01-15", "2006-01-16", "2006-03-28", "2006-03-29", "2006-03-30", "2006-04-24", "2006-05-26", "2006-06-18"]
    fig, ax = plt.subplots(len(bd), 1, figsize=(10, 2.0 * len(bd)), sharex=True)
    for i, d in enumerate(bd):
        a = ax[i]
        for r, ls in (("p5f_2006", "-"), ("p5e_2006", ":")):
            for s, c in (("data", "C0"), ("c1", "C3")):
                h = H[(H.run == r) & (H.date == d) & (H.series == s)].set_index("hour_utc").n.reindex(range(24), fill_value=0)
                a.step(np.arange(24) + 0.5, h.values, where="mid", color=c, ls=ls,
                       label=f"{'frozen' if r.startswith('p5f') else 'era'} {'data' if s == 'data' else 'scrambled control'}")
        b = a.twinx()
        for ss, c in (("CN4", "k"), ("all", "C2")):
            w = E[(E.date == d) & (E.sset == ss) & E.tremor]
            hh = np.bincount(((w.win.astype(np.int64) * 150 + 150) // 3600), minlength=24)[:24]
            b.bar(np.arange(24) + (0.3 if ss == "all" else 0.7), hh, width=0.35, color=c, alpha=0.3,
                  label=f"env tremor windows ({ss})")
        b.set_ylabel("env win/h", fontsize=7)
        a.set_ylabel("cand./h", fontsize=7)
        a.set_title(f"{d} (UTC; PT = UTC-8)", fontsize=9, loc="left")
        if i == 0:
            a.legend(fontsize=6, ncol=4, loc="upper right"); b.legend(fontsize=6, loc="upper left")
    ax[-1].set_xlabel("hour of UTC day")
    ax[-1].set_xlim(0, 24)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_hourly.png", dpi=110); plt.close(fig)

    # ---- figure 2: family positions firing (data minus control counts), along-strike vs down-dip
    fd = ["2006-01-15", "2006-03-28", "2006-03-29", "2006-04-24", "2006-05-26", "2006-06-18"]
    fig, ax = plt.subplots(2, len(fd), figsize=(3.0 * len(fd), 6), sharex=True, sharey=True)
    for j, d in enumerate(fd):
        for i, r in enumerate(("p5f_2006", "p5e_2006")):
            a = ax[i, j]
            a.scatter(F.s_km, F.d_km, s=4, c="0.8")
            g = FC[(FC.run == r) & (FC.date == d)].pivot_table(index="fam", columns="series", values="n", fill_value=0)
            if "c1" not in g: g["c1"] = 0
            if "data" not in g: g["data"] = 0
            ex = (g["data"] - g["c1"])
            p = F.reindex(ex.index)
            pos_ = ex > 0
            a.scatter(p.s_km[pos_], p.d_km[pos_], s=3 * ex[pos_].clip(0, 80), c="C0", alpha=0.6, label="data > control")
            a.scatter(p.s_km[~pos_], p.d_km[~pos_], s=3 * (-ex[~pos_]).clip(0, 80), c="C3", alpha=0.6, label="control > data")
            a.set_title(f"{d} {'frozen' if i == 0 else 'era'}\n{int(g['data'].sum())} data / {int(g['c1'].sum())} ctrl", fontsize=8)
            if j == 0: a.set_ylabel("down-dip km (NE +)")
            if i == 1: a.set_xlabel("along-strike km (NW +)")
    ax[0, 0].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_families.png", dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
