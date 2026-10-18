#!/usr/bin/env python
"""Track D, D1(b) -- APPLY the calibrated sensitivity model to every production family-day:
results/trackD/mc_family_day.parquet with u, v, Mc(f,d) and P_det at a reference magnitude, per tier.

The day terms are formed exactly as in the calibration (`d1_mc_model.py` docstring):
  u = log10( mad_scale(f,d) / mad_ref(f, epoch) ), mad_ref = the median of the 2012 ablation run that
      used the SAME station set as the epoch (so u is a pure day-to-day anomaly);
  v = 0.5*log10( sum_s w(f,s)^2 * e8_s(d)/e8_ref(s, instrument epoch) / sum_s w(f,s)^2 ), over the
      stations that carried a usable 2012 template for f AND contributed that day (from the family-day's
      `epoch` station list, its `sites`, minus `excluded_sta`).  e8_ref is the median over that station's
      INSTRUMENT epoch across the whole record (referee years included), so an instrument gain change
      cancels and only relative noise -- including the seasonal cycle -- survives."""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd
from scipy.stats import norm
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"
OUT = f"{EXT}/results"
NETSTA = dict(PGC="CN.PGC", VGZ="CN.VGZ", NLLB="CN.NLLB", SNB="CN.SNB", LRIV="UW.LRIV", B001="PB.B001",
              B003="PB.B003", B004="PB.B004", B005="PB.B005", B006="PB.B006", B007="PB.B007",
              B009="PB.B009", B010="PB.B010", B011="PB.B011", B926="PB.B926")
SITE_OF = {s: ("PGC" if s in ("PGC", "B009", "B010", "B011") else
               "LRIV" if s in ("LRIV", "B005", "B006", "B007") else s) for s in NETSTA}
SET_REF_TAG = {"S15": "a1_mergeB_L12", "S14": "d1_abl2012_E2", "S13": "d1_abl2012_E3"}
# instrument epochs used to normalise e8 (station -> list of (start, end)); from results/station_epochs.csv
# plus the A-2 confirmed sensor boundaries (SNB 2017-06-08 and VGZ 2017-03-03 are same-band sensor swaps).
INSTR_BOUNDARY = {"PGC": ["2017-08-10"], "NLLB": ["2017-02-23"], "VGZ": ["2017-03-03"],
                  "SNB": ["2011-03-17", "2017-06-08"], "LRIV": ["2019-05-30"]}


def instr_epoch(sta, date):
    b = INSTR_BOUNDARY.get(sta, [])
    return sum(1 for x in b if date >= x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mag-ref", type=float, default=1.8)
    a = ap.parse_args()
    EPOCHS, TIERS, TSRC = d_tiers.load()
    SET_OF = {e["epoch_id"]: e["station_set"] for e in EPOCHS}
    print(f"[tiers] {TSRC}")
    P = pd.read_csv(f"{TD}/mc_model_params.csv")
    w = pd.read_csv(f"{IN}/e5r_weights_snr_2012.csv")
    W = {(int(r.fam), r.sta): float(r.w) for r in w.itertuples()}
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    qc = qc[qc.usable.astype(bool) & qc.template.astype(bool)]
    TPL = qc.groupby("bkey").sta.apply(set).to_dict()

    # ---- e8 reference per (netsta, instrument epoch), over the referee years AND every production month
    dn = [pd.read_parquet(f"{TD}/daily_noise_ref.parquet")]
    dn += [pd.read_parquet(p) for p in sorted(glob.glob(f"{TC}/daily_noise_*.parquet"))]
    cols = ["netsta", "year", "doy", "status", "e8_Z", "e8_N", "e8_E"]
    DN = pd.concat([d[[c for c in cols if c in d.columns]] for d in dn], ignore_index=True)
    DN = DN[DN.status.isin(["cached", "built"])].copy()
    DN["e8"] = DN[["e8_Z", "e8_N", "e8_E"]].mean(axis=1)
    DN["date"] = [str((pd.Timestamp(f"{int(y)}-01-01") + pd.Timedelta(days=int(d) - 1)).date())
                  for y, d in zip(DN.year, DN.doy)]
    S_OF = {v: k for k, v in NETSTA.items()}
    DN["sta"] = DN.netsta.map(S_OF)
    DN = DN[DN.sta.notna()]
    DN["iep"] = [instr_epoch(s, d) for s, d in zip(DN.sta, DN.date)]
    ref = DN.groupby(["sta", "iep"]).e8.median().to_dict()
    DN["r"] = DN.e8 / [ref[(s, i)] for s, i in zip(DN.sta, DN.iep)]
    RMAP = {(r.sta, r.date): float(r.r) for r in DN.itertuples() if np.isfinite(r.r)}

    # ---- mad_ref(f, epoch) from the 2012 run with that station set
    madref = {}
    for ss, tag in SET_REF_TAG.items():
        if not os.path.exists(f"{OUT}/e5_budget_{tag}.csv"):
            continue
        d = pd.read_csv(f"{OUT}/e5_budget_{tag}.csv")
        d = d[d.L == 12.0]
        d = d[d["T"] == sorted(d["T"].unique())[0]]
        m = d.drop_duplicates(["fam", "doy"]).groupby("fam").madn.median()
        for f_, v_ in m.items():
            madref[(ss, int(f_))] = float(v_)

    # ---- production family-days
    fdf = sorted(glob.glob(f"{TC}/family_days_*.parquet"))
    if not fdf:
        raise SystemExit("[d1b-apply] no Track C months on disk")
    FD = pd.concat([pd.read_parquet(p) for p in fdf], ignore_index=True).drop_duplicates(["fam", "date"])
    FD["epoch_id"] = d_tiers.epoch_of(FD.date.to_numpy().astype(str), EPOCHS)
    FD = FD[FD.epoch_id != ""].reset_index(drop=True)
    FD["station_set"] = FD.epoch_id.map(SET_OF).fillna("S15")

    # u
    mr = np.array([madref.get((e, int(f)), np.nan) for e, f in zip(FD.station_set, FD.fam)])
    FD["mad_ref"] = mr
    FD["u"] = np.log10(FD.mad_scale / FD.mad_ref)
    # v
    vs = np.full(len(FD), np.nan)
    ep_sta = {}
    for i, r in enumerate(FD.itertuples()):
        key = r.epoch
        if key not in ep_sta:
            ep_sta[key] = set(x.split(":")[0] for x in str(key).split("|") if x)
        stas = ep_sta[key]
        sites = set(str(r.sites).split(",")) if isinstance(r.sites, str) else set()
        excl = set(str(r.excluded_sta).split(",")) if isinstance(r.excluded_sta, str) and r.excluded_sta else set()
        tpl = TPL.get(int(r.fam), set())
        use = [s for s in stas & tpl if SITE_OF.get(s, s) in sites and s not in excl]
        num = den = 0.0
        for s in use:
            ww = W.get((int(r.fam), s))
            rr = RMAP.get((s, r.date))
            if ww is None or rr is None or not np.isfinite(rr):
                continue
            num += ww * ww * rr; den += ww * ww
        if den > 0:
            vs[i] = 0.5 * np.log10(num / den)
    FD["v"] = vs
    FD["u"] = FD.u.fillna(0.0); FD["v"] = FD.v.fillna(0.0)

    # ---- Mc per tier
    out = FD[["fam", "date", "epoch_id", "station_set", "observable", "n_sites", "n_sta", "quiet_flag",
              "mad_scale", "mad_ref", "u", "v"]].copy().rename(columns={"fam": "family_id"})
    for tier in ("T_ext", "T_core"):
        pr = P[(P.tier == tier)]
        pe = {(r.epoch_id, int(r.family_id)): r.b_family for r in pr.itertuples()}
        g = pr.groupby("epoch_id").first()
        b0 = g.b0.to_dict(); al = g.alpha.to_dict(); c1 = g.c1.to_dict(); c2 = g.c2.to_dict()
        eta0 = np.array([b0.get(e, np.nan) + pe.get((e, int(f)), np.nan)
                         + c1.get(e, np.nan) * u + c2.get(e, np.nan) * v
                         for e, f, u, v in zip(out.epoch_id, out.family_id, out.u, out.v)])
        alp = np.array([al.get(e, np.nan) for e in out.epoch_id])
        out[f"Mc_{tier}"] = a.mag_ref - eta0 / alp
        out[f"Pdet_{tier}_M{a.mag_ref:g}"] = norm.cdf(eta0)
    # the referee windows are late-summer ETS windows: flag family-days whose day terms fall OUTSIDE the
    # 1-99 % range the coefficients were calibrated over, where Mc is an extrapolation.
    g2 = P[P.tier == "T_core"].groupby("epoch_id").first()
    for c, lo, hi in (("u", "u_calib_lo", "u_calib_hi"), ("v", "v_calib_lo", "v_calib_hi")):
        l = np.array([g2[lo].get(e, np.nan) for e in out.epoch_id])
        h = np.array([g2[hi].get(e, np.nan) for e in out.epoch_id])
        out[f"{c}_in_calib_range"] = (out[c].to_numpy() >= l) & (out[c].to_numpy() <= h)
    out["day_terms_in_calib_range"] = out.u_in_calib_range & out.v_in_calib_range
    out.loc[~out.observable, [c for c in out.columns if c.startswith("Mc_") or c.startswith("Pdet_")]] = np.nan
    out["template_less"] = ~out.family_id.isin(TPL.keys())
    out.to_parquet(f"{TD}/mc_family_day.parquet", index=False)
    print(out[["Mc_T_ext", "Mc_T_core", "Pdet_T_ext_M1.8", "Pdet_T_core_M1.8", "u", "v"]]
          .describe().round(3).to_string())
    print(f"[write] {TD}/mc_family_day.parquet  ({len(out)} family-days, "
          f"{out.date.min()}..{out.date.max()})")


if __name__ == "__main__":
    main()
