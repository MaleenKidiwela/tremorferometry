#!/usr/bin/env python
"""P4 quiet flag, step 2: option 3 = CN4 array-energy quiet AND not within +-3 d of a Bostock window.
Rules (written before computing): p4/P4.md, "Quiet flag (option 3)".

  python bw_p4_qflag.py     (idempotent)
    * builds Z(d) from every p4/tables/energy/cn4_hourly_<Y>.parquet present;
    * calibrates q on 2010-01-10..2011-12-31 vs PNSN ONCE -> p4/tables/qflag_calibration.json (frozen; never refit);
    * validates on 2012-2013 ONCE when their energy files and P5F_2012/2013_DONE exist -> qflag_validation.json;
    * writes the flag for every day 2003-01-01..2013-12-31 -> p4/tables/quiet_flag_cn4.csv (re-written each run).
"""
import os, sys, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import bw_p4_common as P

EN = f"{P.TAB}/energy"
CAL = f"{P.TAB}/qflag_calibration.json"
VAL = f"{P.TAB}/qflag_validation.json"
FLAG = f"{P.TAB}/quiet_flag_cn4.csv"
EPOCHS = {"PGC": ["2000-06-28", "2017-08-10"], "NLLB": ["2003-01-01", "2017-02-21"],
          "VGZ": ["2003-02-27", "2007-02-28", "2009-07-31", "2017-03-03"],
          "SNB": ["2000-05-25", "2006-04-20", "2008-12-11", "2011-03-17", "2017-06-08"]}
CAL0, CAL1 = "2010-01-10", "2011-12-31"
REF0, REF1 = "2010-01-01", "2011-12-31"
VAL0, VAL1 = "2012-01-01", "2013-12-31"
PAD = 3


def epoch_id(sta, dates):
    b = pd.to_datetime(EPOCHS[sta]); d = pd.to_datetime(dates)
    return np.searchsorted(b.values, d.values, side="right") - 1


def build_Z():
    fs = sorted(glob.glob(f"{EN}/cn4_hourly_*.parquet"))
    years = sorted(int(f[-12:-8]) for f in fs)
    E = pd.concat([pd.read_parquet(f) for f in fs], ignore_index=True)
    E = E[(E.hour >= 0) & np.isfinite(E.e) & (E.e > 0)].copy()
    E["bnear"] = P.bostock_near(E.date.to_numpy(), PAD)
    E["ep"] = -1
    for s in EPOCHS:
        m = E.sta == s
        E.loc[m, "ep"] = epoch_id(s, E.loc[m, "date"])
    refs, prov = {}, {}
    for (s, ep), g in E.groupby(["sta", "ep"]):
        a, b = EPOCHS[s][ep], EPOCHS[s][ep + 1]
        overlaps_cal = a <= REF1 and b > REF0
        gg = g[~g.bnear]
        if overlaps_cal:
            gg = gg[(gg.date >= REF0) & (gg.date <= REF1)]
            prov[(s, ep)] = not ({2010, 2011} <= set(years))
        else:
            need = set(range(max(int(a[:4]), 2003), min(int(b[:4]), 2013) + 1))
            prov[(s, ep)] = not need <= set(years)
        refs[(s, ep)] = float(np.median(gg.e)) if len(gg) else np.nan
    E["ref"] = [refs[(s, ep)] for s, ep in zip(E.sta, E.ep)]
    E["prov"] = [prov[(s, ep)] for s, ep in zip(E.sta, E.ep)]
    E["r"] = E.e / E.ref
    H = E.groupby(["date", "hour"]).agg(net=("r", "median"), n_sta=("r", "count"), prov=("prov", "max")).reset_index()
    H = H[H.n_sta >= 2]
    Z = H.groupby("date").agg(Z=("net", "median"), n_hours=("net", "count"), n_sta=("n_sta", "median"),
                              ref_provisional=("prov", "max")).reset_index()
    Z.loc[Z.n_hours < 12, "Z"] = np.nan
    R = pd.DataFrame([dict(sta=s, epoch_start=EPOCHS[s][ep], epoch_end=EPOCHS[s][ep + 1], ref=v, provisional=prov[(s, ep)])
                      for (s, ep), v in refs.items()])
    return Z, R, years


def flag(cal, thr):
    """cal: DataFrame over a FULL calendar (date, Z) -> active, energy_quiet, bostock_near, quiet_opt3."""
    D = cal.copy()
    D["active"] = (D.Z >= thr).fillna(False).astype(bool)
    a = D.active.to_numpy().astype(int)
    near = np.convolve(a, np.ones(2 * PAD + 1, int), mode="same") > 0
    D["energy_quiet"] = D.Z.notna() & ~near
    D["bostock_near"] = P.bostock_near(D.date.to_numpy(), PAD)
    D["quiet_opt3"] = D.energy_quiet & ~D.bostock_near
    # the other flags, reported beside option 3 (P4.md: "on failure, report both flags")
    D["quiet_bostock3"] = ~D.bostock_near                                  # > 3 d from every Bostock window
    D["quiet_bostock30"] = ~P.bostock_near(D.date.to_numpy(), 30)           # option 1: > 30 d from every window
    return D


def sens_spec(D, lo, hi):
    m = (D.date >= lo) & (D.date <= hi) & D.pnsn_quiet.notna()
    x = D[m]; pq = x.pnsn_quiet.astype(bool); fq = x.quiet_opt3.astype(bool)
    sens = float((~fq[~pq]).mean()) if (~pq).any() else np.nan
    spec = float(fq[pq].mean()) if pq.any() else np.nan
    return sens, spec, int(pq.sum()), int((~pq).sum()), int((fq & pq).sum()), int((fq & ~pq).sum())


def sigma_p(years):
    """per (date): sum of frozen forward WAVE_same P over data and control candidates, and observable family-days."""
    import joblib
    sys.path.insert(0, f"{P.LFE}/scripts")
    import svi_ref_fit as SF
    M = joblib.load(f"{P.REFINE}/conf/models.joblib")["models"]["WAVE_same"]
    out = []
    for y in years:
        d = f"{P.BW}/mf/feat/p5f_{y}"
        for p in sorted(glob.glob(f"{d}/*.parquet")):
            if p.endswith("_fd.parquet") or os.path.getsize(p) == 0:
                continue
            ds = os.path.basename(p)[:10]
            A = pd.read_parquet(p)
            fd = pd.read_parquet(p.replace(".parquet", "_fd.parquet"), columns=["series", "observable"])
            row = dict(date=ds, n_fd_data=int(((fd.series == "data") & fd.observable.astype(bool)).sum()),
                       n_fd_c1=int(((fd.series == "c1") & fd.observable.astype(bool)).sum()), sp_data=0.0, sp_c1=0.0,
                       n_data=0, n_c1=0)
            if len(A):
                for c in SF.WAVE:
                    if c not in A:
                        A[c] = np.nan
                X = A[SF.WAVE].astype(float).to_numpy()
                pr = M.predict_proba(X)[:, 1]
                for s in ("data", "c1"):
                    m = (A.series == s).to_numpy()
                    row[f"sp_{s}"] = float(pr[m].sum()); row[f"n_{s}"] = int(m.sum())
            out.append(row)
    return pd.DataFrame(out)


def main():
    Z, R, years = build_Z()
    R.to_csv(f"{P.TAB}/qflag_station_refs.csv", index=False)
    cal = pd.DataFrame(dict(date=[str(d.date()) for d in pd.date_range("2003-01-01", "2013-12-31")]))
    cal = cal.merge(Z, on="date", how="left")
    Q = P.pnsn_quiet("2003-01-01", "2013-12-31")[["date", "quiet", "n_tremor_day"]].rename(columns={"quiet": "pnsn_quiet"})
    cal = cal.merge(Q, on="date", how="left")
    # ---- calibration (once, frozen)
    if not os.path.exists(CAL):
        if not {2010, 2011} <= set(years):
            raise SystemExit("[qflag] calibration needs the 2010 and 2011 energy files")
        cm = (cal.date >= CAL0) & (cal.date <= CAL1) & (cal.pnsn_quiet == True) & cal.Z.notna()   # noqa: E712
        zq = cal.loc[cm, "Z"].to_numpy()
        rows = []
        for q in np.round(np.arange(0.50, 0.995, 0.01), 2):
            thr = float(np.quantile(zq, q))
            D = flag(cal, thr)
            se, sp, npq, nnq, fq_pq, fq_nq = sens_spec(D, CAL0, CAL1)
            rows.append(dict(q=q, thr=thr, sens=se, spec=sp, J=se + sp - 1, n_pnsn_quiet=npq, n_pnsn_notquiet=nnq,
                             n_flag_quiet=fq_pq + fq_nq, n_flag_quiet_but_pnsn_active=fq_nq))
        C_ = pd.DataFrame(rows); C_.to_csv(f"{P.TAB}/qflag_calibration_curve.csv", index=False)
        best = C_.sort_values(["J", "q"], ascending=[False, False]).iloc[0]
        json.dump(dict(q=float(best.q), thr=float(best.thr), sens_cal=float(best.sens), spec_cal=float(best.spec),
                       J=float(best.J), cal_days=[CAL0, CAL1], n_cal_pnsn_quiet_days_with_Z=int(cm.sum()),
                       frozen_pt=os.popen("TZ=America/Los_Angeles date").read().strip(),
                       rule="q maximising sens+spec of option 3 vs PNSN quiet on 2010-01-10..2011-12-31; ties -> larger q"),
                  open(CAL, "w"), indent=1)
        P.log(f"[qflag] calibrated q={best.q:.2f} thr={best.thr:.4f} sens {best.sens:.3f} spec {best.spec:.3f}")
    cj = json.load(open(CAL)); thr = cj["thr"]
    D = flag(cal, thr)
    D["year_cached"] = D.date.str[:4].astype(int).isin(years)
    D["ref_provisional"] = D.ref_provisional.astype("boolean")
    D["q"] = cj["q"]; D["thr"] = thr
    D.to_csv(FLAG + ".tmp", index=False); os.replace(FLAG + ".tmp", FLAG)
    P.log(f"[qflag] flag written for {int(D.Z.notna().sum())} days with Z (years {years}); option-3 quiet days "
          f"{int(D.quiet_opt3.sum())}")
    # ---- validation (once) on 2012-13
    if not os.path.exists(VAL) and {2012, 2013} <= set(years) and all(os.path.exists(f"{P.BW}/P5F_{y}_DONE") for y in (2012, 2013)):
        se, sp, npq, nnq, fq_pq, fq_nq = sens_spec(D, VAL0, VAL1)
        S = sigma_p([2012, 2013]); S.to_csv(f"{P.TAB}/qflag_sigmaP_days_2012_2013.csv", index=False)
        V = D[(D.date >= VAL0) & (D.date <= VAL1)][["date", "pnsn_quiet", "quiet_opt3"]].merge(S, on="date", how="inner")
        V["month"] = V.date.str[:7]

        def rate(x):
            n = x.n_fd_data.sum()
            return (x.sp_data.sum() - x.sp_c1.sum()) / n if n else np.nan
        Vp = V[V.pnsn_quiet == True]; Vf = V[V.quiet_opt3 == True]                     # noqa: E712
        r_p, r_f = rate(Vp), rate(Vf)
        months = []
        for m, g in V.groupby("month"):
            gp, gf = g[g.pnsn_quiet == True], g[g.quiet_opt3 == True]                  # noqa: E712
            months.append(dict(month=m, n_days_pnsn=len(gp), n_days_opt3=len(gf), rate_pnsn=rate(gp), rate_opt3=rate(gf),
                               ratio=rate(gf) / rate(gp) if len(gp) and rate(gp) else np.nan))
        Mo = pd.DataFrame(months); Mo.to_csv(f"{P.TAB}/qflag_validation_monthly.csv", index=False)
        el = Mo[(Mo.n_days_pnsn >= 5) & (Mo.n_days_opt3 >= 5)]
        frac_ok = float(((el.ratio >= 0.9) & (el.ratio <= 1.1)).mean()) if len(el) else np.nan
        pooled = r_f / r_p if r_p else np.nan
        ok = dict(spec=sp >= 0.9, sens=se >= 0.8, sigmaP_pooled=bool(0.9 <= pooled <= 1.1),
                  sigmaP_months=bool(np.isfinite(frac_ok) and frac_ok >= 0.8))
        json.dump(dict(sens=se, spec=sp, n_pnsn_quiet=npq, n_pnsn_notquiet=nnq, n_flag_quiet_and_pnsn_quiet=fq_pq,
                       n_flag_quiet_but_pnsn_active=fq_nq, rate_pnsn=r_p, rate_opt3=r_f, pooled_ratio=pooled,
                       n_months_eligible=int(len(el)), frac_months_within_10pct=frac_ok, checks=ok,
                       pass_=bool(all(ok.values())), q=cj["q"], thr=thr,
                       written_pt=os.popen("TZ=America/Los_Angeles date").read().strip()),
                  open(VAL, "w"), indent=1, default=float)
        P.log(f"[qflag] validation 2012-13: sens {se:.3f} spec {sp:.3f} sigmaP ratio {pooled:.3f} months-ok {frac_ok} "
              f"-> {'PASS' if all(ok.values()) else 'FAIL'}")


if __name__ == "__main__":
    main()
