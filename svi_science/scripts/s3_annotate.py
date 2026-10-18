"""S3 descriptive add-ons (post-hoc, labelled as such): front coherence class, start over-spread p, speed of the
activity-return streaks, margin-scale pattern for SVI/NW episodes that overlap in time."""
import json
import numpy as np
import pandas as pd
import importlib.util

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s3"
spec = importlib.util.spec_from_file_location("k", f"{ROOT}/svi_science/scripts/s3_kinematics.py"); KM = importlib.util.module_from_spec(spec); spec.loader.exec_module(KM)


def main():
    K = pd.read_csv(f"{OUT}/episodes_kinematics.csv")
    K["coherent"] = (np.sign(K.lfe_end10 - K.lfe_start10) == np.sign(K.lfe_speed)) & (K.lfe_speed.abs() >= 1.5)
    K["front_class"] = np.where(~K.migrates, "no clear migration", np.where(~K.coherent, "no single front",
                                np.where(K.lfe_speed > 0, "one front toward NW", "one front toward SE")))
    B = pd.read_csv(f"{OUT}/bursts.csv"); B = B[B.series == "LFE"]
    EP = pd.read_csv(f"{OUT}/episodes_all.csv", parse_dates=["w0", "w1"])
    rs = []
    for reg in ["SVI", "NW"]:
        E = pd.read_parquet(f"{OUT}/events_{reg.lower()}.parquet")
        if "selected_v2" not in E:
            E["selected_v2"] = True
        for ep, b in B[B.region == reg].groupby("episode"):
            w0 = EP[(EP.region == reg) & (EP.episode == ep)].w0.iloc[0]
            h = E[(E.ep == ep) & (E.series == "data") & (E.P_same >= 0.75) & E.selected_v2 & E.counted_once]
            d = KM.pt_day(h.time, w0); dr = int(b.return_day_idx.min())
            m = (d >= dr - 1) & (d <= dr + 1)
            t = KM.days(h.time[m], w0); a = h.as_km.to_numpy()[m]
            v = KM.ts_slope(t, a) if m.sum() >= 10 else np.nan
            rs.append(dict(region=reg, episode=ep, n_returns=len(b), first_return_day_PT=b.sort_values("return_day_idx").return_day_PT.iloc[0],
                           return_window_events=int(m.sum()), return_window_speed_km_d=v,
                           return_window_span_km=float(np.ptp(a)) if m.sum() else np.nan))
    R = pd.DataFrame(rs)
    K = K.merge(R[["region", "episode", "return_window_speed_km_d", "return_window_span_km"]], on=["region", "episode"], how="left")
    K.to_csv(f"{OUT}/episodes_kinematics.csv", index=False)
    R.to_csv(f"{OUT}/return_streaks.csv", index=False)
    # starts more spread than uniform (post-hoc direction of the pre-registered test)
    rep = json.load(open(f"{OUT}/repeatability.json"))
    for reg in rep:
        rep[reg]["p_starts_more_spread_than_uniform"] = 1 - rep[reg]["p_uniform"] + 1e-4
        k = K[K.region == reg]
        rep[reg]["n_one_front"] = int(k.front_class.str.startswith("one front").sum())
    json.dump(rep, open(f"{OUT}/repeatability.json", "w"), indent=1)
    # margin-scale pattern
    X = pd.read_csv(f"{OUT}/cross_region_overlaps.csv")
    ks = K[K.region == "SVI"].set_index("episode"); kn = K[K.region == "NW"].set_index("episode")
    pat = []
    for _, r in X.iterrows():
        cs, cn = ks.loc[r.svi_episode, "front_class"], kn.loc[r.nw_episode, "front_class"]
        if cs.endswith("NW") and cn.endswith("NW"):
            p = "one sweep from Washington to Vancouver Island (toward NW)"
        elif cs.endswith("SE") and cn.endswith("SE"):
            p = "one sweep from Vancouver Island to Washington (toward SE)"
        elif cs.endswith("NW") and cn.endswith("SE"):
            p = "started between the two (Strait of Juan de Fuca) and split both ways"
        else:
            p = f"mixed (SVI: {cs}; NW: {cn})"
        pat.append(p)
    X["svi_front"] = [ks.loc[e, "front_class"] for e in X.svi_episode]; X["nw_front"] = [kn.loc[e, "front_class"] for e in X.nw_episode]
    X["margin_pattern"] = pat
    X.to_csv(f"{OUT}/cross_region_overlaps.csv", index=False)
    print(K[["region", "episode", "front_class", "lfe_start10", "lfe_end10", "lfe_speed", "front_speed_onsets", "trem_start10", "trem_speed",
             "n_bursts", "n_bursts_shuffle_mean", "return_window_speed_km_d"]].round(1).to_string())
    print(X.to_string()); print(R.round(1).to_string())


if __name__ == "__main__":
    main()
