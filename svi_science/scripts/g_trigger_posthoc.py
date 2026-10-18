#!/usr/bin/env python
"""g_trigger_posthoc (2026-10-11, POST-HOC, written after the pre-registered tests) -- chance check for the one
earthquake that met the individual 'named' rule (WA, Kermadec 2019-06-15 M7.3): how often does a burst of
LFEs / tremor START after a long quiet near the region reference, and how many such starts would fall inside the
W_inst windows of all trials by chance?

Onset definition (fixed before computing): first item at time t with (i) none in the previous 72 h and (ii) >= 20
items in [t, t + 12 h); items = PNSN epicentres within 100 km of the region reference (2010-01-06 .. 2026-05-27), or
counted-once LFEs (data, P rule) of the region.  Expected number in the trial windows = onset rate x total window
length (overlaps merged); P(>= observed) from the Poisson law.  Also the windows' onset count by status.

    python svi_science/scripts/g_trigger_posthoc.py  -> tables/posthoc_onsets.csv, tables/posthoc_onset_list.csv
"""
import numpy as np, pandas as pd
from scipy import stats

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/g_trigger/tables"
NS = 10 ** 9; H = 3600 * NS
REF = {"SVI": (48.4837, -123.7485), "WA": (47.8432, -123.0282)}
S0, S1 = pd.Timestamp("2010-01-06").value, pd.Timestamp("2026-05-28").value


def onsets(t, quiet=72 * H, nmin=20, win=12 * H):
    t = np.sort(t); out = []
    for i in range(len(t)):
        if (i == 0 or t[i] - t[i - 1] >= quiet) and np.searchsorted(t, t[i] + win) - i >= nmin:
            out.append(t[i])
    return np.array(out, dtype=np.int64)


def merged_len(W):
    W = sorted(W); tot = 0; cs, ce = None, None
    for s, e in W:
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return tot + (ce - cs if cs is not None else 0)


def main():
    tr = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["tns"] = pd.to_datetime(tr.time).astype("datetime64[ns]").astype("int64")
    rows, lst = [], []
    for reg in ("SVI", "WA"):
        la, lo = REF[reg]
        dk = 111.195 * np.hypot(tr.lat - la, (tr.lon - lo) * np.cos(np.radians(la)))
        items = {"tremor100": tr.tns[dk <= 100].to_numpy()}
        ev = pd.read_parquet(f"{T}/{'svi' if reg == 'SVI' else 'wa'}_events.parquet")
        items["lfe"] = np.sort(ev.time_utc[ev.series == "data"].astype("datetime64[ns]").astype("int64").to_numpy())
        TT = pd.read_csv(f"{T}/inst_trials_{reg}.csv")
        TT = TT[(TT.inst_s >= S0) & (TT.inst_s < S1) & TT.status.isin(["between", "ets", "pad"])]
        for kind, t in items.items():
            t = t[(t >= S0) & (t < S1)]
            for nmin in (10, 20):
                on = onsets(t, nmin=nmin)
                rate = len(on) / (S1 - S0)
                for st in ("all", "ets", "between"):
                    q = TT if st == "all" else TT[TT.status == st]
                    W = list(zip(q.inst_s, q.inst_e))
                    L = merged_len(W)
                    hit = [o for o in on if any(s <= o < e for s, e in W)]
                    exp = rate * L
                    rows.append(dict(region=reg, items=kind, nmin=nmin, trials=st, n_trials=len(q), n_onsets_total=len(on),
                                     onsets_per_year=len(on) / ((S1 - S0) / (365.25 * 86400 * NS)),
                                     window_hours=L / H, expected_in_windows=exp, observed_in_windows=len(hit),
                                     p_poisson=stats.poisson.sf(len(hit) - 1, exp) if len(hit) else 1.0))
                    for o in hit:
                        m = q[(q.inst_s <= o) & (q.inst_e > o)].iloc[0]
                        lst.append(dict(region=reg, items=kind, nmin=nmin, trials=st, onset_PT=pd.Timestamp(o).tz_localize("UTC")
                                        .tz_convert("America/Los_Angeles").strftime("%Y-%m-%d %H:%M %Z"),
                                        eq=m.place, eq_PT=m.t0_PT, mag=m.mag, sigma_kpa=m.sigma_kpa,
                                        minutes_after_origin=(o - m.t0_ns) / 60e9))
    R = pd.DataFrame(rows); R.to_csv(f"{T}/posthoc_onsets.csv", index=False)
    L = pd.DataFrame(lst).drop_duplicates(); L.to_csv(f"{T}/posthoc_onset_list.csv", index=False)
    pd.set_option("display.width", 250)
    print(R.to_string()); print(L.to_string())


if __name__ == "__main__":
    main()
