#!/usr/bin/env python
"""DUC_FIN step 3 -- raw-count repeat windows for dv/v around EVERY >= 8 x MAD Ducellier candidate, 2008-01..2012-07.

Uses scripts/cut_windows.py UNCHANGED (imported) with a new run profile 'duc' registered at run time:
  candidates  OUT/catalogue_final.parquet (duc_fin_join.py; all of them, with xmad, cc, nchannel, nlive, mad, ccsum,
              thr_perm, n_tpl, n_new_tpl, scan (= catalogue_final.src) and the pass_his flag)
  offsets     EVERY station in the raw folders gets an S-arrival offset (duc_fin_common.arrival_offsets): his template's
              energy peak at the family's template stations, CVM1.7 1-D moveout + the family's correction elsewhere.
              (The previous cut fell back to the detection time for every station without a template.)
  window      -15 s .. +30 s around that arrival, raw counts, native rate, every channel of the preferred band, on every
              station with a raw day file (traces/ and traces_duc/)
  metadata    channel epochs from the LOCAL StationXML files (resp/<NET.STA>.xml); no FDSN query, nothing downloaded
Output CUT_OUT = results/cut/duc/: windows_<M>.h5, segments_<M>.parquet, index_<M>.parquet (offset_source per row:
template / moveout_cvm17 / detection_time), detections.parquet (det_id), stationdays*.parquet, cut_summary.csv.

  python scripts/duc_fin_cut.py cut [--workers 6] [--months 2009-06] [--max-days 1]
  python scripts/duc_fin_cut.py verify [--n 20]
  python scripts/duc_fin_cut.py stackcheck
"""
import os, sys, json, glob, argparse, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
import cut_windows as CW
from t_common import log, TR

RUN = "duc"


def raw_stations():
    s = set()
    for r in (TR, DC.TR2):
        if os.path.isdir(r): s |= {x for x in os.listdir(r) if "." in x and os.path.isdir(f"{r}/{x}")}
    return sorted(s)


def prof_duc():
    D = pd.read_parquet(f"{F.OUT}/catalogue_final.parquet")
    D["fam"] = D.fam.astype(str); D["t"] = D.t.astype("datetime64[ns]")
    D["xmad"] = D.cc / D.mad; D["pass_his"] = D.pass_his.astype(bool)
    D = D.rename(columns={"src": "scan"})                 # cut_windows' segment table owns the column name 'src
    p = f"{F.OUT}/arrival_offsets_all.csv"
    if not os.path.exists(p):
        F.arrival_offsets(all_stations=raw_stations()).to_csv(p, index=False)
    O = pd.read_csv(p, dtype={"fam": str})
    OFF = pd.DataFrame(dict(fam=O.fam, netsta=O.netsta, offset_s=O.offset_s, offset_source=O.offset_source))
    return dict(det=D, off=OFF, stations=F.stations(), t0=F.T0, t1=F.T1, roots=[TR, DC.TR2], flags=["pass_his"],
                stats=["xmad", "cc", "nchannel", "nlive", "mad", "ccsum", "thr_perm", "n_tpl", "n_new_tpl", "scan"],
                tpl=None, done="DUC_FIN_DONE", failed="DUC_FIN_FAILED", label="ducellier")


def local_channel_epochs(netstas, t0, t1, path):
    """cut_windows.channel_epochs from the local StationXML (resp/<NET.STA>.xml) -- no network."""
    if os.path.exists(path):
        E = pd.read_csv(path, dtype={"loc": str}); E["loc"] = E["loc"].fillna("").astype(str)
        if set(netstas) <= set(E.netsta): return E
    from obspy import read_inventory
    rows = []
    for ns in netstas:
        x = f"{DC.RESP}/{ns}.xml"
        if not os.path.exists(x):
            rows.append(dict(netsta=ns, loc="", cha="", start="", end="", sr=np.nan, sensor="no local StationXML", prov="")); continue
        inv = read_inventory(x, format="STATIONXML")
        for net in inv:
            for st in net:
                for ch in st:
                    if ch.code[:2] not in ("HH", "BH", "EH", "SH"): continue
                    sens = ch.response.instrument_sensitivity if ch.response is not None else None
                    rows.append(dict(netsta=ns, loc=ch.location_code, cha=ch.code, start=str(ch.start_date)[:19],
                                     end=str(ch.end_date)[:19] if ch.end_date else "", sr=ch.sample_rate,
                                     sensor=(ch.sensor.description or ch.sensor.model or "") if ch.sensor else "",
                                     gain=sens.value if sens else np.nan, az=ch.azimuth, dip=ch.dip, prov=f"local:{x}"))
    E = pd.DataFrame(rows); E.to_csv(path, index=False); E["loc"] = E["loc"].fillna("").astype(str)
    return E


def setup():
    CW.PROFILES[RUN] = prof_duc
    CW.channel_epochs = local_channel_epochs


def fix_offset_source(out):
    """cut_windows writes ONE offset_source text for every non-fallback row; put the per-(family, station) source."""
    O = pd.read_csv(f"{F.OUT}/arrival_offsets_all.csv", dtype={"fam": str})
    src = {(r.fam, r.netsta): r.offset_source for r in O.itertuples()}
    n = {}
    for p in sorted(glob.glob(f"{out}/index_*.parquet")):
        I = pd.read_parquet(p)
        if not len(I): continue
        s = np.array([src.get((str(f), ns), "detection_time") if osrc != "detection_time" else "detection_time"
                      for f, ns, osrc in zip(I.fam.astype(str), I.netsta, I.offset_source)], dtype=object)
        I["offset_source"] = s; I.to_parquet(p, index=False)
        for k, v in pd.Series(s).value_counts().items(): n[k] = n.get(k, 0) + int(v)
    json.dump(n, open(f"{out}/offset_source_counts.json", "w"), indent=1)
    log(f"[cut] offset sources over all windows: {n}")


def cutguard(out):
    """station-day read guard for the cut: every Ducellier station-day that duc_fix's plan says has a raw file must
    appear in stationdays with status ok; read errors / no usable band of any station count as errors."""
    SD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{out}/stationdays_*.parquet"))], ignore_index=True)
    months = {os.path.basename(p)[5:12] for p in glob.glob(f"{out}/done_*.json")}
    seen = {(r.netsta, r.day): r.status for r in SD.itertuples()}
    rows = []
    for m in sorted(months):
        p = f"{F.FIX}/data/plan_{m}.csv"
        if not os.path.exists(p): continue
        P = pd.read_csv(p, dtype={"loc": str})
        P = P[P.day <= SD[SD.day.str[:7] == m].day.max()]          # (a smoke run cuts only the first days of a month)
        for r in P.itertuples():
            ns = f"{r.net}.{r.sta}"; st = seen.get((ns, r.day))
            if str(r.fetch).split(":")[0] not in ("have", "data"): rows.append((r.day, ns, "nodata", str(r.fetch)))
            elif st is None: rows.append((r.day, ns, "missing", "plan has a raw file, the cut did not see it"))
            elif st == "ok": rows.append((r.day, ns, "ok", ""))
            else: rows.append((r.day, ns, "error", st))
    duc = {x for _, x, _, _ in rows}
    for r in SD.itertuples():
        if r.netsta in duc: continue
        rows.append((r.day, r.netsta, "ok" if r.status == "ok" else "error", "" if r.status == "ok" else r.status))
    F.write_status("cut", "0", rows)
    F.guard("cut")


# ------------------------------------------------------------------------------------------------ stack check
def his_chain(x, fs, inv, seed, t0):
    """his get_data.py chain on one raw window: detrend, 5 s hann taper, remove_response VEL pre_filt (0.2,0.5,10,15)
    water level 80, bandpass 1.5-9 zero-phase, resample to 20 sps."""
    from obspy import Trace, UTCDateTime
    from scipy.signal import resample_poly
    from fractions import Fraction
    n, s, l, c = seed.split(".")
    tr = Trace(np.asarray(x, np.float64)); tr.stats.update(dict(network=n, station=s, location=l, channel=c, sampling_rate=fs,
                                                                starttime=UTCDateTime(t0)))
    tr.detrend("linear"); tr.taper(type="hann", max_percentage=None, max_length=5.0)
    tr.remove_response(inventory=inv, output="VEL", pre_filt=(0.2, 0.5, 10.0, 15.0), water_level=80.0)
    tr.filter("bandpass", freqmin=1.5, freqmax=9.0, zerophase=True)
    fr = Fraction(20.0 / fs).limit_denominator(1000)
    y = resample_poly(tr.data, fr.numerator, fr.denominator) if fr != 1 else tr.data
    return y


def ncc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float(np.dot(a, b) / (np.sqrt(np.dot(a, a) * np.dot(b, b)) + 1e-30))


def stackcheck(out):
    from obspy import read_inventory
    I = CW.load_index(out); I["fam"] = I.fam.astype(str); I["month"] = I.h5_file.str[8:15]
    J = I[(I.pass_his) & (~I.partial) & (I.n_gap == 0) & (I.offset_source == "template")]
    T = DC.his_templates(); dur = {r.family: float(r.duration) for r in DC.families().itertuples()}
    g = J.groupby(["fam", "netsta", "month"]).det_id.nunique().reset_index(name="n").sort_values("n", ascending=False)
    res = None
    for b in g.head(5).itertuples():
        sta = b.netsta.split(".")[1]; tp = {c: v for (s2, c), v in T.get(b.fam, {}).items() if s2 == sta}
        if not tp: continue
        K = J[(J.fam == b.fam) & (J.netsta == b.netsta) & (J.month == b.month)]
        inv = read_inventory(f"{DC.RESP}/{b.netsta}.xml")
        L = int(round(dur[b.fam] * 20)) + 1
        stk = {}; nst = {}; single = []
        for det, gg in K.groupby("det_id"):
            per = {}
            for r in gg.to_dict("records"):
                cha = r["cha"]
                if cha not in tp: continue
                x, t0, fs = CW.read_window(r, out)
                try: y = his_chain(x, fs, inv, r["seed_id"], t0)
                except Exception: continue
                tdet = pd.Timestamp(r["t_det"]).value / 1e9
                i0 = int(round((tdet - t0) * 20))           # template first sample inside the window
                seg = np.full(L, np.nan); a0, a1 = max(i0, 0), min(i0 + L, y.size)
                if a1 - a0 < L // 2: continue
                seg[a0 - i0:a1 - i0] = y[a0:a1]; per[cha] = seg
            if not per: continue
            nrm = np.sqrt(np.nanmean(np.concatenate([np.nan_to_num(v) ** 2 for v in per.values()]))) + 1e-30
            for c, v in per.items():
                stk[c] = stk.get(c, 0) + np.nan_to_num(v) / nrm; nst[c] = nst.get(c, 0) + 1
                ok = np.isfinite(v)
                single.append(ncc(np.nan_to_num(v)[ok], tp[c][ok]))
        if not stk: continue
        ccs, ctl = [], []
        for c, s_ in stk.items():
            h = tp[c]; ok = s_ != 0
            best = max(ncc(np.roll(s_, k)[ok], h[ok]) for k in range(-10, 11))
            far = max(abs(ncc(np.roll(s_, k)[ok], h[ok])) for k in list(range(-200, -99)) + list(range(100, 201)))
            ccs.append(best); ctl.append(far)
        res = dict(run=RUN, fam=b.fam, netsta=b.netsta, month=b.month, channels=sorted(stk), n_windows=int(max(nst.values())),
                   cc_stack_vs_his_template=round(float(np.mean(ccs)), 3), cc_single_median=round(float(np.nanmedian(single)), 3),
                   cc_control_shift_5_10s=round(float(np.mean(ctl)), 3),
                   note="stack of the cut raw windows through his chain (response removal, 1.5-9 Hz, 20 sps), cut at the "
                        "detection time (= his template's first sample), compared with HIS template (+-0.5 s); control = the "
                        "same stack shifted by 5-10 s")
        res["pass"] = bool(res["cc_stack_vs_his_template"] >= 0.5 and res["cc_stack_vs_his_template"] > res["cc_control_shift_5_10s"] + 0.2)
        break
    if res is None:
        res = dict(run=RUN, error="no family-station-month with template-offset windows", **{"pass": False})
    json.dump(res, open(f"{out}/STACKCHECK.json", "w"), indent=1)
    log(f"[stackcheck] {res}")
    return res["pass"]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["cut", "verify", "stackcheck", "fixsrc", "guard"])
    ap.add_argument("--workers", type=int, default=6); ap.add_argument("--months", default="")
    ap.add_argument("--max-days", type=int, default=0); ap.add_argument("--n", type=int, default=20)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    setup(); out = F.CUT_OUT; os.makedirs(out, exist_ok=True)
    if a.what == "cut":
        CW.cut(RUN, a.workers, [m for m in a.months.split(",") if m] or None, out, a.max_days)
        fix_offset_source(out)
    elif a.what == "fixsrc": fix_offset_source(out)
    elif a.what == "guard": cutguard(out)
    elif a.what == "verify": sys.exit(0 if CW.verify(RUN, a.n, out) else 5)
    elif a.what == "stackcheck": sys.exit(0 if stackcheck(out) else 6)
