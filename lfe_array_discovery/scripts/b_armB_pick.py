#!/usr/bin/env python
"""ARM B stage 2 -- LFEDetect(`jms`) picks over the 2005 POLARIS bed.  Prefix `b_`.

Merlin's Step-1 audit (Q3, `results/merlin_review_step1_2026-09-20.md`): the installed LFEDetect
weights `cascadia`, `cms`, `jcms` and `jcs` are trained on Bostock et al. 2015 Cascadia stacks, so
a north control with them is a TRAINING-SET test.  **Only `jms` is legitimate** (Guerrero, Frank
et al. 2014 + San Andreas, Shelly 2017 + Nankai, JMA -- no Cascadia), and Washington must later
run the same weights.

Data: the RAW 100 Hz 3-C miniSEED of `traces/bed2005_polaris/`, NOT the frozen 20 Hz 1-8 Hz cache.
The model card is sampling_rate 20, filter bandpass 1-8 Hz zerophase, 1200 samples = 60 s,
norm="peak"; `seisbench.annotate` performs that whole chain itself at the native rate and in the
training-time operation order, so feeding it the already-filtered/decimated cache would apply the
band twice and change the normalisation.  PB boreholes are cached/recorded as (EHZ, EH1, EH2) --
the horizontals are NOT rotated to N/E; seisbench's `flexible_horizontal_components` default
(True) accepts Z12 as ZNE, which is the documented behaviour for rotationally-invariant models.

THRESHOLD LEVELS.  `picks_from_annotations` triggers on prob > t and off at t/2, so the pick set
at threshold t is exactly {picks with thr_first <= t <= prob} (e2_pick's nesting argument).  Both
the model's DEFAULT threshold (0.3, `_annotate_args["*_threshold"]`) and ONE LOWER level (0.15)
are picked and stored with `thr_first`, so either level is reconstructed exactly from the table.

Usage
  python scripts/b_armB_pick.py --out results/armB/picks [--days ets|quiet|all] [--smoke STA:DOY]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
STUDY = os.path.dirname(HERE)
TRACES = f"{STUDY}/traces/bed2005_polaris"
MODEL = "jms"
THR_DEFAULT = 0.3
THR_LOW = 0.15
SWEEP = (THR_LOW, THR_DEFAULT)


def bed():
    cfg = json.load(open(f"{STUDY}/results/step1c_bed.json"))
    ets = pd.date_range(cfg["day0"], cfg["day1"])
    quiet = pd.to_datetime(cfg["quiet_days"])
    return cfg, {int(d.dayofyear) for d in ets}, {int(d.dayofyear) for d in quiet}


def det_stations():
    """The 14 DETECTION stations of the frozen 1c bed (withheld JRBC/VGZ/SNB/LZB excluded)."""
    os.environ.setdefault("W_BED", "c")
    import w_common as W
    return list(W.sites()[3]), list(W.WITHHELD)


def trace_dir(sta):
    m = glob.glob(f"{TRACES}/*.{sta}")
    return m[0] if m else ""


def load_model(device="cuda"):
    import torch
    import seisbench.models as sbm
    m = sbm.LFEDetect.from_pretrained(MODEL, version_str="1")
    dev = device if (device != "cuda" or torch.cuda.is_available()) else "cpu"
    m.to(dev)
    m.eval()
    return m, dev


def model_card(m):
    a = m.get_model_args()
    return dict(name=MODEL, citation=str(getattr(m, "weights_docstring", ""))[:1200],
                sampling_rate=a.get("sampling_rate"), component_order=a.get("component_order"),
                filter_args=a.get("filter_args"), filter_kwargs=a.get("filter_kwargs"),
                in_samples=a.get("in_samples"), norm=a.get("norm"), phases=a.get("phases"),
                default_threshold=THR_DEFAULT, low_threshold=THR_LOW,
                default_args={k: v[1] for k, v in m._annotate_args.items()})


def picks_for_day(model, st, batch_size, overlap):
    ann = model.annotate(st, batch_size=batch_size, overlap=overlap)
    rows = {}
    for phase in ("P", "S"):
        sel = ann.select(channel=f"LFEDetect_{phase}")
        if len(sel) == 0:
            continue
        for thr in SWEEP:                      # ascending: thr_first set on first appearance
            for pk in model.picks_from_annotations(sel, thr, phase):
                key = (phase, round(float(pk.peak_time.timestamp) * 20))
                if key in rows:
                    continue
                rows[key] = dict(phase=phase, time=pk.peak_time.datetime,
                                 prob=float(pk.peak_value), thr_first=float(thr))
    return list(rows.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="results/armB/picks")
    ap.add_argument("--days", default="all", choices=["ets", "quiet", "all"])
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--overlap", type=int, default=600)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--torch-threads", type=int, default=4)
    ap.add_argument("--smoke", default="", help="STA:DOY -- one station-day, nothing written")
    ap.add_argument("--stations", default="", help="comma list; default = the 14 detection stations")
    ap.add_argument("--doys", default="", help="comma list of day-of-year; default = the bed days")
    a = ap.parse_args()
    import torch
    torch.set_num_threads(a.torch_threads)
    from obspy import read

    cfg, doy_ets, doy_q = bed()
    det, withheld = det_stations()
    model, dev = load_model(a.device)
    card = model_card(model)
    print("[model] LFEDetect/%s on %s" % (MODEL, dev))
    print("[model] training regions:", card["citation"].replace("\n", " ")[:300])
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    json.dump(card, open(os.path.join(os.path.dirname(a.out), "jms_model_card.json"), "w"),
              indent=1, default=str)

    if a.smoke:
        sta, doy = a.smoke.split(":")
        p = f"{trace_dir(sta)}/2005/{int(doy):03d}.mseed"
        t0 = time.time()
        st = read(p)
        rows = picks_for_day(model, st, a.batch_size, a.overlap)
        D = pd.DataFrame(rows)
        print(f"[smoke] {sta} doy {doy}: {len(D)} picks in {time.time()-t0:.1f} s "
              f"({[str(tr.id) for tr in st]})")
        if len(D):
            print(D.groupby(["phase", "thr_first"]).size().to_dict())
            print("  n(prob>=%.2f) = %d ; n(prob>=%.2f) = %d"
                  % (THR_LOW, int((D.prob >= THR_LOW).sum()),
                     THR_DEFAULT, int((D.prob >= THR_DEFAULT).sum())))
        return

    keep_doy = doy_ets | doy_q if a.days == "all" else (doy_ets if a.days == "ets" else doy_q)
    if a.doys:
        keep_doy &= {int(x) for x in a.doys.split(",")}
    if a.stations:
        det = [s for s in a.stations.split(",")]
    tot = dict(files=0, picks=0, t=time.time())
    for sta in det:
        d = trace_dir(sta)
        if not d:
            print(f"[skip] {sta}: no trace directory")
            continue
        outp = f"{a.out}_{sta}.parquet"
        if os.path.exists(outp):
            print(f"[skip] {sta}: {outp} exists")
            continue
        files = sorted(p for p in glob.glob(f"{d}/2005/*.mseed")
                       if int(os.path.basename(p).split(".")[0]) in keep_doy)
        rows = []
        t0 = time.time()
        for p in files:
            doy = int(os.path.basename(p).split(".")[0])
            try:
                st = read(p)
            except Exception as e:                                          # noqa: BLE001
                print(f"[warn] {sta} {doy} read failed: {type(e).__name__}: {str(e)[:120]}")
                continue
            try:
                r = picks_for_day(model, st, a.batch_size, a.overlap)
            except Exception as e:                                          # noqa: BLE001
                print(f"[warn] {sta} {doy} annotate failed: {type(e).__name__}: {str(e)[:160]}")
                continue
            for x in r:
                x.update(sta=sta, net=os.path.basename(d).split(".")[0], doy=doy,
                         kind=("ets" if doy in doy_ets else "quiet"))
            rows += r
            tot["files"] += 1
        D = pd.DataFrame(rows, columns=["sta", "net", "doy", "kind", "phase", "time", "prob",
                                        "thr_first"])
        if len(D):
            D["time"] = pd.to_datetime(D.time).astype("datetime64[ns]")
        D.to_parquet(outp, index=False)
        tot["picks"] += len(D)
        nd = max(len(files), 1)
        print(f"[pick] {sta}: {len(files)} station-days, {len(D)} picks "
              f"({len(D)/nd:.0f}/station-day), {time.time()-t0:.0f} s -> {outp}", flush=True)
    print(f"[pick] TOTAL {tot['files']} station-days, {tot['picks']} picks, "
          f"{time.time()-tot['t']:.0f} s")


if __name__ == "__main__":
    main()
