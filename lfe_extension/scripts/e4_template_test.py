#!/usr/bin/env python
"""E4 step 1 — THE DECIDING EXPERIMENT: do Bostock's own catalogued LFEs register against his own
family stacks on the 2010 OUTER-RING network, with the components fixed and an EMPIRICAL moveout?

E3 step 4 answered "no" (0.60 MAD) but it (i) fed the templates' N against the data's E, (ii) used a
single common network lag with a Crosson-predicted moveout whose station-to-station error is
0.3-0.5 s, and (iii) used an 8-s template for a signal coherent over ~2 s.  All three are fixed here.

METHOD
  templates  the family's own SeisBench stack at that station, segment S-PRE .. S-PRE+L, L in
             {2, 4, 8} s, PRE = 0.5 s, each COMPONENT L2-normalised separately (so the station cc
             is the mean of 3 per-channel normalised cc and the network cc is the mean over all
             available channels -- the scaling that makes the sqrt(N_channels) check meaningful).
  placement  origin = Bostock DB time + the per-family offset of E3 step 5;
             template start at station s = origin + Crosson ttP(family hypocentre, s)
                                          + (stack S sample - 800)/20  - PRE
             + (pass 2) the EMPIRICAL moveout correction = the median best single-station lag.
  cc series  full-day continuous normalised correlation per channel, GPU rFFT (the data FFT is
             computed once per station-day and reused by every family and every template length).
  background per (family, station, day) and per (family, day) for the network: the MEDIAN and MAD
             of that template's OWN continuous cc series over that day -- not a best-of-129
             distribution, and not a Gaussian assumption.
  station    value at a Bostock event = max over an INDEPENDENT +-3 s lag scan; its argmax is the
             per-event best lag, and the per-(family, station) median of those is the empirical
             moveout written for pass 2.
  network    mean over all available channels of the per-station 3-C cc, each station shifted by its
             own (predicted + empirical) moveout; value at an event = max over a COMMON +-1 s lag.
  controls   (i) the continuous series itself gives the false alarms per template per day: peaks
             above 6 / 8 MAD, declustered by the template length;
             (ii) N random times scored exactly like the events (max over the same lag scan);
             (iii) a per-station TIME-SHIFTED network: each station's cc series is rolled by its own
             distinct offset drawn in [300, 3000] s before the channels are averaged, which destroys
             inter-station coherence while keeping every station's own cc statistics -> the chance
             level of the network statistic.

ACCEPTANCE (Merlin E3-audit prescription 1): "template route open on our station class" if the
median network MAD at Bostock M >= 1.8 is >= 3 with a 2-4 s template AND the recall of those events
at an 8-MAD network threshold is >= 30 % with shifted-control chance <= 10 %.  FAILURE if the median
is < 1.5 MAD after the fix and the empirical moveout.

Outputs (all e4_-prefixed; nothing of E2/E3 is touched)
  results/e4_tt_events_<tag>.parquet     per (event, template length): station and network MAD
  results/e4_tt_station_<tag>.csv        per (family, station, L): MAD summary + median best lag
  results/e4_tt_network_<tag>.csv        per (L, magnitude bin): median/p90 MAD, fractions, recall
  results/e4_tt_fa_<tag>.csv             false alarms per template per day
  results/e4_moveout_<tag>.csv           the empirical moveout table (pass 1 -> pass 2 input)
"""
from __future__ import annotations
import argparse, os, sys, time, glob, json
import numpy as np
import pandas as pd
import h5py
from scipy.signal import butter, sosfiltfilt

ROOT = "/home/jovyan/tremorferometry/lfe_catalog_repair"
RES = f"{ROOT}/results"
CACHE = f"{ROOT}/cache20hz"
CACHES = [CACHE]
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa: E402

FS = 20.0
PRE = 0.5
MAGBINS = [(-9, 1.5, "M<1.5"), (1.5, 1.8, "M1.5-1.8"), (1.8, 2.1, "M1.8-2.1"), (2.1, 9, "M>=2.1")]
R_EARTH_KM_DEG = 111.195


def bp20(x):
    sos = butter(4, [1.0 / (FS / 2), 8.0 / (FS / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6))
    return np.hypot((np.asarray(lon1) - lon2) * R_EARTH_KM_DEG * cs,
                    (np.asarray(lat1) - lat2) * R_EARTH_KM_DEG)


def load_templates(stations, lengths):
    """-> tpl[(fam, sta)][L] = (3, L) float32 per-component L2-normalised; fam hypocentre; S-P."""
    m = pd.read_csv(f"{SD}/metadata.csv")
    m = m[m.station_code.isin(stations) & m.trace_S_arrival_sample.notna()]
    tpl, meta, nrm = {}, {}, {}
    Lmax = int(round(max(lengths) * FS))
    with h5py.File(f"{SD}/waveforms.hdf5", "r") as f:
        for r in m.itertuples():
            b, idx = r.trace_name.split("$")
            p = idx.split(",")
            w = np.asarray(f["data"][b][int(p[0]), :3, :int(p[2].lstrip(":"))], float)
            w = bp20(w)
            s0 = int(round(r.trace_S_arrival_sample - PRE * FS))
            if s0 < 0 or s0 + Lmax > w.shape[1]:
                continue
            k = (float(r.source_family_id), r.station_code)
            ok = True
            d = {}
            for L in lengths:
                n = int(round(L * FS))
                seg = w[:, s0:s0 + n].astype(np.float32)
                nn = np.linalg.norm(seg, axis=1)
                if np.any(nn <= 0):
                    ok = False
                    break
                d[L] = (seg / nn[:, None]).astype(np.float32)
            if not ok:
                continue
            tpl[k] = d
            meta.setdefault(float(r.source_family_id),
                            dict(lat=float(r.source_latitude_deg), lon=float(r.source_longitude_deg),
                                 dep=float(r.source_depth_km)))
            meta[float(r.source_family_id)][f"sp_{r.station_code}"] = \
                float(r.trace_S_arrival_sample - 800.0) / FS
            nrm[k] = float(r.trace_S_snr_db) if np.isfinite(r.trace_S_snr_db) else np.nan
    return tpl, meta, nrm


def decluster_count(x, thr, n, valid):
    """Peaks above `thr` (per row) that are the max of a window of n samples -> counts per row."""
    import torch
    import torch.nn.functional as Fnn
    b = (x >= thr[:, None]) & valid[None, :]
    if not bool(b.any()):
        return np.zeros(x.shape[0], np.float64)
    mx = Fnn.max_pool1d(x.unsqueeze(1), kernel_size=2 * (n // 2) + 1, stride=1,
                        padding=n // 2).squeeze(1)
    return ((x >= mx - 1e-12) & b).float().sum(1).cpu().numpy()


def decluster_idx_multi(x, thrs, n, valid):
    """E5: declustered peak INDICES of a 1-D series for several thresholds (max-pool computed once).
    x (N,) tensor, thrs (nT,) tensor, valid (N,) bool -> list of int numpy arrays, one per threshold."""
    import torch
    import torch.nn.functional as Fnn
    mx = Fnn.max_pool1d(x[None, None], kernel_size=2 * (n // 2) + 1, stride=1,
                        padding=n // 2)[0, 0]
    ispk = (x >= mx - 1e-12) & valid
    return [torch.nonzero(ispk & (x >= t)).squeeze(1).cpu().numpy() for t in thrs]


def min_absdiff(a, b_sorted):
    """E5: for each a, the distance to the nearest element of the sorted array b (inf if b empty)."""
    a = np.asarray(a, float)
    if len(a) == 0:
        return np.zeros(0)
    if len(b_sorted) == 0:
        return np.full(len(a), np.inf)
    j = np.searchsorted(b_sorted, a)
    lo = b_sorted[np.clip(j - 1, 0, len(b_sorted) - 1)]
    hi = b_sorted[np.clip(j, 0, len(b_sorted) - 1)]
    return np.minimum(np.abs(a - lo), np.abs(a - hi))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="PGC,LZB,VGZ,SNB,YOUB,PFB")
    ap.add_argument("--netsta", default="PGC=CN.PGC,LZB=CN.LZB,VGZ=CN.VGZ,SNB=CN.SNB,"
                                       "YOUB=CN.YOUB,PFB=CN.PFB")
    ap.add_argument("--year", type=int, default=2010)
    ap.add_argument("--lengths", default="2,4,6,8,12")
    ap.add_argument("--tag", default="ring2010_p1")
    ap.add_argument("--moveout-in", default="")
    ap.add_argument("--sta-lag", type=float, default=3.0)
    ap.add_argument("--net-lag", type=float, default=1.0)
    ap.add_argument("--n-random", type=int, default=400, help="random control times per day")
    ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--max-days", type=int, default=0)
    ap.add_argument("--doy-from", type=int, default=0)
    ap.add_argument("--doy-to", type=int, default=999)
    # AMENDMENT 2026-09-11 (Merlin vet F2): the Aug-build / Sep-test split is impossible (the 2010
    # ETS is over by Sep 1; top families have 235-327 August events and 1-35 September events).
    # Hold-out = ALTERNATE 2-HOUR BLOCKS with a 60-s guard at every block edge.  Never by parity
    # (35 % of events have a same-family neighbour within 30 s).
    ap.add_argument("--half", default="all", choices=["all", "build", "test"])
    ap.add_argument("--block-s", type=float, default=7200.0)
    ap.add_argument("--guard-s", type=float, default=60.0)
    ap.add_argument("--sta-weight", default="", help="net.sta=nchan overrides, e.g. STW=1")
    ap.add_argument("--extra-tpl", default="", help="npz of extra (fam,sta) templates (step 2)")
    ap.add_argument("--extra-tpl-stations", default="",
                    help="restrict --extra-tpl to these stations (step 2: the ring 2010 re-stacks "
                         "and the transferred stations are the same file but different legs)")
    ap.add_argument("--drop-stations", default="",
                    help="drop these stations' SeisBench templates (step 2 q-control: replace the "
                         "ring's published stacks with our own 2010 re-stacks)")
    # E5 (2026-09-11, Merlin E4-audit F1/F11) -- DETECTION-BUDGET and SNR-WEIGHT extensions, both
    # OFF by default so every E4 run is reproduced bit-for-bit without them.
    #   --budget      score the network series of EVERY family with a template on EVERY day (E4 only
    #                 built the network for families with a Bostock event that day, which biases a
    #                 quiet-day false-alarm count), and record per (day, L, family, threshold) the
    #                 declustered detections, the rolled-control detections, the detections matched
    #                 within +-2 s to a Bostock event of the SAME family / of ANY family, the rolled
    #                 detections matched to the same family (= chance matching), and the 8-MAD extras
    #                 within 60 s of a NEIGHBOURING family's event vs the same for random times.
    #                 Also adds a +-2 s common-lag network score per event (net_mad_2s).
    #   --weights     csv (fam, sta, w): per-(family, station) weights multiplying the channel count
    #                 in the network sum AND in the rolled control (F11: weights from build-half stack
    #                 S-SNR only, never from test-half registration). Missing pairs -> 1.0.
    ap.add_argument("--budget", action="store_true")
    ap.add_argument("--budget-thr", default="5,6,7,8,10")
    ap.add_argument("--weights", default="")
    # E5 STEP 5 (2026-09-11, the 2012 from-scratch replication) -- all OFF by default (E4/E5 runs bit-for-bit unchanged):
    #   --no-seisbench     no SeisBench stack is used as a template; family hypocentres from --families-csv (Bostock's own
    #                      families keyed by bkey = fam_id); every template comes from --extra-tpl
    #   --stations-csv / --gt-parquet / --offsets-csv   the 2012 station table, ground truth and per-family offsets (same columns)
    #   --skip-z           do not build the Z-only cc series (memory: 15 stations x ~110 families x 1.7 M samples); staZ_* = sta_*
    #   --catalog-out      in --budget mode also write EVERY declustered detection at the lowest --budget-thr with its network
    #                      MAD, raw network cc, SNR-weighted relative amplitude (units of each station's quiet-day window rms,
    #                      --amp-ref-csv: sta, L, e3_ref), station count and +-2 s matches -> results/e5r_detections_<tag>.parquet
    ap.add_argument("--no-seisbench", action="store_true")
    ap.add_argument("--families-csv", default="")
    ap.add_argument("--stations-csv", default="")
    ap.add_argument("--gt-parquet", default="")
    ap.add_argument("--offsets-csv", default="")
    ap.add_argument("--skip-z", action="store_true")
    ap.add_argument("--catalog-out", action="store_true")
    ap.add_argument("--amp-ref-csv", default="")
    # E6 (2026-09-11, Merlin E5r-audit prescriptions 1-2; results/merlin_audit_e5r_2026-09-11.md) -- all OFF by default,
    # every E4/E5/E5r run is reproduced bit-for-bit without them (the site shifts use their OWN rng so the main stream is untouched):
    #   --sites        SITE groups 'PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007' (stations within ~1 km); every other station
    #                  is its own site.  Adds the SITE-ROLLED control (ONE shift per site, distinct 300-3000 s) next to the
    #                  station-rolled one (F1: the per-station roll destroys within-site AND between-site coherence): budget
    #                  columns n_det_roll_site / n_det_roll_site_test / n_roll_site_match_own, event column net_mad_shift_site,
    #                  netfam shiftsite_fa6/8; with --catalog-out every detection also carries the per-SITE contribution to the
    #                  weighted network sum (con_<site>, sum = net_cc) and the site's weight share (wsh_<site>).
    #   --site-merge   step 2: within a site the members' cc series are w-averaged into ONE series (same per-member moveout), the
    #                  site enters the sum as one "station" with site weight A = max member w (fully coherent noise) or
    #                  B = sqrt(sum w^2) (partially coherent).  Singleton sites are unchanged.  Needs --sites.
    #   --site-cap     step-2 fallback: 'PGC=0.25' caps that site's share of the total weight (merge only).
    #   --coinc-gate   step-2 fallback: 'nsite,z' keeps a detection (data AND both controls) only if >= nsite sites are each
    #                  >= z MAD of their own site series at the detection sample (merge only).
    #   --only-fams    restrict the templates to these families (step 1b: the 5 worst quiet-day families).
    ap.add_argument("--sites", default="")
    ap.add_argument("--site-merge", default="", choices=["", "A", "B"])
    ap.add_argument("--site-cap", default="")
    ap.add_argument("--coinc-gate", default="")
    ap.add_argument("--only-fams", default="")
    # lfe_extension Track A (2026-09-16).  --out-dir: every output file goes here (lfe_catalog_repair stays READ-ONLY input).
    # --min-sites N (Merlin plan-review F8): a family-day is SCORED only if at least N distinct SITES have a usable template AND
    # data that day; with --sites PGC:PGC,B009,B010,B011;LRIV:LRIV,B005,B006,B007 the four PGC-cluster sensors count as ONE site
    # and the four LRIV-cluster sensors as one, so a family carried by the PGC borehole cluster alone is now UNOBSERVABLE rather
    # than scored off a single travel-time.  n_sites is recorded on every budget row and every detection either way.
    ap.add_argument("--out-dir", default="")
    # lfe_extension: extra 20 Hz cache roots searched BEFORE the lfe_catalog_repair cache (first hit wins), so a station-year
    # cached here (e.g. PB.B009/B010 2010, downloaded by Track B) joins the ones already cached in lfe_catalog_repair.
    ap.add_argument("--cache-dirs", default="")
    ap.add_argument("--min-sites", type=int, default=0)
    a = ap.parse_args()
    E6 = bool(a.sites or a.site_merge)
    global RES, CACHES
    if a.out_dir:
        RES = a.out_dir.rstrip("/")
        os.makedirs(RES, exist_ok=True)
        print(f"[out] every output of this run goes to {RES}", flush=True)
    CACHES = [d.rstrip("/") for d in a.cache_dirs.split(",") if d] + [CACHE]
    if len(CACHES) > 1:
        print(f"[cache] search order {CACHES}", flush=True)
    if a.site_merge and not a.sites:
        print("[E6] WARNING --site-merge without --sites: every site is a single station -> identical to the baseline", flush=True)
    GATE = None
    if a.coinc_gate:
        if not a.site_merge:
            raise SystemExit("--coinc-gate needs --site-merge")
        GATE = (int(a.coinc_gate.split(",")[0]), float(a.coinc_gate.split(",")[1]))
    rng_site = np.random.default_rng(6)              # E6: separate stream for the site shifts
    BTHR = [float(x) for x in a.budget_thr.split(",")]
    n_unobs = 0                                   # F8: family-days skipped by --min-sites (unobservable)

    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.set_num_threads(4)
    print(f"[env] {dev} {torch.cuda.get_device_name(0) if dev=='cuda' else ''}", flush=True)

    STAS = a.stations.split(",")
    NETSTA = dict(kv.split("=") for kv in a.netsta.split(","))
    NCH = {s_: 3 for s_ in STAS}   # overridden per station by --sta-weight (e.g. STW=1, OCP=1)
    if a.sta_weight:
        for kv in a.sta_weight.split(","):
            k_, v_ = kv.split("="); NCH[k_] = int(v_)
    LENS = [float(x) for x in a.lengths.split(",")]
    Lsamp = {L: int(round(L * FS)) for L in LENS}
    Lmax = max(Lsamp.values())
    # E6: station -> site map (default: every station is its own site)
    SITE = {s_: s_ for s_ in STAS}
    if a.sites:
        for grp in a.sites.split(";"):
            q_, mem_ = grp.split(":")
            for s_ in mem_.split(","):
                if s_ in SITE:
                    SITE[s_] = q_
        print(f"[E6] sites: " + "; ".join(f"{q}={[s for s in STAS if SITE[s]==q]}" for q in sorted(set(SITE.values())))
              + (f"; merge variant {a.site_merge}" if a.site_merge else "; no merge (control only)")
              + (f"; cap {a.site_cap}" if a.site_cap else "") + (f"; coincidence gate {GATE}" if GATE else ""), flush=True)

    if a.no_seisbench:
        # E5 step 5: FROM-SCRATCH mode -- meta seeded from Bostock's own 2012 families (bkey), templates only via --extra-tpl
        tpl, meta, snr = {}, {}, {}
        for r in pd.read_csv(a.families_csv).itertuples():
            meta[float(r.bkey)] = dict(lat=float(r.lat), lon=float(r.lon), dep=float(r.dep))
        print(f"[tpl] --no-seisbench: {len(meta)} families from {os.path.basename(a.families_csv)}", flush=True)
    else:
        tpl, meta, snr = load_templates(STAS, LENS)
    if a.budget and not a.no_seisbench:
        # E5 --budget (bugfix run 1, 2026-09-11): the network loop visits every family with ANY template on the
        # selected stations, including families whose only SeisBench stacks sit on stations NOT selected (e.g.
        # LZB/YOUB/PFB for the operating subset) but which have a transferred template -> their hypocentre must
        # be known for the placement.  Seed meta with every SeisBench family's hypocentre (KeyError (6.0,'B011')).
        m_all = pd.read_csv(f"{SD}/metadata.csv")
        m_all = m_all[m_all.trace_S_arrival_sample.notna()].groupby("source_family_id").first()
        for fid, r in m_all.iterrows():
            meta.setdefault(float(fid), dict(lat=float(r.source_latitude_deg), lon=float(r.source_longitude_deg),
                                             dep=float(r.source_depth_km)))
    if a.drop_stations:
        dr = set(a.drop_stations.split(","))
        tpl = {k: v for k, v in tpl.items() if k[1] not in dr}
        print(f"[tpl] dropped the published stacks of {sorted(dr)} -> {len(tpl)} left", flush=True)
    if a.extra_tpl:
        z = np.load(a.extra_tpl, allow_pickle=True)
        ex = z["tpl"].item()          # {(fam,sta): {L: (nc,Lsamp)}}, already normalised
        exsp = z["sp"].item()         # {(fam,sta): the S-P to use in the placement formula}
        keep = set(a.extra_tpl_stations.split(",")) if a.extra_tpl_stations else None
        nadd = 0
        for k, v in ex.items():
            if keep is not None and k[1] not in keep:
                continue
            if not all(L in v for L in LENS):
                continue
            tpl[k] = v
            nadd += 1
            if k[0] in meta:
                meta[k[0]][f"sp_{k[1]}"] = exsp[k]
        print(f"[tpl] + {nadd} templates from {os.path.basename(a.extra_tpl)}"
              + (f" (stations {sorted(keep)})" if keep else ""), flush=True)
    if a.only_fams:
        keepf = {float(x) for x in a.only_fams.split(",")}
        tpl = {k: v for k, v in tpl.items() if k[0] in keepf}
        print(f"[E6] --only-fams {sorted(keepf)}: {len(tpl)} templates kept", flush=True)
    stas_tbl = pd.read_csv(a.stations_csv or f"{RES}/e2_stations.csv").set_index("sta")   # E5 step 5: --stations-csv
    print(f"[tpl] {len(tpl)} (family, station) templates over {len({k[1] for k in tpl})} stations, "
          f"{len(meta)} families", flush=True)

    # ---- Bostock events with an accepted per-family time-base offset
    bo = pd.read_parquet(a.gt_parquet or f"{RES}/e2_gt_bostock.parquet")          # E5 step 5: --gt-parquet
    fo = pd.read_csv(a.offsets_csv or f"{RES}/e3_bostock_family_offsets.csv")     # E5 step 5: --offsets-csv
    acc = fo[fo.accept.astype(bool)]
    omap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.bkey.astype(int).isin(omap)].copy()
    bo["time"] = bo.time.astype("datetime64[ns]")
    bo["t_origin"] = (bo.time.astype("int64") / 1e9 + bo.bkey.astype(int).map(omap).to_numpy())
    bo = bo[bo.fam_id.isin(meta.keys())].copy()
    blk = np.floor(bo.t_origin.to_numpy() / a.block_s).astype(np.int64)
    pos_in = bo.t_origin.to_numpy() - blk * a.block_s
    guard = (pos_in < a.guard_s) | (pos_in > a.block_s - a.guard_s)
    bo["half"] = np.where(blk % 2 == 0, "build", "test")
    bo["guard"] = guard
    print(f"[evt] {len(bo)} Bostock events, {bo.fam_id.nunique()} families with a template "
          f"and an accepted offset; blocks {a.block_s:.0f} s -> "
          f"build {int((bo.half=='build').sum())} / test {int((bo.half=='test').sum())}, "
          f"{int(guard.sum())} in the {a.guard_s:.0f}-s edge guard", flush=True)
    if a.half != "all":
        bo = bo[(bo.half == a.half) & (~bo.guard)].copy()
        print(f"[evt] evaluating on the {a.half} half only: {len(bo)} events", flush=True)
    # E5 budget: EVERY Bostock event (all families, both halves) for the "matched to ANY family"
    # and "extras near a NEIGHBOURING family" counts; families without an accepted offset get the
    # median accepted offset (their times are only used at the +-2 s / 60 s matching level).
    bo_all = pd.read_parquet(a.gt_parquet or f"{RES}/e2_gt_bostock.parquet")      # E5 step 5: --gt-parquet
    bo_all["time"] = bo_all.time.astype("datetime64[ns]")
    med_off = float(np.median(list(omap.values())))
    bo_all_t = (bo_all.time.astype("int64") / 1e9
                + bo_all.bkey.astype(int).map(omap).fillna(med_off).to_numpy()).to_numpy()
    bo_all_f = bo_all.fam_id.to_numpy(float)          # NaN for families outside SeisBench
    WGT = {}
    if a.weights:
        Wt = pd.read_csv(a.weights)
        WGT = {(float(r.fam), str(r.sta)): float(r.w) for r in Wt.itertuples()}
        print(f"[weights] {len(WGT)} (family, station) weights from {os.path.basename(a.weights)}:"
              f" median {np.median(list(WGT.values())):.2f}, max {max(WGT.values()):.2f}, "
              f"min {min(WGT.values()):.2f}; missing pairs -> 1.0", flush=True)
    NLAG2 = int(round(2.0 * FS))
    AMPREF = {}
    if a.amp_ref_csv:
        AMPREF = {(str(r.sta), float(r.L)): float(r.e3_ref) for r in pd.read_csv(a.amp_ref_csv).itertuples()}
    cat_rows = []                                     # E5 step 5 --catalog-out

    # ---- predicted template-start offset from origin, per (family, station)
    off = {}
    for fid, mm in meta.items():
        for s in STAS:
            if (fid, s) not in tpl or f"sp_{s}" not in mm:
                continue
            dx = epi_km(mm["lat"], mm["lon"], float(stas_tbl.loc[s, "lat"]),
                        float(stas_tbl.loc[s, "lon"]))
            off[(fid, s)] = float(v1d.tt(dx, mm["dep"], "P")) + mm[f"sp_{s}"] - PRE
    mv = {}
    if a.moveout_in:
        M = pd.read_csv(a.moveout_in)
        for r in M.itertuples():
            if np.isfinite(r.lag_median):
                mv[(float(r.fam), r.sta, float(r.L))] = float(r.lag_median)
        print(f"[moveout] {len(mv)} (family, station, L) corrections loaded from "
              f"{os.path.basename(a.moveout_in)}: median {np.median(list(mv.values())):+.2f} s, "
              f"IQR {np.percentile(list(mv.values()),25):+.2f}..{np.percentile(list(mv.values()),75):+.2f}",
              flush=True)

    days = sorted({int(os.path.basename(f)[:3])
                   for s in STAS
                   for c in CACHES
                   for f in glob.glob(f"{c}/{NETSTA[s]}/{a.year}/*.npz")})
    days = [d for d in days if a.doy_from <= d <= a.doy_to]
    if a.max_days:
        days = days[:a.max_days]
    print(f"[data] {len(days)} days in {a.year}: {days[0]}..{days[-1]}", flush=True)

    rng = np.random.default_rng(4)
    SLAG = int(round(a.sta_lag * FS))
    NLAG = int(round(a.net_lag * FS))
    ev_rows, sta_rows, fa_rows = [], [], []
    bud_rows = []                                     # E5 --budget
    t_start = time.time()

    for doy in days:
        # ---------- load the day ----------
        D, T0, NPT = {}, {}, {}
        for s in STAS:
            p = next((q for q in (f"{c}/{NETSTA[s]}/{a.year}/{doy:03d}.npz" for c in CACHES) if os.path.exists(q)), None)
            if p is None:
                continue
            z = np.load(p)
            D[s] = torch.from_numpy(z["x"]).to(dev)
            T0[s] = float(z["t0"])
            NPT[s] = D[s].shape[1]
        if not D:
            continue
        Nmin = min(NPT.values())
        NF = int(2 ** np.ceil(np.log2(Nmin + Lmax)))
        DH, EN = {}, {}
        for s in D:
            x = D[s][:, :Nmin]
            DH[s] = torch.fft.rfft(x, n=NF)
            # BUGFIX 2026-09-11 (E4 step 1): the prefix sum of squares MUST be float64.  In float32
            # the running total for a station with counts ~1e4 reaches ~3e14, whose ulp is ~2e7, so
            # differencing two prefix sums for a quiet 2-s window (energy ~1e9, but far less in a
            # quiet minute) suffers catastrophic cancellation and can come out ~0 or negative ->
            # cc = num/1e-10 ~ 1e12 and the MAD scale explodes (seen as mad_scale 7e11 at PGC, 4e11
            # at SNB, 1e10 at YOUB, while LZB/VGZ, which are 1000x quieter, were fine).
            xd = x.double()
            c = torch.cumsum(torch.cat([torch.zeros((x.shape[0], 1), device=dev,
                                                     dtype=torch.float64), xd * xd], 1), 1)
            EN[s] = c                                          # prefix sums of squares, float64
        # ---------- per template length ----------
        for L in LENS:
            n = Lsamp[L]
            NV = Nmin - n                                      # valid template-start indices
            # station 3-C and Z-only cc series, per family
            fam_list = sorted({k[0] for k in tpl if k[1] in D})
            S3, SZ, VALID, E3W = {}, {}, {}, {}
            for s in D:
                fs_ = [f for f in fam_list if (f, s) in tpl]
                if not fs_:
                    continue
                out3 = torch.empty((len(fs_), NV), device=dev)
                outz = out3 if a.skip_z else torch.empty((len(fs_), NV), device=dev)   # E5 step 5 --skip-z: alias, no Z-only series
                ew = (EN[s][:, n:Nmin + 1] - EN[s][:, :NV + 1])[:, :NV]          # (3, NV) float64
                efloor = 1e-10 * ew.median()
                e = ew.clamp_min(efloor).sqrt().float()
                for i0 in range(0, len(fs_), a.fam_chunk):
                    ff = fs_[i0:i0 + a.fam_chunk]
                    W = torch.from_numpy(np.stack([tpl[(f, s)][L] for f in ff])).to(dev)  # (nf,nc,n)
                    WH = torch.fft.rfft(W, n=NF)
                    num = torch.fft.irfft(WH.conj() * DH[s][None], n=NF)[:, :, :NV]
                    cc = num / e[None]
                    out3[i0:i0 + len(ff)] = cc.mean(1)
                    if not a.skip_z:
                        outz[i0:i0 + len(ff)] = cc[:, 0]
                S3[s] = (fs_, out3)
                SZ[s] = outz
                e3 = ew.sum(0)                            # 3-C window energy, (NV,) float64
                VALID.setdefault(s, {})[L] = (e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0)
                if a.catalog_out:
                    E3W[s] = e3.float()                     # E5 step 5: 3-C window energy for the relative amplitude
                    if (s, L) not in AMPREF:
                        AMPREF[(s, L)] = float(e3.median().item())   # fallback: the day's own median (labelled in the ledger)
                del ew, e3
            # ---------- station-level scoring ----------
            evd = bo[(bo.t_origin >= min(T0.values())) &
                     (bo.t_origin <= min(T0[s] + NPT[s] / FS for s in D))]
            evd = evd[(evd.t_origin >= min(T0.values()) + 60) &
                      (evd.t_origin <= min(T0[s] + Nmin / FS for s in D) - 120)]
            # BUGFIX 2026-09-11 (E5 step 5, quiet-day FA): a day with NO catalogued event was skipped entirely, so in --budget
            # mode the 15 pre-SSE quiet days of 2012 (0 Bostock events) produced no detections and no budget rows (E5's 2010
            # quiet days had <= 50 events, so this never surfaced).  In --budget mode the day is processed with an empty event set.
            if not len(evd) and not a.budget:
                continue
            fam_ev = evd.fam_id.to_numpy()
            t_ev = evd.t_origin.to_numpy()
            rt = rng.uniform(min(T0.values()) + 60,
                             min(T0[s] + Nmin / FS for s in D) - 120, a.n_random)
            store = {}                        # (sta) -> dict of per-event arrays
            for s in D:
                if s not in S3:
                    continue
                fs_, cc3 = S3[s]
                fidx = {f: i for i, f in enumerate(fs_)}
                # FIX 2026-09-11 (E4): gaps merged with fill_value=0 give a window energy of ~0, an
                # undefined cc and a MAD of ~0 -> MAD scores of 1e11 (seen at PFB, whose 2010 day
                # files are fragmented into 231 traces).  Mask every sample whose 3-C window energy
                # is below 1e-6 of the day's median before ANY median/MAD/extremum is taken.
                good_t = VALID[s][L]
                cc3m = cc3[:, good_t]
                med = cc3m.median(1).values
                mad = (cc3m - med[:, None]).abs().median(1).values.clamp_min(1e-6)
                szm = SZ[s][:, good_t]
                medz = szm.median(1).values
                madz = (szm - medz[:, None]).abs().median(1).values.clamp_min(1e-6)
                del cc3m, szm
                # event positions
                # dtypes explicit (E5 step 5): on an event-free day these lists are EMPTY and np.array([]) is float64,
                # which made `sel & (...)` raise TypeError in the first --budget re-run over the quiet days.
                sel = np.array([f in fidx for f in fam_ev], dtype=bool)
                idxf = np.array([fidx.get(f, 0) for f in fam_ev], dtype=np.int64)
                o = np.array([off.get((f, s), np.nan) for f in fam_ev], dtype=float)
                o = o + np.array([mv.get((f, s, L), 0.0) for f in fam_ev])
                pos = np.rint((t_ev + o - T0[s]) * FS).astype(np.int64)
                good = sel & np.isfinite(o) & (pos - SLAG >= 0) & (pos + SLAG + 1 <= NV)
                res3 = np.full(len(evd), np.nan); resZ = np.full(len(evd), np.nan)
                lag = np.full(len(evd), np.nan)
                if good.any():
                    gi = np.where(good)[0]
                    ii = torch.from_numpy(idxf[gi]).to(dev)
                    base = torch.from_numpy(pos[gi] - SLAG).to(dev)
                    cols = base[:, None] + torch.arange(2 * SLAG + 1, device=dev)[None]
                    win3 = cc3[ii[:, None], cols]
                    v3, l3 = win3.max(1)
                    winz = SZ[s][ii[:, None], cols]
                    vz, _ = winz.max(1)
                    res3[gi] = ((v3 - med[ii]) / mad[ii]).cpu().numpy()
                    resZ[gi] = ((vz - medz[ii]) / madz[ii]).cpu().numpy()
                    lag[gi] = (l3.cpu().numpy() - SLAG) / FS
                # random control at the same scan width, same families
                rpos = {}
                rr3 = np.full((len(fs_), a.n_random), np.nan, np.float32)
                rrz = np.full((len(fs_), a.n_random), np.nan, np.float32)
                rb = np.rint((rt - T0[s]) * FS).astype(np.int64)
                okr = (rb - SLAG >= 0) & (rb + SLAG + 1 <= NV)
                if okr.any():
                    colsr = torch.from_numpy(rb[okr] - SLAG).to(dev)[:, None] + \
                        torch.arange(2 * SLAG + 1, device=dev)[None]
                    wr = cc3[:, colsr]                        # (nf, nr, 2S+1)
                    vr = wr.max(2).values
                    rr3[:, np.where(okr)[0]] = ((vr - med[:, None]) /
                                                mad[:, None]).cpu().numpy()
                    # the SINGLE-CHANNEL (Z) random-time background, needed for the channel-gain
                    # check: the max-over-scan background of a 1-channel cc series is NOT the same
                    # as that of the 3-C mean, so the Z excess must be referenced to its own.
                    vz_ = SZ[s][:, colsr].max(2).values
                    rrz[:, np.where(okr)[0]] = ((vz_ - medz[:, None]) /
                                                madz[:, None]).cpu().numpy()
                # false alarms from the continuous series: declustered peaks per day
                # false alarms per template per day = DECLUSTERED peaks above the threshold in the
                # day's own continuous series: a sample counts only if it is the maximum of a window
                # one template long centred on it (1-8 Hz cc series oscillate, so a rising-edge count
                # over-counts one excursion several times).  Masked (gap) samples excluded.
                nfa6 = decluster_count(cc3, med + 6 * mad, n, good_t)
                nfa8 = decluster_count(cc3, med + 8 * mad, n, good_t)
                store[s] = dict(res3=res3, resZ=resZ, lag=lag)
                for i, f in enumerate(fs_):
                    ok = np.isfinite(res3) & (fam_ev == f)
                    sta_rows.append(dict(
                        doy=doy, L=L, fam=f, sta=s, n_ev=int(ok.sum()),
                        mad3_median=float(np.nanmedian(res3[ok])) if ok.any() else np.nan,
                        madZ_median=float(np.nanmedian(resZ[ok])) if ok.any() else np.nan,
                        lag_median=float(np.nanmedian(lag[ok])) if ok.any() else np.nan,
                        mad_scale=float(mad[i].item()), med_scale=float(med[i].item()),
                        rand3_median=float(np.nanmedian(rr3[i])),
                        rand3_p99=float(np.nanpercentile(rr3[i], 99)),
                        randZ_median=float(np.nanmedian(rrz[i])),
                        randZ_p99=float(np.nanpercentile(rrz[i], 99)),
                        fa6_per_day=float(nfa6[i] * 86400.0 * FS / Nmin),
                        fa8_per_day=float(nfa8[i] * 86400.0 * FS / Nmin),
                        n_valid_s=float(Nmin / FS)))
            # ---------- network ----------
            nres = np.full(len(evd), np.nan)
            nres_shift = np.full(len(evd), np.nan)
            nres_shift_site = np.full(len(evd), np.nan)   # E6: SITE-rolled network score at the events
            nres2 = np.full(len(evd), np.nan)          # E5: +-2 s common-lag score
            nch = np.zeros(len(evd), int)
            netfam = []
            rndn = []
            shift_s = {s: int(rng.uniform(300, 3000) * FS) for s in D}
            # E6 step 1 (Merlin E5r-audit F1/F4): ONE shift per SITE, shared by all its stations, so this control destroys
            # BETWEEN-site coherence only and keeps the within-site (co-located sensor) coherence that the per-station roll
            # above destroys.  Drawn from a separate rng (the main stream, hence every E4/E5 number, is unchanged).
            shift_site = {q: int(rng_site.uniform(300, 3000) * FS) for q in sorted({SITE[s] for s in D})} if E6 else {}
            # E5 --budget: every family with a template on a station that has data today, not only
            # the families with a catalogued event today (quiet-day FA must count them all).
            fam_iter = sorted(set(fam_ev)) if not a.budget else \
                sorted({k[0] for k in tpl if k[1] in D and k[1] in S3})
            for f in fam_iter:
                ss = [s for s in D if (f, s) in tpl and s in S3 and (f, s) in off]   # E5: needs a placement
                if not ss:
                    continue
                n_sites_f = len({SITE[s] for s in ss})       # F8: distinct SITES with a usable template AND data today
                if a.min_sites and n_sites_f < a.min_sites:
                    n_unobs += 1
                    continue
                o = {s: off[(f, s)] + mv.get((f, s, L), 0.0) for s in ss}
                # A station's cc index j means "template start at sample j", i.e. absolute time
                # T0[s] + j/FS, i.e. ORIGIN time T0[s] + j/FS - o[s].  Different stations can have
                # different day-file start times, so build the network on a COMMON origin axis:
                # network index m <-> origin tref + m/FS, station index = k[s] + m with
                # k[s] = round((tref - T0[s] + o[s]) * FS), shifted up by b so every k[s] >= 0.
                tref0 = min(T0[s] for s in ss)
                k = {s: int(round((tref0 - T0[s] + o[s]) * FS)) for s in ss}
                b = max(0, -min(k.values()))
                k = {s: k[s] + b for s in ss}
                tref = tref0 + b / FS
                NVo = NV - max(k.values())
                if NVo <= 1000:
                    continue
                # AMENDMENT (Merlin vet F4): threshold the NETWORK SUM normalised by the number of
                # CHANNELS contributing that day, never per-station thresholds.  station3c is the
                # mean of that station's channels, so the channel-weighted mean is
                #   sum_s w_s * station3c_s / sum_s w_s   with w_s = that station's channel count.
                acc_n = torch.zeros(NVo, device=dev)
                acc_s = torch.zeros(NVo, device=dev)
                wtot = 0.0
                weff = {}                                 # E6: effective per-station weight in the sum (amplitude / contributions)
                site_ser, site_stats, groups, Wsite = {}, {}, {}, {}
                if not E6:
                    for s in ss:
                        fs_, cc3 = S3[s]
                        i = fs_.index(f)
                        # E5: w = channel count x SNR weight (1.0 unless --weights); the rolled control
                        # below uses the SAME w, so its false-alarm rate is that of the weighted sum.
                        w = float(NCH.get(s, 3)) * float(WGT.get((f, s), 1.0))
                        acc_n += w * cc3[i, k[s]:k[s] + NVo]
                        acc_s += w * torch.roll(cc3[i, k[s]:k[s] + NVo], shift_s[s] % NVo)
                        wtot += w
                    acc_n /= wtot; acc_s /= wtot
                    vmask = torch.ones(NVo, dtype=torch.bool, device=dev)
                    vmask_s = torch.ones(NVo, dtype=torch.bool, device=dev)
                    for s in ss:
                        vmask &= VALID[s][L][k[s]:k[s] + NVo]
                        # BUGFIX 2026-09-11 (E5 step 1, run 1): the ROLLED control must roll each station's gap
                        # mask WITH its series -- with the unrolled mask, a station's gap (cc undefined, huge)
                        # lands on samples the mask calls valid and the "chance" rate explodes (50.8 detections
                        # /template/day at 8 MAD on doy 251, where LRIV has a 0.2 % gap, vs 0.5-0.9 on other quiet
                        # days).  Applied in E5 (--budget / --weights) runs; E4's ring runs had no gaps at the
                        # stations used, leg C (LRIV/STW) may have been mildly inflated.
                        vmask_s &= torch.roll(VALID[s][L][k[s]:k[s] + NVo], shift_s[s] % NVo)
                    if not (a.budget or a.weights):
                        vmask_s = vmask                       # E4 behaviour, bit-for-bit
                else:
                    # E6: the SAME weighted sum built site by site (no merge: bit-for-bit the baseline up to float32 summation
                    # order), plus the SITE-rolled control acc_ss (one shift per site) and, with --site-merge, the within-site
                    # w-weighted average of the members' series entering the sum as ONE station with weight Wsite.
                    acc_ss = torch.zeros(NVo, device=dev)
                    vmask = torch.ones(NVo, dtype=torch.bool, device=dev)
                    vmask_s = torch.ones(NVo, dtype=torch.bool, device=dev)
                    vmask_ss = torch.ones(NVo, dtype=torch.bool, device=dev)
                    for s in ss:
                        groups.setdefault(SITE[s], []).append(s)
                    wm = {s: float(NCH.get(s, 3)) * float(WGT.get((f, s), 1.0)) for s in ss}
                    for q, mem in groups.items():
                        if a.site_merge == "A":
                            Wsite[q] = max(wm[s] for s in mem)                        # fully coherent noise: no sqrt(N) gain
                        elif a.site_merge == "B":
                            Wsite[q] = float(np.sqrt(sum(wm[s] ** 2 for s in mem)))   # partially coherent
                        else:
                            Wsite[q] = sum(wm[s] for s in mem)                        # no merge: plain sum (baseline weights)
                    if a.site_cap and a.site_merge and len(Wsite) > 1:
                        for kv in a.site_cap.split(","):
                            q_, c_ = kv.split("=")[0], float(kv.split("=")[1])
                            if q_ in Wsite:
                                oth = sum(v for q, v in Wsite.items() if q != q_)
                                Wsite[q_] = min(Wsite[q_], c_ / (1.0 - c_) * oth)     # share of the total <= c_
                    for q, mem in groups.items():
                        wsum_m = sum(wm[s] for s in mem)
                        cs = torch.zeros(NVo, device=dev); cs_s = torch.zeros(NVo, device=dev)
                        vm = torch.ones(NVo, dtype=torch.bool, device=dev); vm_s = torch.ones(NVo, dtype=torch.bool, device=dev)
                        for s in mem:
                            fs_, cc3 = S3[s]
                            i = fs_.index(f)
                            seg = cc3[i, k[s]:k[s] + NVo]
                            msk = VALID[s][L][k[s]:k[s] + NVo]
                            cs += (wm[s] / wsum_m) * seg                             # within-site w-weighted MEAN
                            cs_s += (wm[s] / wsum_m) * torch.roll(seg, shift_s[s] % NVo)   # members rolled separately
                            weff[s] = Wsite[q] * wm[s] / wsum_m                      # = wm[s] when Wsite = sum (no merge)
                            vm &= msk; vm_s &= torch.roll(msk, shift_s[s] % NVo)
                        sh = shift_site[q] % NVo
                        acc_n += Wsite[q] * cs; acc_s += Wsite[q] * cs_s; acc_ss += Wsite[q] * torch.roll(cs, sh)
                        vmask &= vm; vmask_s &= vm_s; vmask_ss &= torch.roll(vm, sh)
                        wtot += Wsite[q]
                        if GATE:
                            site_ser[q] = (cs, cs_s, vm, sh)
                    acc_n /= wtot; acc_s /= wtot; acc_ss /= wtot
                    if GATE:
                        for q, (cs, cs_s, vm, sh) in site_ser.items():
                            am_ = cs[vm]; mq = am_.median()
                            site_stats[q] = (mq, (am_ - mq).abs().median().clamp_min(1e-6))
                fa = {}
                ctrls = (("net", acc_n, vmask), ("shift", acc_s, vmask_s)) + ((("shiftsite", acc_ss, vmask_ss),) if E6 else ())
                for name, arr, vm in ctrls:
                    am = arr[vm]
                    med = am.median(); mad = (am - med).abs().median().clamp_min(1e-6)
                    for T in (6, 8):
                        fa[f"{name}_fa{T}"] = float(
                            decluster_count(arr[None, :], (med + T * mad)[None], n, vm)[0]
                            * 86400.0 * FS / Nmin)
                    if name == "net":
                        medn, madn = med, mad
                    elif name == "shift":
                        meds, mads = med, mad
                    else:
                        medss, madss = med, mad
                gi = np.where(fam_ev == f)[0]
                pos = np.rint((t_ev[gi] - tref) * FS).astype(np.int64)
                okp = (pos - NLAG >= 0) & (pos + NLAG + 1 <= NVo)
                if okp.any():
                    cols = torch.from_numpy(pos[okp] - NLAG).to(dev)[:, None] + \
                        torch.arange(2 * NLAG + 1, device=dev)[None]
                    vn = acc_n[cols].max(1).values
                    vs = acc_s[cols].max(1).values
                    nres[gi[okp]] = ((vn - medn) / madn).cpu().numpy()
                    nres_shift[gi[okp]] = ((vs - meds) / mads).cpu().numpy()
                    if E6:
                        nres_shift_site[gi[okp]] = ((acc_ss[cols].max(1).values - medss) / madss).cpu().numpy()
                    nch[gi[okp]] = 3 * len(ss)
                okp2 = (pos - NLAG2 >= 0) & (pos + NLAG2 + 1 <= NVo)
                if okp2.any():
                    cols2 = torch.from_numpy(pos[okp2] - NLAG2).to(dev)[:, None] + \
                        torch.arange(2 * NLAG2 + 1, device=dev)[None]
                    nres2[gi[okp2]] = ((acc_n[cols2].max(1).values - medn) / madn).cpu().numpy()
                if a.budget:
                    tw0, tw1 = tref, tref + NVo / FS
                    t_own = np.sort(t_ev[fam_ev == f])
                    t_own = t_own[(t_own >= tw0) & (t_own <= tw1)]
                    inw = (bo_all_t >= tw0) & (bo_all_t <= tw1)
                    t_all = np.sort(bo_all_t[inw])
                    t_oth = np.sort(bo_all_t[inw & ~(bo_all_f == f)])
                    idx_n = decluster_idx_multi(acc_n, torch.stack([medn + T * madn for T in BTHR]),
                                                n, vmask)
                    idx_s = decluster_idx_multi(acc_s, torch.stack([meds + T * mads for T in BTHR]),
                                                n, vmask_s)
                    idx_ss = decluster_idx_multi(acc_ss, torch.stack([medss + T * madss for T in BTHR]),
                                                 n, vmask_ss) if E6 else [np.zeros(0, np.int64) for _ in BTHR]
                    if GATE:
                        # E6 step-2 fallback: two-site coincidence -- keep a peak only if >= GATE[0] sites are each >= GATE[1] MAD
                        # of their OWN site series at that sample; applied identically to the data and to both rolled controls
                        # (torch.roll(x, sh)[j] = x[(j - sh) % N], so the site-rolled sample is read at (j - sh)).
                        def _gate(idx, which):
                            if not len(idx):
                                return idx
                            ii_ = torch.from_numpy(idx).to(dev)
                            cnt = torch.zeros(len(idx), device=dev)
                            for q_, (cs_, css_, vm_, sh_) in site_ser.items():
                                mq_, sq_ = site_stats[q_]
                                v_ = cs_[ii_] if which == 0 else (css_[ii_] if which == 1 else cs_[(ii_ - sh_) % NVo])
                                cnt += ((v_ - mq_) / sq_ >= GATE[1]).float()
                            return idx[(cnt >= GATE[0]).cpu().numpy()]
                        idx_n = [_gate(x, 0) for x in idx_n]; idx_s = [_gate(x, 1) for x in idx_s]; idx_ss = [_gate(x, 2) for x in idx_ss]
                    if a.catalog_out and len(idx_n[0]):
                        # E5 step 5: every declustered detection at the LOWEST threshold (thresholding this list at T >= that
                        # reproduces the T-level set exactly: the peak set is threshold-independent).  Amplitude = SNR-weighted
                        # mean over the contributing stations of cc3_s x sqrt(E3_s / E3_ref_s) at the detection sample.
                        ii0 = torch.from_numpy(idx_n[0]).to(dev)
                        vals = acc_n[ii0]
                        amp_acc = torch.zeros(len(ii0), device=dev); wsum = 0.0
                        con = {q: torch.zeros(len(ii0), device=dev) for q in groups}     # E6: per-site contribution to the sum
                        for s in ss:
                            fs_, cc3 = S3[s]; i = fs_.index(f)
                            w = weff[s] if E6 else float(NCH.get(s, 3)) * float(WGT.get((f, s), 1.0))
                            js = ii0 + k[s]
                            amp_acc += w * cc3[i, js] * torch.sqrt(E3W[s][js] / AMPREF[(s, L)])
                            if E6:
                                con[SITE[s]] += w * cc3[i, js]
                            wsum += w
                        td0 = tref + idx_n[0] / FS
                        extra = {}
                        if E6:
                            # con_<site> = that site's weighted contribution / total weight (sum over sites = net_cc);
                            # wsh_<site> = the site's share of the total weight -> share of the sum = con / net_cc vs wsh
                            for q, mem in groups.items():
                                extra[f"con_{q}"] = (con[q] / wsum).cpu().numpy()
                                extra[f"wsh_{q}"] = float(sum(weff[s] for s in mem) / wsum)
                            extra["n_site"] = len(groups)
                        cat_rows.append(pd.DataFrame(dict(
                            doy=doy, L=L, fam=f, t=td0, net_mad=((vals - medn) / madn).cpu().numpy(), net_cc=vals.cpu().numpy(),
                            amp=(amp_acc / wsum).cpu().numpy(), n_sta=len(ss), n_sites=n_sites_f, match_own=min_absdiff(td0, t_own) <= 2.0,
                            match_any=min_absdiff(td0, t_all) <= 2.0, medn=float(medn), madn=float(madn), **extra)))
                    rok = (rt >= tw0) & (rt <= tw1)
                    n_rand = int(rok.sum())
                    n_rand60 = int((min_absdiff(rt[rok], t_oth) <= 60.0).sum())
                    for iT, T in enumerate(BTHR):
                        td = tref + idx_n[iT] / FS
                        tsd = tref + idx_s[iT] / FS
                        m_own = min_absdiff(td, t_own) <= 2.0
                        m_any = min_absdiff(td, t_all) <= 2.0
                        r_own = min_absdiff(tsd, t_own) <= 2.0
                        ex = td[~m_own]
                        # held-out split: detections / own events falling in TEST blocks (odd 2-h
                        # blocks, outside the 60-s guard) -- templates and moveout never saw them
                        def _tst(tt):
                            bb = np.floor(tt / a.block_s).astype(np.int64); pp = tt - bb * a.block_s
                            return (bb % 2 == 1) & (pp >= a.guard_s) & (pp <= a.block_s - a.guard_s)
                        tst, tsts, own_tst = _tst(td), _tst(tsd), _tst(t_own)
                        e6cols = {}
                        if E6:
                            tssd = tref + idx_ss[iT] / FS
                            e6cols = dict(n_det_roll_site=int(len(tssd)), n_det_roll_site_test=int(_tst(tssd).sum()),
                                          n_roll_site_match_own=int((min_absdiff(tssd, t_own) <= 2.0).sum()),
                                          n_site=len(groups), w_share_PGC=float(sum(weff[s] for s in groups.get("PGC", [])) / wtot))
                        bud_rows.append(dict(
                            n_det_test=int(tst.sum()), n_match_own_test=int((m_own & tst).sum()),
                            n_match_any_test=int((m_any & tst).sum()),
                            n_det_roll_test=int(tsts.sum()), n_ev_own_test=int(own_tst.sum()), **e6cols,
                            doy=doy, L=L, fam=f, T=T, n_sta=len(ss), n_sites=n_sites_f,
                            n_ch=int(sum(NCH.get(s, 3) for s in ss)),
                            n_det=int(len(td)), n_det_roll=int(len(tsd)),
                            n_match_own=int(m_own.sum()), n_match_any=int(m_any.sum()),
                            n_roll_match_own=int(r_own.sum()),
                            n_extra=int(len(ex)),
                            n_extra_near60_other=int((min_absdiff(ex, t_oth) <= 60.0).sum()),
                            n_rand=n_rand, n_rand_near60_other=n_rand60,
                            n_ev_own=int(len(t_own)), n_ev_all=int(len(t_all)),
                            n_valid_s=float(vmask.sum().item() / FS),
                            medn=float(medn), madn=float(madn)))
                rb = np.rint((rt - tref) * FS).astype(np.int64)
                okr = (rb - NLAG >= 0) & (rb + NLAG + 1 <= NVo)
                if okr.any():
                    colsr = torch.from_numpy(rb[okr] - NLAG).to(dev)[:, None] + \
                        torch.arange(2 * NLAG + 1, device=dev)[None]
                    vr = acc_n[colsr].max(1).values
                    vrs = acc_s[colsr].max(1).values
                    rn = ((vr - medn) / madn).cpu().numpy()
                    rs = ((vrs - meds) / mads).cpu().numpy()
                    rndn.append(rn)
                    netfam.append(dict(doy=doy, L=L, fam=f, n_sta=len(ss), n_sites=n_sites_f,
                                       n_ch=int(sum(NCH.get(s, 3) for s in ss)),
                                       net_rand_med=float(np.median(rn)),
                                       net_rand_p99=float(np.percentile(rn, 99)),
                                       net_rand_frac6=float(np.mean(rn >= 6)),
                                       net_rand_frac8=float(np.mean(rn >= 8)),
                                       shift_rand_med=float(np.median(rs)), **fa))
            # ---------- accumulate ----------
            rows = pd.DataFrame(dict(
                doy=doy, L=L, t_origin=t_ev, fam=fam_ev, mag=evd.mag.to_numpy(),
                half=evd.half.to_numpy(), guard=evd.guard.to_numpy(),
                net_mad=nres, net_mad_shift=nres_shift, net_mad_2s=nres2, n_ch=nch,
                **({"net_mad_shift_site": nres_shift_site} if E6 else {}),
                **{f"sta_{s}": store[s]["res3"] if s in store else np.nan for s in STAS},
                **{f"staZ_{s}": store[s]["resZ"] if s in store else np.nan for s in STAS},
                **{f"lag_{s}": store[s]["lag"] if s in store else np.nan for s in STAS}))
            if len(evd):
                ev_rows.append(rows)
            fa_rows.extend(netfam)
            del S3, SZ, VALID, E3W
        del D, DH, EN
        torch.cuda.empty_cache()
        print(f"  [day {doy}] {time.time()-t_start:6.0f} s elapsed", flush=True)
        if E6 and a.budget:
            # E6 (GPU may go offline mid-pass): per-day PARTIAL checkpoint of the budget rows and detections (overwritten daily;
            # the final files below carry no ".partial").
            pd.DataFrame(bud_rows).to_csv(f"{RES}/e5_budget_{a.tag}.partial.csv", index=False)
            if a.catalog_out and cat_rows:
                pd.concat(cat_rows, ignore_index=True).to_parquet(f"{RES}/e5r_detections_{a.tag}.partial.parquet", index=False)

    # E6: budget rows and detections are written FIRST (they are the deliverable of a --budget pass), then the event tables;
    # a run over event-free days only (step 1b: quiet days) has no event rows -> empty frame instead of a concat crash.
    if a.catalog_out and cat_rows:
        CAT = pd.concat(cat_rows, ignore_index=True)
        CAT.to_parquet(f"{RES}/e5r_detections_{a.tag}.parquet", index=False)
        print(f"[write] results/e5r_detections_{a.tag}.parquet ({len(CAT)} detections at >= {min(BTHR)} MAD)", flush=True)
    if a.budget:
        pd.DataFrame(bud_rows).to_csv(f"{RES}/e5_budget_{a.tag}.csv", index=False)
        print(f"[write] results/e5_budget_{a.tag}.csv ({len(bud_rows)} family-day-threshold rows)",
              flush=True)
        for p_ in (f"{RES}/e5_budget_{a.tag}.partial.csv", f"{RES}/e5r_detections_{a.tag}.partial.parquet"):
            if os.path.exists(p_):
                os.remove(p_)
    EV = pd.concat(ev_rows, ignore_index=True) if ev_rows else pd.DataFrame(columns=["doy", "L", "t_origin", "fam", "mag", "half", "guard", "net_mad"])
    EV.to_parquet(f"{RES}/e4_tt_events_{a.tag}.parquet", index=False)
    ST = pd.DataFrame(sta_rows)
    ST.to_csv(f"{RES}/e4_tt_station_{a.tag}.csv", index=False)
    FA = pd.DataFrame(fa_rows)          # per (day, L, family): network background + FA rates
    FA.to_csv(f"{RES}/e4_tt_netfam_{a.tag}.csv", index=False)
    # empirical moveout for pass 2: per (family, station, L) median best lag over all days
    mvrows = []
    for (f, s, L), g in (ST.groupby(["fam", "sta", "L"]) if len(ST) else []):    # E6: no station rows on event-free-only runs
        w = g.dropna(subset=["lag_median"])
        if not len(w):
            continue
        mvrows.append(dict(fam=f, sta=s, L=L, n_ev=int(g.n_ev.sum()),
                           lag_median=float(np.average(w.lag_median, weights=np.maximum(w.n_ev, 1))),
                           mad3_median=float(np.nanmedian(g.mad3_median))))
    pd.DataFrame(mvrows).to_csv(f"{RES}/e4_moveout_{a.tag}.csv", index=False)
    print(f"[write] {RES}/e4_tt_{{events,station,netfam}}_{a.tag}, e4_moveout_{a.tag}.csv "
          f"({len(EV)} event rows)", flush=True)
    if a.min_sites:
        print(f"[F8] --min-sites {a.min_sites}: {n_unobs} family-day(-L) combinations skipped as UNOBSERVABLE "
              f"(< {a.min_sites} distinct sites with a usable template and data)", flush=True)


if __name__ == "__main__":
    main()
