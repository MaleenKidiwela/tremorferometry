#!/usr/bin/env python
"""Track C (production) — shared configuration and loaders.

THE FROZEN CHAIN (Merlin F6 + Track A step 1, results/frozen_operating_point.json):
  20 Hz, 1-8 Hz Butterworth-4 zero-phase, component order Z / N|1 / E|2 resolved by e4_chanorder and
  NOT rotated; 12-s templates; mergeB site merge; SNR weights; per (family, day) network MAD;
  T_ext = 8.25, T_core = 8.75; HF-ratio veto at PGC (10-19 Hz over 2-8 Hz), cut 3.6384.
Nothing here re-decides any of that: every number is read from the frozen JSON where one exists.

Read-only inputs live in ../lfe_catalog_repair/results (templates, weights, moveout, families,
stations, amplitude reference).  Nothing outside lfe_extension is ever written.
"""
from __future__ import annotations
import os, sys, json, glob
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)                                   # .../lfe_extension
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"                                         # READ-ONLY
RES = f"{EXT}/results"
TRC = f"{RES}/trackC"
CACHE = f"{EXT}/cache20hz"
TRACES = f"{EXT}/traces"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
sys.path.insert(0, HERE)

FS = 20.0
PRE = 0.5
R_EARTH_KM_DEG = 111.195

# ---- the frozen operating point -------------------------------------------------------------
def frozen():
    with open(f"{RES}/frozen_operating_point.json") as f:
        fz = json.load(f)
    # a newer Track A iteration may name its own templates in results/TRACK_A_RESULTS.md; c_run_month
    # reads that file and overrides `templates` when it does (see resolve_templates).
    return fz


def resolve_templates(fz, override=""):
    """-> (path, tag).  The 2012 set unless TRACK_A_RESULTS.md names a newer one or --templates is given."""
    if override:
        return override, os.path.basename(override).replace(".npz", "")
    md = f"{RES}/TRACK_A_RESULTS.md"
    if os.path.exists(md):
        import re
        txt = open(md).read()
        m = re.search(r"^\s*FINAL_TEMPLATES\s*[:=]\s*(\S+\.npz)\s*$", txt, re.M)
        if m:
            p = m.group(1)
            p = p if os.path.isabs(p) else os.path.join(EXT, p)
            if os.path.exists(p):
                return p, os.path.basename(p).replace(".npz", "")
    return fz["inputs"]["templates"]["path"], os.path.basename(fz["inputs"]["templates"]["path"]).replace(".npz", "")


# ---- network ---------------------------------------------------------------------------------
# Station order of the FROZEN run (scripts/a1_run_mergeB_L12.sh --stations): the float32 site sum is
# accumulated in THIS order, so keeping it makes the production series bit-comparable with Track A.
STAS_ORDER = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
              "B009", "B010", "B011", "B926"]


def order_stations(stas):
    k = {s: i for i, s in enumerate(STAS_ORDER)}
    return sorted(stas, key=lambda s: (k.get(s, 999), s))


# Instrument-epoch boundaries AFTER the 2012 templates were built (Merlin F13 + station_epochs.csv +
# the codex triplet-epoch table).  A station-day at or after its boundary is NOT covered by a template
# built on 2012 instruments: it stays `2012_provisional` until Track A's A-2 agent delivers
# results/a4_templates_<STA>_<YYYY-MM-DD>.npz for that epoch (see templates_for_epoch()).
POST_2012_BOUNDARY = {"NLLB": "2017-02-23", "VGZ": "2017-03-03", "SNB": "2017-06-08",
                      "PGC": "2017-08-10", "LRIV": "2019-05-30"}


# ---- NETWORK EPOCHS (coordinator, 2026-09-17) -------------------------------------------------
# The detector's noise floor is a property of the NETWORK, not of the year: at a fixed MAD a smaller
# network has a higher quiet-day false-alarm rate (Track D measured ~3x at the frozen T after SNB was
# lost in 2021-03).  The two tiers are therefore re-derived per epoch by the SAME false-alarm-only
# rule that defined them (see scripts/c_tiers_by_epoch.py).  Month-granularity, first..last inclusive.
NETWORK_EPOCHS = [("E0", "2014-01", "2017-02", "15 stations, pre-CN-upgrade"),
                  ("E1", "2017-03", "2018-12", "post CN sensor upgrade (PGC/VGZ/NLLB/SNB Trillium)"),
                  ("E2", "2019-01", "2021-03", "no PFB; LRIV dropped from 2019-05-30"),
                  ("E3", "2021-04", "2099-12", "no SNB")]


def epoch_of_month(ym):
    for name, a, b, _ in NETWORK_EPOCHS:
        if a <= ym <= b:
            return name
    return "E?"


def tiers_by_epoch():
    """results/trackC/tiers_by_epoch.csv -> {epoch: (T_core, T_ext)}; {} when not derived yet."""
    p = f"{TRC}/tiers_by_epoch.csv"
    if not os.path.exists(p):
        return {}
    T = pd.read_csv(p)
    return {str(r.epoch): (float(r.T_core), float(r.T_ext)) for r in T.itertuples()
            if np.isfinite(r.T_core) and np.isfinite(r.T_ext)}


SITES = {"PGC": ["PGC", "B009", "B010", "B011"], "LRIV": ["LRIV", "B005", "B006", "B007"]}


def station_table():
    s = pd.read_csv(f"{IN}/e5r_stations.csv")
    return s.set_index("sta")


def site_map(stas):
    sm = {s: s for s in stas}
    for q, mem in SITES.items():
        for s in mem:
            if s in sm:
                sm[s] = q
    return sm


def epochs():
    e = pd.read_csv(f"{RES}/station_epochs.csv")
    e["start"] = pd.to_datetime(e.start)
    e["end"] = pd.to_datetime(e.end)
    e["netsta"] = e.net + "." + e.sta
    e["epoch_id"] = e.netsta + ":" + e.band + ":" + e.start.dt.strftime("%Y-%m-%d")
    return e


def alive(ep, day):
    """Stations (and their band/epoch_id) alive on a pandas Timestamp `day`.

    The epoch END is EXCLUSIVE, matching scripts/download_plan.py: a station-day ON the end date was
    queued `no_epoch` and never fetched, so treating it as inclusive here made five months
    (2014-06 CN.LZB 06-17, 2014-09 CN.YOUB 09-03, 2015-01 C8.SOKB 01-18, 2017-02 CN.NLLB 02-21,
    2018-12 CN.PFB 12-20) wait forever on a station-day that will never arrive.  On a band-switch date
    the two epochs abut, so the station stays alive through the NEW epoch and nothing is lost."""
    d = pd.Timestamp(day)
    return ep[(ep.start <= d) & (d < ep.end)]


def manifest():
    m = pd.read_csv(f"{RES}/download_manifest.csv", dtype={"doy": "Int64"})
    m["netsta"] = m.net + "." + m.sta
    return m


def trace_path(netsta, year, doy, note=""):
    """Our own traces first, then the read-only lfe_catalog_repair copy named in the manifest note."""
    p = f"{TRACES}/{netsta}/{year}/{doy:03d}.mseed"
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return p
    if isinstance(note, str) and note.endswith(".mseed") and os.path.exists(note) and os.path.getsize(note) > 0:
        return note
    q = f"{REP}/traces/{netsta}/{year}/{doy:03d}.mseed"
    if os.path.exists(q) and os.path.getsize(q) > 0:
        return q
    return None


def cache_path(netsta, year, doy):
    return f"{CACHE}/{netsta}/{year}/{doy:03d}.npz"


# ---- month helpers ---------------------------------------------------------------------------
def month_days(ym):
    a = pd.Timestamp(ym + "-01")
    b = a + pd.offsets.MonthEnd(1)
    return pd.date_range(a, b, freq="D")


def epi_km(lat1, lon1, lat2, lon2):
    cs = np.cos(np.radians(48.6))
    return np.hypot((np.asarray(lon1) - lon2) * R_EARTH_KM_DEG * cs,
                    (np.asarray(lat1) - lat2) * R_EARTH_KM_DEG)


# ---- frozen template / weight / moveout tables ------------------------------------------------
def load_frozen_tables(fz, tpl_path, L):
    """-> tpl {(fam,sta):{L:(3,n)}}, sp {(fam,sta): S-P}, fam_meta DataFrame (bkey lat lon dep),
    wgt {(fam,sta): w}, mv {(fam,sta,L): lag}, ampref {(sta,L): e3_ref}."""
    z = np.load(tpl_path, allow_pickle=True)
    tpl = z["tpl"].item()
    sp = z["sp"].item()
    tpl = {k: v for k, v in tpl.items() if L in v}
    F = pd.read_csv(fz["inputs"]["families"]["path"])
    W = pd.read_csv(fz["inputs"]["weights"]["path"])
    wgt = {(float(r.fam), str(r.sta)): float(r.w) for r in W.itertuples()}
    M = pd.read_csv(fz["inputs"]["moveout"]["path"])
    mv = {(float(r.fam), str(r.sta), float(r.L)): float(r.lag_median)
          for r in M.itertuples() if np.isfinite(r.lag_median)}
    A = pd.read_csv(f"{IN}/e5r_amp_ref_2012.csv")
    ampref = {(str(r.sta), float(r.L)): float(r.e3_ref) for r in A.itertuples()}
    return tpl, sp, F, wgt, mv, ampref


def placement(tpl, sp, F, mv, stas_tbl, L):
    """off[(fam, sta)] = ttP + (S-P) - PRE + empirical moveout — the FROZEN placement of
    e4_template_test.py (lines 'predicted template-start offset from origin')."""
    import c_vel1d as v1d
    Fi = F.set_index("bkey")
    off = {}
    for (f, s) in tpl:
        if (f, s) not in sp or s not in stas_tbl.index:
            continue
        if int(f) not in Fi.index:
            continue
        m = Fi.loc[int(f)]
        dx = epi_km(float(m.lat), float(m.lon), float(stas_tbl.loc[s, "lat"]), float(stas_tbl.loc[s, "lon"]))
        off[(f, s)] = (float(v1d.tt(dx, float(m.dep), "P")) + float(sp[(f, s)]) - PRE
                       + mv.get((f, s, L), 0.0))
    return off


def family_box(F, pad=0.2):
    return (F.lat.min() - pad, F.lat.max() + pad, F.lon.min() - pad, F.lon.max() + pad)


def share_ok_matrix(stas=None):
    """Frozen cross-family station-sharing matrix (Merlin Step 7 / E6 F8): |A & B| / min(|A|,|B|) >= 0.5
    over the 2012 usable-template station sets.  bkey-indexed bool matrix."""
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    if stas:
        qc = qc[qc.sta.isin(stas)]
    sets = qc[qc.usable & qc.template].groupby("bkey").sta.apply(set).to_dict()
    nb = int(max(sets) + 1) if sets else 1
    S = np.zeros((nb + 1, nb + 1), bool)
    for fa, A in sets.items():
        for fb, B in sets.items():
            S[fa, fb] = len(A & B) / max(min(len(A), len(B)), 1) >= 0.5
    return S


def local_max_keep(t, v, w=2.0):
    """max-network-MAD cross-family assignment (the pre-E6 rule; kept as the alternative column)."""
    o = np.argsort(t, kind="stable"); ts = t[o]; vs = v[o]; n = len(ts)
    lo = np.searchsorted(ts, ts - w, "left"); hi = np.searchsorted(ts, ts + w, "right")
    keep = np.fromiter((vs[i] >= vs[lo[i]:hi[i]].max() for i in range(n)), bool, n)
    out = np.zeros(n, bool); out[o] = keep
    return out


def local_max_keep_cc(t, v, fam, share_ok, w=2.0):
    """FROZEN cross-family assignment (Merlin Step 7 / E6 F8), verbatim from e5r_score.py:
    a detection is dropped iff a detection of a STATION-SHARING family within +-w s has a higher RAW
    weighted cc (ties: the earlier one wins).  Families on disjoint station sets never suppress."""
    o = np.argsort(t, kind="stable"); ts = t[o]; vs = v[o]; fs = fam[o]; n = len(ts)
    lo = np.searchsorted(ts, ts - w, "left"); hi = np.searchsorted(ts, ts + w, "right")
    keep = np.ones(n, bool)
    for i in range(n):
        fi = fs[i]; vi = vs[i]
        for j in range(lo[i], hi[i]):
            if j == i or not share_ok[fi, fs[j]]:
                continue
            if vs[j] > vi or (vs[j] == vi and j < i):
                keep[i] = False
                break
    out = np.zeros(n, bool); out[o] = keep
    return out


def eff_epoch_start(sta, day, ep):
    """The epoch a station-day belongs to for TEMPLATE purposes: the later of its StationXML epoch start
    and its post-2012 instrument boundary (if the day is at or after it).  Returns 'YYYY-MM-DD'."""
    d = pd.Timestamp(day)
    m = ep[(ep.sta == sta) & (ep.start <= d) & (d < ep.end)]        # end EXCLUSIVE, see alive()
    start = m.start.max() if len(m) else pd.Timestamp("1900-01-01")   # a day ON a band-switch date
                                                                     # belongs to the LATER epoch
    b = POST_2012_BOUNDARY.get(sta)
    if b is not None and d >= pd.Timestamp(b) and pd.Timestamp(b) > start:
        start = pd.Timestamp(b)
    return start.strftime("%Y-%m-%d")


def is_post_boundary(sta, day):
    b = POST_2012_BOUNDARY.get(sta)
    return bool(b is not None and pd.Timestamp(day) >= pd.Timestamp(b))


def epoch_policy():
    """results/trackC/templates_by_epoch.csv — the per-(station, effective epoch) TEMPLATE POLICY.
    Columns (only `sta` and `epoch_start` are required):
      sta          station code
      epoch_start  the EFFECTIVE epoch start this row applies to, 'YYYY-MM-DD' (see eff_epoch_start:
                   the later of the StationXML epoch start and the station's post-2012 boundary)
      npz          per-epoch template file to use instead of the frozen 2012 set ('' = keep 2012)
      label        the epoch_template label to stamp ('' = decide automatically)
      drop         1 = EXCLUDE this station from the network for this epoch entirely
      note         free text
    Returns {(sta, epoch_start): row-dict}."""
    m = f"{TRC}/templates_by_epoch.csv"
    if not os.path.exists(m):
        return {}
    T = pd.read_csv(m, dtype=str).fillna("")
    out = {}
    for r in T.itertuples():
        out[(str(r.sta), str(r.epoch_start))] = dict(
            npz=str(getattr(r, "npz", "")), label=str(getattr(r, "label", "")),
            drop=str(getattr(r, "drop", "")).strip() in ("1", "true", "True", "yes"),
            note=str(getattr(r, "note", "")))
    return out


def epoch_template_file(sta, epoch_start, policy=None):
    """Per-epoch template drop-in.  Highest precedence: the `npz` column of the policy above.  Then,
    by convention, results/a4_templates_<STA>_<YYYY-MM-DD>.npz or ..._<YYYYMMDD>.npz.  Returns a path
    or None — no code change is needed when the files appear."""
    pol = (policy if policy is not None else epoch_policy()).get((sta, epoch_start))
    if pol and pol["npz"]:
        q = pol["npz"]
        q = q if os.path.isabs(q) else os.path.join(EXT, q)
        if os.path.exists(q):
            return q
        raise SystemExit(f"[tpl] templates_by_epoch.csv names {pol['npz']} for {sta}@{epoch_start} "
                         f"but it does not exist")
    if pol is not None and not pol["npz"]:
        return None            # policy says explicitly: keep the frozen 2012 set for this epoch
    for cand in (f"{RES}/a4_templates_{sta}_{epoch_start}.npz",
                 f"{RES}/a4_templates_{sta}_{epoch_start.replace('-', '')}.npz"):
        if os.path.exists(cand):
            return cand
    return None


# ---- per-epoch HF-veto cuts ------------------------------------------------------------------
# results/a4_hf_veto_cuts.csv (written by the Track A-2 agent).  Columns: `cut` is required; `sta`,
# `fs`, `epoch_start`, `epoch_end` are optional constraints, and the MOST CONSTRAINED matching row
# wins.  The frozen cut (3.6384, PGC BH 40 Hz) is the fallback, and is used with a warning when no
# row matches — it is NOT epoch-invariant (the 10-19 Hz band sits inside the BH anti-alias roll-off).
def hf_cuts(band="10-19"):
    """results/a4_hf_veto_cuts.csv.  TWO schemas are accepted:
      (a) Track A-2's:  station_epoch ('PGC_BH_2010-08-01..2017-08-10'), band ('10-19'), fs, cut, ...
      (b) the generic:  sta, fs, epoch_start, epoch_end, cut
    `cut` is required; everything else is an optional constraint and the MOST CONSTRAINED matching row
    wins.  `band` MUST match the band the stored `hfr` was computed in (default 10-19 Hz) — a cut
    calibrated in another band is silently wrong, so rows in other bands are DROPPED here."""
    p = f"{RES}/a4_hf_veto_cuts.csv"
    if not os.path.exists(p):
        return []
    T = pd.read_csv(p).fillna("")
    rows = []
    for r in T.to_dict("records"):
        if r.get("cut", "") == "":
            continue
        if "band" in r and str(r["band"]).strip() and str(r["band"]).strip() != band:
            continue                                   # a cut for a band we do not measure
        sta = str(r.get("sta", "") or "")
        start = str(r.get("epoch_start", "") or r.get("start", "") or "")
        end = str(r.get("epoch_end", "") or r.get("end", "") or "")
        se = str(r.get("station_epoch", "") or "")
        if se:
            # '<STA>_<BAND>_<start>..<end>'  (end may be a bare year, e.g. '2026')
            head, _, span = se.rpartition("_")
            sta = sta or head.split("_")[0]
            if ".." in span:
                a_, b_ = span.split("..", 1)
                start = start or a_.strip()
                end = end or (b_.strip() if "-" in b_ else f"{b_.strip()}-12-31")
        rows.append(dict(sta=sta, fs=r.get("fs", ""), start=start, end=end, cut=float(r["cut"]),
                         method=str(r.get("method", "")),
                         src=(se or f"{sta}:{start}..{end}") + " | " + str(r.get("method", ""))[:40]))
    return rows


def hf_cut_for(rows, sta, date, fs, default):
    """-> (cut, source string).  The MOST CONSTRAINED matching row wins.  Among equally constrained
    rows the ANCHOR wins — a row whose `method` says it is the frozen referee cut — so a per-epoch
    recalibration can never silently move the epoch the operating point was frozen on; otherwise the
    last matching row wins.  Ambiguity that the anchor rule does not settle is reported."""
    best, best_n, best_anchor, amb = default, -1, False, []
    for r in rows:
        n = 0
        if r["sta"]:
            if r["sta"] != sta:
                continue
            n += 1
        if r["fs"] not in ("", None) and not (isinstance(r["fs"], float) and np.isnan(r["fs"])):
            # the recorded native rate can be 99.99999 Hz on CN fragments -> compare snapped
            if not (np.isfinite(fs) and abs(round(float(r["fs"])) - round(float(fs))) < 1e-6):
                continue
            n += 1
        if r["start"]:
            if str(date) < r["start"]:
                continue
            n += 1
        if r["end"]:
            if str(date) >= r["end"]:
                continue
            n += 1
        anchor = "FROZEN" in r.get("method", "").upper()
        if n > best_n or (n == best_n and (anchor or not best_anchor)):
            if n == best_n and best != r["cut"] and not (anchor or best_anchor):
                amb.append(r["src"])
            best, best_n, best_anchor = r["cut"], n, anchor or (n == best_n and best_anchor)
    src = "frozen_default"
    for r in rows:
        if r["cut"] == best:
            src = r["src"]
            break
    if amb:
        src += f" [ambiguous, also matched: {'; '.join(amb)}]"
    return best, src
