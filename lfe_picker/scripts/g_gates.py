#!/usr/bin/env python
"""g_gates — the reusable, detector-agnostic FAMILY-ACCEPTANCE pipeline (Gates 0-3).

Input  : a PROPOSAL = a set of detection times on the ORIGIN axis + a proposed interface cell
         (lat, lon, dep), from ANY detector (template matched filter now, the P2 model later), plus
         a region / station configuration.
Output : per-gate statistics with every null, and PASS / FAIL / UNTESTABLE per gate and overall.

Spec: lfe_picker/results/merlin_review_plan_2026-09-17.md, section "The gates, quantified".

  Gate 0 (identity / class) -- is this a NEW source, and is it LFE-class?
      distance to every accepted family; multi-site stack cc to the NEAREST accepted family against
      that station's cross-family cc p95 (results/gates/g_crossfam_cc.csv); HF ratio (10-19 / 2-8 Hz)
      against the LFE core distribution (median ~0.25, veto >= 3.5).
      PASS = (min distance >= 5 km OR cc below the per-station chance level at >= 2 sites)
             AND the HF class test does not veto.

  Gate 1 (withheld-station stack) -- is there a SOURCE there at all?
      >= 2 stations at >= 2 distinct sites excluded from detection for the whole run; N >= 20 members;
      statistic 10 log10 E(S-1..S+5) / E(S-12..S-6) on the 3-C 20 Hz stack at the predicted S
      (cell + Crosson 1-D + a5 station correction).
      Nulls: 200 matched random-time stacks (need excess > p99 AND >= +3 dB); 100 moveout-shuffled
      stacks at another cell whose S offset differs by >= 5 s (need the true-moveout excess > p95).
      Positive control: a known Bostock family at MATCHED N and comparable distance must pass at the
      same station, else the verdict is UNTESTABLE rather than FAIL.

  Gate 2 (location with controls)
      S picks on the per-site family stacks, >= 4 sites within 60 km, azimuthal gap <= 180 deg,
      the a5 estimator on the 2-D Slab2 interface grid + a free-depth sanity run.
      Controls: shuffled-pick null (observed misfit < 0.5x null), synthetic floor at that geometry
      (0.3 s S noise, +-0.5 s statics; expect <= 3 km), leave-one-site-out spread <= 5 km, free depth
      within +-10 km of Slab2 and not at a grid bound.

  Gate 3 (tremor timing, SOFT)
      E_local = (fraction of members on days with PNSN tremor within 25 km of the cell) / (calendar
      fraction of such days); E_far the same for 75-150 km.  Need E_local >= 3, E_local/E_far >= 1.5,
      p < 0.01 from 1,000 circular shifts of >= 30 d.  A family that fails Gate 3 but passes Gates
      1-2 strongly is reported as an 'inter-ETS candidate', never rejected.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402
import g_locate as L                                             # noqa: E402

PASS, FAIL, UNTESTABLE = "PASS", "FAIL", "UNTESTABLE"


# ================================================================================ proposal & config
@dataclass
class Proposal:
    pid: str
    times: np.ndarray                 # absolute epoch seconds, ORIGIN axis
    lat: float
    lon: float
    dep: float = np.nan               # NaN -> Slab2 at (lat, lon)
    offset_search_s: float = 0.0      # >0: search a CONSTANT time offset (the same search is applied
                                      # to every null, so the extra freedom is paid for)
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.times = np.sort(np.asarray(self.times, float))
        if not np.isfinite(self.dep):
            self.dep = float(G.slab2_depth([self.lat], [self.lon])[0])


@dataclass
class GateConfig:
    # region
    sites: tuple = ("PGC", "VGZ", "SNB", "NLLB", "LRIV", "B926", "B001", "B003", "B004")
    withheld: tuple = ("VGZ", "SNB", "B926")     # excluded from detection for the whole run
    grid_lat0: float = 48.50
    grid_lon0: float = -123.75
    grid_halfwidth: float = 62.0
    # gate 0
    g0_min_km: float = 5.0
    g0_min_sites_below: int = 2
    g0_reg_db: float = 2.0            # a site "registers" when the stack S-window excess >= this
    g0_max_frac_above: float = 0.5    # alias if >= half the registering sites are at/above chance
    g0_hf_veto: float = 3.5
    g0_hf_station: tuple = ("PGC", "B011")
    # gate 1
    g1_min_members: int = 20
    g1_n_time_null: int = 200
    g1_n_moveout_null: int = 100
    g1_min_excess_db: float = 3.0
    g1_moveout_min_dt: float = 5.0
    g1_max_dist_km: float = 60.0                  # blind-spot flag, not a bar
    # gate 2
    g2_min_sites_60km: int = 4
    g2_max_gap_deg: float = 180.0
    g2_min_pick_snr: float = 3.0
    g2_align: bool = True                         # cc-align the pick stacks (e5r align_stack)
    g2_n_pick_null: int = 200
    g2_n_synth: int = 50
    g2_synth_bar_km: float = 3.0
    g2_loso_bar_km: float = 5.0
    g2_depth_bar_km: float = 10.0
    g2_null_ratio: float = 0.5
    # gate 3
    g3_local_km: float = 25.0
    g3_far_km: tuple = (75.0, 150.0)
    g3_min_elocal: float = 3.0
    g3_min_ratio: float = 1.5
    g3_p: float = 0.01
    g3_n_shift: int = 1000
    g3_min_shift_d: float = 30.0
    seed: int = 20260917


# ================================================================================ context
class GateContext:
    def __init__(self, cfg: GateConfig = None, day_store=None, accepted=None, verbose=True):
        self.cfg = cfg or GateConfig()
        self.verbose = verbose
        self.S = G.stations()
        self.corr = G.station_corrections()
        self.store = day_store or G.DayStore()
        self.nstore = G.NativeStore()
        self.site_ll = {q: (float(self.S.loc[mem[0], "lat"]), float(self.S.loc[mem[0], "lon"]))
                        for q, mem in G.SITE_MEMBERS.items()}
        self.site_members = {q: [m for m in mem if m in self.S.index] for q, mem in G.SITE_MEMBERS.items()}
        # accepted families (Gate 0 reference).  Default = the 130 Bostock families.
        self.accepted = (accepted if accepted is not None
                         else pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv")[
                             ["bkey", "lat", "lon", "dep"]])
        Z = np.load(f"{G.REP}/results/e5r_templates_2012.npz", allow_pickle=True)
        self._tplz = Z
        self.TPL = Z["tpl"].item()
        cf = pd.read_csv(f"{G.GATES}/g_crossfam_cc.csv").set_index("sta")
        self.ccp95 = {s: float(cf.loc[s, "cc_p95"]) for s in cf.index}
        self._grid = None
        self._tremor = None
        self._controls = None

    # ----- lazy heavy objects
    @property
    def grid(self):
        if self._grid is None:
            t = time.time()
            self._grid = L.Grid(self.cfg.sites, self.site_ll, self.cfg.grid_lat0,
                                self.cfg.grid_lon0, self.cfg.grid_halfwidth)
            if self.verbose:
                print(f"[ctx] location grid {len(self._grid.E)} nodes x {len(self._grid.gz)} depths "
                      f"in {time.time() - t:.1f} s", flush=True)
        return self._grid

    @property
    def tremor(self):
        if self._tremor is None:
            T = pd.read_csv(f"{G.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv",
                            usecols=["time", "lat", "lon"])
            T["time"] = pd.to_datetime(T.time).astype("datetime64[ns]")
            T["day"] = T.time.values.astype("datetime64[D]")
            self._tremor = T
        return self._tremor

    @property
    def controls(self):
        if self._controls is None:
            p = f"{G.GATES}/g_gate1_controls.csv"
            self._controls = pd.read_csv(p) if os.path.exists(p) else None
        return self._controls

    # ----- helpers
    def netsta(self, sta):
        return str(self.S.loc[sta, "netsta"])

    def pred_S(self, sta, lat, lon, dep):
        """Predicted OBSERVED S travel time from the cell to the station (Crosson + a5 correction)."""
        d = float(G.epi_km(lat, lon, float(self.S.loc[sta, "lat"]), float(self.S.loc[sta, "lon"])))
        ttS, ttP = G.ttSP(d, dep)
        cS, cP = self.corr.get(sta, (0.0, 0.0))
        return d, ttS + cS, ttP + cP


# ================================================================================ Gate 1 machinery
def _stack_sdb(store, netsta, s_times, offset_grid=None):
    """Stack at the predicted S and return (sdb, n_used, stack, best_offset).
    With an offset grid the stack is cut wide ONCE and the read-out window is slid: a CONSTANT
    per-family offset is exactly a shift of the measurement window, so the search costs nothing."""
    if offset_grid is None:
        st, n = store.stack(netsta, s_times)
        if st is None or n == 0:
            return np.nan, 0, None, 0.0
        return G.sdb(st, G.IDX_S), n, st, 0.0
    pre = 50.0
    nw = int(round(100.0 * G.FS))
    idx0 = int(round(pre * G.FS))
    W, keep = store.cut(netsta, np.asarray(s_times, float) - pre, n=nw)
    if W is None or len(keep) == 0:
        return np.nan, 0, None, 0.0
    st = W.mean(0)
    best, boff = -np.inf, 0.0
    for off in offset_grid:
        i = idx0 + int(round(off * G.FS))
        v = G.sdb(st, i)
        if v > best:
            best, boff = v, float(off)
    return best, len(keep), st, boff


def _gate1_one_station(ctx, prop, sta, rng, n_time=None, n_mo=None):
    cfg = ctx.cfg
    n_time = cfg.g1_n_time_null if n_time is None else n_time
    n_mo = cfg.g1_n_moveout_null if n_mo is None else n_mo
    netsta = ctx.netsta(sta)
    dist, ttS, ttP = ctx.pred_S(sta, prop.lat, prop.lon, prop.dep)
    og = (np.arange(-prop.offset_search_s, prop.offset_search_s + 1e-9, 0.25)
          if prop.offset_search_s > 0 else None)
    s_times = prop.times + ttS
    sdb_true, n_used, stack, boff = _stack_sdb(ctx.store, netsta, s_times, og)
    out = dict(sta=sta, site=G.STA2SITE.get(sta, sta), dist_km=dist, ttS=ttS,
               n_members=len(prop.times), n_used=int(n_used), sdb_db=sdb_true,
               offset_fit_s=boff)
    if n_used < cfg.g1_min_members:
        out["verdict"] = UNTESTABLE
        out["reason"] = f"n_used {n_used} < {cfg.g1_min_members}"
        return out
    out["peak_off_s"] = G.peak_offset_s(stack, G.IDX_S if og is None else int(round(50 * G.FS)) + int(round(boff * G.FS)))
    # ---- random-time null: same station, same days, same N
    ts = pd.to_datetime(prop.times, unit="s")
    doys = sorted(set(zip(ts.year.to_numpy().tolist(), ts.dayofyear.to_numpy().tolist())))
    spans = []
    for (y, d) in doys:
        v = ctx.store.get(netsta, int(y), int(d))
        if v is not None:
            spans.append((v[0], v[0] + v[1].shape[1] / G.FS))
    if not spans:
        out["verdict"] = UNTESTABLE
        out["reason"] = "no data on member days"
        return out
    lo = np.array([a for a, _ in spans]); hi = np.array([b for _, b in spans])
    nul = []
    for _ in range(n_time):
        k = rng.integers(0, len(spans), n_used)
        t = lo[k] + rng.random(n_used) * (hi[k] - lo[k] - 120.0) + 60.0
        v, m, _s, _o = _stack_sdb(ctx.store, netsta, t, og)
        if m >= cfg.g1_min_members:
            nul.append(v)
    nul = np.asarray(nul, float)
    if nul.size < 20:
        out["verdict"] = UNTESTABLE
        out["reason"] = "random-time null failed"
        return out
    out.update(null_mean_db=float(nul.mean()), null_sd_db=float(nul.std()),
               null_p99_db=float(np.quantile(nul, .99)), n_null=int(nul.size))
    out["excess_db"] = float(sdb_true - nul.mean())
    out["p_time"] = float(((nul >= sdb_true).sum() + 1) / (nul.size + 1))
    pass_time = bool(sdb_true > out["null_p99_db"] and out["excess_db"] >= cfg.g1_min_excess_db)
    # ---- moveout-shuffled null: SAME times, another cell whose S offset differs by >= 5 s
    gr = ctx.grid
    j = np.arange(len(gr.lat))
    dd = np.hypot((gr.lat - float(ctx.S.loc[sta, "lat"])) * 111.2,
                  (gr.lon - float(ctx.S.loc[sta, "lon"])) * 111.2 * np.cos(np.radians(48.6)))
    xs, zs, TP, TS = G.tt_table()
    ttS_grid = G._interp2(xs, zs, TS, dd, gr.z_if) + ctx.corr.get(sta, (0.0, 0.0))[0]
    ok = j[np.abs(ttS_grid - ttS) >= cfg.g1_moveout_min_dt]
    mo = []
    if ok.size:
        pick = rng.choice(ok, size=min(n_mo, ok.size), replace=ok.size < n_mo)
        for k in pick:
            v, m, _s, _o = _stack_sdb(ctx.store, netsta, prop.times + ttS_grid[k], og)
            if m >= cfg.g1_min_members:
                mo.append(v)
    mo = np.asarray(mo, float)
    if mo.size < 20:
        out["verdict"] = UNTESTABLE
        out["reason"] = "moveout null failed"
        return out
    mo_ex = mo - nul.mean()
    out.update(n_moveout=int(mo.size), moveout_mean_db=float(mo.mean()),
               moveout_p95_excess_db=float(np.quantile(mo_ex, .95)),
               p_moveout=float(((mo >= sdb_true).sum() + 1) / (mo.size + 1)))
    pass_mo = bool(out["excess_db"] > out["moveout_p95_excess_db"])
    out["pass_time_null"] = pass_time
    out["pass_moveout_null"] = pass_mo
    out["verdict"] = PASS if (pass_time and pass_mo) else FAIL
    out["far_station_flag"] = bool(dist > cfg.g1_max_dist_km)
    return out


def gate1(ctx, prop, rng=None, controls=True):
    cfg = ctx.cfg
    rng = rng or np.random.default_rng(cfg.seed)
    per = [_gate1_one_station(ctx, prop, s, rng) for s in cfg.withheld]
    sites = {r["site"] for r in per if r["verdict"] == PASS}
    npass = len(sites)
    testable = [r for r in per if r["verdict"] in (PASS, FAIL)]
    res = dict(gate="gate1", per_station=per, n_sites_pass=npass,
               n_stations_testable=len(testable),
               all_withheld_far=bool(all(r.get("far_station_flag", True) for r in per)))
    ctl = _gate1_control(ctx, prop) if controls else None
    res["positive_control"] = ctl
    if npass >= 1:
        res["verdict"] = PASS
    elif not testable:
        res["verdict"] = UNTESTABLE
        res["reason"] = "; ".join(r.get("reason", "") for r in per)
    elif ctl is not None and not ctl.get("any_pass", False):
        res["verdict"] = UNTESTABLE
        res["reason"] = "positive control fails at matched N / distance -> the test has no power here"
    else:
        res["verdict"] = FAIL
    return res


def _gate1_control(ctx, prop):
    """A known Bostock family at matched N and comparable distance must pass at the same withheld
    station, else the gate has no power for this proposal and the verdict becomes UNTESTABLE."""
    C = ctx.controls
    if C is None or len(C) == 0:
        return dict(available=False)
    N = len(prop.times)
    rows = []
    for sta in ctx.cfg.withheld:
        d = float(G.epi_km(prop.lat, prop.lon, float(ctx.S.loc[sta, "lat"]),
                           float(ctx.S.loc[sta, "lon"])))
        c = C[(C.sta == sta)]
        if len(c) == 0:
            continue
        c = c.assign(dn=(c.N - N).abs(), dd=(c.dist_km - d).abs())
        nb = c.N.unique()
        nsel = nb[np.argmin(np.abs(nb - N))]
        c = c[c.N == nsel].sort_values("dd")
        if len(c) == 0:
            continue
        sub = c.head(3)
        rows.append(dict(sta=sta, N_matched=int(nsel), dist_km=d,
                         n_control=int(len(sub)), n_pass=int((sub.verdict == PASS).sum()),
                         ctl_dist_km=float(sub.dist_km.iloc[0]),
                         ctl_excess_db=float(sub.excess_db.median())))
    if not rows:
        return dict(available=False)
    return dict(available=True, per_station=rows,
                any_pass=bool(any(r["n_pass"] > 0 for r in rows)))


# ================================================================================ Gate 0
def gate0(ctx, prop, rng=None):
    cfg = ctx.cfg
    A = ctx.accepted
    # leave-one-out: a proposal that IS a known family must be judged against the OTHER accepted
    # families, otherwise Gate 0 can only ever say "you are yourself".  The P3 hidden-family protocol
    # is the same operation.  drop_bkey = -1 (the alias control) keeps the full accepted list.
    db = int(prop.meta.get("drop_bkey", -1)) if prop.meta else -1
    if db >= 0:
        A = A[A.bkey != db]
    d = G.epi_km(A.lat.values, A.lon.values, prop.lat, prop.lon)
    i = int(np.argmin(d))
    out = dict(gate="gate0", min_dist_km=float(d[i]), nearest_bkey=int(A.bkey.values[i]),
               nearest_lat=float(A.lat.values[i]), nearest_lon=float(A.lon.values[i]),
               n_accepted_within_5km=int((d < cfg.g0_min_km).sum()))
    nb = int(A.bkey.values[i])
    # ---- multi-site stack cc to the NEAREST accepted family
    per = []
    for site, members in ctx.site_members.items():
        if site not in cfg.sites:
            continue
        for sta in members:
            key = (float(nb), sta)
            if key not in ctx.TPL or 12.0 not in ctx.TPL[key]:
                continue
            tpl = np.asarray(ctx.TPL[key][12.0], float)
            dist, ttS, _ = ctx.pred_S(sta, prop.lat, prop.lon, prop.dep)
            st, n = ctx.store.stack(ctx.netsta(sta), prop.times + ttS)
            if st is None or n < 10:
                continue
            cc = G.cc_template(st, tpl)
            p95 = ctx.ccp95.get(sta, 0.45)
            reg = G.sdb(st, G.IDX_S)
            per.append(dict(sta=sta, site=site, n_used=int(n), cc=float(cc), cc_p95=p95,
                            sdb_db=float(reg), registers=bool(reg >= cfg.g0_reg_db),
                            ratio=float(abs(cc) / p95), below=bool(abs(cc) < p95)))
    out["cc_per_station"] = per
    bysite, bysite_reg = {}, {}
    for r in per:
        bysite.setdefault(r["site"], []).append(r["ratio"])
        if r["registers"]:
            bysite_reg.setdefault(r["site"], []).append(r["ratio"])
    site_ratio = {s: float(max(v)) for s, v in bysite.items()}
    site_ratio_reg = {s: float(max(v)) for s, v in bysite_reg.items()}
    out["site_ratio"] = site_ratio
    out["n_sites_cc"] = len(site_ratio)
    out["n_sites_registering"] = len(site_ratio_reg)
    out["n_sites_below_chance"] = int(sum(1 for v in site_ratio.values() if v < 1.0))
    out["n_reg_above_chance"] = int(sum(1 for v in site_ratio_reg.values() if v >= 1.0))
    out["frac_reg_above_chance"] = (float(out["n_reg_above_chance"]) / len(site_ratio_reg)
                                    if site_ratio_reg else np.nan)
    out["max_site_ratio"] = float(max(site_ratio.values())) if site_ratio else np.nan
    distinct_by_distance = out["min_dist_km"] >= cfg.g0_min_km
    # Merlin's literal rule -- reported, NOT used (validation: it false-passes 100 % of the alias
    # control, because at low-SNR sites the cc is trivially below chance)
    out["distinct_by_waveform_literal"] = bool(
        out["n_sites_below_chance"] >= cfg.g0_min_sites_below)
    # adopted rule: among the sites where the proposal's own stack REGISTERS an arrival, the waveform
    # must be below the cross-family chance level at most of them
    distinct_by_cc = bool(len(site_ratio_reg) >= cfg.g0_min_sites_below
                          and out["frac_reg_above_chance"] < cfg.g0_max_frac_above)
    out["distinct_by_distance"] = bool(distinct_by_distance)
    out["distinct_by_waveform"] = bool(distinct_by_cc)
    # ---- HF class test
    hf = _hf_class(ctx, prop)
    out.update(hf)
    if out["n_sites_registering"] < cfg.g0_min_sites_below and not distinct_by_distance:
        out["verdict"] = UNTESTABLE
        out["reason"] = (f"only {out['n_sites_registering']} site(s) where the proposal stack "
                         f"registers -- the identity comparison has no power")
        return out
    if hf.get("hf_median") is not None and np.isfinite(hf["hf_median"]) and hf["hf_median"] >= cfg.g0_hf_veto:
        out["verdict"] = FAIL
        out["reason"] = f"HF ratio {hf['hf_median']:.2f} >= veto {cfg.g0_hf_veto} (not LFE class)"
        return out
    out["verdict"] = PASS if (distinct_by_distance or distinct_by_cc) else FAIL
    if out["verdict"] == FAIL:
        out["reason"] = (f"{out['min_dist_km']:.2f} km from accepted family {nb} and its waveform is "
                         f"at/above the cross-family chance level at {out['n_sites_cc'] - out['n_sites_below_chance']}"
                         f"/{out['n_sites_cc']} sites")
    return out


def _hf_class(ctx, prop, max_members=200):
    """Per-member E(10-19 Hz)/E(2-8 Hz) at the nearest HF-capable station, on the NATIVE trace.
    The COHERENT stack cannot carry 10-19 Hz (per-member alignment error >> half a period), so the
    class statistic is the MEDIAN over members -- exactly the statistic Track A calibrated
    (LFE core median 0.247 at PGC, veto 3.5)."""
    t = prop.times
    if len(t) > max_members:
        t = t[np.linspace(0, len(t) - 1, max_members).astype(int)]
    for sta in ctx.cfg.g0_hf_station:
        if sta not in ctx.S.index:
            continue
        _d, ttS, _ = ctx.pred_S(sta, prop.lat, prop.lon, prop.dep)
        s_times = t + ttS
        ts = pd.to_datetime(s_times, unit="s")
        by = {}
        for k, tt in zip(zip(ts.year.to_numpy().tolist(), ts.dayofyear.to_numpy().tolist()), s_times):
            by.setdefault(k, []).append(float(tt))
        r = G.hf_ratio_members(ctx.nstore, ctx.netsta(sta), by)
        if r.size >= 10:
            return dict(hf_station=sta, hf_n=int(r.size), hf_median=float(np.median(r)),
                        hf_p75=float(np.quantile(r, .75)))
    return dict(hf_station=None, hf_n=0, hf_median=np.nan, hf_p75=np.nan)


# ================================================================================ Gate 2
def _pick_stack(st, idx_s, ttS, ttP, min_snr, fs=G.FS):
    """S on sqrt(H1^2+H2^2) within +-3 s of the predicted S, P on Z within +-2.5 s of the predicted
    P; envelope-onset + AIC (a5_eq_corrections.env_onset recipe), accepted at peak/noise-rms >= min_snr."""
    from a5_env import env_onset
    H = np.sqrt(st[1] ** 2 + st[2] ** 2)
    iP0 = idx_s - int(round((ttS - ttP) * fs))
    nnoise = max(60, iP0 - int(3.0 * fs))
    iS, snrS = env_onset(H, idx_s - int(3 * fs), idx_s + int(3 * fs), nnoise, fs)
    iP, snrP = env_onset(st[0], iP0 - int(2.5 * fs), iP0 + int(2.5 * fs), nnoise, fs)
    dS = (iS - idx_s) / fs if (iS is not None and snrS >= min_snr) else np.nan
    dP = (iP - iP0) / fs if (iP is not None and snrP >= min_snr) else np.nan
    return dS, snrS, dP, snrP


def gate2(ctx, prop, rng=None):
    cfg = ctx.cfg
    rng = rng or np.random.default_rng(cfg.seed + 1)
    gr = ctx.grid
    picks = []
    for site in cfg.sites:
        for sta in ctx.site_members.get(site, []):
            dist, ttS, ttP = ctx.pred_S(sta, prop.lat, prop.lon, prop.dep)
            W, keep = ctx.store.cut(ctx.netsta(sta), prop.times + ttS - G.HALF_W)
            if W is None or len(keep) < 10:
                continue
            n = len(keep)
            # cc-aligned stack (the e5r recipe a5 picked on).  Gate 1 deliberately does NOT align.
            st = G.align_stack(W)[0] if cfg.g2_align else W.mean(0)
            cS, cP = ctx.corr.get(sta, (0.0, 0.0))
            dS, snrS, dP, snrP = _pick_stack(st, G.IDX_S, ttS, ttP, cfg.g2_min_pick_snr)
            picks.append(dict(sta=sta, site=site, dist_km=dist, n_used=int(n),
                              # travel time with the station correction REMOVED (a5 convention)
                              t_S=(ttS + dS - cS) if np.isfinite(dS) else np.nan,
                              t_P=(ttP + dP - cP) if np.isfinite(dP) else np.nan,
                              snrS=snrS, snrP=snrP))
    P = pd.DataFrame(picks)
    out = dict(gate="gate2", n_station_stacks=int(len(P)),
               n_station_S=int(P.t_S.notna().sum()) if len(P) else 0)
    if len(P) == 0:
        out["verdict"] = UNTESTABLE
        out["reason"] = "no stacks"
        return out
    Sg = P.groupby("site").agg(t_S=("t_S", "median"), t_P=("t_P", "median"),
                               dist_km=("dist_km", "median"), snrS=("snrS", "max"),
                               n_mem_S=("t_S", "count")).reset_index()
    sites = list(gr.sites)
    tS = np.array([float(Sg.set_index("site").t_S.get(s, np.nan)) for s in sites])
    tP = np.array([float(Sg.set_index("site").t_P.get(s, np.nan)) for s in sites])
    sel = [i for i, s in enumerate(sites) if np.isfinite(tS[i])]
    picked = [sites[i] for i in sel]
    near = [s for s in picked
            if G.epi_km(prop.lat, prop.lon, *ctx.site_ll[s]) <= 60.0]
    gap = L.azimuth_gap(prop.lat, prop.lon, picked, ctx.site_ll)
    out.update(n_sites_picked=len(sel), n_sites_60km=len(near), gap_deg=gap,
               sites=",".join(picked))
    if len(near) < cfg.g2_min_sites_60km or gap > cfg.g2_max_gap_deg or len(sel) < 3:
        out["verdict"] = UNTESTABLE
        out["reason"] = (f"geometry: {len(near)} picked sites within 60 km (need "
                         f"{cfg.g2_min_sites_60km}), gap {gap:.0f} deg")
        return out
    r_if = L.locate_interface(gr, tS, tP, sel)
    r_fr = L.locate_free(gr, tS, tP, sel)
    if r_if is None:
        out["verdict"] = UNTESTABLE
        out["reason"] = "locate failed"
        return out
    out["interface"] = r_if
    out["free_depth"] = r_fr
    out["dist_to_proposal_km"] = float(G.epi_km(r_if["lat"], r_if["lon"], prop.lat, prop.lon))
    z_slab = float(G.slab2_depth([r_if["lat"]], [r_if["lon"]])[0])
    out["z_slab_km"] = z_slab
    out["free_dz_km"] = float(r_fr["z"] - z_slab) if r_fr else np.nan
    # ---- shuffled-pick null (truth-free: the misfit, not a distance)
    nul = []
    vals = tS[sel]
    for _ in range(cfg.g2_n_pick_null):
        p = rng.permutation(len(sel))
        if np.all(p == np.arange(len(sel))):
            continue
        t2 = tS.copy(); t3 = tP.copy()
        for a, i in enumerate(sel):
            t2[i] = vals[p[a]]
            t3[i] = np.nan
        r = L.locate_interface(gr, t2, t3, sel)
        if r is not None:
            nul.append(r["rms_s"])
    nul = np.asarray(nul, float)
    r_if_S = L.locate_interface(gr, tS, np.full_like(tP, np.nan), sel)
    out["rms_obs_s"] = float(r_if_S["rms_s"])
    out["rms_null_med_s"] = float(np.median(nul)) if nul.size else np.nan
    out["null_ratio"] = float(out["rms_obs_s"] / out["rms_null_med_s"]) if nul.size else np.nan
    # ---- synthetic floor at THIS geometry
    tPm = [np.isfinite(tP[i]) for i in sel]
    med, p90, _ = L.synthetic_floor(gr, ctx.site_ll, picked, r_if["lat"], r_if["lon"], r_if["z"],
                                    tPm, n_real=cfg.g2_n_synth, seed=cfg.seed)
    out["synth_med_km"], out["synth_p90_km"] = med, p90
    # ---- leave-one-site-out
    loso = []
    for drop in range(len(sel)):
        s2 = [x for a, x in enumerate(sel) if a != drop]
        if len(s2) < 3:
            continue
        r = L.locate_interface(gr, tS, tP, s2)
        if r is not None:
            loso.append((sites[sel[drop]], G.epi_km(r["lat"], r["lon"], r_if["lat"], r_if["lon"])))
    out["loso"] = [dict(dropped=a, dkm=float(b)) for a, b in loso]
    out["loso_max_km"] = float(max(b for _a, b in loso)) if loso else np.nan
    # ---- verdict
    checks = dict(
        null=bool(np.isfinite(out["null_ratio"]) and out["null_ratio"] < cfg.g2_null_ratio),
        synth=bool(np.isfinite(med) and med <= cfg.g2_synth_bar_km),
        loso=bool(np.isfinite(out["loso_max_km"]) and out["loso_max_km"] <= cfg.g2_loso_bar_km),
        depth=bool(r_fr is not None and abs(out["free_dz_km"]) <= cfg.g2_depth_bar_km
                   and not r_fr["at_bound"]),
        bounds=bool(not r_if["at_bound"]))
    out["checks"] = checks
    out["verdict"] = PASS if all(checks.values()) else FAIL
    if out["verdict"] == FAIL:
        out["reason"] = "failed: " + ",".join(k for k, v in checks.items() if not v)
    return out


# ================================================================================ Gate 3
def gate3(ctx, prop, evaluated_days=None, rng=None):
    cfg = ctx.cfg
    rng = rng or np.random.default_rng(cfg.seed + 2)
    T = ctx.tremor
    d = G.epi_km(T.lat.values, T.lon.values, prop.lat, prop.lon)
    tday = T.day.values.astype("datetime64[D]").astype(np.int64)     # integer day numbers throughout
    loc_days = np.unique(tday[d <= cfg.g3_local_km])
    far_days = np.unique(tday[(d >= cfg.g3_far_km[0]) & (d <= cfg.g3_far_km[1])])
    mem_days = (pd.to_datetime(prop.times, unit="s").values
                .astype("datetime64[D]").astype(np.int64))
    if evaluated_days is None:
        evaluated_days = np.arange(mem_days.min(), mem_days.max() + 1)
    else:
        evaluated_days = np.asarray(evaluated_days, dtype="datetime64[D]").astype(np.int64)
    span_d = float(len(evaluated_days))
    is_loc = np.isin(evaluated_days, loc_days)
    is_far = np.isin(evaluated_days, far_days)
    out = dict(gate="gate3", n_members=len(prop.times), span_days=span_d,
               n_local_tremor_days=int(is_loc.sum()), n_far_tremor_days=int(is_far.sum()))
    f_loc_cal = out["n_local_tremor_days"] / span_d if span_d else np.nan
    f_far_cal = out["n_far_tremor_days"] / span_d if span_d else np.nan
    f_loc_mem = float(np.isin(mem_days, loc_days).mean())
    f_far_mem = float(np.isin(mem_days, far_days).mean())
    out.update(frac_cal_local=f_loc_cal, frac_cal_far=f_far_cal,
               frac_mem_local=f_loc_mem, frac_mem_far=f_far_mem)
    out["E_local"] = float(f_loc_mem / f_loc_cal) if f_loc_cal > 0 else np.nan
    out["E_far"] = float(f_far_mem / f_far_cal) if f_far_cal > 0 else np.nan
    out["ratio"] = float(out["E_local"] / out["E_far"]) if (out["E_far"] and out["E_far"] > 0) else np.nan
    # ---- circular-shift null (>= 30 d) of the member time series on the evaluated calendar
    min_sh = cfg.g3_min_shift_d
    if span_d < 2 * min_sh + 1:
        out["verdict"] = UNTESTABLE
        out["reason"] = (f"evaluated span {span_d:.0f} d < 2 x {min_sh:.0f} d: a circular shift of "
                         f">= {min_sh:.0f} d does not exist on this calendar")
        return out
    d0 = int(evaluated_days[0])
    off = (mem_days - d0).astype(int) % int(span_d)
    E = []
    for _ in range(cfg.g3_n_shift):
        sh = int(rng.integers(min_sh, span_d - min_sh))
        sd = evaluated_days[(off + sh) % int(span_d)]
        E.append(float(np.isin(sd, loc_days).mean()) / f_loc_cal if f_loc_cal > 0 else np.nan)
    E = np.asarray(E, float)
    out["p_shift"] = float(((E >= out["E_local"]).sum() + 1) / (E.size + 1))
    out["null_E_local_med"] = float(np.median(E))
    ok = (np.isfinite(out["E_local"]) and out["E_local"] >= cfg.g3_min_elocal
          and np.isfinite(out["ratio"]) and out["ratio"] >= cfg.g3_min_ratio
          and out["p_shift"] < cfg.g3_p)
    out["verdict"] = PASS if ok else FAIL
    return out


# ================================================================================ driver
def evaluate(ctx, prop, gates=(0, 1, 2, 3), evaluated_days=None, rng=None):
    rng = rng or np.random.default_rng(ctx.cfg.seed)
    t0 = time.time()
    res = dict(pid=prop.pid, lat=prop.lat, lon=prop.lon, dep=prop.dep,
               n_members=int(len(prop.times)), meta=dict(prop.meta))
    if 0 in gates:
        res["gate0"] = gate0(ctx, prop, rng)
    if 1 in gates:
        res["gate1"] = gate1(ctx, prop, rng)
    if 2 in gates:
        res["gate2"] = gate2(ctx, prop, rng)
    if 3 in gates:
        res["gate3"] = gate3(ctx, prop, evaluated_days, rng)
    hard = [res[f"gate{i}"]["verdict"] for i in (0, 1, 2) if f"gate{i}" in res]
    if FAIL in hard:
        v = FAIL
    elif UNTESTABLE in hard:
        v = UNTESTABLE
    else:
        g3 = res.get("gate3", {}).get("verdict")
        v = "PASS_inter_ETS_candidate" if g3 == FAIL else PASS
    res["verdict"] = v
    res["runtime_s"] = time.time() - t0
    return res


def flatten(r):
    """One row per proposal for the results tables."""
    o = dict(pid=r["pid"], lat=r["lat"], lon=r["lon"], dep=r["dep"], n_members=r["n_members"],
             verdict=r["verdict"], runtime_s=r.get("runtime_s", np.nan))
    o.update({f"meta_{k}": v for k, v in r.get("meta", {}).items()})
    g = r.get("gate0")
    if g:
        o.update(g0=g["verdict"], g0_min_dist_km=g.get("min_dist_km"),
                 g0_nearest=g.get("nearest_bkey"), g0_n_sites_cc=g.get("n_sites_cc"),
                 g0_n_reg=g.get("n_sites_registering"),
                 g0_n_below=g.get("n_sites_below_chance"),
                 g0_frac_above=g.get("frac_reg_above_chance"),
                 g0_literal=g.get("distinct_by_waveform_literal"),
                 g0_max_ratio=g.get("max_site_ratio"),
                 g0_hf_median=g.get("hf_median"), g0_hf_n=g.get("hf_n"))
    g = r.get("gate1")
    if g:
        best = None
        for p in g["per_station"]:
            if best is None or (p.get("excess_db", -99) or -99) > (best.get("excess_db", -99) or -99):
                best = p
        o.update(g1=g["verdict"], g1_sites_pass=g.get("n_sites_pass"),
                 g1_best_sta=best.get("sta") if best else None,
                 g1_best_excess_db=best.get("excess_db") if best else np.nan,
                 g1_best_n=best.get("n_used") if best else 0,
                 g1_best_dist_km=best.get("dist_km") if best else np.nan,
                 g1_best_p_time=best.get("p_time") if best else np.nan,
                 g1_best_mo_p95=best.get("moveout_p95_excess_db") if best else np.nan,
                 g1_ctl_pass=(g.get("positive_control") or {}).get("any_pass"))
    g = r.get("gate2")
    if g:
        o.update(g2=g["verdict"], g2_n_sites_60km=g.get("n_sites_60km"), g2_gap=g.get("gap_deg"),
                 g2_dist_km=g.get("dist_to_proposal_km"), g2_null_ratio=g.get("null_ratio"),
                 g2_synth_km=g.get("synth_med_km"), g2_loso_km=g.get("loso_max_km"),
                 g2_free_dz=g.get("free_dz_km"))
    g = r.get("gate3")
    if g:
        o.update(g3=g["verdict"], g3_E_local=g.get("E_local"), g3_E_far=g.get("E_far"),
                 g3_ratio=g.get("ratio"), g3_p=g.get("p_shift"), g3_span_d=g.get("span_days"))
    return o


def load_proposals(path):
    """Proposal file: parquet/csv with columns pid, t (epoch s), lat, lon [, dep, offset_search_s]."""
    D = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path)
    out = []
    for pid, g in D.groupby("pid", sort=False):
        r0 = g.iloc[0]
        meta = {c: (int(r0[c]) if c in ("truth_bkey", "drop_bkey") else r0[c])
                for c in g.columns if c.startswith("meta_")
                or c in ("kind", "truth_bkey", "drop_bkey")}
        out.append(Proposal(pid=str(pid), times=g.t.to_numpy(float), lat=float(r0.lat),
                            lon=float(r0.lon),
                            dep=float(r0.dep) if "dep" in g.columns and np.isfinite(r0.dep) else np.nan,
                            offset_search_s=float(r0.offset_search_s) if "offset_search_s" in g.columns else 0.0,
                            meta=meta))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("proposals")
    ap.add_argument("--out", default=f"{G.GATES}/g_gates_out")
    ap.add_argument("--gates", default="0,1,2,3")
    ap.add_argument("--withheld", default="VGZ,SNB,B926")
    ap.add_argument("--n-time-null", type=int, default=200)
    ap.add_argument("--n-moveout-null", type=int, default=100)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-controls", action="store_true")
    a = ap.parse_args()
    cfg = GateConfig(withheld=tuple(a.withheld.split(",")), g1_n_time_null=a.n_time_null,
                     g1_n_moveout_null=a.n_moveout_null)
    ctx = GateContext(cfg)
    props = load_proposals(a.proposals)
    if a.limit:
        props = props[:a.limit]
    gates = tuple(int(x) for x in a.gates.split(","))
    rows, full = [], []
    t0 = time.time()
    for k, p in enumerate(props):
        r = evaluate(ctx, p, gates=gates)
        rows.append(flatten(r))
        full.append(r)
        print(f"[{k+1}/{len(props)}] {p.pid} N={len(p.times)} -> {r['verdict']} "
              f"({r['runtime_s']:.1f} s)", flush=True)
    R = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    R.to_csv(a.out + ".csv", index=False)
    json.dump(full, open(a.out + ".json", "w"), indent=1, default=str)
    print(f"[done] {len(R)} proposals in {time.time()-t0:.0f} s -> {a.out}.csv")
    print(R.verdict.value_counts().to_string())


if __name__ == "__main__":
    main()
