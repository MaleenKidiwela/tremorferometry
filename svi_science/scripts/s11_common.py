"""s11_common -- shared paths/constants/helpers for S11 (deep vs shallow dv/v with kernels).
Rules: svi_science/results/s11_depth/PLAN.md. eps = -dv/v. Wraps s10_* machinery (never edits it)."""
from __future__ import annotations
import os, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s11_depth"
TB, FIG, LOGD = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs"
MULTI = f"{ROOT}/svi_science/results/s10_dvv/multi"
MTB, MWIN = f"{MULTI}/tables", f"{MULTI}/win"
STATIONS = ["PB.B011", "PB.B009", "PB.B010", "PB.B926", "CN.PGC", "CN.VGZ", "CN.SNB"]
SITE = {"PB.B009": "PGC", "PB.B010": "PGC", "PB.B011": "PGC", "CN.PGC": "PGC", "PB.B926": "B926",
        "CN.VGZ": "VGZ", "CN.SNB": "SNB"}
SITES = ["PGC", "B926", "VGZ", "SNB"]
COLOC = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
# windows past the picked S (s); "25" = headline 2-5 s
WINS = {"25": (2.0, 5.0), "36": (3.0, 6.0), "2035": (2.0, 3.5), "3560": (3.5, 6.0)}
HEAD = "25"
RAIL = 0.0099
N_HEAD, CC_HEAD, N_COL = 50, 0.8, 20
SEED = 20261012
VEL = f"{ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv"
LAYERS = np.array([0.0, 0.1, 0.3, 1.0, 3.0, 5.0, 10.0, 25.0, 35.0, 45.0, np.inf])
LAYER_NAMES = ["0-0.1", "0.1-0.3", "0.3-1", "1-3", "3-5", "5-10", "10-25", "25-35", "35-45", ">45"]
ZONES = ["src", "sta", "rest"]
CATS = [f"{l}|{z}" for l in LAYER_NAMES for z in ZONES]       # layer-major order
GROUPS = {"s": [0, 1, 2, 3, 4], "m": [5, 6], "d": [7, 8], "b": [9]}   # layer indices


def group_coef(a, grp):
    """a (..., 30 categories) -> sum over the layers of a group (all zones)"""
    a = np.asarray(a).reshape(a.shape[:-1] + (len(LAYER_NAMES), len(ZONES)))
    return a[..., GROUPS[grp], :].sum((-1, -2))


def episodes():
    """episode units of S10 multi with a mid date"""
    u = pd.read_csv(f"{MTB}/episode_units.csv")
    pn = f"{TB}/ext/episode_units_new.csv"            # Revision 2: bursts formed only by the added families' days
    if os.path.exists(pn):
        u = pd.concat([u, pd.read_csv(pn)], ignore_index=True)
    u["mid"] = pd.to_datetime(u.start) + (pd.to_datetime(u.end) - pd.to_datetime(u.start)) / 2
    return u.set_index("episode")


def annual_fit(t, y, s):
    """WLS y = c + A cos + B sin (annual); t datetime; -> amp, sig_amp, peak_doy, chi2red"""
    t = pd.to_datetime(pd.Series(t)); ph = 2 * np.pi * (t.dt.dayofyear.to_numpy() - 1) / 365.25
    ok = np.isfinite(y) & np.isfinite(s) & (s > 0)
    if ok.sum() < 5:
        return dict(amp=np.nan, amp_sig=np.nan, peak_doy=np.nan, chi2red=np.nan, n=int(ok.sum()))
    G = np.column_stack([np.ones(ok.sum()), np.cos(ph[ok]), np.sin(ph[ok])]); w = 1 / s[ok] ** 2
    N = G.T @ (G * w[:, None]); m = np.linalg.solve(N, G.T @ (w * y[ok])); C = np.linalg.inv(N)
    r = y[ok] - G @ m; chi2 = float(np.sum(w * r ** 2) / max(ok.sum() - 3, 1))
    C = C * max(chi2, 1.0)
    A = float(np.hypot(m[1], m[2]))
    J = np.array([m[1], m[2]]) / max(A, 1e-12)
    sA = float(np.sqrt(J @ C[1:, 1:] @ J))
    pk = (np.degrees(np.arctan2(m[2], m[1])) % 360) / 360 * 365.25 + 1
    return dict(amp=A, amp_sig=sA, peak_doy=float(pk), chi2red=chi2, n=int(ok.sum()), c=float(m[0]), a_cos=float(m[1]),
                a_sin=float(m[2]))


class Log:
    def __init__(self, path):
        self.f = open(path, "a")

    def __call__(self, m):
        s = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
        print(s, flush=True); self.f.write(s + "\n"); self.f.flush()


def ensure():
    for d in (TB, FIG, LOGD):
        os.makedirs(d, exist_ok=True)


def pt_now():
    import subprocess
    return subprocess.run(["date", "+%Y-%m-%d %-I:%M %p PT"], env={"TZ": "America/Los_Angeles"}, capture_output=True,
                          text=True).stdout.strip()


# ------------------------------------------------------------------------------ region switch (PLAN §5, step 7)
# S11_REGION unset or "SVI" -> everything above, unchanged. "WA_NW" / "WA_SW" -> the same pipeline on Washington: every
# input/output under OUT/<wa_nw|wa_sw>/ (tables, logs, figures, DONE markers); the WA inventory (s11_wa_inventory.py) writes
# the region tables that replace the S10-multi ones (MTB/MWIN) and the station/site table read here.
REGION = os.environ.get("S11_REGION", "SVI").upper()
DEC1 = ("2004-01-01", "2014-12-31")                    # Q2 first decade (SVI)
MAJ_LIST = None                                          # Q3 major list (SVI: s11_q3.MAJ15)
if REGION != "SVI":
    _R = {"WA_NW": dict(sub="wa_nw", vel="royer_NW_WA", cat="NW", famreg="Royer_NW_WA"),
          "WA_SW": dict(sub="wa_sw", vel="royer_SW_WA", cat="SW", famreg="Royer_SW_WA")}[REGION]
    REGION_CFG = _R
    OUT = f"{OUT}/{_R['sub']}"
    TB, FIG, LOGD = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs"
    MULTI = OUT
    MTB, MWIN = f"{TB}/inv", f"{TB}/inv/win"
    VEL = f"{ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/{_R['vel']}.csv"
    DEC1 = ("2007-01-01", "2014-12-31")                  # WA repeats start in 2007 (PLAN §5)
    N_HEAD = 20                                          # WA headline cell rule n >= 20 (PLAN §5)
    STATIONS, SITE, SITES, COLOC = [], {}, [], []
    _st = f"{MTB}/stations.csv"
    if os.path.exists(_st):
        _S = pd.read_csv(_st)
        _S = _S[_S.use.astype(bool)] if "use" in _S else _S
        STATIONS = list(_S.netsta)
        SITE = dict(zip(_S.netsta, _S.site))
        SITES = sorted(set(SITE.values()))
    _mj = f"{TB}/maj_list.csv"
    if os.path.exists(_mj):
        MAJ_LIST = list(pd.read_csv(_mj).episode.astype(str))
