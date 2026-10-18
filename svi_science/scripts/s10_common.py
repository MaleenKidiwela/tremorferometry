"""s10_common -- shared constants for the S10 dv/v pilot (rules: svi_science/results/s10_dvv/pilot/PLAN.md)."""
from __future__ import annotations
import os, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s10_dvv/pilot"
TB, FIG, WIN, LOGD = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/win", f"{OUT}/logs"
CUT = f"{ROOT}/lfe_extension/results/cut/svi_ext"
RAW_BW = f"{ROOT}/lfe_extension/traces_backward"
RAW_FW = f"{ROOT}/lfe_extension/traces"
FEAT_FW = f"{ROOT}/lfe_extension/results/refine/feat"
FEAT_BW = f"{ROOT}/svi_science/results/backward/mf/feat"
FAMS = (10, 87, 17)
FS = 40.0                       # stack sampling rate
T0W, T1W = -14.0, 10.0          # stored window about expected S (s)
NW = int(round((T1W - T0W) * FS)) + 1
TAU = T0W + np.arange(NW) / FS
CUT_PRE, CUT_POST = 15.0, 30.0  # cut-store span about expected S (raw reads use the same span -> identical filtering)
FS_DET = 20.0                   # detector frame rate (svi_ref_common.FS)
SITE = {"B009": "PGC", "B010": "PGC", "B011": "PGC", "B005": "LRIV", "B006": "LRIV", "B007": "LRIV"}
POLARIS = ["C8.KHVB", "C8.MGCB", "C8.SHDB", "C8.SHVB", "C8.TWBB", "PO.KELB", "PO.KLNB", "PO.SILB", "PO.SSIB", "PO.TSJB",
           "PO.TWGB", "PO.TWKB", "CN.LZB"]
CORE = ["PB.B011", "PB.B926", "PB.B009", "PB.B010", "CN.PGC", "CN.VGZ", "CN.SNB"]
POL_END = "2007-01-01"          # POLARIS stations used 2003-2006 only

# response epochs, end exclusive (PLAN §3 and §9)
EPOCHS = {
    "PB.B009": [("E1", "2005-11-02", "2100-01-01")], "PB.B010": [("E1", "2005-11-02", "2100-01-01")],
    "PB.B011": [("E1", "2005-11-02", "2100-01-01")], "PB.B926": [("E1", "2007-11-13", "2100-01-01")],
    "CN.PGC": [("E1", "2000-06-28", "2015-01-12"), ("E2", "2015-01-12", "2017-08-10"), ("E3", "2017-08-10", "2100-01-01")],
    "CN.VGZ": [("E1", "2003-02-27", "2007-02-28"), ("E2", "2007-02-28", "2009-07-31"), ("E3", "2009-07-31", "2017-03-03"),
               ("E4", "2017-03-03", "2100-01-01")],
    "CN.SNB": [("E1", "2000-05-25", "2006-04-20"), ("E2", "2006-04-20", "2008-12-11"), ("E3", "2008-12-11", "2011-03-17"),
               ("E4", "2011-03-17", "2017-06-08"), ("E5", "2017-06-08", "2021-03-28")],
    "C8.SHVB": [("E1", "2003-02-27", "2004-09-30")],
    "C8.MGCB": [("E0", "2003-06-27", "2003-08-11"), ("E1", "2003-08-11", "2006-05-13")],
    "C8.TWBB": [("E0", "2002-11-27", "2003-07-30"), ("E1", "2003-07-30", "2006-06-20")],
    "C8.SHDB": [("E1", "2003-02-28", "2005-01-11")],
    "C8.KHVB": [("E1", "2003-02-28", "2004-08-05")],
    "PO.TWKB": [("E0", "2002-11-20", "2003-08-01"), ("E1", "2003-08-01", "2006-05-14")],
    "PO.TSJB": [("E0", "2003-03-04", "2003-07-31"), ("E1", "2003-07-31", "2006-06-19")],
    "PO.KELB": [("Ea", "2002-09-12", "2003-07-31"), ("Eb", "2003-07-31", "2004-06-04")],
    "PO.KLNB": [("E1", "2004-06-04", "2006-06-22")],
    "PO.SILB": [("E0", "2003-02-22", "2003-08-03"), ("E1", "2003-08-03", "2006-07-04")],
    "PO.SSIB": [("E0", "2002-12-04", "2003-07-25"), ("E1", "2003-07-25", "2006-06-16")],
    "PO.TWGB": [("E0", "2002-11-25", "2003-08-02"), ("E1", "2003-08-02", "2005-05-10")],
    "CN.LZB": [("E1", "2001-01-10", "2006-07-01")],
}
# instrument / network changes drawn on every panel (gate f)
CHANGES = [("2014-06-17", "LZB ends"), ("2015-01-12", "PGC response sub-change"), ("2017-02-23", "NLLB BH->HH"),
           ("2017-03-03", "VGZ sensor swap"), ("2017-06-08", "SNB sensor swap"), ("2017-08-10", "PGC BH->HH"),
           ("2019-05-30", "LRIV BH->HH"), ("2021-03-28", "SNB ends"), ("2006-07-01", "era->frozen detector"),
           ("2014-01-01", "backward->forward catalogue"), ("2007-02-28", "VGZ CMG3ESP->CMG3T"),
           ("2006-04-20", "SNB CMG40T->CMG3ESP"), ("2008-12-11", "SNB digitiser"), ("2009-07-31", "VGZ BH->HH"),
           ("2011-03-17", "SNB BH->HH"), ("2003-07-31", "KELB gain change")]
CHANGES = sorted(CHANGES)


def epoch_of(netsta, dates):
    d = pd.Series(np.asarray(dates, dtype=str))
    out = np.full(len(d), "", dtype=object)
    for eid, a, b in EPOCHS.get(netsta, []):
        out[((d >= a) & (d < b)).to_numpy()] = eid
    if netsta in POLARIS:
        out[(d >= POL_END).to_numpy()] = ""
    return out


class Log:
    def __init__(self, path):
        self.f = open(path, "a")

    def __call__(self, m):
        s = f"[{time.strftime('%H:%M:%S')}] {m}"
        print(s, flush=True); self.f.write(s + "\n"); self.f.flush()


def ensure():
    for d in (TB, FIG, WIN, LOGD):
        os.makedirs(d, exist_ok=True)
