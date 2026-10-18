#!/usr/bin/env python
"""Assemble the FAMILY-DAY table for the daily fixed-effects joint inversion (Merlin-audited spec,
2026-08-22 rounds 1+2: family-level offsets fleet-wide — see notes/2026-08-22_daily_joint_inversion_design.md).

Rows = one per (family, day): NO aggregation, NO gap filling — a day with no stack has no row (gaps
stay gaps, exactly the series the HTML map shows). Certified families only, cc_max > 0.70.
Each family maps to its 0.1-deg-snapped cell, so (tag, cell) still indexes the SAME G.npz kernel rows
(families within a cell share the cell's kernel; the offset c is per family x era in the solver).

Gates (B3, all row-DROPS, counts logged):
  rail        |eps| >= 0.99 x the file's search edge (per file: orig 0.02, unif05 0.05, eps05 recovered)
  cycle-skip  |eps| > 0.02 AND cc < 0.80 (wide-range files only — the map's gate)
  outlier     |eps - 90d-rolling-median| > 8 x robust MAD of the family's detrended series
Era id (I6, exact asymmetric contamination): a daily value at t blends two instruments iff
t in [boundary, boundary+30 d) (trailing stack) -> those rows DROPPED; era = #boundaries <= t.
Family-era minimum: family-eras with < 30 surviving rows dropped (a loose offset would soak field signal).

Writes (RESDIR): family_days.parquet [tag, fam, cell, date, dvv(%), cc, n_stack, era, src]
                 families_daily.csv  per-(tag,fam,era) summary: n, span, sigma (robust MAD), railed_frac
CSV preference per station: _unif05 > _eps05 > original (DAILY_FILES=original forces originals).
Run with the same RESDIR/CELL_DEG as the target grid.
"""
import os, glob
import numpy as np, pandas as pd

CC_MIN = 0.70
CELL_DEG = float(os.environ.get("CELL_DEG", "0.10"))
OUT = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
FILES = os.environ.get("DAILY_FILES", "prefer-unif05")
STEP = 2e-4                      # grid step (0.02/201 = 0.05/501)
MIN_FAM_ERA = 30                 # Merlin round-2: min rows per family-era, else drop its rows
MAD_K = 8.0                      # outlier gate, x robust MAD of detrended family series


def cert_fams(stem):
    f1 = f"data/{stem}_causality_cert.csv"
    out = set()
    if os.path.exists(f1):
        c = pd.read_csv(f1); out = set(c[c.reliable].fam)
    else:
        f2 = f"data/{stem}_fwd_vs_rev_coda.csv"
        if os.path.exists(f2):
            c = pd.read_csv(f2); out = set(c[c.ratio > 1.5].fam)
    bp = "data/family_user_blocklist.csv"
    if out and os.path.exists(bp):
        b = pd.read_csv(bp)
        out -= set(b[b.stem == stem].fam)          # user-vetoed families (2026-08-31 mechanism)
    return out


def famll(pid):
    a = pid.split("__")[0].split("_")
    return float(a[0]), float(a[1])


def snap(v):
    return round(round(v / CELL_DEG) * CELL_DEG, 2)


def pick_file(tag):
    u = f"data/daily_dvv_{tag}_Z_2to4_unif05.csv"
    e = f"data/daily_dvv_{tag}_Z_2to4_eps05.csv"
    o = f"data/daily_dvv_{tag}_Z_2to4.csv"
    if FILES != "original":
        if os.path.exists(u): return u, "unif05", 0.05
        if os.path.exists(e): return e, "eps05", None
    return o, "orig", 0.02


# era boundaries per tag (used-band physical response changes).
# 2026-08-31 USER RULING: era-split ONLY stations where a visible cross-era dv/v change was
# confirmed by visual review (data/era_split_whitelist.txt). All other stations are treated as
# single-era — their metadata changes were judged not to produce a measurable step, so no
# straddle drops and one offset c per family.
et = pd.read_csv(f"{OUT}/era_table.csv")
_wl = set()
if os.path.exists("data/era_split_whitelist.txt"):
    _wl = {x.strip().upper() for x in open("data/era_split_whitelist.txt") if x.strip()}
BND = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan")
       for _, r in et.iterrows()
       if str(r.tag).upper().replace("P90F40", "") in _wl}

rows, fam_rows = [], []
drops = dict(rail=0, cycleskip=0, outlier=0, fam_era_min=0, straddle=0)
tags = sorted({os.path.basename(f).split("daily_dvv_")[1].split("_Z_2to4")[0]
               for f in glob.glob("data/daily_dvv_*_Z_2to4*.csv")})
n_sta, src_count = 0, {}
for tag in tags:
    if "rev" in tag.lower() or "MIRROR" in tag:
        continue
    stem = tag.lower().replace("p90f40", "")
    fams = cert_fams(stem)
    if not fams:
        continue
    fp, src, emax = pick_file(tag)
    if not os.path.exists(fp):
        continue
    d = pd.read_csv(fp, usecols=["patch", "date", "dvv", "cc_max", "n_stack"])
    d = d[(d.cc_max > CC_MIN) & (d.patch.isin(fams))]
    if not len(d):
        continue
    if emax is None:
        emax = float(np.ceil(d.dvv.abs().max() / STEP) * STEP)   # eps05: recover per-station range
    d["date"] = pd.to_datetime(d["date"], errors="coerce")
    d = d.dropna(subset=["date"]).sort_values(["patch", "date"])
    # user date-range veto (2026-08-31 mechanism): drop this station's rows inside vetoed windows
    dbp = "data/station_date_blocklist.csv"
    if os.path.exists(dbp):
        for _, br in pd.read_csv(dbp).iterrows():
            if br.stem == stem:
                k = (d.date >= pd.Timestamp(br.t0)) & (d.date <= pd.Timestamp(br.t1))
                drops["date_veto"] = drops.get("date_veto", 0) + int(k.sum())
                d = d[~k]
    # --- B3a rail drop (per-file edge) ---
    k = d.dvv.abs() >= 0.99 * emax
    drops["rail"] += int(k.sum()); d = d[~k]
    # --- B3b cycle-skip gate (wide-range files: freed values with weak cc are skips, not velocity) ---
    if emax > 0.021:
        k = (d.dvv.abs() > 0.02) & (d.cc_max < 0.80)
        drops["cycleskip"] += int(k.sum()); d = d[~k]
    # --- era id + asymmetric straddle drop: value at t blended iff t in [b, b+30d) ---
    bs = BND.get(tag, [])
    if bs:
        t = d.date.values
        era = np.zeros(len(d), np.int16); bad = np.zeros(len(d), bool)
        for b in bs:
            b64 = np.datetime64(b)
            era += (t >= b64).astype(np.int16)
            bad |= (t >= b64) & (t < b64 + np.timedelta64(30, "D"))
        drops["straddle"] += int(bad.sum())
        d = d[~bad].copy(); d["era"] = era[~bad]
    else:
        d["era"] = np.int16(0)
    if not len(d):
        continue
    d["dvv"] = d.dvv * 100.0                                     # percent (eps convention unchanged)
    # --- B3c robust outlier gate on the family's detrended series ---
    keep_parts = []
    for fam, g in d.groupby("patch", sort=False):
        ser = g.set_index("date").dvv
        det = (ser - ser.rolling("90D").median()).values          # trailing 90-d median (gap-aware)
        mad = 1.4826 * np.nanmedian(np.abs(det - np.nanmedian(det)))
        if np.isfinite(mad) and mad > 0:
            k = np.abs(det - np.nanmedian(det)) > MAD_K * mad
            drops["outlier"] += int(np.nansum(k)); g = g[~k]
        sig = float(mad) if np.isfinite(mad) and mad > 0 else float(np.nanstd(ser.values)) or 0.2
        for era, ge in g.groupby("era"):
            if len(ge) < MIN_FAM_ERA:
                drops["fam_era_min"] += len(ge); continue
            keep_parts.append(ge)
            fam_rows.append(dict(tag=tag, fam=fam, era=int(era), src=src, n=len(ge),
                                 t0=str(ge.date.min().date()), t1=str(ge.date.max().date()),
                                 sigma=round(sig, 4)))
    if not keep_parts:
        continue
    g = pd.concat(keep_parts)
    la, lo = zip(*g.patch.map(famll))
    g["cell"] = [f"{snap(a):.2f}_{snap(b):.2f}" for a, b in zip(la, lo)]
    g = g.rename(columns={"patch": "fam", "cc_max": "cc"})
    g["tag"] = tag; g["src"] = src
    g["dvv"] = g.dvv.round(4); g["cc"] = g.cc.round(3)
    rows.append(g[["tag", "fam", "cell", "date", "dvv", "cc", "n_stack", "era", "src"]])
    n_sta += 1; src_count[src] = src_count.get(src, 0) + 1

fd = pd.concat(rows, ignore_index=True)
fams = pd.DataFrame(fam_rows)
os.makedirs(OUT, exist_ok=True)
fd.to_parquet(f"{OUT}/family_days.parquet", index=False)
fams.to_csv(f"{OUT}/families_daily.csv", index=False)
print(f"stations {n_sta} (sources {src_count}) | family-eras {len(fams)} | family-days {len(fd):,} | "
      f"span {fd.date.min().date()} -> {fd.date.max().date()}")
print(f"drops: {drops}")
print(f"sigma (family robust MAD) median {fams.sigma.median():.3f}% | families {fd.fam.nunique()} | cells {fd.cell.nunique()}")
print(f"wrote {OUT}/family_days.parquet, families_daily.csv")
