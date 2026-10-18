#!/usr/bin/env python
"""PNET check, step 2: selection-free 'outside the 12-s template window' score per candidate and real-share estimates.

out_cc  = correlation of a single candidate window OUTSIDE the template window ([-8, 0) s: the P wave; [12, 24) s: coda)
          with an independent reference stack of the same family, station and template key; coef-weighted over the
          candidate's stations (identical to quietwave/scripts/qw_extra.py).  Also the P part and the coda part alone.
Reference stacks (built here from the new windows, both models must agree so neither model is favoured):
          NW: tremor-day data with P_own_WAVE >= 0.85 AND P_own_WAVE_PNET >= 0.85;  SW / Arm B: both >= 0.5.
          Two-fold by day-of-month parity: tremor-day rows are scored against the stack of the OTHER parity; quiet and
          unknown rows (other days) against the stack of all tremor days.  >= 5 windows per stack.
Real share (rough, family-matched, 90 % day-block bootstrap; random times = 0):
  same_rule   set mean / tremor-day mean of rows passing the SAME rule (the quiet-day test's definition: 1 = as real as
              tremor-day detections selected the same way).  Families with >= 10 (NW) / >= 5 (SW, Arm B) tremor rows.
  vs_both     set mean / tremor-day mean of rows BOTH models pass at 0.75 (a common yardstick; low-strength real events
              score below 1 on it, so it is a floor for weak sets)."""
import os, sys, glob, json
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_extension_wa/scripts")
import wa_prod_common as P
W = P.W
OD = f"{P.OUT}/quietwave/pnet_check"; WK = f"{OD}/work"; TB = f"{OD}/tables"
os.makedirs(TB, exist_ok=True)
n = W.NT; FS = W.FS
PRE = int(20 * FS)
SEGS = {"P": slice(PRE - int(8 * FS), PRE), "coda": slice(PRE + n, PRE + n + int(12 * FS))}
S = pd.read_parquet(f"{OD}/rows.parquet")
S["par"] = np.where(S.date.str[-2:].astype(int) % 2 == 0, "even", "odd")
S["w"] = S.P_own_WAVE.fillna(0); S["p"] = S.P_own_WAVE_PNET.fillna(0)
thr = np.where(S.reg == "NW", 0.85, 0.5)
S["is_ref"] = (S.series == "data") & S.tremor30 & (S.w >= thr) & (S.p >= thr)
REFROWS = set(S.row_id[S.is_ref])
PAR = dict(zip(S.row_id, S.par)); FID = dict(zip(S.row_id, S.fid))


def tkeys(ds):
    return json.load(open(f"{P.OUT}/mf/prod/run_{ds}.json"))["tkeys"]


def pass1(ds):
    """partial reference sums per (fid, sta, tkey, parity)."""
    a = pd.read_parquet(f"{WK}/{ds}_win.parquet")
    a = a[a.row_id.isin(REFROWS) & (a.widx >= 0)]
    if not len(a):
        return {}
    wav = np.load(f"{WK}/{ds}.npz")["wav"]; tk = tkeys(ds); out = {}
    for r in a.itertuples():
        x = wav[r.widx].astype(np.float32)
        if not np.isfinite(x).all():
            continue
        k = (FID[r.row_id], r.sta, tk.get(r.sta), PAR[r.row_id])
        s, c = out.get(k, (0, 0)); out[k] = (s + x, c + 1)
    return out


REF = None


def setref(R):
    global REF
    REF = R


def pass2(ds):
    a = pd.read_parquet(f"{WK}/{ds}_win.parquet")
    a = a[a.widx >= 0]
    if not len(a):
        return None
    wav = np.load(f"{WK}/{ds}.npz")["wav"]; tk = tkeys(ds)
    rec = []
    for r in a.itertuples():
        f = FID[r.row_id]; t = tk.get(r.sta)
        par = PAR[r.row_id] if r.row_id in TREM else "all"
        ref = REF.get((f, r.sta, t, {"even": "odd", "odd": "even", "all": "all"}[par]))
        if ref is None:
            continue
        x = wav[r.widx].astype(np.float32)
        if not np.isfinite(x).all():
            continue
        ch = (np.abs(x).sum(1) > 0) & (np.abs(ref).sum(1) > 0)
        if not ch.any():
            continue
        v = {}
        for lab, segs in (("out_cc", ("P", "coda")), ("out_P", ("P",)), ("out_coda", ("coda",))):
            A = np.concatenate([x[ch][:, SEGS[s]].ravel() for s in segs]); B = np.concatenate([ref[ch][:, SEGS[s]].ravel() for s in segs])
            den = np.sqrt(np.dot(A, A) * np.dot(B, B)); v[lab] = np.dot(A, B) / den if den > 0 else np.nan
        rec.append((r.row_id, r.coef, v["out_cc"], v["out_P"], v["out_coda"], r.rms_sig / r.rms_pre if r.rms_pre > 0 else np.nan))
    return pd.DataFrame(rec, columns=["row_id", "coef", "out_cc", "out_P", "out_coda", "esnr"])


TREM = set(S.row_id[S.tremor30])

if __name__ == "__main__":
    days = sorted(S.date.unique())
    # replication check
    ROW = pd.concat([pd.read_parquet(p) for p in glob.glob(f"{WK}/*_row.parquet")], ignore_index=True)
    S = S.merge(ROW, on="row_id", how="left")
    chk = S[(S.series != "random") & (S.n_sta_used > 0)]
    d = np.abs(chk.net_re - chk.net_cc)
    rep = dict(rows=len(S), rows_cut=int((S.n_sta_used > 0).sum()), checked=len(chk), max_abs_diff=float(d.max()),
               frac_match_1e3=float((d < 1e-3).mean()))
    print(rep, flush=True); json.dump(rep, open(f"{TB}/replication_check.json", "w"), indent=1)
    with Pool(16) as pool:
        parts = pool.map(pass1, days, chunksize=8)
    acc = {}
    for d_ in parts:
        for k, (s, c) in d_.items():
            s0, c0 = acc.get(k, (0, 0)); acc[k] = (s0 + s, c0 + c)
    R = {}
    for (f, s, t, par), (sm, c) in acc.items():
        a0, c0 = R.get((f, s, t, "all"), (0, 0)); R[(f, s, t, "all")] = (a0 + sm, c0 + c)
        R[(f, s, t, par)] = (sm, c)
    R = {k: sm / c for k, (sm, c) in R.items() if c >= 5}
    print("reference stacks", len(R), flush=True)
    del parts, acc
    with Pool(16, initializer=setref, initargs=(R,)) as pool:
        out = pool.map(pass2, days, chunksize=8)
    Wn = pd.concat([o for o in out if o is not None], ignore_index=True)
    C = Wn.groupby("row_id").apply(lambda x: pd.Series(dict(
        out_cc=np.average(x.out_cc, weights=x.coef), out_P=np.average(x.out_P, weights=x.coef),
        out_coda=np.average(x.out_coda, weights=x.coef), esnr=np.nanmedian(x.esnr), n_ref_sta=len(x)))).reset_index()
    C = C.merge(S[["row_id", "date", "par", "reg", "cls", "series", "sel", "fid", "w", "p", "net_cc", "net_mad"]], on="row_id")
    C["e_snr"] = C.esnr.clip(0, 3) - 1
    C.to_parquet(f"{TB}/per_candidate.parquet", index=False)
    print(C.groupby(["reg", "cls", "series", "sel"]).size().to_string(), flush=True)
