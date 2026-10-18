#!/usr/bin/env python
"""DUC_FIN step 2 -- re-template the NEW station epochs that start inside the data on disk (<= 2012-07), with the
procedure checked on the known answer (results/DUC_EXT_FIXED_RESULTS.md: re-template from our pass-1 detections made
WITHOUT that station; known-answer check median cc 0.694 with his template on the same channel).

It runs scripts/duc_fix_retemplate.py UNCHANGED (imported), finalising every new epoch that starts before 2012-08 with
the detections of 2008-01..2012-07 (the 24-month template window is cut at 2012-07: no later data on disk), into its own
folder OUT/retpl/ (results/duc_fix/mf1 and station_epochs.csv are symlinked READ-ONLY; duc_fix's own templates_new /
retemplate_state are not touched).  Keep / drop per (family, channel) = duc_fix_retemplate's rule: N >= 20 windows and
the stack coherence above all of >= 10 (of 20) scrambled controls (the same detection times on days shifted by 3-40 days).

Also writes OUT/retpl/epoch_summary.csv (kept / dropped per epoch and station, reasons), OUT/retpl/template_vs_his.csv
(each kept new template against HIS template of the old channel of that station: max cc and max |cc| within +-1 s and the
lag -> polarity / timing changes at the analogue -> digital switch) and OUT/status/retpl_0.csv (station-day read check:
every station-day of a new epoch's window that duc_fix's plan says is on disk must have a readable cache file and a prep
entry; guard > 3 %).

  python scripts/duc_fin_retpl.py [--workers 6] [--upto 2012-07]
"""
import os, sys, json, argparse, shutil
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
from t_common import log


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); ap.add_argument("--upto", default="2012-07")
    ap.add_argument("--epochs", default="", help="smoke: only these epochs")
    a = ap.parse_args()
    D = f"{F.OUT}/retpl"; os.makedirs(D, exist_ok=True)
    for nm in ["mf1", "station_epochs.csv"]:
        if not os.path.exists(f"{D}/{nm}"): os.symlink(f"{F.FIX}/{nm}", f"{D}/{nm}")
    import duc_fix_retemplate as R
    import duc_fix_clean as CL
    assert os.path.abspath(CL.PREP) == os.path.abspath(f"{F.FIX}/prep")
    E = DC.epochs()
    fp = DC.families(); used = set(",".join(fp.stations).split(","))
    mend = pd.Timestamp(a.upto + "-01") + pd.offsets.MonthBegin(1)
    newE = E[(E.status == "new") & E.sta.isin(used)].copy()
    newE["on_disk"] = newE.start < mend
    if a.epochs:
        keep = set(a.epochs.split(","))
        newE = newE[newE.epoch.isin(keep) | ~newE.on_disk]
    # ---- station-day read check (light: cache file opens, has the epoch's channels, prep lists the station)
    plan = F.plan_status(); rows = []
    for r in newE[newE.on_disk].itertuples():
        ns = f"{r.net}.{r.sta}"
        for d in pd.date_range(r.start, min(r.end, mend), freq="D", inclusive="left"):
            ds = str(d.date()); st = plan.get((ns, ds))
            if st not in ("built", "cached"):
                rows.append((ds, ns, "nodata", str(st))); continue
            p = DC.cache_path(r.net, r.sta, d)
            if not os.path.exists(p): rows.append((ds, ns, "missing", "no cache file")); continue
            try:
                z = np.load(p); ch = [str(c) for c in z["chans"]]
                if not set(r.chans.split(",")) <= set(ch): rows.append((ds, ns, "error", f"chans {ch}")); continue
            except Exception as e:
                rows.append((ds, ns, "error", str(e)[:80])); continue
            pp = f"{CL.PREP}/{ds}.npz"
            try:
                ok = os.path.exists(pp) and ns in set(str(s) for s in np.load(pp)["stas"])
            except Exception as e:
                ok = False
            rows.append((ds, ns, "ok" if ok else "missing", "" if ok else "station not in prep (short/zero day or prep gap)"))
    F.write_status("retpl", "0", rows)
    # ---- run duc_fix_retemplate on the redirected folder (forced finalisation at --upto)
    for f in ["templates_new.csv", "templates_new.npz", "retemplate_state.csv", "rerun_months.txt"]:
        if os.path.exists(f"{D}/{f}"): shutil.move(f"{D}/{f}", f"{D}/{f}.prev")
    R.FIX = D; R.STATE = f"{D}/retemplate_state.csv"; DC.FIX = D
    if a.epochs:
        _E = DC.epochs
        DC.epochs = lambda: _E()[lambda x: (x.status != "new") | x.epoch.isin(set(a.epochs.split(",")))]
        R.epochs = DC.epochs
    sys.argv = ["x", "--upto", a.upto, "--final", "--workers", str(a.workers)]
    R.main()
    DC.FIX = F.FIX
    # ---- summaries
    T = pd.read_csv(f"{D}/templates_new.csv", dtype={"fam": str})
    T["reason"] = np.where(T.keep == 1, "kept", np.where(T.N < R.NMIN, f"N<{R.NMIN} windows",
                   np.where(T.n_ctrl.fillna(0) < 10, "<10 controls", "stack not above every scrambled control")))
    T.to_csv(f"{D}/templates_new.csv", index=False)
    HT = DC.his_templates(); Z = np.load(f"{D}/templates_new.npz")
    cmp_ = []
    for r in T[T.keep == 1].itertuples():
        old = [(c, v) for (s2, c), v in HT.get(r.fam, {}).items() if s2 == r.sta and c[-1] == r.cha[-1]]
        if not old or r.key not in Z.files: continue
        c0, h = old[0]; st = Z[r.key]
        if h.size != st.size: continue
        a_ = (st - st.mean()) / (st.std() + 1e-30); b_ = (h - h.mean()) / (h.std() + 1e-30)
        cc = np.array([np.dot(np.roll(a_, s), b_) / a_.size for s in range(-20, 21)])
        cmp_.append(dict(fam=r.fam, epoch=r.epoch, sta=r.sta, cha=r.cha, his_cha=c0, N=r.N, cc_max=cc.max(),
                         lag_max_s=(np.argmax(cc) - 20) / 20.0, abscc_max=np.abs(cc).max(), lag_abs_s=(np.argmax(np.abs(cc)) - 20) / 20.0,
                         sign_at_abs=float(np.sign(cc[np.argmax(np.abs(cc))]))))
    Cm = pd.DataFrame(cmp_); Cm.to_csv(f"{D}/template_vs_his.csv", index=False)
    es = []
    for r in newE.itertuples():
        t = T[T.epoch == r.epoch]
        es.append(dict(epoch=r.epoch, sta=r.sta, start=str(r.start.date()), chans=r.chans, sr=r.sr,
                       change=r.change, on_disk=bool(r.on_disk), family_channels=len(t), kept=int(t.keep.sum()),
                       dropped=int(len(t) - t.keep.sum()), families_kept=int(t[t.keep == 1].fam.nunique()),
                       families=int(t.fam.nunique()), median_N=float(t.N.median()) if len(t) else np.nan,
                       median_cc_his=float(Cm[Cm.epoch == r.epoch].cc_max.median()) if len(Cm) and (Cm.epoch == r.epoch).any() else np.nan,
                       median_abscc_his=float(Cm[Cm.epoch == r.epoch].abscc_max.median()) if len(Cm) and (Cm.epoch == r.epoch).any() else np.nan,
                       decision=("RE-TEMPLATED" if t.keep.sum() > 0 else "DROPPED (no family-channel passes the scrambled control)")
                       if r.on_disk else "PENDING: starts after 2012-07, no data on disk (no downloads in this run)",
                       reasons=json.dumps(t.reason.value_counts().to_dict()) if len(t) else ""))
    ES = pd.DataFrame(es); ES.to_csv(f"{D}/epoch_summary.csv", index=False)
    od = ES[ES.on_disk]
    log(f"[retpl] epochs on disk {len(od)}: re-templated {int((od.kept > 0).sum())}, dropped {int((od.kept == 0).sum())}; "
        f"family-channels kept {int(od.kept.sum())}/{int(od.family_channels.sum())}; pending (no data) {int((~ES.on_disk).sum())}; "
        f"kept templates vs his old channel: median cc {Cm.cc_max.median() if len(Cm) else np.nan:.3f}, median |cc| "
        f"{Cm.abscc_max.median() if len(Cm) else np.nan:.3f}, polarity reversed (|cc| at negative cc) {int((Cm.sign_at_abs < 0).sum()) if len(Cm) else 0}/{len(Cm)}")
    F.guard("retpl")


if __name__ == "__main__":
    main()
