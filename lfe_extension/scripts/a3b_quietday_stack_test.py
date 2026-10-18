#!/usr/bin/env python
"""Track A step 3b (user amendment 2026-09-16) -- are the quiet-day detections false alarms, or REAL inter-ETS repeats?

Quiet-day detections are an UPPER BOUND on the false-alarm rate, not a measurement: families do fire between ETS
episodes, and Bostock's catalog (which has no events outside ETS windows) cannot say otherwise.  The test, per family:
stack the 3-C 20 Hz waveform at a station that was EXCLUDED from that family's network sum, at the family's quiet-day
detection times, aligned on the predicted S (origin + ttS(Crosson) + the frozen template's S-pick refinement at that
station), and compare the S-window energy and the cc against the family's template with a MATCHED RANDOM-TIME stack
(same station, same quiet days, same N, --draws realisations).  ETS-day own-matched detections are the positive control.
Withheld stations: VGZ, SNB, B926 (113 of the 116 templated families have a usable template at one of them), taken in
order of the frozen S-SNR; the detections come from a run with all three EXCLUDED from the network sum, so the withheld
station never contributed to the detection it is asked to confirm.
Verdicts: PASS = S-SNR excess >= 3 dB AND |envelope peak - predicted S| <= 2 s AND null p < 0.01 (real repeats);
FAIL = excess <= 1 dB; INSUFFICIENT = fewer than --min-n detections.  Output: results/a3b_quietday_stack_test.csv."""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e5r_templates as T  # noqa: E402
import d3_vel1d as v1d  # noqa: E402

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
FS = T.FS
QUIET = list(range(228, 243))
WITHHELD = ["VGZ", "SNB", "B926"]


def snr_db(st, s_idx, half_s=2.0):
    n = int(round(2 * half_s * FS)); a = s_idx - int(round(half_s * FS))
    sig = st[:, a:a + n]; noi = st[:, 0:n]
    return float(10 * np.log10(max((sig ** 2).mean(), 1e-30) / max((noi ** 2).mean(), 1e-30)))


def peak_off(st, s_idx, search_s=6.0):
    env = T.envelope(st); w = int(round(search_s * FS)); a, b = max(s_idx - w, 0), min(s_idx + w, len(env))
    return (a + int(np.argmax(env[a:b])) - s_idx) / FS


def cc_tpl(st, s_idx, tpl, lag_s=1.0):
    """max 3-C cc of the stack against the family's own template (L2-normalised per component) over +-lag_s."""
    nL = tpl.shape[1]; i0 = s_idx - int(round(1.0 * FS)); lag = int(round(lag_s * FS))
    best = -9.0
    for k in range(-lag, lag + 1):
        a = i0 + k
        if a < 0 or a + nL > st.shape[1]:
            continue
        seg = st[:, a:a + nL]
        nn = np.linalg.norm(seg, axis=1)
        if np.any(nn <= 0):
            continue
        best = max(best, float((seg / nn[:, None] * tpl).sum() / tpl.shape[0]))
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="a3b_noVSB")
    ap.add_argument("--thresholds", default="8.25,8.75")
    ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--draws", type=int, default=120)
    ap.add_argument("--min-n", type=int, default=20)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--hf-parquet", default=""); ap.add_argument("--hf-cut", type=float, default=3.6384)
    ap.add_argument("--out", default=f"{OUT}/a3b_quietday_stack_test.csv")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    S = pd.read_csv(f"{IN}/e5r_stations.csv").set_index("sta")
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv"); qc = qc[qc.usable & qc.template]
    z = np.load(f"{IN}/e5r_templates_2012.npz", allow_pickle=True); TPL = z["tpl"].item()
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet"); gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    pd_ = gt.groupby("doy").size(); ETS = set(pd_[pd_ >= 100].index)
    det = pd.read_parquet(f"{OUT}/e5r_detections_{a.tag}.parquet"); det = det[det.L == a.L]
    if a.hf_parquet:
        H = pd.read_parquet(a.hf_parquet); H = H[H.L == a.L][["doy", "fam", "t", "hfr_PGC"]]
        det = det.merge(H, on=["doy", "fam", "t"], how="left")
        v = det.hfr_PGC.to_numpy(); det = det[~(np.isfinite(v) & (v > a.hf_cut))]
    THR = [float(x) for x in a.thresholds.split(",")]
    # per family: the withheld station with the best frozen S-SNR among VGZ/SNB/B926
    pick_sta = {}
    for f, g in qc[qc.sta.isin(WITHHELD)].groupby("bkey"):
        r = g.sort_values("snr_db", ascending=False).iloc[0]
        pick_sta[int(f)] = (r.sta, float(r.s_pick_s))
    n = int(round(2 * T.HALF_W * FS)); idx = int(round(T.HALF_W * FS))
    DATA = {s: T.load_days(S.loc[s, "netsta"], 2012) for s in WITHHELD}
    rows = []; STK = {}
    for s in WITHHELD:
        D = DATA[s]
        qd = [d for d in QUIET if d in D]
        t_lo = min(D[d][0] for d in qd); t_hi = max(D[d][0] + D[d][1].shape[1] / FS for d in qd)
        fams = [f for f, (ss, _) in pick_sta.items() if ss == s]
        print(f"[{s}] {len(fams)} families, quiet days {len(qd)}", flush=True)
        for f in fams:
            fr = F.loc[f]; pick = pick_sta[f][1]
            dx = T.epi_km(fr.lat, fr.lon, float(S.loc[s, "lat"]), float(S.loc[s, "lon"]))
            ttS = float(v1d.tt(dx, fr.dep, "S"))
            tp = TPL.get((float(f), s), {}).get(8.0)
            for Tt in THR:
                dq = det[(det.fam == float(f)) & det.doy.isin(QUIET) & (det.net_mad >= Tt)]
                de = det[(det.fam == float(f)) & det.doy.isin(ETS) & det.match_own & (det.net_mad >= Tt)]
                for name, tt in (("quiet", np.sort(dq.t.to_numpy())), ("ets_control", np.sort(de.t.to_numpy()))):
                    base = dict(fam=f, sta=s, T=Tt, set=name, n_det=len(tt))
                    if len(tt) < a.min_n:
                        rows.append(dict(**base, verdict="INSUFFICIENT")); continue
                    W, _ = T.cut_windows(D, tt + ttS + pick - T.HALF_W, n)
                    if W is None or len(W) < a.min_n:
                        rows.append(dict(**base, verdict="INSUFFICIENT (no data)")); continue
                    st = W.mean(0); sn = snr_db(st, idx); pk = peak_off(st, idx)
                    cc = cc_tpl(st, idx, tp) if tp is not None else np.nan
                    nul, ncc = [], []
                    for _ in range(a.draws):
                        Wn, _ = T.cut_windows(D, rng.uniform(t_lo, t_hi, len(W)) + ttS + pick - T.HALF_W, n)
                        if Wn is not None and len(Wn) >= a.min_n:
                            stn = Wn.mean(0); nul.append(snr_db(stn, idx))
                            if tp is not None:
                                ncc.append(cc_tpl(stn, idx, tp))
                    nul = np.array(nul); ncc = np.array(ncc)
                    exc = sn - nul.mean() if len(nul) else np.nan
                    p = float((nul >= sn).sum() + 1) / (len(nul) + 1) if len(nul) else np.nan
                    pcc = float((ncc >= cc).sum() + 1) / (len(ncc) + 1) if len(ncc) else np.nan
                    if name == "quiet":
                        vd = ("PASS" if (exc >= 3.0 and abs(pk) <= 2.0 and p < 0.01) else
                              "FAIL" if exc <= 1.0 else "inconclusive")
                    else:
                        vd = "control: registers" if (exc >= 3.0 and abs(pk) <= 2.0) else "control: weak"
                    rows.append(dict(**base, n_stacked=len(W), snr_db=sn, null_mean_db=float(nul.mean()) if len(nul) else np.nan,
                                     null_sd_db=float(nul.std()) if len(nul) else np.nan, excess_db=exc, p_null=p,
                                     cc_template=cc, cc_null_mean=float(ncc.mean()) if len(ncc) else np.nan, p_cc=pcc,
                                     peak_offset_s=pk, n_null=len(nul), verdict=vd))
                    STK[f"{s}|{int(f)}|{Tt:g}|{name}"] = st.astype(np.float32)
                    if tp is not None:
                        STK[f"tpl|{s}|{int(f)}"] = tp.astype(np.float32)
    R = pd.DataFrame(rows); R.to_csv(a.out, index=False)
    np.savez_compressed(a.out.replace(".csv", "_stacks.npz"), **STK)
    pd.set_option("display.width", 300)
    for Tt in THR:
        q = R[(R.set == "quiet") & (R["T"] == Tt)]
        c = R[(R.set == "ets_control") & (R["T"] == Tt)]
        tested = q[q.verdict.isin(["PASS", "FAIL", "inconclusive"])]
        npass = int((tested.verdict == "PASS").sum())
        det_pass = tested[tested.verdict == "PASS"].n_det.sum(); det_all = tested.n_det.sum()
        print(f"\n[T = {Tt}] families testable {len(tested)} (of {q.fam.nunique()} with quiet-day detections): "
              f"PASS {npass}, FAIL {int((tested.verdict=='FAIL').sum())}, inconclusive {int((tested.verdict=='inconclusive').sum())}; "
              f"quiet-day detections in PASS families {int(det_pass)} / {int(det_all)} = {det_pass/max(det_all,1):.1%}")
        print(f"          ETS positive control: registers {int((c.verdict=='control: registers').sum())} / "
              f"{int(c.verdict.str.startswith('control').sum())}")
    print(f"[write] {a.out}")


if __name__ == "__main__":
    main()
