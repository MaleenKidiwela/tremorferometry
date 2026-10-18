#!/usr/bin/env python
"""E6 step 1b -- POSITIVE LFE-vs-noise STACK TEST at the EXCLUDED stations (Merlin E5r-audit F2, prescription 1b).
The network was re-run with CN.VGZ and CN.SNB EXCLUDED from the sum (e4_template_test.py --only-fams 103,32,75,128,74, 13 stations,
tag e6_1b_noVGZSNB).  At its >= 8-MAD QUIET-day detection times (doy 228-242) the VGZ and SNB 3-C 20-Hz data are stacked with the
family moveout (predicted S = origin + ttS(Crosson) + the E5r template's S-pick refinement at that station; E5r cut_windows: rms-
normalised windows, linear stack, NO alignment refinement).  Stack S-SNR (dB) = energy S-2..S+2 over the 4-s pre-window at the
start of the +-20-s cut (S-20..S-16, before P at every distance here), against a same-N random-time stack over the same quiet days
(>= 20 draws).  POSITIVE CONTROL: the same families' ETS-day OWN-MATCHED detections (same excluded network), same N where possible.
Verdicts (pre-registered): NOISE if S-SNR <= null mean + 1 dB at N >= 300; REAL inter-ETS LFEs if >= null + 3 dB with the envelope
peak within +-2 s of the predicted S; otherwise INCONCLUSIVE.  Per family, per station, and pooled (all five families aligned on
their own predicted S).  Also S-SNR of nested sub-stacks (N = 30, 100, 300, all) for figures/e6_stack_test.png.
Writes results/e6_stack_test_1b.csv, results/e6_stack_test_1b_stacks.npz."""
from __future__ import annotations
import argparse, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e5r_templates as T  # noqa: E402
import d3_vel1d as v1d  # noqa: E402
RES = T.RES; FS = T.FS
FAMS = [103, 32, 75, 128, 74]
STAS = ["VGZ", "SNB"]
QUIET = list(range(228, 243))
NSUB = (30, 100, 300)


def snr_db(st, s_idx, half_s=2.0):
    n = int(round(2 * half_s * FS)); a = s_idx - int(round(half_s * FS))
    sig = st[:, a:a + n]; noi = st[:, 0:n]                       # pre-window = the first 4 s of the cut (S-20..S-16)
    return float(10 * np.log10((sig ** 2).mean() / (noi ** 2).mean()))


def peak_offset_s(st, s_idx, search_s=6.0):
    env = T.envelope(st); w = int(round(search_s * FS)); a, b = max(s_idx - w, 0), min(s_idx + w, len(env))
    return (a + int(np.argmax(env[a:b])) - s_idx) / FS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="e6_1b_noVGZSNB"); ap.add_argument("--thr", type=float, default=8.0)
    ap.add_argument("--draws", type=int, default=24); ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    S = pd.read_csv(f"{RES}/e5r_stations.csv").set_index("sta")
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv").set_index("bkey")
    qc = pd.read_csv(f"{RES}/e5r_template_qc_2012.csv")
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet"); gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size(); ETS = set(per_day[per_day >= 100].index)
    det = pd.read_parquet(f"{RES}/e5r_detections_{a.tag}.parquet"); det = det[(det.L == 8.0) & (det.net_mad >= a.thr)]
    print(f"[det] {a.tag}: {len(det)} detections >= {a.thr} MAD; quiet {int(det.doy.isin(QUIET).sum())}, ETS own-matched "
          f"{int((det.doy.isin(ETS) & det.match_own).sum())}", flush=True)
    n = int(round(2 * T.HALF_W * FS)); idx_nom = int(round(T.HALF_W * FS))
    rows, stacks = [], {}
    pooled = {s: {"quiet": [], "ets": [], "null": [[] for _ in range(a.draws)]} for s in STAS}
    for s in STAS:
        D = T.load_days(S.loc[s, "netsta"], 2012)
        qd = [d for d in QUIET if d in D]
        t_lo = min(D[d][0] for d in qd); t_hi = max(D[d][0] + D[d][1].shape[1] / FS for d in qd)
        for f in FAMS:
            fr = F.loc[f]
            dx = T.epi_km(fr.lat, fr.lon, float(S.loc[s, "lat"]), float(S.loc[s, "lon"]))
            ttS = float(v1d.tt(dx, fr.dep, "S"))
            q = qc[(qc.sta == s) & (qc.bkey == f) & qc.usable & qc.template]
            pick = float(q.s_pick_s.iloc[0]) if len(q) else 0.0      # the E5r template's S refinement at this station
            dq = det[(det.fam == float(f)) & det.doy.isin(QUIET)]
            de = det[(det.fam == float(f)) & det.doy.isin(ETS) & det.match_own]
            sets = {"quiet": np.sort(dq.t.to_numpy()), "ets_control": np.sort(de.t.to_numpy())}
            for name, tt in sets.items():
                if len(tt) < 10:
                    rows.append(dict(fam=f, sta=s, set=name, N=len(tt), note="too few")); continue
                W, keep = T.cut_windows(D, tt + ttS + pick - T.HALF_W, n)
                if W is None or len(W) < 10:
                    rows.append(dict(fam=f, sta=s, set=name, N=0, note="no data")); continue
                st = W.mean(0); snr = snr_db(st, idx_nom); pk = peak_offset_s(st, idx_nom)
                sub = {f"snr_N{k}": (snr_db(W[rng.choice(len(W), k, replace=False)].mean(0), idx_nom) if k < len(W) else np.nan) for k in NSUB}
                # null: same-N random QUIET-day times, same station, same offset (>= 20 draws)
                nul, nul_pk = [], []
                for _ in range(a.draws):
                    Wn, _ = T.cut_windows(D, rng.uniform(t_lo, t_hi, len(W)) + ttS + pick - T.HALF_W, n)
                    if Wn is not None and len(Wn) >= 10:
                        stn = Wn.mean(0); nul.append(snr_db(stn, idx_nom)); nul_pk.append(peak_offset_s(stn, idx_nom))
                nul = np.array(nul)
                excess = snr - nul.mean()
                if name == "quiet":
                    verdict = ("NOISE" if (excess <= 1.0 and len(W) >= 300) else
                               "REAL" if (excess >= 3.0 and abs(pk) <= 2.0) else
                               "inconclusive (N<300)" if len(W) < 300 else "inconclusive")
                else:
                    verdict = "control: registers" if (excess >= 3.0 and abs(pk) <= 2.0) else "control: weak"
                rows.append(dict(fam=f, sta=s, set=name, N=len(W), ttS_s=ttS, pick_s=pick, has_template=bool(len(q)), snr_db=snr,
                                 null_mean_db=float(nul.mean()), null_sd_db=float(nul.std()), null_max_db=float(nul.max()), excess_db=excess,
                                 peak_offset_s=pk, null_peak_abs_median_s=float(np.median(np.abs(nul_pk))), verdict=verdict, **sub))
                stacks[f"{s}|{f}|{name}"] = st.astype(np.float32)
                key = "quiet" if name == "quiet" else "ets"
                pooled[s][key].append(W)
                if name == "quiet":
                    for i in range(min(a.draws, len(nul))):
                        Wn, _ = T.cut_windows(D, rng.uniform(t_lo, t_hi, len(W)) + ttS + pick - T.HALF_W, n)
                        if Wn is not None:
                            pooled[s]["null"][i].append(Wn)
                print(f"  [{s} fam {f} {name}] N {len(W)} S-SNR {snr:+.2f} dB, null {nul.mean():+.2f} +- {nul.std():.2f} (max {nul.max():+.2f}), "
                      f"excess {excess:+.2f} dB, envelope peak {pk:+.2f} s -> {verdict}", flush=True)
        # pooled over the five families (each aligned on its own predicted S)
        for key in ("quiet", "ets"):
            if not pooled[s][key]:
                continue
            W = np.concatenate(pooled[s][key]); st = W.mean(0); snr = snr_db(st, idx_nom); pk = peak_offset_s(st, idx_nom)
            sub = {f"snr_N{k}": (snr_db(W[rng.choice(len(W), k, replace=False)].mean(0), idx_nom) if k < len(W) else np.nan) for k in NSUB}
            if key == "quiet":
                nul = np.array([snr_db(np.concatenate(ws).mean(0), idx_nom) for ws in pooled[s]["null"] if ws])
                excess = snr - nul.mean()
                verdict = ("NOISE" if (excess <= 1.0 and len(W) >= 300) else "REAL" if (excess >= 3.0 and abs(pk) <= 2.0) else "inconclusive")
                rows.append(dict(fam=-1, sta=s, set="pooled_quiet", N=len(W), snr_db=snr, null_mean_db=float(nul.mean()), null_sd_db=float(nul.std()),
                                 null_max_db=float(nul.max()), excess_db=excess, peak_offset_s=pk, verdict=verdict, **sub))
                print(f"  [{s} POOLED quiet] N {len(W)} S-SNR {snr:+.2f} dB vs null {nul.mean():+.2f} +- {nul.std():.2f}: excess {excess:+.2f} dB, peak {pk:+.2f} s -> {verdict}", flush=True)
            else:
                rows.append(dict(fam=-1, sta=s, set="pooled_ets_control", N=len(W), snr_db=snr, peak_offset_s=pk, **sub))
                print(f"  [{s} POOLED ETS own-matched control] N {len(W)} S-SNR {snr:+.2f} dB, peak {pk:+.2f} s", flush=True)
            stacks[f"{s}|pooled|{key}"] = st.astype(np.float32)
    R = pd.DataFrame(rows); R.to_csv(f"{RES}/e6_stack_test_1b.csv", index=False)
    np.savez_compressed(f"{RES}/e6_stack_test_1b_stacks.npz", **stacks)
    pd.set_option("display.width", 300)
    print(R.round(2).to_string(index=False))
    print(f"[write] results/e6_stack_test_1b.csv, e6_stack_test_1b_stacks.npz")


if __name__ == "__main__":
    main()
