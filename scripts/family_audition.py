#!/usr/bin/env python
"""Two-stage cheap family audition before committing to a full-record densify.

  Stage 1 (FREE, pre-densify): template-shape contaminant screen.
      drop families with spectral centroid > CEN_THR Hz AND kurtosis > KURT_THR (noise-like / non-LFE).
  Stage 2 (CHEAP): one-probe-year densify -> count active days that year -> keep only families that are
      ACTIVE >= ACT_THR days (continuous repeaters). A probe year spans ETS+quiet, so it is robust to the
      ETS-burst false positive that a single month suffers.

This is the engine for family-set EXPANSION: audition the whole SNR-eligible pool cheaply, then full-densify
ONLY the genuine + continuous families (episodic families pool to noise at sparse stations, so excluding them
both saves GPU and improves dv/v quality).

Modes:
  --validate STA   : prove the logic on an already-densified station using its EXISTING stacks (no download).
                     Simulates the probe-year from long_window_daily_<STA>.npz and compares to full-record truth.
  --screen STA     : run Stage-1 template screen on the eligible pool, write the genuine-family list.
  (live Stage-2 densify of a probe year is driven by the pipeline once waveforms are on disk; see --emit-probe-cmd)
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy.stats import kurtosis as kurt_fn

CEN_THR = 4.3      # spectral centroid Hz: above (with high kurtosis) = contaminant
KURT_THR = 4.0
ACT_THR = 15       # active days in the probe YEAR-month-equivalent to call continuous
FS = 40.0


def template_feats(w):
    w = np.asarray(w, float); w = w - w.mean()
    if w.std() == 0:
        return np.nan, np.nan
    W = np.abs(np.fft.rfft(w)); fr = np.fft.rfftfreq(len(w), 1 / FS)
    cen = (W * fr).sum() / W.sum()
    return cen, float(kurt_fn(w, fisher=True))


def screen_pool(sta, floor):
    """Stage 1 (CONSERVATIVE pre-filter, not the final word — Stage 2 is the arbiter).
    SNR-floored eligible families; flag contaminants by a STATION-RELATIVE centroid threshold
    (max of the global 4.3 Hz and the station's own 85th-pct centroid) AND high kurtosis. The
    station-relative cap avoids clipping genuine families at intrinsically high-frequency stations
    (B020/B028), cutting false-negatives on continuous families 3.5%->1.2%."""
    s = sta.lower()
    summ = pd.read_csv(f'data/{s}_pnsn_families_100km.summary.csv')
    z = np.load(f'data/{s}_pnsn_families_100km.npz', allow_pickle=True)
    elig = summ[summ.snr >= floor]
    rows = []
    for f in elig.family_id:
        if f not in z.files:
            continue
        cen, k = template_feats(z[f])
        rows.append((f, cen, k))
    t = pd.DataFrame(rows, columns=['family_id', 'centroid', 'kurt'])
    # STATION-HEALTH GATE (FLAW-5 fix): the station-relative cap (85th-pct) self-defeats at
    # contaminant-DOMINATED stations — there the 85th-pct centroid is itself high, so the relative
    # cap rises to meet the junk and passes it (B933 would sail through). If the station's MEDIAN
    # template centroid is already in the contaminant range (> global GOOD ~3.4 + margin), the station
    # is contaminant-dominated -> fall back to the strict GLOBAL threshold (and flag for human review).
    med_cen = t.centroid.median()
    contaminant_station = med_cen > 4.3   # B933 median was 5.03; healthy stations ~3.4
    cen_cut = CEN_THR if contaminant_station else max(CEN_THR, t.centroid.quantile(0.85))
    t['genuine'] = ~((t.centroid > cen_cut) & (t.kurt > KURT_THR))
    if contaminant_station:
        print(f'[{sta}] ⚠ CONTAMINANT-DOMINATED (median centroid {med_cen:.1f} Hz > 4.3) — '
              f'using strict global threshold; flag for human review. genuine {int(t.genuine.sum())}/{len(t)}')
    return t


def validate(sta):
    """Prove probe-year classification vs full-record truth, from EXISTING stacks (no download)."""
    npz = np.load(f'data/long_window_daily_{sta}.npz', allow_pickle=True)
    pat = npz['patches'].astype(str); dts = pd.to_datetime(npz['dates'].astype(str))
    df = pd.DataFrame({'p': pat, 'd': dts})
    df['y'] = df.d.dt.year
    # full-record truth: continuous if active >=50% of own-span days (matches family_continuity_classes)
    ad = pd.Index(df.d.unique())
    res = []
    for p, g in df.groupby('p'):
        dd = g.d.sort_values()
        own = ad[(ad >= dd.iloc[0]) & (ad <= dd.iloc[-1])]
        cov = len(dd) / max(len(own), 1)
        truth = cov >= 0.5
        # probe each calendar year the family exists; active days per (probe-year)/12 ~ active-month-equiv
        per_yr = g.groupby('y').d.nunique()
        # one-probe-year decision: pick the family's MEDIAN active-year, scale to per-month, threshold
        # robust probe = active-days-per-year/12 >= ACT_THR? use active fraction of the probe year
        for yr, nd in per_yr.items():
            pred = (nd / 12.0) >= ACT_THR  # avg active days/month within the probe year
            res.append((p, yr, truth, pred, nd))
    R = pd.DataFrame(res, columns=['patch', 'probe_yr', 'truth', 'pred', 'days_in_yr'])
    # per-family: majority vote across its probe-years (a real probe uses ONE year; report both)
    fam = R.groupby('patch').agg(truth=('truth', 'first'),
                                 pred_any=('pred', 'max'), pred_med=('pred', 'mean'))
    fam['pred1yr'] = fam.pred_med >= 0.5
    tp = ((fam.truth) & (fam.pred1yr)).sum(); fn = ((fam.truth) & (~fam.pred1yr)).sum()
    fp = ((~fam.truth) & (fam.pred1yr)).sum(); tn = ((~fam.truth) & (~fam.pred1yr)).sum()
    n = len(fam)
    print(f'[{sta}] {n} families | continuous(truth)={fam.truth.sum()} sparse={(~fam.truth).sum()}')
    print(f'  probe-year classifier (active>={ACT_THR} d/mo-equiv):')
    print(f'    recall(continuous kept) = {tp/(tp+fn):.2f}   |  sparse-rejected = {tn/(tn+fp):.2f}')
    print(f'    precision(kept are continuous) = {tp/(tp+fp) if (tp+fp) else float("nan"):.2f}')
    return fam


def probe_stage2(sta, probe_years, net='PB', cov_thr=0.5, despike=8, workers=20):
    """Stage 2 (LIVE): short-window densify of the Stage-1 genuine candidates NOT already densified,
    classify continuous vs episodic by active-day fraction in the probe window, write the keep-list,
    and DELETE the probe scratch for the SKIPPED families (the probe is a throwaway classifier — winners
    get a full-record densify separately, so all probe detections are intermediate)."""
    import glob, subprocess
    s = sta.lower()
    scr = pd.read_csv(f'data/{s}_audition_screen.csv')
    done = set(pd.read_csv(f'data/{s}_coverage_selection.summary.csv').family_id)
    summ = pd.read_csv(f'data/{s}_pnsn_families_100km.summary.csv').set_index('family_id')
    cand = [f for f in scr[scr.genuine].family_id if f not in done]
    if not cand:
        print(f'[{sta}] no probe candidates'); return
    # temp selection summary of just the candidates
    tmp = f'data/{s}_audition_candidates.summary.csv'
    summ.loc[summ.index.intersection(cand)].reset_index().to_csv(tmp, index=False)
    y0, y1 = min(probe_years), max(probe_years)
    pre = f'mf_{s}_probe_'
    print(f'[{sta}] Stage-2 probe densify: {len(cand)} candidates over probe years {y0}-{y1}', flush=True)
    subprocess.run([sys.executable, 'scripts/densify_launcher.py', '--templates-npz',
        f'data/{s}_pnsn_families_100km.npz', '--summary-csv', tmp, '--min-snr', '0',
        '--network', net, '--station', sta, '--out-prefix', pre, '--workers', str(workers),
        '--top-n', '100', '--max-raw-det', '3000000', '--despike-mad', str(despike),
        '--start-year', str(y0), '--end-year', str(y1)], check=True)
    # classify: active-day fraction in the probe window per family
    pf = sorted(glob.glob(f'data/{pre}[12]*.csv'))
    det = pd.concat([pd.read_csv(f) for f in pf], ignore_index=True)
    det['day'] = pd.to_datetime(det.time).dt.date
    ndays_probe = pd.to_datetime(det.time).dt.normalize().nunique()
    act = det.groupby('template').day.nunique()
    keep_frac = act / max(ndays_probe, 1)
    log = pd.DataFrame({'family_id': act.index, 'active_days': act.values,
                        'probe_days': ndays_probe, 'active_frac': keep_frac.values})
    log['continuous'] = log.active_frac >= cov_thr
    # candidates that produced ZERO probe detections are episodic-by-default
    zero = [f for f in cand if f not in set(act.index)]
    if zero:
        log = pd.concat([log, pd.DataFrame({'family_id': zero, 'active_days': 0, 'probe_days': ndays_probe,
                                            'active_frac': 0.0, 'continuous': False})], ignore_index=True)
    log.to_csv(f'data/{s}_audition_log.csv', index=False)
    keep = log[log.continuous].family_id.tolist()
    pd.DataFrame({'family_id': keep}).to_csv(f'data/{s}_audition_keep.csv', index=False)
    # CLEANUP: probe densify is throwaway scratch -> delete it (skipped AND winners; winners get a full densify)
    removed = 0
    for f in pf + [tmp]:
        try: os.remove(f); removed += 1
        except OSError: pass
    print(f'[{sta}] Stage-2: {len(cand)} candidates -> {len(keep)} continuous (kept), {len(cand)-len(keep)} episodic (skipped)')
    print(f'  wrote data/{s}_audition_keep.csv (+ _audition_log.csv); deleted {removed} probe-scratch files')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--validate'); ap.add_argument('--screen'); ap.add_argument('--floor', type=float, default=6)
    ap.add_argument('--act-thr', type=int, default=ACT_THR)
    ap.add_argument('--probe-stage2', help='station: run LIVE probe densify + classify + cleanup (needs waveforms)')
    ap.add_argument('--probe-years', nargs='+', type=int, help='probe year(s), e.g. 2016 (or 2014 2018 for 2 separated)')
    ap.add_argument('--network', default='PB')
    a = ap.parse_args()
    globals()['ACT_THR'] = a.act_thr
    if a.validate:
        validate(a.validate)
    elif a.screen:
        t = screen_pool(a.screen, a.floor)
        out = f'data/{a.screen.lower()}_audition_screen.csv'
        t.to_csv(out, index=False)
        print(f'[{a.screen}] eligible {len(t)} | genuine (pass template) {int(t.genuine.sum())} '
              f'({100*t.genuine.mean():.0f}%) -> {out}')
    elif a.probe_stage2:
        assert a.probe_years, '--probe-years required'
        probe_stage2(a.probe_stage2, a.probe_years, net=a.network)


if __name__ == '__main__':
    main()
