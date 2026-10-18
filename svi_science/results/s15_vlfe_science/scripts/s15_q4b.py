"""S15 Amendment 5 parts without injections: per-episode VLFE vs concurrent-LFE along-strike centroids; VLFEs in the
regional slip-rate history (column); largest Mw and completeness-corrected count >= Mw 3.6 (b-value substitute)."""
import json
import numpy as np, pandas as pd
from scipy.stats import chi2
import s15_common as C
from s15_q3 import compl_table


def pois_ci(k):
    lo = 0.5 * chi2.ppf(0.025, 2 * k) if k > 0 else 0.0
    return lo, 0.5 * chi2.ppf(0.975, 2 * k + 2)


def main():
    Q1 = pd.read_csv(f'{C.OUT_T}/q1_events.csv')
    E = C.episodes().set_index('episode')
    comp = compl_table()
    with_lfe = Q1.n_lfe_10m.values >= 1
    # 1. centroids
    rows = []
    for e, g in Q1[with_lfe].groupby('episode'):
        if len(g) < 5:
            continue
        rows.append(dict(episode=e, direction=g.direction.iloc[0], n=len(g), vlfe_along_median=float(g.vlfe_along.median()),
                         lfe_anchor_along_median=float(g.anchor_along.median()),
                         diff_km=float(g.vlfe_along.median() - g.anchor_along.median()),
                         diff_per_event_median=float((g.vlfe_along - g.anchor_along).median()),
                         vlfe_lat_median=float(g.lat_best.median()),
                         lfe_anchor_lat_median=float(C.LAT0 + g.anchor_y.median() / C.KM)))
    CE = pd.DataFrame(rows); CE.to_csv(f'{C.OUT_T}/q4_centroids.csv', index=False)
    allr = dict(n=int(with_lfe.sum()), diff_per_event_median=float((Q1.vlfe_along - Q1.anchor_along)[with_lfe].median()),
                by_direction={d: float((Q1.vlfe_along - Q1.anchor_along)[with_lfe & (Q1.direction == d).values].median()) for d in ('NW-ward', 'SE-ward')},
                vlfe_along_p10_p90=[float(np.percentile(Q1.vlfe_along, q)) for q in (10, 90)],
                lfe_anchor_along_p10_p90=[float(np.nanpercentile(Q1.anchor_along[with_lfe], q)) for q in (10, 90)])
    # 2. slip-rate column
    SL = pd.read_parquet(f'{C.ROOT}/lfe_extension/results/slip_compare/slip_svi_family_daily.parquet', columns=['fid', 'date', 'slip_rate_mm_per_day'])
    reg = SL.groupby('date').slip_rate_mm_per_day.mean(); reg.index = pd.to_datetime(reg.index)
    H = pd.read_parquet(f'{C.OUT_T}/q2c_hourly.parquet')          # zone LFE counts per hour (episodes with >= 5 VLFEs)
    sr = []
    for e in H.episode.unique():
        h = H[H.episode == e]
        lday = pd.to_datetime(h.t_mid.values, unit='s').normalize()
        w = h.n_zone.values
        lr = reg.reindex(lday).values
        ok = np.isfinite(lr) & (w > 0)
        lr_w = np.repeat(lr[ok], w[ok].astype(int))                   # one value per zone LFE
        g = Q1[Q1.episode == e]
        vday = pd.to_datetime(g.epoch.values, unit='s').normalize()
        vr = reg.reindex(vday).values
        if len(lr_w) < 10 or not np.isfinite(vr).any():
            continue
        pct = np.array([100 * np.mean(lr_w < x) for x in vr[np.isfinite(vr)]])
        win = reg[(reg.index >= pd.Timestamp(E.loc[e, 'start'])) & (reg.index <= pd.Timestamp(E.loc[e, 'end']) + pd.Timedelta(days=1))]
        peak = win.idxmax() if len(win) else pd.NaT
        sr.append(dict(episode=e, n_vlfe=int(np.isfinite(vr).sum()), n_zone_lfe=len(lr_w), vlfe_rate_median=float(np.nanmedian(vr)),
                       lfe_rate_median=float(np.median(lr_w)), vlfe_pct_in_lfe_median=float(np.median(pct)), vlfe_pct_in_lfe_mean=float(np.mean(pct)),
                       peak_rate_day=str(peak.date()) if pd.notna(peak) else None,
                       vlfe_minus_peak_days_median=float(np.median((vday[np.isfinite(vr)] - peak).days)) if pd.notna(peak) else np.nan,
                       lfe_minus_peak_days_median=float(np.median(np.repeat((lday[ok] - peak).days, w[ok].astype(int)))) if pd.notna(peak) else np.nan))
    SR = pd.DataFrame(sr); SR.to_csv(f'{C.OUT_T}/q5orig_sliprate.csv', index=False)
    # 3. b-value substitute
    bv = []
    for e, g in Q1.groupby('episode'):
        if len(g) < 5:
            continue
        k = int((g.Mw_scan >= 3.6).sum()); c36 = comp.get(e, {}).get('compl_3.6', np.nan)
        lo, hi = pois_ci(k)
        bv.append(dict(episode=e, n=len(g), Mw_max=float(g.Mw_scan.max()), Mw_median=float(g.Mw_scan.median()), n_ge36=k, compl_3_6=c36,
                       n_ge36_corr=k / c36 if c36 else np.nan, n_ge36_corr_lo=lo / c36 if c36 else np.nan, n_ge36_corr_hi=hi / c36 if c36 else np.nan,
                       n_ge38=int((g.Mw_scan >= 3.8).sum())))
    BV = pd.DataFrame(bv); BV.to_csv(f'{C.OUT_T}/q_bsubstitute.csv', index=False)
    json.dump(allr, open(f'{C.OUT_T}/q4_centroids_all.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(CE.round(2).to_string()); print(json.dumps(allr, indent=1)); print(SR.round(2).to_string()); print(BV.round(2).to_string())


if __name__ == '__main__':
    main()
