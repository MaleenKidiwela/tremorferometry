"""S5 setup, step 4: combine body + ocean-load stress, Coulomb, constituent split, validation figures.

Output: results/s5_setup/tidal_stress_svi_hourly.parquet
  fid, time (UTC, hourly 2014-01-01 00:00 .. 2026-09-30 23:00),
  tau_Pa, sigma_n_Pa, coulomb_0p1, coulomb_0p4               (total = body + ocean)
  body_tau_Pa, body_sigma_n_Pa, ocean_tau_Pa, ocean_sigma_n_Pa
  tau_updip_Pa, tau_strike_Pa (+ body_/ocean_ parts)          (alternative slip directions)
Conventions: sigma_n tension (unclamping) positive; tau positive = promotes hanging-wall slip
in the 236 deg convergence direction projected on the Slab2 plane; Coulomb = tau + mu' sigma_n.
Also: tidal_constituents_per_family.csv, variant_sensitivity.csv, figures/s5_*.png
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from s5_tides_common import OUT, ROOT, MAJOR, design, MODEL_CONST

FIG = ROOT / 'svi_science/figures'
FIG.mkdir(exist_ok=True)


def main():
    body = pd.read_parquet(OUT / 'body_stress_hourly.parquet')
    ocean = pd.read_parquet(OUT / 'ocean_stress_hourly.parquet')
    df = body[['fid', 'time', 'body_tau_Pa', 'body_sigma_n_Pa', 'body_tau_updip_Pa', 'body_tau_strike_Pa']].merge(
        ocean, on=['fid', 'time'], how='inner', validate='one_to_one')
    assert len(df) == len(body) == len(ocean)
    for c in ('tau', 'sigma_n', 'tau_updip', 'tau_strike'):
        df[f'{c}_Pa'] = (df[f'body_{c}_Pa'] + df[f'ocean_{c}_Pa']).astype(np.float32)
    df['coulomb_0p1'] = (df.tau_Pa + 0.1 * df.sigma_n_Pa).astype(np.float32)
    df['coulomb_0p4'] = (df.tau_Pa + 0.4 * df.sigma_n_Pa).astype(np.float32)
    cols = ['fid', 'time', 'tau_Pa', 'sigma_n_Pa', 'coulomb_0p1', 'coulomb_0p4', 'body_tau_Pa', 'body_sigma_n_Pa',
            'ocean_tau_Pa', 'ocean_sigma_n_Pa', 'tau_updip_Pa', 'tau_strike_Pa', 'body_tau_updip_Pa',
            'body_tau_strike_Pa', 'ocean_tau_updip_Pa', 'ocean_tau_strike_Pa']
    df = df[cols].sort_values(['fid', 'time']).reset_index(drop=True)
    df.to_parquet(OUT / 'tidal_stress_svi_hourly.parquet', index=False)
    print('rows', len(df), 'families', df.fid.nunique(), df.time.min(), df.time.max(), flush=True)

    # ---- constituent split per family (harmonic LSQ over the full record)
    fids = np.sort(df.fid.unique()); t = df.time[df.fid == fids[0]].values
    consts = MAJOR + ['m4', 'j1', 'oo1', 'mu2', 'nu2', 'l2', 't2', '2q1', 'rho1', 'm1', 'mf', 'mm']
    X = np.c_[design(t, consts), np.ones(len(t))]
    pinv = np.linalg.pinv(X)
    ser = ['tau_Pa', 'sigma_n_Pa', 'coulomb_0p4', 'body_tau_Pa', 'body_sigma_n_Pa', 'ocean_tau_Pa', 'ocean_sigma_n_Pa']
    rows = []
    for f in fids:
        sub = df[df.fid == f]
        for s in ser:
            c = pinv @ sub[s].values.astype(float)
            n = len(consts)
            for i, k in enumerate(consts[:10]):
                hc = c[i] + 1j * c[n + i]
                rows.append(dict(fid=f, series=s, const=k, amp_Pa=abs(hc), phase_deg=(-np.degrees(np.angle(hc))) % 360))
            rows.append(dict(fid=f, series=s, const='std', amp_Pa=sub[s].std(), phase_deg=np.nan))
            rows.append(dict(fid=f, series=s, const='p2p_99', amp_Pa=np.subtract(*np.percentile(sub[s], [99.5, 0.5])), phase_deg=np.nan))
    C = pd.DataFrame(rows); C.to_csv(OUT / 'tidal_constituents_per_family.csv', index=False)
    summ = C.pivot_table(index=['series'], columns='const', values='amp_Pa', aggfunc='median')[
        ['m2', 's2', 'n2', 'k1', 'o1', 'p1', 'std', 'p2p_99']].round(0)
    summ.to_csv(OUT / 'tidal_constituents_median_summary.csv')
    print(summ.to_string(), flush=True)

    # ---- variant sensitivity of the ocean part (M2, K1, O1 vector differences)
    T = pd.read_csv(OUT / 'ocean_constituent_stress.csv')
    T['hc'] = T.re + 1j * T.im
    def var(model, fill, depth, zones):
        q = T[(T.model == model) & (T.fill == fill) & (T.depth == depth) & (T.zone.isin(zones))]
        return q.groupby(['fid', 'comp', 'const']).hc.sum()
    ref = var('EOT20', 'gauge', 'cat', ['inner', 'outer'])
    vs = {'EOT20 nearest-fill (no gauges)': var('EOT20', 'nearest', 'cat', ['inner', 'outer']),
          'GOT5.5 gauge-fill': var('GOT5.5', 'gauge', 'cat', ['inner', 'outer']),
          'Slab2 depth instead of catalogue depth': var('EOT20', 'gauge', 'slab2', ['inner', 'outer']),
          'inner 330 km only (no 330-700 km annulus)': var('EOT20', 'gauge', 'cat', ['inner'])}
    rows = []
    for name, v in vs.items():
        d = (v - ref).abs() / ref.abs()
        r = v.abs() / ref.abs()
        for comp in ('tau', 'sigma_n'):
            for k in ('m2', 'k1', 'o1'):
                idx = (slice(None), comp, k)
                rows.append(dict(variant=name, comp=comp, const=k, median_amp_ratio=r.loc[idx].median(),
                                 median_rel_vector_diff=d.loc[idx].median(), p90_rel_vector_diff=d.loc[idx].quantile(0.9)))
    S = pd.DataFrame(rows); S.to_csv(OUT / 'variant_sensitivity.csv', index=False)
    print(S.round(3).to_string(), flush=True)

    # ---- figures: one family, one month + spectrum
    fam = pd.read_csv(ROOT / 'svi_science/results/s0/families.csv')
    f0 = int(fam.sort_values('n075_data_dedup_prod', ascending=False).fid.iloc[0])
    sub = df[df.fid == f0].set_index('time')
    m = sub.loc['2016-02-01':'2016-02-29']
    fig, ax = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for a, s, lab in [(ax[0], 'tau', 'shear stress τ (236° slip dir.)'), (ax[1], 'sigma_n', 'normal stress σn (tension +)')]:
        a.plot(m.index, m[f'body_{s}_Pa'] / 1e3, lw=0.8, color='#7a7a7a', label='solid-earth (body) tide')
        a.plot(m.index, m[f'ocean_{s}_Pa'] / 1e3, lw=0.8, color='#2a7ab0', label='ocean load')
        a.plot(m.index, m[f'{s}_Pa'] / 1e3, lw=1.2, color='#c0392b', label='total')
        a.set_ylabel('kPa'); a.set_title(lab, fontsize=10, loc='left'); a.grid(alpha=0.3)
    ax[0].legend(ncol=3, fontsize=8, loc='upper right')
    ax[2].plot(m.index, m.coulomb_0p1 / 1e3, lw=1, color='#555', label="Coulomb μ'=0.1")
    ax[2].plot(m.index, m.coulomb_0p4 / 1e3, lw=1, color='#d35400', label="Coulomb μ'=0.4")
    ax[2].set_ylabel('kPa'); ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)
    ax[2].set_xlabel('day of 2016-02 (UTC)')
    import matplotlib.dates as mdates
    ax[2].xaxis.set_major_formatter(mdates.DateFormatter('%d'))
    r = fam[fam.fid == f0].iloc[0]
    fig.suptitle(f'Tidal stress on the plate interface, SVI family {f0} ({r.lat:.2f}N {abs(r.lon):.2f}W, {r.dep:.0f} km)', fontsize=11)
    fig.tight_layout(); fig.savefig(FIG / 's5_tidal_stress_month.png', dpi=130); plt.close(fig)

    fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    win = np.hanning(len(sub)); freqs = np.fft.rfftfreq(len(sub), d=1 / 24.0)  # cycles/day
    for ax, q, lab in [(axs[0], 'tau', 'τ (236° slip direction)'), (axs[1], 'sigma_n', 'σn (tension +)')]:
        for pre, col, lw, nm in [('', '#c0392b', 1.4, 'total'), ('ocean_', '#2a7ab0', 0.8, 'ocean load'), ('body_', '#555555', 0.8, 'body tide')]:
            x = sub[f'{pre}{q}_Pa'].values
            A = np.abs(np.fft.rfft((x - x.mean()) * win)) * 2 / win.sum()
            ax.semilogy(freqs, A, lw=lw, color=col, label=nm)
        for nm, f in [('O1', 0.9295), ('K1', 1.0027), ('N2', 1.8960), ('M2', 1.9323), ('S2', 2.0)]:
            ax.axvline(f, color='k', lw=0.4, ls=':'); ax.text(f + 0.005, 2.5e3, nm, fontsize=8)
        ax.set_xlim(0.8, 2.1); ax.set_ylim(1, 6e3); ax.set_ylabel('amplitude (Pa)')
        ax.set_title(lab, fontsize=10, loc='left'); ax.grid(alpha=0.3); ax.legend(fontsize=8, loc='center')
    axs[1].set_xlabel('frequency (cycles/day)')
    fig.suptitle(f'Amplitude spectrum, family {f0}, hourly 2014-01..2026-09 (P1 sits next to K1)', fontsize=10)
    fig.tight_layout(); fig.savefig(FIG / 's5_tidal_stress_spectrum.png', dpi=130); plt.close(fig)

    # ---- map: M2 / K1 amplitude of tau and sigma_n, ocean/body ratio
    geo = pd.read_csv(OUT / 'family_fault_geometry.csv')
    P = C[C.const.isin(['m2', 'k1'])].pivot_table(index='fid', columns=['series', 'const'], values='amp_Pa')
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.2))
    for a, (s, k) in zip(ax, [('tau_Pa', 'm2'), ('tau_Pa', 'k1'), ('sigma_n_Pa', 'm2'), ('sigma_n_Pa', 'k1')]):
        v = P[(s, k)].reindex(geo.fid).values
        sc = a.scatter(geo.lon, geo.lat, c=v, s=18, cmap='viridis')
        plt.colorbar(sc, ax=a, label='Pa'); a.set_title(f'{s.replace("_Pa","")} {k.upper()} amplitude (total)', fontsize=9)
        a.set_aspect(1 / np.cos(np.radians(48.6)))
    fig.tight_layout(); fig.savefig(FIG / 's5_tidal_amplitude_map.png', dpi=120); plt.close(fig)
    print('fig family', f0)


if __name__ == '__main__':
    main()
