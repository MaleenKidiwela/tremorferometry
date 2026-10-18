"""S15 raw figures (before verdicts): Q1 Δt / rank / minutes; Q2 phase histograms, stress deciles, detection vs tide."""
import numpy as np, pandas as pd, json
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import s15_common as C

BLUE, ORANGE, AQUA, GRAY, INK, MUTED = '#2a78d6', '#eb6834', '#1baf7a', '#9a9993', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': '#e6e5e0', 'grid.linewidth': 0.6})
T, FG = C.OUT_T, C.OUT_F


def q1():
    V = pd.read_csv(f'{T}/q1_events.csv')
    u = pd.read_csv(f'{T}/q1b_u_rank.csv').u_rank_reanchored.values
    ctrlR = np.load(f'{T}/q1_ctrl_dt_reanchored.npy')
    ctrlA = np.load(f'{T}/q1_ctrl_dt_concurrent.npy')
    ok = V.n_lfe_10m.values >= 1
    okR = ok & np.isfinite(u)
    obs = V.dt_d.values[okR]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    b = np.arange(-2, 16.5, 0.5)
    a = ax[0, 0]
    a.hist(np.clip(ctrlR.ravel(), b[0], b[-1] - 0.01), b, density=True, color=GRAY, alpha=0.6, label='random-draw control (re-anchored pseudo-VLFEs)')
    a.hist(np.clip(obs, b[0], b[-1] - 0.01), b, density=True, histtype='step', lw=2, color=BLUE, label=f'VLFEs with LFEs (n {okR.sum()})')
    a.axvline(np.median(obs), color=BLUE, lw=1.5, ls='--'); a.axvline(np.median(ctrlR), color=MUTED, lw=1.5, ls='--')
    a.set_xlabel('Δt = VLFE time − first LFE of concurrent families (days)'); a.set_ylabel('density')
    a.set_title('Day scale: time since local LFE onset', loc='left'); a.legend(frameon=False, fontsize=8)
    a = ax[0, 1]
    a.hist(u[okR], np.linspace(0, 1, 11), color=BLUE, alpha=0.85, rwidth=0.95)
    a.axhline(okR.sum() / 10, color=MUTED, ls='--', lw=1.2); a.text(0.02, okR.sum() / 10 * 1.03, 'uniform = VLFE is a random LFE', color=MUTED, fontsize=8)
    a.set_xlabel('rank among pseudo-VLFE times, same families/episode (0 first, 1 last)')
    a.set_ylabel('VLFEs'); a.set_title(f'Rank in the local sequence: mean {np.nanmean(u[okR]):.2f}', loc='left')
    a = ax[1, 0]
    mR = np.median(ctrlR, 0); mA = np.median(ctrlA, 0)
    a.hist(mR, 40, color=GRAY, alpha=0.7, label='1,000 control medians (re-anchored)')
    a.hist(mA, 40, histtype='step', color=MUTED, lw=1.2, label='1,000 control medians (fixed-weight, pre-registered)')
    a.axvline(np.median(obs), color=BLUE, lw=2, label=f'observed median {np.median(obs):.2f} d')
    a.axvline(np.percentile(mR, 99), color=INK, lw=0.8, ls=':', label='99th percentile (re-anchored)')
    a.set_xlabel('median Δt (days)'); a.set_ylabel('draws'); a.legend(frameon=False, fontsize=8); a.set_title('Pooled test (Merlin rule)', loc='left')
    a = ax[1, 1]
    E = pd.read_csv(f'{T}/q1b_day_per_episode.csv'); E = E[(E.label == 'R_reanchored_excl_self') & E.episode.isin(C.BIG9)]
    E = E.set_index('episode').loc[C.BIG9].reset_index()
    col = [BLUE if d == 'NW-ward' else ORANGE for d in E.direction]
    a.bar(range(len(E)), E.diff_d, color=col, width=0.7)
    for i, r in E.iterrows():
        a.text(i, r.diff_d + 0.04, f'n{r.n}', ha='center', fontsize=7, color=MUTED)
    a.set_xticks(range(len(E))); a.set_xticklabels(E.episode, rotation=40, ha='right'); a.axhline(0, color=INK, lw=0.8)
    a.set_ylabel('median Δt − control median (days)')
    a.set_title('Per big episode (blue NW-ward front, orange SE-ward)', loc='left')
    fig.tight_layout(); fig.savefig(f'{FG}/q1_day_scale.png', dpi=150); plt.close(fig)

    # minutes
    z = np.load(f'{T}/q1_minutes_hist.npz'); f = z['fine']; c = (f[:-1] + f[1:]) / 2 / 60
    fig, a = plt.subplots(figsize=(8, 3.6))
    hc = z['h_ctrl'] * z['h_obs'][np.abs(c) > 2].sum() / z['h_ctrl'][np.abs(c) > 2].sum()
    a.step(c, hc, where='mid', color=GRAY, lw=1.5, label='around random LFEs of the same families (scaled to the 2–10 min level)')
    a.step(c, z['h_obs'], where='mid', color=BLUE, lw=2, label='around the 431 VLFEs with LFEs')
    for x in (-2, 0, 2):
        a.axvline(x, color=MUTED, lw=0.6, ls=':')
    a.set_xlabel('LFE time − VLFE time (min), 20-s bins'); a.set_ylabel('LFEs per bin'); a.legend(frameon=False, fontsize=8)
    a.set_title('Minutes scale: pool LFEs around VLFE times', loc='left')
    fig.tight_layout(); fig.savefig(f'{FG}/q1_minutes.png', dpi=150); plt.close(fig)

    # variants per direction
    D = pd.concat([pd.read_csv(f'{T}/q1_day_per_direction.csv'), pd.read_csv(f'{T}/q1b_day_per_direction.csv')])
    labs = [('R_reanchored_excl_self', 'concurrent LFEs\n(re-anchored control)'), ('A_concurrent_431', 'concurrent LFEs\n(fixed control)'),
            ('B_location_431', 'VLFE location ±20 km\n(431 with LFEs)'), ('C_location_109', 'VLFE location ±20 km\n(no LFEs, 83 with pool)')]
    fig, a = plt.subplots(figsize=(8, 3.6))
    for j, (l, nm) in enumerate(labs):
        for k, (d, cc) in enumerate((('NW-ward', BLUE), ('SE-ward', ORANGE))):
            r = D[(D.label == l) & (D.direction == d)]
            if len(r):
                a.bar(j * 3 + k, r.diff_d.iloc[0], color=cc, width=0.85, label=d if j == 0 else None)
    a.set_xticks([j * 3 + 0.5 for j in range(len(labs))]); a.set_xticklabels([n for _, n in labs], fontsize=8)
    a.axhline(0, color=INK, lw=0.8); a.set_ylabel('median Δt − control (days)'); a.legend(frameon=False, title='front direction', fontsize=8)
    a.set_title('Anchor variants by front direction (location variant carries the ~11.5 km NW location bias)', loc='left', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FG}/q1_variants_direction.png', dpi=150); plt.close(fig)


def q2():
    Z = np.load(f'{T}/q2_phases.npz'); ok = Z['okc']
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    b = np.linspace(-180, 180, 13); cb = (b[:-1] + b[1:]) / 2
    for a, key, name in ((ax[0, 0], '', 'M2'), (ax[0, 1], '_k1', 'K1')):
        v = np.degrees(np.angle(np.exp(1j * Z['vlfe' + key][ok])))
        l = np.degrees(np.angle(np.exp(1j * Z['lfe' + key][ok].ravel())))
        r = np.degrees(np.angle(np.exp(1j * Z['rand' + key][ok].ravel())))
        n = ok.sum()
        hv = np.histogram(v, b)[0] / n; hl = np.histogram(l[np.isfinite(l)], b)[0] / np.isfinite(l).sum(); hr = np.histogram(r, b)[0] / len(r)
        a.bar(cb, hv / hr, width=27, color=BLUE, alpha=0.85, label=f'VLFEs (n {n})')
        a.plot(cb, hl / hr, 'o-', color=ORANGE, lw=2, ms=5, label='matched LFE draws (same families, ±24 h)')
        a.axhline(1, color=MUTED, lw=1, ls='--')
        a.set_xlabel(f'{name} phase of up-dip shear (deg; 0 = peak, −90 = peak rate)')
        a.set_ylabel('rate / random times'); a.set_xticks([-180, -90, 0, 90, 180]); a.legend(frameon=False, fontsize=8)
        a.set_title(f'{name} phase histogram (relative to random times ±24 h)', loc='left')
    a = ax[1, 0]
    S = pd.read_csv(f'{T}/q2_stress_deciles.csv')
    a.plot(S.decile, S.vlfe_over_random, 'o-', color=BLUE, lw=2, ms=6, label='VLFEs')
    a.plot(S.decile, S.lfe_over_random, 's-', color=ORANGE, lw=2, ms=5, label='matched LFE draws')
    a.axhline(1, color=MUTED, lw=1, ls='--'); a.set_yscale('log'); a.minorticks_off(); a.set_yticks([0.25, 0.5, 1, 2]); a.set_yticklabels(['0.25', '0.5', '1', '2'])
    a.set_xlabel('decile of total tidal up-dip shear (1 = most negative)'); a.set_ylabel('rate / random times')
    a.legend(frameon=False, fontsize=8); a.set_title('Stress level (total tide, not one constituent)', loc='left')
    a = ax[1, 1]
    P = pd.read_csv(f'{T}/q2_per_episode.csv'); P = P[P.episode != 'all'].reset_index(drop=True)
    x = np.arange(len(P))
    a.vlines(x, P.a_lfe_p01, P.a_lfe_p99, color=ORANGE, lw=6, alpha=0.35, label='matched LFEs, 1st–99th pct')
    a.plot(x, P.a_lfe_median, 's', color=ORANGE, ms=6, label='matched LFEs, median')
    a.plot(x, P.a_vlfe, 'o', color=BLUE, ms=8, label='VLFEs')
    for i, r in P.iterrows():
        a.text(i, r.a_vlfe + 0.06, f'n{r.n}', ha='center', fontsize=7, color=MUTED)
    a.set_xticks(x); a.set_xticklabels(P.episode, rotation=40, ha='right'); a.set_ylabel('M2 modulation amplitude a = 2R')
    a.legend(frameon=False, fontsize=8, loc='upper left'); a.set_title('Per big episode', loc='left')
    fig.tight_layout(); fig.savefig(f'{FG}/q2_tides.png', dpi=150); plt.close(fig)

    J = pd.read_csv(f'{T}/q2_injections_phase.csv')
    fig, a = plt.subplots(figsize=(7, 3.4))
    b = np.linspace(-np.pi, np.pi, 13); cb = np.degrees((b[:-1] + b[1:]) / 2)
    for mw, cc in ((3.6, AQUA), (3.8, BLUE)):
        s = J[J.Mw == mw]; ph = np.angle(np.exp(1j * s.th_m2.values))
        r = [s.rec.values[(ph >= b[i]) & (ph < b[i + 1])].mean() for i in range(12)]
        a.plot(cb, r, 'o-', color=cc, lw=2, ms=5, label=f'injected Mw {mw} (n {len(s):,})')
    a.set_ylim(0.4, 0.9); a.set_xticks([-180, -90, 0, 90, 180])
    a.set_xlabel('M2 phase of up-dip shear at the injection time (deg)'); a.set_ylabel('recovered fraction')
    a.legend(frameon=False, fontsize=8); a.set_title('Detection check: injection recovery vs M2 phase (C0/C1, 17 runs)', loc='left')
    fig.tight_layout(); fig.savefig(f'{FG}/q2_detection_vs_tide.png', dpi=150); plt.close(fig)

    # detection vs local day
    Dd = pd.read_csv(f'{T}/detect_recovery_vs_local_day.csv')
    fig, a = plt.subplots(figsize=(6, 3.2))
    for k, (s, cc) in enumerate((('Mw3.6', AQUA), ('Mw3.8', BLUE))):
        d = Dd[(Dd.set == s) & Dd.bin.isin(['<0', '0-2', '2-5', '>=5'])]
        x = np.arange(4) + (k - 0.5) * 0.25
        a.errorbar(x, d.recovery, yerr=[d.recovery - d.ci_lo, d.ci_hi - d.recovery], fmt='o', color=cc, ms=6, capsize=0, label=f'injected {s}')
    a.set_xticks(range(4)); a.set_xticklabels(['before', '0–2', '2–5', '≥5']); a.set_ylim(0.4, 0.9)
    a.set_xlabel('days since local LFE onset (families ≤ 20 km of the injection source)'); a.set_ylabel('recovered fraction')
    a.legend(frameon=False, fontsize=8); a.set_title('Detection check: recovery vs day in the local sequence', loc='left', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{FG}/q1_detection_vs_local_day.png', dpi=150); plt.close(fig)


if __name__ == '__main__':
    q1(); q2()
