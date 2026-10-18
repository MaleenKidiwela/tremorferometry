"""S15 raw figures for Q3 (partition) and Q5 (migration rate)."""
import json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import s15_common as C
from s15_figs import BLUE, ORANGE, AQUA, GRAY, INK, MUTED

YEL, VIO = '#eda100', '#4a3aa7'
T, FG = C.OUT_T, C.OUT_F


def q3():
    P = pd.read_csv(f'{T}/q3_partition.csv')
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    a = ax[0]; x = np.arange(len(P))
    for col, c, m, lab in (('lfe_self_m0_svi', GRAY, 's', 'LFE self-moment (family mag_med)'),
                           ('cont_m0', BLUE, 'o', 'continuous 20–50 s (S13)'),
                           ('disc_m0_raw', ORANGE, '^', 'discrete VLFEs, raw'),
                           ('disc_m0_div_c34', ORANGE, 'v', 'discrete ÷ compl_3.4'),
                           ('gnss_tierB_M0', VIO, 'D', 'GNSS Tier B (2014+)')):
        y = P[col].where(P[col] > 0)
        a.plot(x, y, m, color=c, ms=7 if m != 'v' else 6, mfc=c if m != 'v' else 'none', label=lab, ls='none')
    a.set_yscale('log'); a.set_xticks(x); a.set_xticklabels(P.episode, rotation=50, ha='right', fontsize=8)
    a.set_ylabel('moment (N m)'); a.set_ylim(5e13, 3e20); a.legend(frameon=False, fontsize=7, loc='upper center', ncol=3); a.set_title('Per-episode moments', loc='left')
    for i, s in enumerate(P.subset):
        if s == 'matched':
            a.axvspan(i - 0.4, i + 0.4, color='#e6e5e0', zorder=0)
    a = ax[1]
    M = P[P.n_vlfe > 0]
    cols = [BLUE if s == 'matched' else (ORANGE if s == 'flagged' else GRAY) for s in M.subset]
    a.scatter(M.reg_slip_mm, M.vlfe_per_100_zone_lfe_corr, c=cols, s=55, zorder=3)
    for r in M.itertuples():
        a.annotate(r.episode, (r.reg_slip_mm, r.vlfe_per_100_zone_lfe_corr), fontsize=7, color=MUTED, xytext=(4, 2), textcoords='offset points')
    a.set_yscale('log'); a.set_xlabel('regional slip increment (vbICA, mm)'); a.set_ylabel('VLFEs per 100 zone LFEs ÷ compl_3.4')
    a.set_title('Size test (blue = matched-completeness subset)', loc='left')
    fig.tight_layout(); fig.savefig(f'{FG}/q3_partition.png', dpi=150); plt.close(fig)


def q5():
    Z = np.load(f'{T}/q5_draws.npz'); E = pd.read_csv(f'{T}/q5_events.csv'); ok = Z['ok']
    edges = Z['edges']
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    a = ax[0]; b = np.linspace(0, 1.0, 26)
    lp = np.clip(np.maximum(Z['lfe_pace'], 0), 0, 0.999); lp = lp[np.isfinite(lp)]
    dp = np.clip(np.maximum(Z['pace_draw'][ok].ravel(), 0), 0, 0.999); dp = dp[np.isfinite(dp)]
    vp = np.clip(np.maximum(E.pace_local.values[ok], 0), 0, 0.999)
    a.hist(lp, b, density=True, color=GRAY, alpha=0.45, label='all pool LFEs (tercile edges)')
    a.hist(dp, b, density=True, histtype='step', color=MUTED, lw=1.5, label='matched pseudo-VLFEs (±24 h)')
    a.hist(vp, b, density=True, histtype='step', color=BLUE, lw=2, label=f'VLFEs (n {ok.sum()})')
    for e_ in edges:
        a.axvline(e_, color=INK, lw=0.7, ls=':')
    a.set_xlabel('local passage pace (d/km; 0 = fast or simultaneous, clipped at 1)'); a.set_ylabel('density')
    a.legend(frameon=False, fontsize=7.5); a.set_title('Local passage pace', loc='left')
    S = json.load(open(f'{T}/q5_summary.json'))
    a = ax[1]
    for k, (key, c, lab) in enumerate((('main_matched_24h', BLUE, 'local pace vs ±24 h draws (verdict)'),
                                       ('variant_whole_episode', GRAY, 'local pace vs whole-episode draws'),
                                       ('daily_speed_tercile_matched_24h', ORANGE, 'daily front speed vs ±24 h draws'))):
        r = S[key]['per_lfe_ratio']
        if key == 'daily_speed_tercile_matched_24h':
            r = r                                      # index 0 = fast, 2 = slow (reordered in s15_q5)
        a.plot([0, 1, 2], r, 'o-', color=c, lw=2, ms=6, label=lab)
    a.axhline(1, color=MUTED, ls='--', lw=1); a.set_xticks([0, 1, 2]); a.set_xticklabels(['fast', 'medium', 'slow'])
    a.set_ylabel('VLFE share ÷ control share'); a.legend(frameon=False, fontsize=7.5); a.set_title('VLFE per LFE by speed tercile', loc='left')
    a = ax[2]
    TI = pd.read_csv(f'{T}/q5_tide_by_tercile.csv'); TI = TI[TI.set == 'declustered_1h']
    a.plot([0, 1, 2], TI.a_vlfe, 'o-', color=BLUE, lw=2, ms=7, label='VLFEs (declustered)')
    a.plot([0, 1, 2], TI.a_lfe_median, 's-', color=ORANGE, lw=2, ms=6, label='matched LFEs, median')
    for i, r in enumerate(TI.itertuples()):
        a.text(i, r.a_vlfe + 0.04, f'n{r.n}', ha='center', fontsize=7, color=MUTED)
    a.set_xticks([0, 1, 2]); a.set_xticklabels(['fast', 'medium', 'slow']); a.set_ylabel('M2 modulation a')
    a.legend(frameon=False, fontsize=7.5); a.set_title('Tidal modulation by local pace', loc='left')
    fig.tight_layout(); fig.savefig(f'{FG}/q5_migration.png', dpi=150); plt.close(fig)


if __name__ == '__main__':
    q3(); q5()
