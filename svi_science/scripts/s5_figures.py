"""S5 figures: modulation vs tidal stress, sensitivity vs days since front/onset, per-era comparison."""
import json, numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

R = Path('/home/jovyan/tremorferometry/svi_science/results/s5')
FIG = Path('/home/jovyan/tremorferometry/svi_science/figures')
C_B, C_O, C_O1, C_O2, C_N = '#2a78d6', '#eb6834', '#eb6834', '#1baf7a', '#8a8984'
DT = ['<1d', '1-2d', '2-3d', '3-4d', '4-6d', '6-10d', '>=10d']
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
                     'grid.color': '#e4e3df', 'grid.linewidth': 0.6})

T2 = pd.read_csv(R / 'step2_modulation_table.csv')
H2 = np.load(R / 'step2_hists.npz')
has3 = (R / 'step3_table_data_075_once.csv').exists()
if has3:
    T3 = pd.concat([pd.read_csv(f) for f in sorted(R.glob('step3_table_*.csv'))], ignore_index=True); H3 = np.load(R / 'step3_hists.npz')


def ratio(H, k):
    o, s = H[f'obs_{k}'], H[f'sur_{k}']
    return o / np.where(s > 0, s, np.nan), o


# ---------------- 1. modulation vs tidal stress
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
zc = (np.linspace(-2.5, 2.5, 11)[:-1] + 0.25)
pc = np.degrees(np.linspace(-np.pi, np.pi, 13)[:-1] + np.pi / 12)
sets = [('Bostock 2003–2013', H2, C_B, 'o')] + ([('ours 2014–2026 (P≥0.75, once)', H3, C_O, 's')] if has3 else [])
for lab, H, col, mk in sets:
    r, o = ratio(H, 'tauz_all'); r, o = r[0], o[0]
    e = r / np.sqrt(np.maximum(o, 1))
    ax[0].errorbar(zc, r, yerr=e, color=col, marker=mk, ms=5, lw=2, capsize=0, label=lab)
    r, o = ratio(H, 'phM2_all'); r, o = r[0], o[0]
    ax[1].errorbar(pc, r, yerr=r / np.sqrt(np.maximum(o, 1)), color=col, marker=mk, ms=5, lw=2, label=lab)
for a in ax:
    a.axhline(1, color='#52514e', lw=0.8)
ax[0].set_xlabel('tidal shear τ at the LFE time (family std; + promotes thrust)')
ax[0].set_ylabel('LFEs observed / expected')
ax[0].set_title('LFE rate vs tidal shear stress (236° slip direction)', loc='left')
ax[1].set_xlabel('M2 phase of τ (deg; 0 = peak shear, −90 = peak shear rate)')
ax[1].set_title('LFE rate vs M2 tidal phase', loc='left'); ax[1].set_xticks([-180, -90, 0, 90, 180])
ax[0].legend(frameon=False, fontsize=8)
fig.text(0.01, 0.01, 'expected = whole-day-shift surrogates; bars = Poisson (understate bursts)', fontsize=7, color='#52514e')
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIG / 's5_modulation_vs_stress.png', dpi=140); plt.close(fig)


# ---------------- 2. sensitivity vs days since front / onset
def series(T, setname, label):
    d = T[(T.label == label)] if setname is None else T[(T.set == setname) & (T.label == label)]
    d = d.set_index('name').reindex(DT)
    return d


fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
x = np.arange(7)
curves = [('Bostock 2003–13, days since family onset', series(T2, None, 'dtbin'), C_B, 'o', -0.15)]
if has3:
    curves += [('ours 2014–26, days since PNSN front', series(T3, 'data_075_once', 'dtbin'), C_O, 's', 0.0),
               ('ours 2014–26, days since family onset (majors)', series(T3, 'data_075_once', 'dtown_major'), '#4a3aa7', '^', 0.15)]
for lab, d, col, mk, off in curves:
    ok = d.n >= 30
    ax[0].errorbar(x[ok] + off, d.a_m2_pct[ok], yerr=200 * d.Csurr_sd_m2[ok], color=col, marker=mk, ms=5, lw=2, label=lab)
    ax[1].plot(x[ok] + off, d.theta0_m2_deg[ok], color=col, marker=mk, ms=5, lw=2, label=lab)
ax[0].set_ylabel('M2 modulation a (%)')
ax[1].set_ylabel('preferred M2 phase θ0 (deg)')
ax[1].axhline(0, color='#52514e', lw=0.8); ax[1].axhline(-90, color='#52514e', lw=0.8, ls=':')
ax[1].text(6.3, -85, 'peak stress rate', fontsize=7, ha='right'); ax[1].text(6.3, 5, 'peak stress', fontsize=7, ha='right')
for a in ax:
    a.set_xticks(x); a.set_xticklabels(DT); a.set_xlabel('time since front / onset at the family')
ax[0].set_title('Tidal sensitivity within the episode', loc='left'); ax[1].set_title('Phase of the response', loc='left')
ax[0].legend(frameon=False, fontsize=7)
fig.tight_layout(); fig.savefig(FIG / 's5_sensitivity_vs_days.png', dpi=140); plt.close(fig)

# ---------------- 3. per-era comparison
if has3:
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    b = T2[T2.label == 'all'].iloc[0]
    h = T3[(T3.set == 'data_075_once') & (T3.label == 'half')].set_index('name')
    rows = [('Bostock\n2003–13', b, C_B), ('ours\n2014–19', h.loc['2014-2019'], C_O1), ('ours\n2020–26', h.loc['2020-2026'], C_O2)]
    sc = T3[(T3.label == 'all')].set_index('set')
    for nm in [x for x in ('c1_050_once', 'c1_all_raw') if x in sc.index]:
        rows.append((f'scrambled\n{"P≥0.5" if "050" in nm else "all P"}', sc.loc[nm], C_N))
    for i, (nm, r, col) in enumerate(rows):
        ax[0].bar(i, r.a_m2_pct, color=col, width=0.6)
        ax[0].errorbar(i, r.a_m2_pct, yerr=200 * r.Csurr_sd_m2, color='#0b0b0b', lw=1)
        ax[1].errorbar(i, r.theta0_m2_deg, yerr=np.degrees(r.Csurr_sd_m2 / max(abs(r.a_m2_pct) / 200, 1e-3)), color=col, marker='o', ms=7, lw=2)
    for a in ax[:2]:
        a.set_xticks(range(len(rows))); a.set_xticklabels([r[0] for r in rows], fontsize=8)
    ax[0].set_ylabel('M2 modulation amplitude a (%)'); ax[0].set_title('Size of the tidal response', loc='left')
    ax[1].set_ylabel('θ0 (deg)'); ax[1].axhline(0, color='#52514e', lw=0.8); ax[1].set_ylim(-180, 180)
    ax[1].set_title('Phase (0 = peak shear)', loc='left')
    # normalised rise curves: response per bin / mean response over days 2-6
    def norm(d):
        v = d.Cex_m2.values.astype(float); n = d.n.values.astype(float)
        v[~(n >= 30)] = np.nan
        return v / np.nanmean(v[2:5])
    fh = lambda p: T3[(T3.set == 'data_075_once') & (T3.label == 'dtfront_half') & T3.name.str.endswith(p)].assign(
        name=lambda q: q.name.str.split('|').str[0]).set_index('name').reindex(DT)
    for lab, d, col, mk, ls in [('Bostock 2003–13, since family onset', series(T2, None, 'dtbin'), C_B, 'o', '-'),
                                ('ours 2014–26, since family onset', series(T3, 'data_075_once', 'dtown_major'), '#4a3aa7', '^', '-'),
                                ('ours 2014–19, since PNSN front', fh('2014-19'), C_O1, 's', '--'),
                                ('ours 2020–26, since PNSN front', fh('2020-26'), C_O2, 'D', '--')]:
        ax[2].plot(x, norm(d), color=col, marker=mk, lw=2, ms=5, ls=ls, label=lab)
    ax[2].set_ylabel('M2 response / mean of days 2–6')
    ax[2].set_xticks(x); ax[2].set_xticklabels(DT); ax[2].axhline(1, color='#52514e', lw=0.8)
    ax[2].set_title('Shape of the within-episode rise', loc='left')
    ax[2].legend(frameon=False, fontsize=7)
    fig.text(0.01, 0.01, 'error bars: tidal-sampling SD from surrogates only; episodes scatter ~1.4–3× more (see RESULTS_S5.md)', fontsize=7, color='#52514e')
    fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIG / 's5_era_comparison.png', dpi=140); plt.close(fig)
print('figures done')
