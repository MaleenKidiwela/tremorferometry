#!/usr/bin/env python
"""S13 scale-up summary (PLAN.md Amendment 1 tests A1-A6) + figures. Table choice scalef_ > refit_ > scale_;
Amendment 4 v2 channel-guard check on tables not run with v2 (units v2 would change -> NaN, listed in
tables/scale_guard_excluded.csv). SVI only: WA withdrawn 2026-10-14 pending the origin-time fix.
Writes tables/scale_{svi,wa}_episodes.csv, tables/scale_{svi,wa}_families.csv, tables/scale_between.csv,
tables/scale_verdict.json, figures/s13_scale_*.png."""
import os, sys, glob, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C

RES = f'{C.ROOT}/svi_science/results'
RNG = np.random.default_rng(C.SEED)
GREY = '#9a9a9a'


GUARD_LOG = []


def pick(reg, b):
    """Table preference (PLAN.md Amendments 3-4): scalef_ (gap fill) > refit_ (guarded re-run) > scale_ (first pass)."""
    for pre in ('scalef', 'refit', 'scale'):
        f = f'{C.TB}/{pre}_{reg}_{b}_units.csv'
        if os.path.isfile(f):
            return f, pre
    return None, None


def guard(U, f, pre):
    """Amendment 4 v2 check: units from a table run with the v2 guard (guard_rule == 'A4v2') are used as they are. Units
    from any other table (first-pass scale_, or v1 cap-rule) are used only if v2 would drop nothing in them and no cap
    was applied (then the numbers are identical); otherwise set to NaN and listed in tables/scale_guard_excluded.csv."""
    if 'guard_rule' in U.columns and (U.guard_rule == C.GUARD_RULE).all():
        return U
    S = pd.read_csv(f.replace('_units.csv', '_stacomp.csv'))
    bad = []
    for u, g in S.groupby('unit', sort=False):
        gg = {s: x.set_index('comp').reindex(list('NEZ')) for s, x in g.groupby('sta', sort=False)}
        _, drop, _ = C.guard2({s: x.sigma_ctrl_nm.values for s, x in gg.items()}, {s: x.n_used.values for s, x in gg.items()})
        if any(d.any() for d in drop.values()):
            bad.append(u)
    if 'guard_capped' in U.columns:
        bad += list(U.unit[U.guard_capped.astype(bool)])
    bad = sorted(set(bad))
    for u in bad:
        GUARD_LOG.append(dict(table=os.path.basename(f), unit=u))
    cols = [c for c in U.columns if c.split('_')[0] in ('thrust', 'strikeslip', 'steepdip') or c in ('vr_free', 'mt_free')]
    U = U.copy(); U.loc[U.unit.isin(bad), cols] = np.nan
    return U


def perm_spearman(x, y, n=5000):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[m], np.asarray(y)[m]
    if len(x) < 5:
        return dict(n=int(len(x)), rho=np.nan, p=np.nan)
    r = spearmanr(x, y).statistic
    rr = np.array([spearmanr(x, RNG.permutation(y)).statistic for _ in range(n)])
    return dict(n=int(len(x)), rho=float(r), p=float((1 + (np.abs(rr) >= abs(r)).sum()) / (n + 1)))


def region(reg):
    E = pd.read_csv(f'{C.TB}/{reg}_episodes.csv')
    U = []
    for e in E[E.built].episode:
        f, pre = pick(reg, e)
        if f:
            U.append(guard(pd.read_csv(f), f, pre).assign(episode=e, table=pre))
    U = pd.concat(U, ignore_index=True)
    P = U[U.unit == 'pooled'].set_index('episode')
    pc = ['n', 'n_sta', 'thrust_m0', 'thrust_z', 'thrust_ctrl_sd', 'thrust_vr', 'thrust_m0_shift_p1', 'thrust_m0_shift_m1',
          'n_nb10', 'guard_dropped', 'max_share', 'max_share_ch', 'loso_m0_min', 'loso_m0_max', 'loso_min_without',
          'loso_max_without']
    ep = E.set_index('episode').join(P[[c for c in pc if c in P.columns]], rsuffix='_u')
    for h in ('pooled_odd', 'pooled_even'):
        hh = U[U.unit == h].set_index('episode')
        ep[f'z_{h}'] = hh.thrust_z; ep[f'm0_{h}'] = hh.thrust_m0
    ep['table'] = U[U.unit == 'pooled'].set_index('episode').table
    ep['m0_total'] = ep.thrust_m0 * ep.n
    ep['detected'] = ep.thrust_z >= 3
    sz = pd.read_csv(f'{RES}/s11_depth/round6/{"svi" if reg == "svi" else "wa_nw"}/tables/r6_sizes.csv').set_index('episode')
    ep = ep.join(sz[[c for c in ('G', 'G_WA', 'M0', 'T', 'Du', 'L') if c in sz.columns]], rsuffix='_s')
    dd = pd.read_csv(f'{RES}/s11_depth/round7/{"svi" if reg == "svi" else "wa_nw"}/tables/r7_hint_dD.csv').set_index('episode')
    ep = ep.join(dd[['dD', 'sig']].rename(columns={'dD': 'dvv_dD', 'sig': 'dvv_sig'}))
    ep.to_csv(f'{C.TB}/scale_{reg}_episodes.csv')
    # families
    Fm = U[U.unit.str.startswith('fam') & U.thrust_m0.notna() & (U.thrust_ctrl_sd > 0)].copy()
    Fm['w'] = 1 / Fm.thrust_ctrl_sd ** 2
    fam = Fm.groupby('unit').apply(lambda g: pd.Series(dict(
        n_ep=len(g), n_lfe=g.n.sum(), m0_wmean=(g.thrust_m0 * g.w).sum() / g.w.sum(), m0_sd=1 / np.sqrt(g.w.sum()),
        chi2_dof=((g.thrust_m0 - (g.thrust_m0 * g.w).sum() / g.w.sum()) ** 2 * g.w).sum() / max(len(g) - 1, 1),
        n_z3=(g.thrust_z >= 3).sum(), lat=g.src_lat.mean(), lon=g.src_lon.mean(), dep=g.src_dep.mean())))
    fam['z_wmean'] = fam.m0_wmean / fam.m0_sd
    fam.to_csv(f'{C.TB}/scale_{reg}_families.csv')
    Z = U[U.unit.str.startswith('zone')]
    Z.to_csv(f'{C.TB}/scale_{reg}_zones.csv', index=False)
    # tests
    v = {}
    built = ep[ep.thrust_m0.notna()]
    v['A1_n_episodes'] = int(len(built)); v['A1_n_detected_z3'] = int(built.detected.sum())
    v['A1_detected'] = list(built.index[built.detected])
    w = 1 / built.thrust_ctrl_sd ** 2
    mu = (built.thrust_m0 * w).sum() / w.sum()
    chi2 = ((built.thrust_m0 - mu) ** 2 * w).sum()
    from scipy.stats import chi2 as CH
    v['A2_m0_per_lfe_wmean'] = float(mu); v['A2_m0_per_lfe_sd'] = float(1 / np.sqrt(w.sum()))
    v['A2_chi2_dof'] = float(chi2 / max(len(built) - 1, 1)); v['A2_p_const'] = float(CH.sf(chi2, max(len(built) - 1, 1)))
    gcol = 'G' if reg == 'svi' else 'G_WA'
    for nm, col in (('gnss', gcol), ('tremor', 'T'), ('duration', 'Du'), ('lfecount', 'L')):
        if col in built.columns:
            v[f'A3_total_vs_{nm}'] = perm_spearman(built.m0_total.values, built[col].values)
            v[f'A3_perlfe_vs_{nm}'] = perm_spearman(built.thrust_m0.values, built[col].values)
    v['A4_families'] = int(len(fam)); v['A4_fam_chi2_dof_median'] = float(fam.chi2_dof.median())
    fw = 1 / fam.m0_sd ** 2; fmu = (fam.m0_wmean * fw).sum() / fw.sum()
    v['A4_between_family_chi2_dof'] = float(((fam.m0_wmean - fmu) ** 2 * fw).sum() / max(len(fam) - 1, 1))
    v['A6_total_vs_dvv'] = perm_spearman(built.m0_total.values, built.dvv_dD.values)
    v['A6_perlfe_vs_dvv'] = perm_spearman(built.thrust_m0.values, built.dvv_dD.values)
    return ep, fam, v


def figs(reg, ep):
    b = ep[ep.thrust_m0.notna()].copy()
    b['mid'] = pd.to_datetime(b.start) + (pd.to_datetime(b.end) - pd.to_datetime(b.start)) / 2
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True)
    ax[0].plot(b.mid, b.thrust_m0, color=GREY, lw=0.8)
    sc = ax[0].scatter(b.mid, b.thrust_m0, c=b.thrust_z.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8,
                       s=12 + 60 * np.sqrt(b.n / b.n.max()), zorder=3)
    ax[0].axhline(0, color=GREY, lw=0.6)
    ax[0].set_ylabel('VLF moment per LFE (N m)')
    ax[1].plot(b.mid, b.m0_total, color=GREY, lw=0.8)
    ax[1].scatter(b.mid, b.m0_total, c=b.thrust_z.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8, s=25, zorder=3)
    ax[1].set_yscale('symlog', linthresh=1e14); ax[1].set_ylabel('time-locked VLF moment (N m)\n= per LFE × n LFEs')
    fig.colorbar(sc, ax=ax, label='z vs scrambled control')
    ax[0].set_title(f'S13 {reg.upper()}: 20–50 s moment time-locked to LFEs, per major episode (size ∝ √n LFEs)', fontsize=9)
    fig.savefig(f'{C.FIG}/s13_scale_{reg}_episodes.png', dpi=130); plt.close(fig)
    gcol = 'G' if reg == 'svi' else 'T'
    if gcol in b.columns and b[gcol].notna().sum() >= 5:
        fig, ax = plt.subplots(figsize=(5.5, 4.2))
        sc = ax.scatter(b[gcol], b.m0_total, c=b.thrust_z.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8, s=30)
        for e, r in b.iterrows():
            if np.isfinite(r[gcol]):
                ax.annotate(e, (r[gcol], r.m0_total), fontsize=6)
        ax.set_xlabel('GNSS offset D (mm, S1)' if gcol == 'G' else 'PNSN tremor count (per 1,000)')
        ax.set_ylabel('time-locked VLF moment (N m)'); fig.colorbar(sc, ax=ax, label='z')
        fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_scale_{reg}_vs_size.png', dpi=130); plt.close(fig)


def main():
    V = {}
    REGS = ('svi',)   # WA withdrawn 2026-10-14 pending the WA origin-time fix (Merlin F1: wa_shift ~16 s early)
    V['wa'] = 'WITHDRAWN pending WA origin-time fix (Merlin F1, 2026-10-14); WA tables not summarised'
    for reg in REGS:
        if not glob.glob(f'{C.TB}/scale_{reg}_*_units.csv'):
            continue
        ep, fam, v = region(reg)
        figs(reg, ep)
        V[reg] = v
    B = []
    for reg in REGS:
        f, pre = pick(reg, 'between')
        if f:
            B.append(guard(pd.read_csv(f), f, pre).assign(table=pre))
    if B:
        B = pd.concat(B, ignore_index=True)
        B.to_csv(f'{C.TB}/scale_between.csv', index=False)
        V['A5_between'] = {r.unit: dict(n=int(r.n), m0=float(r.thrust_m0), z=float(r.thrust_z), vr=float(r.thrust_vr), table=r.table)
                           for r in B.itertuples()}
    pd.DataFrame(GUARD_LOG, columns=['table', 'unit']).to_csv(f'{C.TB}/scale_guard_excluded.csv', index=False)
    V['guard_excluded_units'] = len(GUARD_LOG)   # Amendment 4: unguarded units breaking the rule, no guarded re-run yet
    json.dump(V, open(f'{C.TB}/scale_verdict.json', 'w'), indent=1, default=float)
    print(json.dumps(V, indent=1, default=float))


if __name__ == '__main__':
    main()
