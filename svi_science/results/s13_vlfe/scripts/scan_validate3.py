#!/usr/bin/env python
"""VLFE scan pilot 3 (PLAN_SCAN.md Amendment S8): SVI 2012 with the ratio veto R1 fixed on 2011.
Run with SCAN_YEAR=2012. Inputs: data/scan_2012/runs/{real, ctrlA_s*, inj2, inj3}/cand2.parquet, inject lists, masks2_2012.
Writes tables/scan3_filter_validation.csv, scan3_bw_known_list.csv, scan3_veto_validation.csv, scan3_veto_points.csv,
scan3_threshold_grid{,_noveto,_r1fixed6,_r1nf}.csv, scan3_injection{..}.csv, scan3_injection_summary{..}.csv,
scan3_detections_2012{..}.csv, scan3_verdict.json, figures/scan3_*.png."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
from scipy.stats import poisson
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
import scan_masks2 as M2
from scan_validate import between_mask, ets_vec, live_years
from scan_validate2 import crit, ok_all, match, CRIT, COLS_KEEP

assert S.YEAR == 2012, 'run with SCAN_YEAR=2012'
R_ = f'{S.SD}/runs'
SHIFTS = [3, 5, 7, 11, 13]
ZT = [5, 6, 7, 8, 10, 12, 15, 20]
ZL = [3, 4, 5, 6, 7, 8, 10]
THR = float(np.sqrt(0.4552790 * 0.6677760))       # S8.1: geometric midpoint, 2011 H&G 4 vs smallest 2011 body-wave ratio
FIX6 = ['GNW', 'COR', 'NEW', 'HAWA', 'OMAK', 'I02D']
KEEP = COLS_KEEP + [f'ratio_{k}' for k in FIX6]
V = dict(threshold_R1=THR)


def fixed6(D):
    X = D[[f'ratio_{k}' for k in FIX6 if f'ratio_{k}' in D]].astype(float)
    return X.notna().sum(axis=1).values, X.median(axis=1, skipna=True).values


def r1(D):
    return ((D.n_ref.fillna(0).values >= 2) & (D.ref_ratio.fillna(0).values > THR))


def r1_fixed6(D):
    n, r = fixed6(D)
    return (n >= 2) & (np.nan_to_num(r) > THR)


def r1_nf(D):
    return ((D.n_ref_nf.fillna(0).values >= 2) & (D.ref_ratio_nf.fillna(0).values > THR))


def noveto(D):
    return np.zeros(len(D), bool)


def load():
    v2 = M2.mask('v2'); v1 = M2.mask('v1')
    real = pd.read_parquet(f'{R_}/real/cand2.parquet')
    nav = np.load(f'{R_}/real/navail.npy')
    A = {s: pd.read_parquet(f'{R_}/ctrlA_s{s}/cand2.parquet') for s in SHIFTS}
    navA = {s: np.load(f'{R_}/ctrlA_s{s}/navail.npy') for s in SHIFTS}
    J2 = pd.read_parquet(f'{S.SD}/inject_{S.YEAR}.parquet'); J3 = pd.read_parquet(f'{S.SD}/inject3_{S.YEAR}.parquet')
    I2 = pd.read_parquet(f'{R_}/inj2/cand2.parquet'); I3 = pd.read_parquet(f'{R_}/inj3/cand2.parquet')
    I2 = I2[I2.variant == 0]; I3 = I3[I3.variant == 0]
    J2['cand'] = match(J2, I2); J3['cand'] = match(J3, I3)

    def attach(J, I):
        X = J.copy()
        for c in KEEP:
            X[c + '_det'] = np.where(X.cand >= 0, I.reindex(X.cand.clip(lower=0)).get(c).values if c in I else np.nan, np.nan)
        return X
    JJ = pd.concat([attach(J2, I2).assign(set='inj2'), attach(J3, I3).assign(set='inj3')], ignore_index=True)
    JJ['found'] = JJ.cand >= 0
    JJ['in_v2_true'] = v2[JJ.tsec.astype(int).values]; JJ['in_v1_true'] = v1[JJ.tsec.astype(int).values]
    JJ['loc_err_km'] = S.dist_km(JJ.lat_det, JJ.lon_det, JJ.lat, JJ.lon)
    JD = JJ.rename(columns={c + '_det': c for c in KEEP if c not in ('lat', 'lon', 'm0')})
    for c in ('fp_interior', 'hf_loud'):
        JD[c] = JD[c].astype(float) > 0.5
    return real, nav, A, navA, JD


def main():
    v2 = M2.mask('v2'); v1 = M2.mask('v1')
    btw = between_mask(); ets = ets_vec()
    real, nav, A, navA, JD = load()
    V['mask_v2_frac'] = float(v2.mean()); V['mask_v1_frac'] = float(v1.mean())
    # ---- S8.6 filter re-validation (S7.7 rule) ----
    vs = JD[(JD.gf == 'standard') & (JD.mw >= 3.4) & JD.found]
    Pv = crit(vs)
    use, fv = [], []
    for c in CRIT:
        f = float(Pv[c].mean()); row = dict(criterion=c, pass_frac=f, n=len(vs), kept=f >= 0.85)
        if c == 'max_share' and f < 0.85:
            f2 = float(Pv['share_excess'].mean()); row.update(fallback='share_excess<=0.25', fallback_pass=f2, kept=False)
            if f2 >= 0.85:
                use.append('share_excess')
        elif f >= 0.85:
            use.append(c)
        fv.append(row)
    fam = JD[(JD.gf == 'standard') & (JD.mw >= 3.6) & JD.found & (JD.loc_kind == 'family')]
    fpk = bool(fam.fp_interior.mean() >= 0.85)
    fv.append(dict(criterion='fp_interior (family injections Mw>=3.6)', pass_frac=float(fam.fp_interior.mean()), n=len(fam), kept=fpk))
    pd.DataFrame(fv).to_csv(f'{C.TB}/scan3_filter_validation.csv', index=False)
    V['filter'] = use; V['fp_interior_kept'] = fpk
    print(pd.DataFrame(fv).round(3).to_string())
    # ---- S8.4 known body-wave list 2012 ----
    r0 = real[real.variant == 0].copy()
    pop = r0[(r0.z >= 7) & (r0.z_loso.fillna(-1) >= 7) & ok_all(r0, use) & ~v1[r0.tsec.values.astype(int)]].copy()
    E = pd.read_csv(S.tab('scan2_eq_events.csv'))
    E = E[E.bw5_t0.notna()].drop_duplicates(subset=['t', 'mag'])
    y0 = S.Y0.value / 1e9
    lab = []
    for r in pop.itertuples():
        t = y0 + r.tsec
        m = E[(E.bw5_t0 <= t) & (E.bw5_t1 >= t)]
        lab.append(f'M{m.sort_values("mag").mag.iloc[-1]} {m.sort_values("mag").time.iloc[-1]} {m.sort_values("mag").dist_km.iloc[-1]:.0f} km' if len(m) else '')
    pop['bw_event'] = lab; pop['known_body_wave'] = pop.bw_event != ''
    pop.to_csv(f'{C.TB}/scan3_survivors.csv', index=False)
    bwc = pop[pop.known_body_wave].copy()
    bwc.to_csv(f'{C.TB}/scan3_bw_known_list.csv', index=False)
    print(len(pop), 'survivors;', len(bwc), 'known body-wave clusters')
    # ---- S8.5 R1 validation ----
    injv = JD[(JD.gf == 'standard') & (JD.mw >= 3.6) & JD.found & (JD.loc_kind == 'family')]
    rows = []
    for name, f in (('R1', r1), ('R1_fixed6', r1_fixed6), ('R1_no_far', r1_nf)):
        nb = int(f(bwc).sum()) if len(bwc) else 0
        need = int(np.ceil(0.8 * len(bwc) - 1e-9))
        ii = float((~f(injv)).mean()) if len(injv) else np.nan
        rows.append(dict(veto=name, bw_n=len(bwc), bw_flagged=nb, bw_needed=need, bw_pass=bool(len(bwc) > 0 and nb >= need),
                         inj_n=len(injv), inj_not_flagged=ii, inj_pass=bool(ii >= 0.85),
                         validated=bool(len(bwc) > 0 and nb >= need and ii >= 0.85),
                         bw_ratio_min=float(bwc.ref_ratio.min()) if len(bwc) else np.nan,
                         bw_ratio_med=float(bwc.ref_ratio.median()) if len(bwc) else np.nan,
                         inj_ratio_med=float(injv.ref_ratio.median()) if len(injv) else np.nan))
    VV = pd.DataFrame(rows); VV.to_csv(f'{C.TB}/scan3_veto_validation.csv', index=False)
    print(VV.round(3).to_string())
    V['R1_validated'] = bool(VV.validated.iloc[0]); V['veto_validation'] = VV.to_dict('records')
    # noise: R1 flag rate among control (a) clusters and real clusters (z >= 7, z_loso >= 7, outside v2, filter)
    def flagrate(D):
        d = D[(D.variant == 0) & (D.z >= 7) & (D.z_loso.fillna(-1) >= 7) & ~v2[D.tsec.values.astype(int)]]
        d = d[ok_all(d, use)]
        return int(len(d)), float(r1(d).mean()) if len(d) else np.nan
    V['R1_flag_rate'] = dict(real=flagrate(real), ctrlA=[flagrate(A[s]) for s in SHIFTS],
                             ctrlB=flagrate(real[real.variant > 0].assign(variant=0)))
    pts = pd.concat([r0.assign(kind='real 2012 clusters'), bwc.assign(kind='known body-wave 2012'),
                     injv.assign(kind='injection Mw>=3.6 (family)'),
                     pd.concat([A[s][A[s].variant == 0] for s in SHIFTS]).assign(kind='control (a) clusters')], ignore_index=True)
    pts[['kind', 'tsec', 'z', 'z_loso', 'ref_excess', 'ref_ratio', 'n_ref']].to_csv(f'{C.TB}/scan3_veto_points.csv', index=False)

    def analyse(veto, validated, sfx, VO, JD):
        JD = JD.copy()

        def hl(D, zt, zl, mask=v2):  # noqa
            if len(D) == 0:
                return np.zeros(0, bool)
            ok = ~mask[D.tsec.values.astype(int)] & ok_all(D, use) & (D.fp_interior.fillna(False).values.astype(bool) | (not fpk))
            ok &= ~np.asarray(veto(D), bool)
            ok &= (D.z.values >= zt) & (D.z_loso.fillna(-np.inf).values >= zl)
            return ok
        r0 = real[real.variant == 0]; rb = real[real.variant > 0]; nd = int(real.variant.max())
        lyb = live_years(nav, v2, btw); ly = live_years(nav, v2); lye = live_years(nav, v2, ets)
        fam38 = JD[(JD.gf == 'standard') & (JD.loc_kind == 'family') & (JD.mw == 3.8) & ~JD.in_v2_true]
        grid = []
        for zt in ZT:
            for zl in ZL:
                a = float(np.mean([hl(A[s][A[s].variant == 0], zt, zl).sum() / live_years(navA[s], v2) for s in SHIFTS]))
                bsel = hl(rb, zt, zl) & btw[rb.tsec.values]
                b = bsel.sum() / nd / lyb
                okr = fam38.found.values.copy(); okr[okr] = hl(fam38[okr], zt, zl)
                p0 = hl(r0, zt, zl)
                grid.append(dict(z_thr=zt, zl_thr=zl, ctrlA_per_yr=a, ctrlB_between_per_yr=b, ab=a + b, rec_mw38=float(okr.mean()),
                                 n_mw38=len(fam38), real_between_per_yr=(p0 & btw[r0.tsec.values]).sum() / lyb,
                                 real_between=int((p0 & btw[r0.tsec.values]).sum()), real_ets=int((p0 & ets[r0.tsec.values]).sum()),
                                 real_ets_per_yr=(p0 & ets[r0.tsec.values]).sum() / lye if lye > 0 else np.nan, real_all=int(p0.sum())))
        Gd = pd.DataFrame(grid); Gd.to_csv(f'{C.TB}/scan3_threshold_grid{sfx}.csv', index=False)
        q = Gd[(Gd.ctrlA_per_yr <= 0.3) & (Gd.ctrlB_between_per_yr <= 0.7)]
        if len(q):
            best = q.sort_values(['rec_mw38', 'ab', 'z_thr'], ascending=[False, True, False]).iloc[0]; VO['pair_qualifies'] = True
        else:
            q2 = Gd[Gd.rec_mw38 >= 0.5]
            best = (q2 if len(q2) else Gd).sort_values(['ab', 'rec_mw38'], ascending=[True, False]).iloc[0]; VO['pair_qualifies'] = False
        zt, zl = float(best.z_thr), float(best.zl_thr)
        VO['pair'] = dict(z_thr=zt, zl_thr=zl, **{k: float(best[k]) for k in ('ctrlA_per_yr', 'ctrlB_between_per_yr', 'ab', 'rec_mw38')})
        JD['recovered'] = False
        f_ = JD.found.values
        JD.loc[f_, 'recovered'] = hl(JD[f_], zt, zl)
        JD['recovered'] &= ~JD.in_v2_true
        JD.to_csv(f'{C.TB}/scan3_injection{sfx}.csv', index=False)
        rows = []
        for (gf, lk), g in JD[~JD.in_v2_true].groupby(['gf', 'loc_kind']):
            for mw, gg in g.groupby('mw'):
                rr = gg[gg.recovered]
                rows.append(dict(gf=gf, loc_kind=lk, mw=mw, n=len(gg), rec=gg.recovered.mean(), loc_med_km=rr.loc_err_km.median(),
                                 loc_p95_km=rr.loc_err_km.quantile(.95)))
        IS = pd.DataFrame(rows); IS.to_csv(f'{C.TB}/scan3_injection_summary{sfx}.csv', index=False)
        f38 = JD[(JD.gf == 'standard') & (JD.loc_kind == 'family') & (JD.mw == 3.8) & ~JD.in_v2_true]
        VO['rec_mw38'] = float(f38.recovered.mean()); VO['n_mw38'] = int(len(f38))
        VO['loc_med_mw38_km'] = float(f38[f38.recovered].loc_err_km.median()) if f38.recovered.any() else np.nan
        f36 = JD[(JD.gf == 'standard') & (JD.loc_kind == 'family') & (JD.mw == 3.6) & ~JD.in_v2_true]
        VO['rec_mw36'] = float(f36.recovered.mean())
        p0 = hl(r0, zt, zl)
        det = r0[p0].copy(); det['between'] = btw[det.tsec.values]; det['ets'] = ets[det.tsec.values]
        det['r1_flag'] = r1(det) if len(det) else []
        det.to_csv(f'{C.TB}/scan3_detections_2012{sfx}.csv', index=False)
        n_ets = int(det.ets.sum()); n_btw = int(det.between.sum())
        r_btw = max(n_btw, 0.5) / lyb
        lam = r_btw * lye
        VO.update(live_yr=ly, live_yr_between=lyb, live_days_ets=lye * 365.25, n_det=int(len(det)), n_ets=n_ets, n_between=n_btw,
                  real_between_per_yr=n_btw / lyb, ets_per_yr=n_ets / lye if lye > 0 else np.nan, poisson_lambda_ets=lam,
                  poisson_p=float(poisson.sf(n_ets - 1, lam)) if n_ets > 0 else 1.0,
                  detections=[dict(time_utc=str(r.time_utc), time_pt=r.time_pt, z=float(r.z), z_loso=float(r.z_loso), mw=float(r.mw),
                                   lat=float(r.lat), lon=float(r.lon), ets=bool(r.ets), between=bool(r.between),
                                   ref_ratio=float(r.ref_ratio), n_ref=int(r.n_ref) if np.isfinite(r.n_ref) else 0,
                                   hf_loud=bool(r.hf_loud)) for r in det.itertuples()])
        ab = VO['pair']['ab']
        VO['GO_items'] = dict(controls_ab_le_1=bool(VO['pair_qualifies'] and ab <= 1), strong=bool(ab <= 0.5),
                              rec_mw38_ge_80=bool(VO['rec_mw38'] >= 0.8), loc_le_15=bool(VO['loc_med_mw38_km'] <= 15),
                              ets_ge3_p001=bool(n_ets >= 3 and VO['poisson_p'] < 0.01), veto_validated=validated, hg='n/a (none listed for 2012)')
        VO['GO'] = bool(all(v for k, v in VO['GO_items'].items() if k not in ('strong', 'hg')))
        print(sfx or 'R1', json.dumps({k: VO[k] for k in ('pair', 'pair_qualifies', 'rec_mw38', 'loc_med_mw38_km', 'n_ets', 'n_between',
                                                         'poisson_p', 'GO_items', 'GO')}, default=float))
        return Gd, JD, det
    Gd, JD1, det = analyse(r1, V['R1_validated'], '', V, JD)
    V['no_veto'] = {}; Gn, JDn, detn = analyse(noveto, False, '_noveto', V['no_veto'], JD)
    V['r1_fixed6'] = {}; analyse(r1_fixed6, bool(VV.validated.iloc[1]), '_r1fixed6', V['r1_fixed6'], JD)
    V['r1_no_far'] = {}; analyse(r1_nf, bool(VV.validated.iloc[2]), '_r1nf', V['r1_no_far'], JD)
    json.dump(V, open(f'{C.TB}/scan3_verdict.json', 'w'), indent=1, default=float)
    figures(real, Gd, Gn, JD1, det, pts, V)


def figures(real, Gd, Gn, JD, det, pts, V):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    v2 = M2.mask('v2')
    zt, zl = V['pair']['z_thr'], V['pair']['zl_thr']
    zmax = np.load(f'{R_}/real/zmax.npy')
    n = 600
    zr = np.nanmax(np.where(np.isfinite(zmax), zmax, -np.inf).reshape(-1, n), axis=1)
    tt = S.Y0 + pd.to_timedelta(np.arange(len(zr)) * n, unit='s')
    mk = v2.reshape(-1, n).mean(axis=1) > 0.5
    r0 = real[(real.variant == 0) & (real.z >= 7)]
    fig, ax = plt.subplots(2, 1, figsize=(13, 7.5))
    sc = None
    for a_, (t0, t1), ttl in ((ax[0], (S.Y0, S.Y0 + pd.Timedelta(days=366)), '2012'),
                              (ax[1], (pd.Timestamp('2012-08-25'), pd.Timestamp('2012-10-05')), '2012-08 ETS')):
        sel = (tt >= t0) & (tt < t1)
        a_.plot(tt[sel], np.clip(zr[sel], 0, 40), color='0.6', lw=0.4, label='max z over grid (10-min max)')
        a_.fill_between(tt[sel], 0, 40, where=mk[sel], color='tab:orange', alpha=0.15, step='mid', label='mask v2')
        a_.axvspan(S.ETS[0], S.ETS[1], color='tab:blue', alpha=0.08, label='SVI ETS 2012-08')
        tr = S.Y0 + pd.to_timedelta(r0.tsec, unit='s'); s_ = (tr >= t0) & (tr < t1)
        fl = r1(r0)
        a_.scatter(tr[s_ & fl], np.clip(r0.z[s_ & fl], 0, 40), marker='x', color='tab:red', s=14, lw=0.8, label='z ≥ 7 clusters flagged by R1')
        tdet = S.Y0 + pd.to_timedelta(det.tsec, unit='s'); sd = (tdet >= t0) & (tdet < t1)
        sc = a_.scatter(tdet[sd], np.clip(det.z[sd], 0, 40), c=det.z[sd].clip(0, 40), cmap='viridis', vmin=0, vmax=40, s=40,
                        zorder=3, edgecolor='k', lw=0.4, label=f'headline detections (z ≥ {zt:g}, z_loso ≥ {zl:g}, R1)')
        a_.set_ylim(0, 40); a_.set_ylabel('z (clipped at 40)'); a_.set_title(f'VLFE scan pilot 3 (R1 ratio veto), {ttl}', fontsize=9)
        a_.legend(fontsize=7, loc='upper left', ncol=2)
    fig.colorbar(sc, ax=ax, label='z')
    fig.savefig(f'{C.FIG}/scan3_z_2012.png', dpi=130); plt.close(fig)
    # separation with / without R1, at the chosen zl and at zl = 7, 10
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3), sharey=True)
    for a_, zz in zip(ax, sorted({zl, 5.0, 7.0, 10.0})[-3:]):
        for G_, ls, tag in ((Gd, '-', 'with R1'), (Gn, ':', 'no veto')):
            g = G_[G_.zl_thr == zz]
            a_.semilogy(g.z_thr, g.real_between_per_yr.clip(lower=0.05), 'o' + ls, color='tab:red', label=f'real between ETS, {tag}')
            a_.semilogy(g.z_thr, g.real_ets_per_yr.clip(lower=0.05), 'D' + ls, color='tab:purple', ms=4, label=f'real in ETS (per ETS live yr), {tag}')
            a_.semilogy(g.z_thr, g.ctrlA_per_yr.clip(lower=0.05), '^' + ls, color='0.3', label=f'(a) shifted, {tag}')
            a_.semilogy(g.z_thr, g.ctrlB_between_per_yr.clip(lower=0.05), 'v' + ls, color='tab:green', label=f'(b) polarity, {tag}')
        a_.axhline(0.3, color='0.3', lw=0.5, ls='--'); a_.axhline(0.7, color='tab:green', lw=0.5, ls='--')
        a_.set_xlabel('z threshold'); a_.set_title(f'2012, z_loso ≥ {zz:g}', fontsize=9)
    ax[0].set_ylabel('per live year (floor 0.05)'); ax[0].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan3_separation.png', dpi=130); plt.close(fig)
    # veto plane: ref_ratio vs z, coloured by z
    fig, a_ = plt.subplots(figsize=(7.5, 5.5))
    for k, mkr, s_ in (('control (a) clusters', '.', 4), ('real 2012 clusters', '.', 6)):
        q = pts[pts.kind == k]
        sc = a_.scatter(q.ref_ratio, q.z, c=q.z.clip(5, 30), cmap='Greys' if 'control' in k else 'viridis', marker=mkr, s=s_,
                        alpha=0.5, label=k)
    for k, col, mkr in (('injection Mw>=3.6 (family)', 'tab:green', 'o'), ('known body-wave 2012', 'tab:red', 's')):
        q = pts[pts.kind == k]
        a_.scatter(q.ref_ratio, q.z, color=col, marker=mkr, s=28, edgecolor='k', lw=0.3, label=f'{k} (n={len(q)})')
    a_.axvline(THR, color='k', lw=0.8, ls='--', label=f'R1 threshold {THR:.3f} (fixed on 2011)')
    a_.set_xscale('log'); a_.set_yscale('log'); a_.set_xlabel('ref_ratio (median reference RMS / near RMS)'); a_.set_ylabel('z')
    a_.legend(fontsize=7); fig.colorbar(sc, ax=a_, label='z (real clusters)')
    a_.set_title('2012: ratio veto R1, out of sample', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan3_veto.png', dpi=130); plt.close(fig)
    # injections
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for lk, col in (('family', 'tab:red'), ('off_family', 'tab:blue')):
        g = JD[(JD.gf == 'standard') & (JD.loc_kind == lk) & ~JD.in_v2_true].groupby('mw').recovered.mean()
        ax[0].plot(g.index, g.values, 'o-', color=col, label=lk)
    ax[0].axhline(0.8, color='k', lw=0.5, ls=':'); ax[0].set_xlabel('injected Mw'); ax[0].set_ylabel('recovered (all headline rules, R1)')
    ax[0].legend(fontsize=7)
    r = JD[JD.recovered]
    sc = ax[1].scatter(r.mw + np.random.default_rng(1).uniform(-0.05, 0.05, len(r)), r.loc_err_km, c=r.z.clip(0, 30), s=14, cmap='viridis')
    ax[1].set_xlabel('injected Mw'); ax[1].set_ylabel('location error (km)'); fig.colorbar(sc, ax=ax[1], label='z')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan3_injection.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
