#!/usr/bin/env python
"""VLFE scan pilot 2: validation, threshold pair, controls, GO/NO-GO, classifications, figures (PLAN_SCAN.md S7.2-S7.10).
Inputs: data/scan/runs/{real_c, ctrlA_s*_c, inj2_c, inj3_c}/cand2.parquet (scan_columns2.py), inject lists, masks2,
tables/scan2_bw_known_list.csv, tables/scan2_decay_verdict.json, pilot-1 tables for the cleaner comparison.
Writes tables/scan2_cleaner.csv, scan2_filter_validation.csv, scan2_veto_validation.csv, scan2_veto_points.csv,
scan2_threshold_grid.csv, scan2_injection.csv, scan2_injection_summary.csv, scan2_detections_2011.csv, scan2_hg.csv,
scan2_classify.csv, scan2_ext_grid.csv, scan2_verdict.json, figures/scan2_*.png."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
from scipy.stats import poisson
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
import scan_masks2 as M2
import scan_columns2 as SC2
from scan_validate import between_mask, ets_vec, live_years

R_ = f'{S.SD}/runs'
SHIFTS = [3, 5, 7, 11, 13]
ZT = [5, 6, 7, 8, 10, 12, 15, 20]
ZL = [3, 4, 5, 6, 7, 8, 10]
P1_FILTER = ['n_sta', 'n_Z_agree', 'share_excess', 'vr_ratio', 'vr_gt_ss', 'vr_gt_steep', 'isolation', 'not_hf_loud']
CRIT = ['n_sta', 'n_Z_agree', 'max_share', 'vr_ratio', 'vr_gt_ss', 'vr_gt_steep', 'isolation']
def _veto(k, nf=False):
    sfx = '_nf' if nf else ''
    def f(D):
        n, e, r = D['n_ref' + sfx], D['ref_excess' + sfx], D['ref_ratio' + sfx]
        base = (n >= 1)
        return {'V1': base & (e >= 3) & (r >= 0.5), 'V2': base & (e >= 3), 'V3': base & (e >= 2) & (r >= 0.3)}[k].fillna(False).values.astype(bool)
    return f


VETOS = {k: _veto(k) for k in ('V1', 'V2', 'V3')}
VETOS_NF = {k: _veto(k, True) for k in ('V1', 'V2', 'V3')}
CLASS_T = [('08-18 21:19', '2011-08-18 21:19:46'), ('08-22 08:21', '2011-08-22 08:21:49'),
           ('08-22 08:25', '2011-08-22 08:25:25'), ('08-22 08:28', '2011-08-22 08:28:22')]
V = {}


def crit(D):
    P = pd.DataFrame(index=D.index)
    P['n_sta'] = D.n_sta >= 3
    P['n_Z_agree'] = D.n_Z_agree >= 2
    P['max_share'] = D.max_share <= 0.5
    P['share_excess'] = D.share_excess <= 0.25
    P['vr_ratio'] = D.vr_ratio >= 0.8
    P['vr_gt_ss'] = D.vr_thrust > D.vr_strikeslip
    P['vr_gt_steep'] = D.vr_thrust > D.vr_steepdip
    P['isolation'] = D.isolation < 2
    P['not_hf_loud'] = ~D.hf_loud.fillna(False).astype(bool)
    return P


def ok_all(D, use):
    P = crit(D)
    ok = np.ones(len(D), bool)
    for c in use:
        ok &= P[c].fillna(False).values.astype(bool)
    return ok


def ts(s):
    return int((pd.Timestamp(s) - S.Y0).total_seconds())


def match(J, D, win=60):
    """per injection: index of the largest-z positive cluster within +-win s (or -1)."""
    D = D.sort_values('tsec'); t = D.tsec.values; z = D.z.values; idx = D.index.values
    out = []
    for x in J.tsec.values:
        lo, hi = np.searchsorted(t, x - win), np.searchsorted(t, x + win, 'right')
        out.append(idx[lo + int(np.argmax(z[lo:hi]))] if hi > lo else -1)
    return np.array(out)


def peak_rows(run, times, win=60):
    """largest positive peak within +-win s of each time in run's cand.parquet (variant 0)."""
    Cd = pd.read_parquet(f'{R_}/{run}/cand.parquet')
    Cd = Cd[(Cd.variant == 0) & (Cd.sign > 0)]
    rows = []
    for lab, t in times:
        m = Cd[(Cd.tsec - t).abs() <= win]
        if len(m):
            r = m.loc[m.z.idxmax()].copy(); r['label'] = lab; rows.append(r)
        else:
            rows.append(pd.Series(dict(variant=0, sign=1, tsec=t, node=-1, z=np.nan, m0=np.nan, label=lab)))
    return pd.DataFrame(rows).reset_index(drop=True)


CLEAN = dict(real='real_c', ctrlA='ctrlA_s{s}_c', inj2='inj2_c', inj3='inj3_c')
FALLBACK = dict(real='real', ctrlA='ctrlA_s{s}', inj2='inj2', inj3='inj3_u')          # S7.12
RUNS = dict(CLEAN)
COLS_KEEP = ['z', 'z_loso', 'node', 'lat', 'lon', 'm0', 'n_sta', 'n_Z_agree', 'max_share', 'share_excess', 'vr_ratio',
             'vr_thrust', 'vr_strikeslip', 'vr_steepdip', 'isolation', 'hf_loud', 'fp_interior', 'grid_edge', 'n_ref',
             'ref_excess', 'ref_ratio', 'n_ref_nf', 'ref_excess_nf', 'ref_ratio_nf', 'dose_excess', 'in_v1', 'in_v2', 'cc']


def load(runs):
    v2 = M2.mask('v2'); v1 = M2.mask('v1')
    real = pd.read_parquet(f'{R_}/{runs["real"]}/cand2.parquet')
    nav = np.load(f'{R_}/{runs["real"]}/navail.npy')
    A = {s: pd.read_parquet(f'{R_}/{runs["ctrlA"].format(s=s)}/cand2.parquet') for s in SHIFTS}
    navA = {s: np.load(f'{R_}/{runs["ctrlA"].format(s=s)}/navail.npy') for s in SHIFTS}
    J2 = pd.read_parquet(f'{S.SD}/inject_2011.parquet'); J3 = pd.read_parquet(f'{S.SD}/inject3_2011.parquet')
    I2 = pd.read_parquet(f'{R_}/{runs["inj2"]}/cand2.parquet'); I3 = pd.read_parquet(f'{R_}/{runs["inj3"]}/cand2.parquet')
    I2 = I2[I2.variant == 0]; I3 = I3[I3.variant == 0]
    J2['cand'] = match(J2, I2); J3['cand'] = match(J3, I3)

    def attach(J, I):
        X = J.copy()
        for c in COLS_KEEP:
            X[c + '_det'] = np.where(X.cand >= 0, I.reindex(X.cand.clip(lower=0)).get(c).values if c in I else np.nan, np.nan)
        return X
    JJ = pd.concat([attach(J2, I2).assign(set='inj2'), attach(J3, I3).assign(set='inj3')], ignore_index=True)
    JJ['found'] = JJ.cand >= 0
    JJ['in_v2_true'] = v2[JJ.tsec.astype(int).values]; JJ['in_v1_true'] = v1[JJ.tsec.astype(int).values]
    JJ['loc_err_km'] = S.dist_km(JJ.lat_det, JJ.lon_det, JJ.lat, JJ.lon)
    JD = JJ.rename(columns={c + '_det': c for c in COLS_KEEP if c not in ('lat', 'lon', 'm0')})
    for c in ('fp_interior', 'hf_loud'):
        JD[c] = JD[c].astype(float) > 0.5
    return real, nav, A, navA, JD


def main():
    v2 = M2.mask('v2'); v1 = M2.mask('v1')
    btw = between_mask(); ets = ets_vec()
    real, nav, A, navA, JD = load(CLEAN)

    # ---------- S7.5 cleaner acceptance ----------
    P1 = pd.read_csv(f'{C.TB}/scan_injection.csv')
    p1_rec = float(P1[(P1.gf == 'standard') & (P1.mw >= 3.6)].recovered.mean())
    s2 = JD[(JD.set == 'inj2') & (JD.gf == 'standard') & (JD.mw >= 3.6)]
    okf = np.zeros(len(s2), bool); hz = s2.found.values
    okf[hz] = ok_all(s2[hz], P1_FILTER)
    rec2 = hz & okf & (s2.z.fillna(0).values >= 7) & ~s2.in_v1_true.values
    p2_rec = float(rec2.mean())
    a_rates = []
    for s in SHIFTS:
        d = A[s][(A[s].variant == 0) & ~A[s].in_v1 & (A[s].z >= 7)]
        a_rates.append(ok_all(d, P1_FILTER).sum() / live_years(navA[s], v1))
    fl = pd.concat([pd.read_csv(f'{R_}/{r}/clean_flags.csv').assign(run=r) for r in
                    ['real_c'] + [f'ctrlA_s{s}_c' for s in SHIFTS] + ['inj2_c', 'inj3_c']])
    fl.to_csv(f'{C.TB}/scan2_cleaner.csv', index=False)
    V['cleaner'] = dict(ctrlA_z7_per_yr_after=float(np.mean(a_rates)), ctrlA_each=[float(x) for x in a_rates],
                        ctrlA_z7_per_yr_pilot1=714.0, inj_rec_mw36_pilot1=p1_rec, inj_rec_mw36_after=p2_rec,
                        keep_frac=p2_rec / p1_rec if p1_rec > 0 else np.nan,
                        cleaned_frac_real={r.sta: float(r.cleaned_s / (365 * 86400)) for r in fl[fl.run == 'real_c'].itertuples()})
    V['cleaner']['inj_pass'] = bool(V['cleaner']['keep_frac'] >= 0.85)
    V['cleaner']['ctrlA_pass'] = bool(V['cleaner']['ctrlA_z7_per_yr_after'] <= 20)
    V['cleaner']['used'] = V['cleaner']['inj_pass']
    print('cleaner', json.dumps(V['cleaner'], indent=1, default=float))
    if not V['cleaner']['used']:
        RUNS.update(FALLBACK)
        real, nav, A, navA, JD = load(RUNS)
    V['runs'] = dict(RUNS)

    # ---------- S7.7 filter re-validation (pilot-2 injections, standard GF, Mw >= 3.4, found z >= 5) ----------
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
    fv.append(dict(criterion='fp_interior (family injections Mw>=3.6)', pass_frac=float(fam.fp_interior.mean()), n=len(fam),
                   kept=bool(fam.fp_interior.mean() >= 0.85)))
    fv.append(dict(criterion='hf_loud (class only)', pass_frac=float((~vs.hf_loud.fillna(False).astype(bool)).mean()), n=len(vs), kept=False))
    pd.DataFrame(fv).to_csv(f'{C.TB}/scan2_filter_validation.csv', index=False)
    V['filter'] = use; V['fp_interior_kept'] = bool(fam.fp_interior.mean() >= 0.85)
    print(pd.DataFrame(fv).round(3).to_string())

    # ---------- S7.3 veto validation ----------
    hg = peak_rows(RUNS['real'], [(f'H&G {k}', ts(t)) for k, t, *_ in S.HG])
    bw = pd.read_csv(f'{C.TB}/scan2_bw_known_list.csv')
    p1 = pd.read_parquet(f'{R_}/real/loso.parquet')
    bw = bw.merge(p1[p1.variant == 0][['tsec', 'm0']], on='tsec', how='left')
    bwr = pd.DataFrame(dict(variant=0, sign=1, tsec=bw.tsec.values, node=bw.node.values, m0=bw.m0.values, z=bw.z.values,
                            label=['bw ' + str(t)[:16] for t in bw.time_utc]))
    clr = peak_rows(RUNS['real'], [(lab, ts(t)) for lab, t in CLASS_T])
    p1c = []
    for lab, t in CLASS_T:
        m = p1[(p1.variant == 0) & ((p1.tsec - ts(t)).abs() <= 5)]
        if len(m):
            p1c.append(dict(variant=0, sign=1, tsec=int(m.tsec.iloc[0]), node=int(m.node.iloc[0]), m0=float(m.m0.iloc[0]),
                            z=float(m.z.iloc[0]), label=lab + ' (pilot-1 time/node)'))
    EV = pd.concat([hg, bwr, clr, pd.DataFrame(p1c)], ignore_index=True)
    EVok = EV[EV.node >= 0].copy(); EVok['node'] = EVok.node.astype(int); EVok['tsec'] = EVok.tsec.astype(int)
    EVc = SC2.eval_rows(RUNS['real'], f'{S.SD}/bp_2011g', f'{S.SD}/bp2_ref', EVok[['variant', 'sign', 'tsec', 'node', 'm0', 'z', 'label']].copy(), 8)
    EVc = SC2.add_meta(EVc)
    EVc.to_csv(f'{C.TB}/scan2_validation_events.csv', index=False)
    hgc = EVc[EVc.label.str.startswith('H&G')]; bwc = EVc[EVc.label.str.startswith('bw ')]
    injv = JD[(JD.gf == 'standard') & (JD.mw >= 3.6) & JD.found & (JD.loc_kind == 'family')]
    vrows = []
    chosen = {}
    for setname, VS in (('all', VETOS), ('no_far', VETOS_NF)):
        chosen[setname] = None
        for k, f in VS.items():
            h24 = hgc[hgc.label.isin(['H&G 2', 'H&G 4'])]
            i_ok = bool(len(h24) == 2 and not f(h24).any())
            nflag = int(f(bwc).sum()); need = int(np.ceil(0.8 * len(bw) - 1e-9))
            ii = float((~f(injv)).mean()) if len(injv) else np.nan
            ok = i_ok and nflag >= need and ii >= 0.85
            vrows.append(dict(ref_set=setname, variant=k, hg2_hg4_not_flagged=i_ok, bw_flagged=nflag, bw_n=len(bw),
                              bw_needed=need, inj_not_flagged=ii, inj_n=len(injv), passes=ok,
                              bw_flagged_list=';'.join(bwc.label[f(bwc)].tolist())))
            if ok and chosen[setname] is None:
                chosen[setname] = k
    VV = pd.DataFrame(vrows); VV.to_csv(f'{C.TB}/scan2_veto_validation.csv', index=False)
    print(VV.drop(columns='bw_flagged_list').to_string())
    V['veto_chosen'] = chosen['all']; V['veto_validated'] = chosen['all'] is not None
    V['veto_chosen_no_far'] = chosen['no_far']
    none = (lambda D: np.zeros(len(D), bool))
    veto = VETOS[chosen['all']] if chosen['all'] else none
    veto_nf = VETOS_NF[chosen['no_far']] if chosen['no_far'] else none
    chosen = chosen['all']
    pts = pd.concat([real[real.variant == 0].assign(kind='real 2011 clusters'), bwc.assign(kind='known body-wave'),
                     hgc.assign(kind='H&G'), injv.assign(kind='injection Mw>=3.6 (family)')], ignore_index=True)
    pts[['kind', 'tsec', 'z', 'ref_excess', 'ref_ratio', 'n_ref', 'label']].to_csv(f'{C.TB}/scan2_veto_points.csv', index=False)

    def analyse(veto, validated, sfx, VO, JD):
        JD = JD.copy()
        # ---------- headline + threshold pair ----------
        def hl(D, zt, zl, mask=v2):  # noqa
            if len(D) == 0:
                return np.zeros(0, bool)
            ok = ~mask[D.tsec.values.astype(int)] & ok_all(D, use) & (D.fp_interior.fillna(False).values.astype(bool) | (not V['fp_interior_kept']))
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
                                 real_all=int(p0.sum())))
        Gd = pd.DataFrame(grid); Gd.to_csv(f'{C.TB}/scan2_threshold_grid{sfx}.csv', index=False)
        q = Gd[(Gd.ctrlA_per_yr <= 0.3) & (Gd.ctrlB_between_per_yr <= 0.7)]
        if len(q):
            best = q.sort_values(['rec_mw38', 'ab', 'z_thr'], ascending=[False, True, False]).iloc[0]; VO['pair_qualifies'] = True
        else:
            q2 = Gd[Gd.rec_mw38 >= 0.5]
            best = (q2 if len(q2) else Gd).sort_values(['ab', 'rec_mw38'], ascending=[True, False]).iloc[0]; VO['pair_qualifies'] = False
        zt, zl = float(best.z_thr), float(best.zl_thr)
        VO['pair'] = dict(z_thr=zt, zl_thr=zl, **{k: float(best[k]) for k in ('ctrlA_per_yr', 'ctrlB_between_per_yr', 'ab', 'rec_mw38')})
        print(Gd.round(3).to_string()); print('pair', VO['pair'])

        # ---------- at the chosen pair ----------
        JD['recovered'] = False
        f_ = JD.found.values
        JD.loc[f_, 'recovered'] = hl(JD[f_], zt, zl)
        JD['recovered'] &= ~JD.in_v2_true
        JD.to_csv(f'{C.TB}/scan2_injection{sfx}.csv', index=False)
        rows = []
        for (gf, lk), g in JD[~JD.in_v2_true].groupby(['gf', 'loc_kind']):
            for mw, gg in g.groupby('mw'):
                rr = gg[gg.recovered]
                rows.append(dict(gf=gf, loc_kind=lk, mw=mw, n=len(gg), rec=gg.recovered.mean(), loc_med_km=rr.loc_err_km.median(),
                                 loc_p95_km=rr.loc_err_km.quantile(.95)))
        IS = pd.DataFrame(rows); IS.to_csv(f'{C.TB}/scan2_injection_summary{sfx}.csv', index=False)
        print(IS.round(3).to_string())
        f38 = JD[(JD.gf == 'standard') & (JD.loc_kind == 'family') & (JD.mw == 3.8) & ~JD.in_v2_true]
        VO['rec_mw38'] = float(f38.recovered.mean()); VO['n_mw38'] = int(len(f38))
        VO['loc_med_mw38_km'] = float(f38[f38.recovered].loc_err_km.median()) if f38.recovered.any() else np.nan
        f36 = JD[(JD.gf == 'standard') & (JD.loc_kind == 'family') & (JD.mw == 3.6) & ~JD.in_v2_true]
        VO['rec_mw36'] = float(f36.recovered.mean())
        p0 = hl(r0, zt, zl)
        det = r0[p0].copy(); det['between'] = btw[det.tsec.values]; det['ets'] = ets[det.tsec.values]
        det.to_csv(f'{C.TB}/scan2_detections_2011{sfx}.csv', index=False)
        n_ets = int(det.ets.sum()); n_btw = int(det.between.sum())
        r_btw = max(n_btw, 0.5) / lyb
        lam = r_btw * lye
        VO.update(live_yr=ly, live_yr_between=lyb, live_days_ets=lye * 365.25, n_det=int(len(det)), n_ets=n_ets, n_between=n_btw,
                 real_between_per_yr=n_btw / lyb, ets_per_yr=n_ets / lye if lye > 0 else np.nan, poisson_lambda_ets=lam,
                 poisson_p=float(poisson.sf(n_ets - 1, lam)) if n_ets > 0 else 1.0)
        hgc2 = hgc.copy(); hgc2['headline'] = hl(hgc2, zt, zl)
        hgc2['veto_flag'] = np.asarray(veto(hgc2), bool); hgc2['in_v2'] = v2[hgc2.tsec.values.astype(int)]
        hgc2.to_csv(f'{C.TB}/scan2_hg{sfx}.csv', index=False)
        h24 = hgc2[hgc2.label.isin(['H&G 2', 'H&G 4'])]
        VO['hg24_headline'] = bool(len(h24) == 2 and h24.headline.all())
        VO['hg'] = {r.label: dict(z=float(r.z), z_loso=float(r.z_loso), headline=bool(r.headline), in_v2=bool(r.in_v2),
                                 veto=bool(r.veto_flag), fp=bool(r.fp_interior), ref_excess=float(r.ref_excess)) for r in hgc2.itertuples()}
        ab = VO['pair']['ab']
        VO['GO_items'] = dict(controls_ab_le_1=bool(VO['pair_qualifies'] and ab <= 1), strong=bool(ab <= 0.5),
                             rec_mw38_ge_80=bool(VO['rec_mw38'] >= 0.8), loc_le_15=bool(VO['loc_med_mw38_km'] <= 15),
                             hg24=VO['hg24_headline'], ets_ge3_p001=bool(n_ets >= 3 and VO['poisson_p'] < 0.01),
                             veto_validated=validated)
        VO['GO'] = bool(all(v for k, v in VO['GO_items'].items() if k != 'strong'))

        return hl, Gd, JD, det, zt, zl

    hl, Gd, JD, det, zt, zl = analyse(veto, V['veto_validated'], '', V, JD)
    V['no_far'] = {}
    analyse(veto_nf, V['veto_chosen_no_far'] is not None, '_nf', V['no_far'], JD)
    # ---------- S7.10 classification ----------
    V['classify'] = classify(EVc, clr, veto, use, zt, zl, hl)
    json.dump(V, open(f'{C.TB}/scan2_verdict.json', 'w'), indent=1, default=float)
    print(json.dumps(V, indent=1, default=float))
    figures(real, A, Gd, JD, det, pts, V, zt, zl, chosen)


def classify(EVc, clr, veto, use, zt, zl, hl):
    out = {}
    fiv = f'{R_}/{RUNS["real"]}/clean_iv.parquet'
    iv = pd.read_parquet(fiv) if os.path.isfile(fiv) else pd.DataFrame(columns=['sta', 't0', 't1'])
    E2 = pd.read_csv(f'{C.TB}/scan2_eq_events.csv')
    rows = []
    for lab, t in CLASS_T:
        tt = ts(t)
        r = EVc[EVc.label == lab]
        p1r = EVc[EVc.label == lab + ' (pilot-1 time/node)']
        cl = iv[(iv.t0 <= tt + 120) & (iv.t1 >= tt - 60)]
        rec = dict(label=lab, time_utc=t, time_pt=pd.Timestamp(t).tz_localize('UTC').tz_convert('America/Los_Angeles').strftime('%Y-%m-%d %H:%M:%S %Z'),
                   cleaned_stations=','.join(sorted(cl.sta.unique())))
        for pre, X in (('c_', r), ('p1_', p1r)):
            if len(X):
                x = X.iloc[0]
                for k in ('z', 'z_loso', 'loso_dropped', 'node', 'lat', 'lon', 'mw', 'n_sta', 'n_ref', 'ref_excess', 'ref_ratio',
                          'dose_excess', 'pw_n', 'pw_vr', 'pw_slow', 'pw_baz', 'fp_interior', 'grid_edge', 'hf_loud', 'in_head2',
                          'in_bw5', 'in_reg', 'in_reg3', 'in_m45', 'ref_list', 'cc'):
                    rec[pre + k] = x.get(k, np.nan)
                rec[pre + 'veto'] = bool(np.asarray(veto(X.iloc[:1]))[0])
        # catalogue events whose any window covers t
        y = S.Y0.value / 1e9 + tt
        hits = []
        for tier in ('head2', 'bw5', 'reg', 'reg3', 'm45'):
            m = E2[(E2[f'{tier}_t0'] <= y) & (E2[f'{tier}_t1'] >= y)]
            for e in m.drop_duplicates(subset=['t']).itertuples():
                hits.append(f'{tier}: {e.time} M{e.mag} {e.magType} {e.dist_km:.0f} km')
        rec['catalogue_hits'] = ' | '.join(hits[:6])
        # rule order (S7.10)
        z_c = rec.get('c_z', np.nan); zl_c = rec.get('c_z_loso', np.nan)
        if not np.isfinite(z_c) or z_c < 7 or (np.isfinite(zl_c) and zl_c < 5) or not np.isfinite(zl_c):
            cls = 'instrument (removed by the step cleaner or single-station: z_loso < 5)'
        elif any(bool(rec.get('c_in_' + k, False)) for k in ('head2', 'bw5', 'reg', 'reg3')):
            cls = 'catalogued earthquake'
        elif rec.get('c_veto', False):
            cls = 'travelling wave (veto)'
        elif rec.get('c_n_sta', 0) >= 3 and np.isfinite(rec.get('c_ref_excess', np.nan)) and rec['c_ref_excess'] < 3:
            cls = 'local candidate'
        else:
            cls = 'unresolved'
        rec['class'] = cls
        rows.append(rec)
    T = pd.DataFrame(rows)
    # extended grid for 08-18 if local candidate
    r18 = T[T.label == '08-18 21:19'].iloc[0]
    if r18['class'] == 'local candidate':
        T.loc[T.label == '08-18 21:19', 'ext_grid'] = json.dumps(ext_grid(ts(CLASS_T[0][1])), default=float)
    T.to_csv(f'{C.TB}/scan2_classify.csv', index=False)
    print(T.T.to_string())
    for _, r in T.iterrows():
        out[r['label']] = dict(cls=r['class'], cleaned=r['cleaned_stations'], z=r.get('c_z'), z_loso=r.get('c_z_loso'),
                               ref_excess=r.get('c_ref_excess'), ref_ratio=r.get('c_ref_ratio'), n_ref=r.get('c_n_ref'),
                               refs=r.get('c_ref_list'), catalogue=r['catalogue_hits'], ext_grid=r.get('ext_grid'))
    return out


def ext_grid(tt):
    """S7.10: z on a grid extended to 49.6 N at time tt (cleaned real_c data, same templates/gains/weights; daily MAD of
    the statistic per node computed on this day)."""
    from scipy.signal import fftconvolve
    from matplotlib.path import Path  # noqa
    import s13_green as G
    kx = 111.19 * np.cos(np.radians(S.LAT0))
    ys = np.arange((47.8 - S.LAT0) * 111.19, (49.6 - S.LAT0) * 111.19 + 1e-6, 10.0)
    xs = np.arange((-124.8 - S.LON0) * kx, (-122.6 - S.LON0) * kx + 1e-6, 10.0)
    X, Y = np.meshgrid(xs, ys)
    lat = S.LAT0 + Y.ravel() / 111.19; lon = S.LON0 + X.ravel() / kx
    dep = S.slab_depth(lat, lon); ok = np.isfinite(dep)
    lat, lon, dep = lat[ok], lon[ok], dep[ok]
    SC2.init(RUNS['real'], f'{S.SD}/bp_2011g', f'{S.SD}/bp2_ref')
    W = SC2.W
    stas = W['stas']; d, u = divmod(tt, 86400)
    cc = W['cc']
    Mt = G.mt_ned(*S.MECHS['thrust'], 1.0)
    g = W['g'][d]
    data = {}
    for i, s in enumerate(stas):
        x = SC2.cleaned_window(i, s, d * 86400, d * 86400, 0, 86400)
        sig = W['DS'][s]['sig'][d] * 1e9
        for c in range(3):
            if W['won'][d, i, c] and g[i, c] > 0 and np.isfinite(sig[c]) and sig[c] > 0:
                data[(i, c)] = (x[c], 1 / sig[c] ** 2)
    res = []
    for la, lo, de in zip(lat, lon, dep):
        num = np.zeros(86400); den = np.zeros(86400)
        tcache = {}
        for (i, c), (x, w) in data.items():
            if i not in tcache:
                tcache[i] = S.template(Mt, (la, lo, de), (cc.loc[stas[i], 'lat'], cc.loc[stas[i], 'lon'])) * 1e21
            T = tcache[i][c] * g[i, c]
            fin = np.isfinite(x)
            xc = np.where(fin, x, 0.0)
            cor = fftconvolve(np.r_[xc, np.zeros(S.NT)], T[::-1], mode='full')[S.NT - 1:S.NT - 1 + 86400]
            # window [u-60, u+120] -> correlation index u-60 (template origin at index 60)
            cor = np.r_[np.full(60, np.nan), cor[:86400 - 60]]
            badc = np.r_[0, np.cumsum(~fin)]
            uu = np.arange(86400); a0 = np.clip(uu - 60, 0, 86400); a1 = np.clip(uu + 121, 0, 86400)
            av = ((badc[a1] - badc[a0]) == 0) & (uu - 60 >= 0) & (uu + 121 <= 86400)
            num += np.where(av, w * cor, 0); den += np.where(av, w * (T ** 2).sum(), 0)
        st = num / np.sqrt(np.where(den > 0, den, np.nan))
        mk = M2.mask('head')[d * 86400:(d + 1) * 86400]
        v = st[::4][~mk[::4] & np.isfinite(st[::4])]
        mad = 1.4826 * np.median(np.abs(v - np.median(v))) if v.size > 1000 else np.nan
        zz = st[max(u - 20, 0):u + 21] / mad
        res.append(dict(lat=la, lon=lo, dep=de, z=float(np.nanmax(zz)) if np.isfinite(zz).any() else np.nan))
    E = pd.DataFrame(res); E.to_csv(f'{C.TB}/scan2_ext_grid.csv', index=False)
    b = E.loc[E.z.idxmax()]
    top = E.lat.max()
    return dict(best_lat=b.lat, best_lon=b.lon, best_z=b.z, on_north_edge=bool(abs(b.lat - top) < 0.01),
                on_ew_edge=bool(abs(b.lon - E.lon.min()) < 0.01 or abs(b.lon - E.lon.max()) < 0.01), n_nodes=len(E))


def figures(real, A, Gd, JD, det, pts, V, zt, zl, chosen):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    v2 = M2.mask('v2')
    zmax = np.load(f'{R_}/{RUNS["real"]}/zmax.npy')
    n = 600
    zr = np.nanmax(np.where(np.isfinite(zmax), zmax, -np.inf).reshape(-1, n), axis=1)
    tt = S.Y0 + pd.to_timedelta(np.arange(len(zr)) * n, unit='s')
    mk = v2.reshape(-1, n).mean(axis=1) > 0.5
    fig, ax = plt.subplots(2, 1, figsize=(13, 7.5))
    for a_, (t0, t1), ttl in ((ax[0], (S.Y0, S.Y0 + pd.Timedelta(days=365)), '2011'),
                              (ax[1], (pd.Timestamp('2011-08-10'), pd.Timestamp('2011-09-06')), '2011-08 ETS')):
        sel = (tt >= t0) & (tt < t1)
        a_.plot(tt[sel], np.clip(zr[sel], 0, 40), color='0.6', lw=0.4, label='max z over grid, cleaned (10-min max)')
        a_.fill_between(tt[sel], 0, 40, where=mk[sel], color='tab:orange', alpha=0.15, step='mid', label='mask v2')
        a_.axvspan(S.ETS[0], S.ETS[1], color='tab:blue', alpha=0.08, label='SVI ETS')
        tdet = S.Y0 + pd.to_timedelta(det.tsec, unit='s'); sd = (tdet >= t0) & (tdet < t1)
        sc = a_.scatter(tdet[sd], np.clip(det.z[sd], 0, 40), c=det.cc[sd], cmap='viridis', vmin=0, vmax=1, s=30, zorder=3,
                        edgecolor='k', lw=0.3, label=f'headline detections (z ≥ {zt:g}, z_loso ≥ {zl:g})')
        if ttl != '2011':
            for k, ts_, *_ in S.HG:
                a_.axvline(pd.Timestamp(ts_), color='r', lw=0.8, ymax=0.15)
                a_.annotate(f'H&G {k}', (pd.Timestamp(ts_), 36), color='r', fontsize=7, rotation=90)
            for lab, t in CLASS_T:
                a_.axvline(pd.Timestamp(t), color='m', lw=0.8, ymin=0.85)
            a_.annotate('08-18 21:19', (pd.Timestamp(CLASS_T[0][1]), 30), color='m', fontsize=7, rotation=90)
            a_.annotate('08-22 triple', (pd.Timestamp(CLASS_T[1][1]), 30), color='m', fontsize=7, rotation=90)
        a_.set_ylim(0, 40); a_.set_ylabel('z (clipped at 40)'); a_.set_title(f'VLFE scan pilot 2 (step-cleaned), {ttl}', fontsize=9)
        a_.legend(fontsize=7, loc='upper left', ncol=2)
    fig.colorbar(sc, ax=ax, label='network CC')
    fig.savefig(f'{C.FIG}/scan2_z_2011.png', dpi=130); plt.close(fig)
    # separation at the chosen zl
    g = Gd[Gd.zl_thr == zl]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].semilogy(g.z_thr, g.real_between_per_yr.clip(lower=0.05), 'o-', color='tab:red', label='real, between ETS')
    ax[0].semilogy(g.z_thr, g.ctrlA_per_yr.clip(lower=0.05), '^-', color='0.3', label='(a) station-day shifted, mean of 5')
    ax[0].semilogy(g.z_thr, g.ctrlB_between_per_yr.clip(lower=0.05), 'v-', color='tab:green', label='(b) polarity scrambled, between ETS')
    ax[0].axhline(0.3, color='0.3', lw=0.5, ls=':'); ax[0].axhline(0.7, color='tab:green', lw=0.5, ls=':')
    ax[0].set_xlabel('z threshold'); ax[0].set_ylabel('per live year (floor 0.05)'); ax[0].legend(fontsize=7)
    ax[0].set_title(f'all pilot-2 headline rules, z_loso ≥ {zl:g}', fontsize=9)
    sc = ax[1].scatter(Gd.ab.clip(lower=0.01), Gd.rec_mw38, c=Gd.z_thr, cmap='viridis', s=12 + 4 * Gd.zl_thr)
    ax[1].set_xscale('log'); ax[1].axvline(1, color='k', lw=0.5, ls=':')
    ax[1].scatter([max(V['pair']['ab'], 0.01)], [V['pair']['rec_mw38']], marker='*', s=200, color='r', label='chosen pair')
    ax[1].set_xlabel('(a) + (b) per year (floor 0.01)'); ax[1].set_ylabel('recovery at Mw 3.8 (family injections)')
    fig.colorbar(sc, ax=ax[1], label='z_thr (size = z_loso thr)'); ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan2_separation.png', dpi=130); plt.close(fig)
    # veto validation
    fig, a_ = plt.subplots(figsize=(7, 5.5))
    p = pts[pts.kind == 'real 2011 clusters']
    sc = a_.scatter(p.ref_ratio, p.ref_excess, c=p.z.clip(5, 30), cmap='Greys', s=4, alpha=0.5, label='real 2011 clusters (z ≥ 5)')
    for k, col, mkr in (('injection Mw>=3.6 (family)', 'tab:green', 'o'), ('known body-wave', 'tab:red', 's'), ('H&G', 'tab:blue', '*')):
        q = pts[pts.kind == k]
        a_.scatter(q.ref_ratio, q.ref_excess, color=col, marker=mkr, s=60 if k == 'H&G' else 25, edgecolor='k', lw=0.3, label=k)
    a_.axhline(3, color='k', lw=0.6, ls='--'); a_.axvline(0.5, color='k', lw=0.6, ls='--')
    a_.set_xscale('log'); a_.set_yscale('symlog'); a_.set_xlabel('ref_ratio (median reference RMS / near RMS)')
    a_.set_ylabel('ref_excess (median, noise units)'); a_.legend(fontsize=7)
    a_.set_title(f'veto validation; chosen variant: {chosen}', fontsize=9); fig.colorbar(sc, ax=a_, label='z')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan2_veto.png', dpi=130); plt.close(fig)
    # injections
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for lk, col in (('family', 'tab:red'), ('off_family', 'tab:blue')):
        g = JD[(JD.gf == 'standard') & (JD.loc_kind == lk) & ~JD.in_v2_true].groupby('mw').recovered.mean()
        ax[0].plot(g.index, g.values, 'o-', color=col, label=lk)
    ax[0].axhline(0.8, color='k', lw=0.5, ls=':'); ax[0].set_xlabel('injected Mw'); ax[0].set_ylabel('recovered (all headline rules)')
    ax[0].legend(fontsize=7)
    r = JD[JD.recovered]
    sc = ax[1].scatter(r.mw + np.random.default_rng(1).uniform(-0.05, 0.05, len(r)), r.loc_err_km, c=r.z.clip(0, 30), s=14, cmap='viridis')
    ax[1].set_xlabel('injected Mw'); ax[1].set_ylabel('location error (km)'); fig.colorbar(sc, ax=ax[1], label='z')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan2_injection.png', dpi=130); plt.close(fig)
    record_sections()


def record_sections():
    """Raw band-passed Z (cleaned near stations + references) around the 08-18 cluster and the 08-22 triple."""
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    SC2.init(RUNS['real'], f'{S.SD}/bp_2011g', f'{S.SD}/bp2_ref')
    W = SC2.W
    fig, ax = plt.subplots(1, 2, figsize=(13, 8))
    for a_, (lab, t, w0, w1) in zip(ax, (('2011-08-18 21:19:46 UTC', CLASS_T[0][1], -600, 900),
                                          ('2011-08-22 08:21-08:28 UTC', CLASS_T[1][1], -600, 1200))):
        tt = ts(t); d, u = divmod(tt, 86400)
        tr = []
        for i, s in enumerate(W['stas']):
            x = SC2.cleaned_window(i, s, tt, tt, w0, w1)[2]
            dist = float(S.dist_km(49.24, -123.18, W['sll'][i, 0], W['sll'][i, 1]))
            tr.append((dist, s, x, 'k'))
        for ref in W['refs']:
            x = np.asarray(ref['X'][2, tt + w0:tt + w1], float) * 1e9
            tr.append((float(S.dist_km(49.24, -123.18, ref['lat'], ref['lon'])), ref['sta'], x, 'tab:blue'))
        tr.sort(key=lambda r: r[0])
        sc = np.nanmax([np.nanmax(np.abs(x)) for _, _, x, _ in tr if np.isfinite(x).any()] + [1])
        tax = np.arange(w0, w1)
        for k, (dist, s, x, col) in enumerate(tr):
            a_.plot(tax, x / sc * 2.5 - k, color=col, lw=0.6)
            a_.text(w0, -k + 0.3, f'{s} {dist:.0f} km  peak {np.nanmax(np.abs(x)) if np.isfinite(x).any() else np.nan:.0f} nm', fontsize=6)
        a_.axvline(0, color='r', lw=0.6)
        a_.set_title(f'{lab}: band-passed Z (0.02-0.05 Hz), common scale; black = pilot stations (cleaned),\n'
                     'blue = references; sorted by distance from 49.24 N -123.18', fontsize=8)
        a_.set_xlabel('s after the pilot-1 cluster time'); a_.set_yticks([])
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan2_classify_sections.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
