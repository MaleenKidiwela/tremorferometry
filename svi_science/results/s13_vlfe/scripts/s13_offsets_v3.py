#!/usr/bin/env python
"""S13 Amendment 6 steps 1-2 (PLAN.md, written 2026-10-14 ~8:50 pm PT before any number): WA per-family origin offsets
from the catalogue alone + a-priori location-suspect flags and duplicate clusters. Reads only the WA catalogue,
families.csv, tpl_library.csv, stations.csv, the CVM v1.7 NW profile and Slab2.

Writes tables/wa_family_offsets.csv, tables/wa3_pairs.csv, tables/wa3_location_flags.csv, tables/wa3_epoch_drift.csv,
tables/wa3_offsets_verdict.json.

Model: t_phys = time_f + d_f; a co-detection of the same event by A and B gives t_B - t_A = d_A - d_B.
offset_f = d_f + c (c from clear S picks: delta = s_off - 0.5 - r/Vs, physical origin relative to time + 0.5),
so to_f = time + 0.5 + offset_f."""
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_units_v2 import slab2_depth

WA = f'{C.ROOT}/lfe_extension_wa/results/wa_prod'
WIN = 12.0
NMIN = 8
FRAC_SUS, DIST_SUS, DEP_MIN = 0.30, 15.0, 20.0
VS = C.VS_KMS
EPOCHS = [('2007-01-01', '2011-12-31', '2007_2011'), ('2012-01-01', '2017-12-31', '2012_2017'),
          ('2018-01-01', '2024-12-31', '2018_2024')]
CVM = f'{C.ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/royer_NW_WA_S.crh'


def load_rows():
    c = pd.read_parquet(f'{WA}/post_v2/catalogue_v2.parquet',
                        columns=['det_id', 'time', 'date', 'family', 'region', 'P_own_pnet', 'selected_v2', 'net_mad'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    return pd.DataFrame(dict(det_id=c.det_id.values, fam=c.family.astype(int).values, P=c.P_own_pnet.values,
                             t=c.time.astype('datetime64[ns]').astype('int64').values / 1e9,
                             date=c.date.astype(str).values, net_mad=c.net_mad.values))


def families():
    F = pd.read_csv(f'{WA}/setup/families.csv').rename(columns={'fid': 'fam', 'lat_use': 'lat', 'lon_use': 'lon'})
    return F.set_index('fam')


def dist_km(la1, lo1, z1, la2, lo2, z2):
    y = (la2 - la1) * 111.19
    x = (lo2 - lo1) * 111.19 * np.cos(np.radians(0.5 * (la1 + la2)))
    return np.sqrt(x ** 2 + y ** 2 + (np.nan_to_num(z2) - np.nan_to_num(z1)) ** 2)


def cvm_slowness():
    L = []
    for ln in open(CVM).read().splitlines()[1:]:
        L.append((float(ln[:5]), float(ln[5:])))
    vs = np.array([a for a, b in L]); zt = np.array([b for a, b in L])

    def avg_slow(z):
        z = max(float(z), 1e-3)
        edges = np.r_[zt, np.inf]
        tt = 0.0
        for i in range(len(vs)):
            a, b = edges[i], min(edges[i + 1], z)
            if b > a:
                tt += (b - a) / vs[i]
        return tt / z
    return avg_slow


def co_detect(R, fams):
    """ordered pairs (A, B): each A row -> nearest B row within +-WIN. Returns DataFrame A, B, n, med, iqr, frac_A."""
    T = {f: np.sort(R.t.values[R.fam.values == f]) for f in fams}
    rows = []
    for A in fams:
        ta = T[A]
        for B in fams:
            if A == B:
                continue
            tb = T[B]
            if len(tb) == 0:
                continue
            k = np.searchsorted(tb, ta)
            lo = np.clip(k - 1, 0, len(tb) - 1); hi = np.clip(k, 0, len(tb) - 1)
            dl = tb[lo] - ta; dh = tb[hi] - ta
            d = np.where(np.abs(dl) <= np.abs(dh), dl, dh)
            m = np.abs(d) <= WIN
            if m.sum() == 0:
                continue
            x = d[m]
            rows.append(dict(A=A, B=B, n=int(m.sum()), med=float(np.median(x)),
                             iqr=float(np.subtract(*np.percentile(x, [75, 25]))), frac_A=float(m.sum() / len(ta)),
                             nA=len(ta), nB=len(tb)))
    return pd.DataFrame(rows)


def components(nodes, edges):
    par = {n: n for n in nodes}

    def f(x):
        while par[x] != x:
            par[x] = par[par[x]]; x = par[x]
        return x
    for a, b in edges:
        ra, rb = f(a), f(b)
        if ra != rb:
            par[ra] = rb
    out = {}
    for n in nodes:
        out.setdefault(f(n), []).append(n)
    return sorted(out.values(), key=len, reverse=True)


def solve(P):
    """P: unordered pairs with columns A, B (A has fewer rows), n, dAB (= d_A - d_B). WLS weight sqrt(n) per pair,
    per connected component (gauge: mean d = 0). Returns d (Series), comp (Series), residual table."""
    nodes = sorted(set(P.A) | set(P.B))
    comps = components(nodes, list(zip(P.A, P.B)))
    d, cid, res = {}, {}, []
    for ci, cn in enumerate(comps):
        idx = {f: i for i, f in enumerate(cn)}
        Q = P[P.A.isin(idx) & P.B.isin(idx)]
        if len(cn) < 2 or len(Q) == 0:
            continue
        G = np.zeros((len(Q) + 1, len(cn))); y = np.zeros(len(Q) + 1); w = np.sqrt(np.sqrt(Q.n.values.astype(float)))
        for k, r in enumerate(Q.itertuples()):
            G[k, idx[r.A]] = w[k]; G[k, idx[r.B]] = -w[k]; y[k] = w[k] * r.dAB
        G[-1, :] = 1.0
        x, *_ = np.linalg.lstsq(G, y, rcond=None)
        for f, i in idx.items():
            d[f] = x[i]; cid[f] = ci
        pred = x[[idx[a] for a in Q.A]] - x[[idx[b] for b in Q.B]]
        res.append(Q.assign(comp=ci, pred=pred, resid=Q.dAB.values - pred))
    return pd.Series(d), pd.Series(cid), (pd.concat(res) if res else pd.DataFrame())


def unordered(D, sel=None):
    """one row per unordered pair from the ordered table: direction = the family with fewer rows matched into the other."""
    X = D[D.nA <= D.nB].copy()
    X = X[~((X.nA == X.nB) & (X.A > X.B))]
    X['dAB'] = X.med               # med = median(t_B - t_A) = d_A - d_B
    X = X[X.n >= NMIN]
    if sel is not None:
        X = X[sel(X)]
    return X


def anchor(d, comp, F, delta, nclear, suspect):
    """per component c = median over anchor families of (delta_f - d_f)."""
    ok = (F.located.astype(bool) & F.location_quality.isin(['A', 'B']) & (F.depth_km >= DEP_MIN))
    anc = [f for f in d.index if ok.get(f, False) and nclear.get(f, 0) >= 1 and not suspect.get(f, False)
           and np.isfinite(delta.get(f, np.nan))]
    cc = {}
    for ci in sorted(set(comp)):
        a = [f for f in anc if comp[f] == ci]
        if a:
            cc[ci] = float(np.median([delta[f] - d[f] for f in a]))
    return cc, anc


def main():
    R = load_rows()
    F = families()
    fams = sorted(set(R.fam) & set(F.index))
    R = R[R.fam.isin(fams)]
    nrow = R.groupby('fam').size()
    v = dict(n_rows=int(len(R)), n_fams=len(fams))
    # ---- co-detections (raw catalogue times)
    D = co_detect(R, fams)
    la, lo, de = F.lat, F.lon, F.depth_km
    D['dist_km'] = dist_km(la.reindex(D.A).values, lo.reindex(D.A).values, de.reindex(D.A).values,
                           la.reindex(D.B).values, lo.reindex(D.B).values, de.reindex(D.B).values)
    D.to_csv(f'{C.TB}/wa3_pairs.csv', index=False)
    # ---- step 2: suspect flags + duplicate clusters
    sus_i = D[(D.frac_A >= FRAC_SUS) & (D.dist_km > DIST_SUS)].groupby('A').apply(
        lambda g: ';'.join(f'{int(b)}:{fr:.2f}@{dk:.0f}km' for b, fr, dk in zip(g.B, g.frac_A, g.dist_km)))
    L = pd.DataFrame(index=pd.Index(fams, name='fam'))
    L['n_rows'] = nrow
    L['family_uid'] = F.family_uid.reindex(fams)
    L['lat'], L['lon'], L['depth_cat'] = la.reindex(fams), lo.reindex(fams), de.reindex(fams)
    L['location_quality'] = F.location_quality.reindex(fams)
    L['sus_codetect'] = sus_i.reindex(fams).fillna('')
    L['sus_depth'] = L.depth_cat.isna() | (L.depth_cat < DEP_MIN)
    L['suspect'] = (L.sus_codetect != '') | L.sus_depth
    E = D[D.frac_A >= FRAC_SUS]
    clus = components(fams, list(zip(E.A, E.B)))
    L['cluster'] = -1; L['cluster_size'] = 1; L['cluster_span_km'] = 0.0
    L['gf_lat'], L['gf_lon'], L['gf_dep'], L['gf_source'] = L.lat, L.lon, L.depth_cat, 'catalogue'
    slab = pd.Series(slab2_depth(L.lat.values, L.lon.values), index=L.index)
    L['slab2_km'] = slab
    qrank = {'A': 0, 'B': 1, 'C': 2}
    for ci, cn in enumerate([c for c in clus if len(c) > 1]):
        M = L.loc[cn]
        span = max(dist_km(M.lat[a], M.lon[a], M.depth_cat[a], M.lat[b], M.lon[b], M.depth_cat[b]) for a in cn for b in cn)
        L.loc[cn, 'cluster'] = ci; L.loc[cn, 'cluster_size'] = len(cn); L.loc[cn, 'cluster_span_km'] = span
        good = M[~M.sus_depth]          # Addendum 6a (b): flag (i) is symmetric, so only the depth flag excludes
        if len(good):
            best = good.assign(q=good.location_quality.map(qrank).fillna(3)).sort_values(['q', 'n_rows'],
                                                                                        ascending=[True, False]).index[0]
            L.loc[cn, 'gf_lat'] = L.lat[best]; L.loc[cn, 'gf_lon'] = L.lon[best]; L.loc[cn, 'gf_dep'] = L.depth_cat[best]
            L.loc[cn, 'gf_source'] = f'cluster_best_{best}'
            L.loc[best, 'gf_source'] = 'catalogue(cluster_best)'
    # suspect-by-depth families keeping their own location -> Slab2 depth
    own = L.gf_source.isin(['catalogue']) & L.sus_depth
    L.loc[own, 'gf_dep'] = L.loc[own, 'slab2_km']; L.loc[own, 'gf_source'] = 'catalogue+slab2_depth'
    # ---- clear-pick deltas
    T = pd.read_csv(f'{WA}/setup/tpl_library.csv')
    S = pd.read_csv(f'{WA}/setup/stations.csv').set_index('netsta')
    Tc = T[T.clear & T.fid.isin(F.index)].copy()
    Tc = Tc[Tc.netsta.isin(S.index)]
    slow = cvm_slowness()
    r = dist_km(F.lat.reindex(Tc.fid).values, F.lon.reindex(Tc.fid).values, F.depth_km.reindex(Tc.fid).values,
                S.lat.reindex(Tc.netsta).values, S.lon.reindex(Tc.netsta).values, np.zeros(len(Tc)))
    Tc['r_km'] = r
    Tc['delta'] = Tc.s_off - 0.5 - r / VS
    Tc['delta_cvm'] = Tc.s_off - 0.5 - r * np.array([slow(z) for z in F.depth_km.reindex(Tc.fid).fillna(35).values])
    Tc = Tc[np.isfinite(Tc.delta)]
    delta = Tc.groupby('fid').delta.median(); delta_cvm = Tc.groupby('fid').delta_cvm.median()
    nclear = Tc.groupby('fid').size()
    Tc.to_csv(f'{C.TB}/wa3_clear_picks.csv', index=False)
    v['clear_picks_used'] = int(len(Tc)); v['clear_fams'] = int(Tc.fid.nunique())
    v['clear_fams_NW'] = int(Tc[Tc.fid.isin(fams)].fid.nunique())
    sus = L.suspect.to_dict()

    def network(Pu, label):
        d, comp, res = solve(Pu)
        cc, anc = anchor(d, comp, F, delta.to_dict(), nclear.to_dict(), sus)
        out = pd.DataFrame(dict(d_f=d, comp=comp))
        out['c'] = out.comp.map(cc)
        out['offset_net'] = out.d_f + out.c
        rms = float(np.sqrt(np.mean(res.resid ** 2))) if len(res) else np.nan
        wrms = float(np.sqrt(np.sum(np.sqrt(res.n) * res.resid ** 2) / np.sum(np.sqrt(res.n)))) if len(res) else np.nan
        a = [f for f in anc if np.isfinite(out.offset_net.get(f, np.nan))]
        dev = np.array([out.offset_net[f] - delta[f] for f in a])
        allc = [f for f in out.index if np.isfinite(out.offset_net[f]) and f in delta.index]
        dev_all = np.array([out.offset_net[f] - delta[f] for f in allc])
        info = dict(label=label, n_pairs=int(len(Pu)), n_fams_net=int(len(d)), n_comp=int(comp.nunique()),
                    comp_sizes=[int(x) for x in comp.value_counts().values], anchored_comps={int(k): round(v_, 2) for k, v_ in cc.items()},
                    n_anchor=len(a), resid_rms=rms, resid_wrms=wrms, resid_absmax=float(res.resid.abs().max()) if len(res) else np.nan,
                    anchor_dev_sd=float(np.std(dev, ddof=1)) if len(dev) > 1 else np.nan,
                    anchor_dev_median=float(np.median(dev)) if len(dev) else np.nan,
                    allclear_dev_sd=float(np.std(dev_all, ddof=1)) if len(dev_all) > 1 else np.nan, n_allclear=len(dev_all))
        # leave-one-out (Addendum 6a): anchors in components with >= 2 anchors, constant from the other anchors
        loo = []
        for f in a:
            oth = [g for g in a if g != f and out.comp[g] == out.comp[f]]
            if oth:
                cl = float(np.median([delta[g] - out.d_f[g] for g in oth]))
                loo.append(dict(fam=f, comp=int(out.comp[f]), dev=float(out.d_f[f] + cl - delta[f])))
        LO = pd.DataFrame(loo)
        info['loo_n'] = int(len(LO)); info['loo_dev_sd'] = float(np.sqrt(np.mean(LO.dev ** 2))) if len(LO) > 1 else np.nan
        info['loo_dev_median'] = float(LO.dev.median()) if len(LO) else np.nan
        info['loo_devs'] = {int(r.fam): round(r.dev, 2) for r in LO.itertuples()} if len(LO) else {}
        info['accept'] = bool(rms < 2 and info['loo_dev_sd'] < 3)
        info['fail_stop'] = bool(info['loo_dev_sd'] > 4)
        return out, res, info, a

    # Addendum 6a (a): primary network from all NW rows with P_own_pnet >= 0.5, pairs n >= 8 and IQR <= 0.5 s
    c = pd.read_parquet(f'{WA}/post_v2/catalogue_v2.parquet', columns=['time', 'family', 'region', 'P_own_pnet'])
    c = c[(c.region == 'NW') & (c.P_own_pnet >= 0.5)]
    Rw = pd.DataFrame(dict(fam=c.family.astype(int).values, t=c.time.astype('datetime64[ns]').astype('int64').values / 1e9))
    Rw = Rw[Rw.fam.isin(fams)]
    Dw = co_detect(Rw, fams)
    Dw.to_csv(f'{C.TB}/wa3_pairs_p05.csv', index=False)
    v['n_rows_p05'] = int(len(Rw))
    Pall = unordered(D)
    net_p, res_p, info_p, anc_p = network(unordered(Dw, lambda X: X.iqr <= 0.5), 'primary_p05_iqr05')
    net_s, res_s, info_s, _ = network(Pall, 'selected_only')
    v['net_primary'] = info_p; v['net_selected_only'] = info_s
    j = net_p.offset_net.dropna().index.intersection(net_s.offset_net.dropna().index)
    dd_ = (net_p.offset_net[j] - net_s.offset_net[j])
    v['selected_vs_primary'] = dict(n=int(len(j)), median_diff=float(dd_.median()) if len(j) else np.nan,
                                    absmax_diff=float(dd_.abs().max()) if len(j) else np.nan,
                                    n_diff_gt1s=int((dd_.abs() > 1).sum()))
    use_variant = False
    v['used'] = 'primary_p05_iqr05'
    net, res = net_p, res_p
    net_s.rename(columns=lambda x: x + '_selonly').to_csv(f'{C.TB}/wa3_network_selected_only.csv')
    res.to_csv(f'{C.TB}/wa3_pair_residuals.csv', index=False)
    # ---- epoch test fam38-fam75 (+ 10 largest pairs, column)
    ep_rows = []
    big = Pall.sort_values('n', ascending=False).head(10)[['A', 'B']].values.tolist()
    if not any({a, b} == {38, 75} for a, b in big):
        big = [[38, 75]] + big
    for A, B in big:
        rec = dict(A=int(A), B=int(B))
        meds = []
        for a0, a1, lab in EPOCHS:
            Re = R[(R.date >= a0) & (R.date <= a1)]
            ta = np.sort(Re.t.values[Re.fam.values == A]); tb = np.sort(Re.t.values[Re.fam.values == B])
            if len(ta) == 0 or len(tb) == 0:
                rec[f'n_{lab}'] = 0; rec[f'med_{lab}'] = np.nan; continue
            k = np.searchsorted(tb, ta); lo_ = np.clip(k - 1, 0, len(tb) - 1); hi_ = np.clip(k, 0, len(tb) - 1)
            dd = np.where(np.abs(tb[lo_] - ta) <= np.abs(tb[hi_] - ta), tb[lo_] - ta, tb[hi_] - ta)
            m = np.abs(dd) <= WIN
            rec[f'n_{lab}'] = int(m.sum()); rec[f'med_{lab}'] = float(np.median(dd[m])) if m.sum() >= NMIN else np.nan
            if m.sum() >= NMIN:
                meds.append(rec[f'med_{lab}'])
        rec['drift'] = float(max(meds) - min(meds)) if len(meds) >= 2 else np.nan
        ep_rows.append(rec)
    EP = pd.DataFrame(ep_rows); EP.to_csv(f'{C.TB}/wa3_epoch_drift.csv', index=False)
    d3875 = EP[((EP.A == 38) & (EP.B == 75)) | ((EP.A == 75) & (EP.B == 38))]
    v['epoch_38_75'] = d3875.iloc[0].to_dict() if len(d3875) else None
    v['epoch_drift_38_75'] = float(d3875.drift.iloc[0]) if len(d3875) else np.nan
    v['per_family_epoch'] = bool(v['epoch_drift_38_75'] > 1.5)
    # ---- offsets table
    O = L[['n_rows', 'family_uid', 'location_quality', 'depth_cat', 'suspect']].copy()
    O = O.join(net[['comp', 'd_f', 'c', 'offset_net']])
    O['delta_f'] = delta.reindex(O.index); O['delta_cvm'] = delta_cvm.reindex(O.index)
    O['n_clear'] = nclear.reindex(O.index).fillna(0).astype(int)
    O['anchor'] = O.index.isin(anc_p if not use_variant else anc_p)
    netmed = float(O.offset_net.dropna().median())
    O['flag'] = np.where(O.offset_net.notna(), 'network',
                         np.where(O.delta_f.notna() & ~O.suspect, 'clear_only', 'default'))
    O['offset_f'] = np.where(O.flag == 'network', O.offset_net, np.where(O.flag == 'clear_only', O.delta_f, netmed))
    # CVM column: same relative offsets, constant from the CVM deltas of the same anchor families
    an = O[O.anchor & O.offset_net.notna()]
    cvm_shift = float(np.median(an.delta_cvm - an.delta_f)) if len(an) else np.nan
    O['offset_cvm'] = O.offset_f + cvm_shift
    v['cvm_constant_shift_s'] = cvm_shift
    v['default_offset'] = netmed
    v['flag_counts'] = O.flag.value_counts().to_dict()
    v['offset_sd_network'] = float(O[O.flag == 'network'].offset_f.std())
    v['offset_range_network'] = [float(O[O.flag == 'network'].offset_f.min()), float(O[O.flag == 'network'].offset_f.max())]
    v['rows_covered_network'] = float(O[O.flag == 'network'].n_rows.sum() / O.n_rows.sum())
    # pair checks quoted by Merlin
    for a, b in ((38, 75), (26, 117), (58, 77)):
        x = D[(D.A == a) & (D.B == b)]
        if len(x):
            v[f'pair_{a}_{b}'] = dict(n=int(x.n.iloc[0]), med=float(x.med.iloc[0]), iqr=float(x.iqr.iloc[0]),
                                      dist_km=float(x.dist_km.iloc[0]))
    # per-family epoch offsets (only if the drift rule fires)
    if v['per_family_epoch']:
        for a0, a1, lab in EPOCHS:
            Re = R[(R.date >= a0) & (R.date <= a1)]
            De = co_detect(Re, sorted(Re.fam.unique()))
            ne, _, ie, _ = network(unordered(De), lab)
            v[f'net_{lab}'] = ie
            O[f'offset_{lab}'] = ne.offset_net.reindex(O.index)
            O[f'offset_{lab}'] = O[f'offset_{lab}'].fillna(O.offset_f)
    O.to_csv(f'{C.TB}/wa_family_offsets.csv')
    L.to_csv(f'{C.TB}/wa3_location_flags.csv')
    v['suspect_n'] = int(L.suspect.sum()); v['suspect_codetect'] = L.index[L.sus_codetect != ''].tolist()
    v['suspect_depth_n'] = int(L.sus_depth.sum())
    v['clusters'] = {int(ci): dict(members=g.index.tolist(), span_km=round(float(g.cluster_span_km.iloc[0]), 1),
                                   gf=g.gf_source.iloc[0]) for ci, g in L[L.cluster >= 0].groupby('cluster')}
    inf = v['net_primary']
    v['step1_accept'] = inf['accept']; v['step1_fail_stop'] = inf['fail_stop']
    json.dump(v, open(f'{C.TB}/wa3_offsets_verdict.json', 'w'), indent=1,
              default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(v, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))


if __name__ == '__main__':
    main()
