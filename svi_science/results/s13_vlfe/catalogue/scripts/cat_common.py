"""Catalogue stage-1 shared code (rules: ../PLAN_CAT.md). Pass rules a-f, e', S5, families/footprint, eq-window
helpers. Read-only use of hg_repro and s13_vlfe scripts/data."""
import os, sys
import numpy as np, pandas as pd

ROOT = '/home/jovyan/tremorferometry'
S13 = f'{ROOT}/svi_science/results/s13_vlfe'
HGR = f'{S13}/hg_repro'
CAT = f'{S13}/catalogue'
TB, FIG, DATA, LOGD = f'{CAT}/tables', f'{CAT}/figures', f'{CAT}/data', f'{CAT}/logs'
for p in (f'{HGR}/scripts', f'{S13}/scripts'):
    if p not in sys.path:
        sys.path.insert(0, p)
FAR5 = ['COR', 'NEW', 'HAWA', 'OMAK', 'I02D', 'MSO', 'WVOR']   # Q0: MSO/WVOR added (absent in 2011-12 -> no change there)
REF_ALL = ['GNW'] + FAR5
FAMREG = ['S_Vancouver_Is', 'Royer_NW_WA', 'Royer_SW_WA', 'Royer_N_Vancouver_Is', 'Olympic_2010', 'Olympic_2010_armB']
PERIODS = {'ETS': [('2011-07-23', '2011-09-06'), ('2012-08-30', '2012-10-11')], 'between': [('2011-09-08', '2012-08-28')]}
BETWEEN_YR = (pd.Timestamp('2012-08-28') - pd.Timestamp('2011-09-08')).days / 365.25


def log(msg, name='cat'):
    line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT {name} {msg}'
    print(line, flush=True)
    open(f'{LOGD}/{name}.log', 'a').write(line + '\n')


def period_of(epoch):
    t = pd.to_datetime(np.asarray(epoch, float), unit='s')
    out = np.full(len(t), 'other', dtype=object)
    for name, wins in PERIODS.items():
        for a, b in wins:
            out[(t >= pd.Timestamp(a)) & (t < pd.Timestamp(b))] = name
    return out


def families():
    F = pd.read_csv(f'{ROOT}/data/catalogue/families_v1_3.csv', low_memory=False)
    F = F[F.region.isin(FAMREG)].copy()
    for a, b in (('pos_lat', 'pos_lon'), ('approx_lat', 'approx_lon')):
        m = F.lat.isna() & F[a].notna()
        F.loc[m, 'lat'] = F.loc[m, a]; F.loc[m, 'lon'] = F.loc[m, b]
    return F[F.lat.notna()][['region', 'family_id', 'lat', 'lon']].reset_index(drop=True)


def dkm(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(la1) - la2) * 111.19, (np.asarray(lo1) - lo2) * 111.19 * np.cos(np.radians(48.3)))


def node_table():
    """per scan grid node: lat, lon, dep, grid_edge, fp_interior, in_footprint_poly, dist_family_km, nearest_family."""
    f = f'{DATA}/node_table.parquet'
    if os.path.isfile(f):
        return pd.read_parquet(f)
    from matplotlib.path import Path
    TP = np.load(f'{S13}/data/scan/templates_grid.npz')
    N = pd.DataFrame(dict(lat=TP['lat'], lon=TP['lon'], dep=TP['dep']))
    import scan_common as S
    la, lo = N.lat.values[:-1], N.lon.values[:-1]
    x, y = S.kmxy(la, lo)
    edge = np.zeros(len(N), bool)
    for i in range(len(la)):
        nb = 0
        for dx, dy in ((10, 0), (-10, 0), (0, 10), (0, -10)):
            nb += np.any((np.abs(x - x[i] - dx) < 1) & (np.abs(y - y[i] - dy) < 1))
        edge[i] = nb < 4
    N['grid_edge'] = edge
    P = pd.read_csv(f'{S13}/tables/scan_footprint_polygon.csv')
    poly = Path(np.c_[P.lon, P.lat])
    N['fp_interior'] = poly.contains_points(np.c_[N.lon, N.lat])
    F = families()
    D = dkm(N.lat.values[:, None], N.lon.values[:, None], F.lat.values[None, :], F.lon.values[None, :])
    k = np.argmin(D, 1)
    N['dist_family_km'] = D[np.arange(len(N)), k]
    N['nearest_family'] = [f'{F.region.values[j]}:{F.family_id.values[j]}' for j in k]
    N.to_parquet(f)
    return N


def eq_table(year):
    name = {2011: 'scan2_eq_events.csv', 2012: 'scan3_eq_events.csv'}.get(year, f'scan{year}_eq_events.csv')
    return pd.read_csv(f'{S13}/tables/{name}')


def head2_maxmag(epochs, years):
    """max ComCat mag of head2 windows overlapping [t0-60, t0+120] (NaN if none)."""
    epochs = np.asarray(epochs, float); years = np.asarray(years)
    out = np.full(len(epochs), np.nan)
    for y in np.unique(years):
        A = eq_table(int(y))
        A = A[A.head2_t0.notna()]
        a0, a1, mg = A.head2_t0.values.astype(float), A.head2_t1.values.astype(float), A.mag.values.astype(float)
        o = np.argsort(a0); a0, a1, mg = a0[o], a1[o], mg[o]
        maxlen = np.max(a1 - a0) if len(a0) else 0
        for k in np.where(years == y)[0]:
            ep = epochs[k]
            lo = np.searchsorted(a0, ep - 60 - maxlen); hi = np.searchsorted(a0, ep + 120)
            sel = slice(lo, hi)
            hit = (a1[sel] > ep - 60) & (a0[sel] < ep + 120)
            if hit.any():
                out[k] = float(mg[sel][hit].max())
    return out


def f_rule(z_loso, z_best):
    """Amendment Y14 (Y11 adopt rule met; Y12 condition met): f' = (z_loso >= 3 and z_loso >= 0.5 z_best) or z_loso >= 5."""
    zl = pd.to_numeric(pd.Series(z_loso), errors='coerce').values; zb = pd.to_numeric(pd.Series(z_best), errors='coerce').values
    return ((zl >= 3) & (zl >= 0.5 * zb)) | (zl >= 5)


def rescore_inj(J):
    """Amendment Y14: injections re-scored from their stored z_loso / z_best (no re-run); f (pre-Y14) kept as columns."""
    if 'z_loso' not in J or 'pass_f' not in J:
        return J
    J = J.copy()
    J['pass_f_preY14'] = J.pass_f; J['battery_e2_preY14'] = J.battery_e2
    J['pass_f'] = f_rule(J.z_loso, J.z_best)
    J['battery_e2'] = (J.pass_a.fillna(False).astype(bool) & J.pass_b.fillna(False).astype(bool) & J.pass_c.fillna(False).astype(bool)
                       & J.pass_e2.fillna(False).astype(bool) & J.pass_f)
    return J


def apply_rules(D):
    """D needs: n_far, ratio_far, travelling, z_best, grid_edge, fp_interior/in_footprint_poly, dist_family_km,
    in_head2, in_bw5, in_reg, in_reg3, head2_maxmag, ex_<FAR5>, z_loso, and S5 inputs (optional)."""
    D = D.copy()
    D['a_testable'] = D.n_far >= 2
    D['pass_a'] = D.a_testable & (D.ratio_far <= 0.551)
    D['pass_a6'] = (D.n6 >= 2) & (D.ratio6 <= 0.551) if 'n6' in D else np.nan
    D['pass_b'] = ~D.travelling.fillna(False).astype(bool)
    D['locatable'] = D.z_best >= 5
    fp = D['in_footprint_poly'] if 'in_footprint_poly' in D else D.fp_interior
    D['pass_c'] = D.locatable & ~D.grid_edge.fillna(True).astype(bool) & (fp.fillna(False).astype(bool) | (D.dist_family_km <= 20))
    for c in ('in_head2', 'in_bw5', 'in_reg', 'in_reg3'):
        D[c] = D[c].fillna(True).astype(bool)
    D['pass_e'] = ~(D.in_head2 | D.in_bw5 | D.in_reg | D.in_reg3)
    exc = np.zeros(len(D), int)
    for t in FAR5:
        if f'ex_{t}' in D:
            exc += (D[f'ex_{t}'] >= 3).fillna(False).astype(int).values
    D['n_far_excited'] = exc
    D['e2_time_fail'] = D.in_bw5 | D.in_reg | D.in_reg3 | (D.in_head2 & (D.head2_maxmag >= 5.5).fillna(False))
    D['pass_e2'] = ~(D.e2_time_fail | (D.in_head2 & (D.n_far_excited >= 2)))
    D['pass_f_preY14'] = (D.z_loso >= 3) & (D.z_loso >= 0.5 * D.z_best)            # f (pre-Y14), kept as a column
    D['pass_f'] = f_rule(D.z_loso, D.z_best)                                        # Amendment Y14: f'
    D['battery'] = D.pass_a & D.pass_b & D.pass_c & D.pass_e & D.pass_f
    D['battery_e2'] = D.pass_a & D.pass_b & D.pass_c & D.pass_e2 & D.pass_f
    D['battery_e2_preY14'] = D.pass_a & D.pass_b & D.pass_c & D.pass_e2 & D.pass_f_preY14
    if 'n_Z_agree' in D:
        vr_ratio = D.vr_thrust / D.vr_free
        comp = np.fmax(D.vr_strikeslip.astype(float), D.vr_steepdip.astype(float))
        D['pass_S5'] = ((D.n_sta >= 3) & (D.n_Z_agree >= 2) & (vr_ratio >= 0.8) & (D.vr_thrust > comp)
                        & (D.isolation < 2) & (D.share_excess <= 0.25)).fillna(False)
    return D
