"""S5 setup, step 2: solid-earth (body) tide stress on the plate interface at each SVI family.

Method (primary, "IERS-FD"):
  * pyTMD.compute.SET_displacements (IERS 2010 conventions: degree 2+3, frequency-dependent
    diurnal/long-period Love numbers, latitude dependence, Montenbruck ephemerides) gives the
    N/E/R surface displacement at the family epicentre and at +-DLAT, +-DLON.
  * Surface strain tensor from first derivatives on a sphere of radius A:
        e_tt = (du_t/dt + u_r)/A,   e_ll = (du_l/dl)/(A sin t) + (u_t cot t + u_r)/A,
        e_tl = [ (du_t/dl)/sin t + du_l/dt - u_l cot t ] / (2A)
    (t = colatitude, l = longitude, u_t = pyTMD 'N' (which is southward, see note), u_l = u_E); ENU: e_ee=e_ll, e_nn=e_tt, e_en=-e_tl.
  * Strain assumed unchanged to 30-45 km depth (depth/A ~ 0.006); stress by PLANE STRESS
    (sigma_zz = sigma_xz = sigma_yz = 0, free-surface condition) with homogeneous moduli
    MU, NU taken from the USGS CVM v1.7 south_VI profile averaged over 29.5-42 km.
  * Traction on the Slab2 plane: sigma_n = n.S.n (tension positive), tau = s.S.n (positive =
    promotes hanging-wall slip in s; s = 236 deg convergence projection, also up-dip & strike).
Check ("analytic"): independent Love-number formulation with nominal h2,l2,h3,l3 from the
lunar/solar degree-2/3 potential, strains by finite differences of W; compared for family 0.
Output: results/s5_setup/body_stress_hourly.parquet  (fid, time, tau/sigma_n components, Pa)
"""
import numpy as np, pandas as pd, time as _time
from pathlib import Path
import pyTMD.compute as C
import pyTMD.astro as AST

ROOT = Path('/home/jovyan/tremorferometry')
OUT = ROOT / 'svi_science/results/s5_setup'
A = 6371.0e3
D = 0.02  # deg stencil half-width (~2 km)

# ---- elastic moduli from USGS CVM v1.7 south_VI, 29.5-42 km --------------------------
prof = pd.read_csv(ROOT / 'data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv')
sel = prof[(prof.ztop_km >= 29.5) & (prof.zbot_km <= 42.0)]
w = (sel.zbot_km - sel.ztop_km).values
vp = np.average(sel.vp, weights=w); vs = np.average(sel.vs, weights=w)
rho = 1.6612 * vp - 0.4721 * vp**2 + 0.0671 * vp**3 - 0.0043 * vp**4 + 0.000106 * vp**5  # Brocher 2005, g/cc
MU = rho * 1e3 * (vs * 1e3) ** 2
NU = (vp**2 - 2 * vs**2) / (2 * (vp**2 - vs**2))


def stress_from_strain(ee, nn, en):
    f = 2 * MU / (1 - NU)
    return f * (ee + NU * nn), f * (nn + NU * ee), 2 * MU * en  # S_ee, S_nn, S_en


def project(See, Snn, Sen, n, s):
    # S has only horizontal components (plane stress)
    Sn_e = See * n[0] + Sen * n[1]
    Sn_n = Sen * n[0] + Snn * n[1]
    return Sn_e * s[0] + Sn_n * s[1]


def surface_strain(uN, uE, uR, lat, dlat, dlon):
    """uX arrays shape (5, nt): [centre, +lat, -lat, +lon, -lon]."""
    th = np.radians(90 - lat)
    dth = np.radians(dlat); dl = np.radians(dlon)
    # NOTE pyTMD 3.0.9 SET_displacements 'N' is projected on the +colatitude unit vector
    # (cos l cos t, sin l cos t, -sin t), i.e. it is SOUTHWARD; verified against the analytic
    # Love-number displacement (corr -0.99995 with -u_theta). So u_theta = +'N'.
    ut = uN; ul = uE
    dut_dth = (ut[2] - ut[1]) / (2 * dth)          # +colat = -lat
    dul_dth = (ul[2] - ul[1]) / (2 * dth)
    dut_dl = (ut[3] - ut[4]) / (2 * dl)
    dul_dl = (ul[3] - ul[4]) / (2 * dl)
    ett = (dut_dth + uR[0]) / A
    ell = dul_dl / (A * np.sin(th)) + (ut[0] / np.tan(th) + uR[0]) / A
    etl = (dut_dl / np.sin(th) + dul_dth - ul[0] / np.tan(th)) / (2 * A)
    return ell, ett, -etl  # e_ee, e_nn, e_en


def analytic_strain(lat, lon, times):
    """Nominal-Love-number body-tide strain from the degree-2/3 potential (check only)."""
    G_M_moon = 4.9028e12; G_M_sun = 1.32712440018e20; g = 9.80665
    h = {2: 0.6078, 3: 0.292}; l = {2: 0.0847, 3: 0.015}
    mjd = (times - np.datetime64('1858-11-17')) / np.timedelta64(1, 'D')
    bodies = [(np.array(AST.lunar_ecef(mjd)), G_M_moon), (np.array(AST.solar_ecef(mjd)), G_M_sun)]
    dd = np.radians(D)
    th0, la0 = np.radians(90 - lat), np.radians(lon)
    offs = {'c': (0, 0), 'tp': (dd, 0), 'tm': (-dd, 0), 'lp': (0, dd), 'lm': (0, -dd),
            'pp': (dd, dd), 'pm': (dd, -dd), 'mp': (-dd, dd), 'mm': (-dd, -dd)}
    out = {}
    for deg in (2, 3):
        Wd = {}
        for k, (a1, a2) in offs.items():
            th, la = th0 + a1, la0 + a2
            r = np.array([np.sin(th) * np.cos(la), np.sin(th) * np.sin(la), np.cos(th)])
            W = 0
            for XYZ, GM in bodies:
                R = np.linalg.norm(XYZ, axis=0)
                mu = (r[:, None] * XYZ).sum(0) / R
                P = 0.5 * (3 * mu**2 - 1) if deg == 2 else 0.5 * (5 * mu**3 - 3 * mu)
                W = W + GM * A**deg / R ** (deg + 1) * P
            Wd[k] = W
        Wt = (Wd['tp'] - Wd['tm']) / (2 * dd); Wl = (Wd['lp'] - Wd['lm']) / (2 * dd)
        Wtt = (Wd['tp'] - 2 * Wd['c'] + Wd['tm']) / dd**2
        Wll = (Wd['lp'] - 2 * Wd['c'] + Wd['lm']) / dd**2
        Wtl = (Wd['pp'] - Wd['pm'] - Wd['mp'] + Wd['mm']) / (4 * dd**2)
        ett = (h[deg] * Wd['c'] + l[deg] * Wtt) / (g * A)
        ell = (h[deg] * Wd['c'] + l[deg] * (Wt / np.tan(th0) + Wll / np.sin(th0) ** 2)) / (g * A)
        etl = l[deg] * (Wtl - Wl / np.tan(th0)) / (g * A * np.sin(th0))
        out[deg] = (ell, ett, -etl)
    return tuple(out[2][i] + out[3][i] for i in range(3)), h[2] * 0  # (ee,nn,en)


def main():
    geo = pd.read_csv(OUT / 'family_fault_geometry.csv')
    t_all = np.arange('2014-01-01', '2026-10-01', np.timedelta64(1, 'h'), dtype='datetime64[ns]')
    print(f'MU={MU/1e9:.1f} GPa NU={NU:.3f} rho={rho:.2f} vp={vp:.2f} vs={vs:.2f}; {len(t_all)} hours', flush=True)
    lat = geo.lat.values; lon = geo.lon.values; nf = len(geo)
    # stencil: [c, +lat, -lat, +lon, -lon]
    slat = np.concatenate([lat, lat + D, lat - D, lat, lat])
    slon = np.concatenate([lon, lon, lon, lon + D, lon - D])
    parts = []
    years = pd.date_range('2014-01-01', '2027-01-01', freq='YS').values
    for y0, y1 in zip(years[:-1], years[1:]):
        t = t_all[(t_all >= y0) & (t_all < y1)]
        if len(t) == 0:
            continue
        t0 = _time.time()
        ds = C.SET_displacements(slon, slat, t, type='time series', standard='datetime',
                                 variable=['N', 'E', 'R'])
        uN = ds['N'].values.reshape(5, nf, -1); uE = ds['E'].values.reshape(5, nf, -1)
        uR = ds['R'].values.reshape(5, nf, -1)
        for i in range(nf):
            ee, nn, en = surface_strain(uN[:, i], uE[:, i], uR[:, i], lat[i], D, D)
            See, Snn, Sen = stress_from_strain(ee, nn, en)
            r = geo.iloc[i]
            n = np.array([r.n_e, r.n_n, r.n_u])
            df = pd.DataFrame({
                'fid': np.int16(r.fid), 'time': t,
                'body_sigma_n_Pa': project(See, Snn, Sen, n, n).astype(np.float32),
                'body_tau_Pa': project(See, Snn, Sen, n, [r.s_e, r.s_n, r.s_u]).astype(np.float32),
                'body_tau_updip_Pa': project(See, Snn, Sen, n, [r.sup_e, r.sup_n, r.sup_u]).astype(np.float32),
                'body_tau_strike_Pa': project(See, Snn, Sen, n, [r.sstr_e, r.sstr_n, r.sstr_u]).astype(np.float32),
                'body_areal_strain': (ee + nn).astype(np.float32),
                'body_S_ee_Pa': See.astype(np.float32), 'body_S_nn_Pa': Snn.astype(np.float32),
                'body_S_en_Pa': Sen.astype(np.float32)})
            parts.append(df)
        print(f'{str(y0)[:4]}: {len(t)} h, {_time.time()-t0:.1f}s', flush=True)
    body = pd.concat(parts, ignore_index=True).sort_values(['fid', 'time'])
    body.to_parquet(OUT / 'body_stress_hourly.parquet', index=False)

    # ---- check vs independent analytic formulation (family 0, 2014-01..2014-03) ----------
    tt = t_all[t_all < np.datetime64('2014-04-01')]
    i = 0
    (ee_a, nn_a, en_a), _ = analytic_strain(lat[i], lon[i], tt)
    b = body[(body.fid == geo.fid[i]) & (body.time < np.datetime64('2014-04-01'))]
    ds = C.SET_displacements(slon[i::nf], slat[i::nf], tt, type='time series', standard='datetime',
                             variable=['N', 'E', 'R'])
    ee, nn, en = surface_strain(ds['N'].values, ds['E'].values, ds['R'].values, lat[i], D, D)
    rows = []
    for name, x, y in [('e_ee', ee, ee_a), ('e_nn', nn, nn_a), ('e_en', en, en_a), ('areal', ee + nn, ee_a + nn_a)]:
        rows.append(dict(component=name, std_iers=np.std(x), std_analytic=np.std(y),
                         ratio=np.std(x) / np.std(y), corr=np.corrcoef(x, y)[0, 1],
                         rms_diff_frac=np.std(x - y) / np.std(y)))
    chk = pd.DataFrame(rows)
    chk.to_csv(OUT / 'body_check_iers_vs_analytic.csv', index=False)
    print(chk.round(4))
    # vertical displacement check (pyTMD R vs h2 W/g implied) is embedded in the strain check
    s = body.groupby('fid')[['body_tau_Pa', 'body_sigma_n_Pa', 'body_S_ee_Pa', 'body_S_nn_Pa']].std()
    print('per-family std (Pa):'); print(s.describe().round(1))


if __name__ == '__main__':
    main()
