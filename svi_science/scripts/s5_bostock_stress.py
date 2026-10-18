"""S5 step 1: tidal stress 2003-01-01 .. 2014-01-31 (hourly UTC) at all 130 Bostock family locations.

Re-uses the S5 setup machinery unchanged:
  geometry  = s5_tides_geometry.py logic (Slab2 plane fit within 10 km, 236 deg slip direction)
  body tide = s5_tides_body.py (pyTMD IERS2010 SET displacements -> surface strain -> plane stress)
  ocean     = s5_tides_ocean.py Boussinesq kernel on the saved EOT20+gauge cell fields (inner + outer)
Acceptance (PLAN_S5.md section 1): ocean constituents vs the setup csv for the shared 115 families, and the
hourly series for 2014-01 vs the setup parquet.
Outputs: results/s5/bostock_fault_geometry.csv, bostock_ocean_constituents.csv,
         tidal_stress_bostock_hourly.parquet, step1_acceptance.json
"""
import json, numpy as np, pandas as pd, xarray as xr, time as _time
from pathlib import Path
import pyTMD.compute as C
from s5_tides_common import ROOT, MODEL_CONST, RHO_W, G0, predict, gc_dist_az
import s5_tides_body as B
from s5_tides_ocean import boussinesq_enu, traction, COMPS

SETUP = ROOT / 'svi_science/results/s5_setup'
OUT = ROOT / 'svi_science/results/s5'
CONV_AZ = 236.0
T0, T1 = '2003-01-01', '2014-02-01'


def geometry():
    e = pd.read_csv(ROOT / 'lfe_catalog_repair/results/e5r_families_2012.csv')
    slab = xr.open_dataset(ROOT / 'data/slab_models/slab2/cas_slab2_dep_02.24.18.grd')
    X, Y = np.meshgrid(slab.x.values - 360.0, slab.y.values)
    Z = -slab.z.values.astype(float)
    rows = []
    for _, f in e.iterrows():
        lat0, lon0 = f.lat, f.lon
        dx = (X - lon0) * 111.195 * np.cos(np.radians(lat0)); dy = (Y - lat0) * 111.195
        m = (np.hypot(dx, dy) <= 10.0) & np.isfinite(Z)
        A = np.c_[np.ones(m.sum()), dx[m], dy[m]]
        c, *_ = np.linalg.lstsq(A, Z[m], rcond=None)
        gx, gy = c[1], c[2]
        dip = np.degrees(np.arctan(np.hypot(gx, gy))); az_dd = np.degrees(np.arctan2(gx, gy)) % 360
        strike = (az_dd - 90) % 360
        d, a = np.radians(dip), np.radians(az_dd)
        n = np.array([np.sin(d) * np.sin(a), np.sin(d) * np.cos(a), np.cos(d)])
        v_dd = np.array([np.cos(d) * np.sin(a), np.cos(d) * np.cos(a), -np.sin(d)])
        s_up = -v_dd
        s_str = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.0])
        h = np.array([np.sin(np.radians(CONV_AZ)), np.cos(np.radians(CONV_AZ)), 0.0])
        s_conv = h - h.dot(n) * n; s_conv /= np.linalg.norm(s_conv)
        rake = np.degrees(np.arctan2(s_conv.dot(s_up), s_conv.dot(s_str)))
        rows.append(dict(fid=int(f.bkey), lat=lat0, lon=lon0, dep_cat_km=f.dep, strike=strike, dip=dip,
                         az_downdip=az_dd, rake_conv236=rake, n_e=n[0], n_n=n[1], n_u=n[2],
                         s_e=s_conv[0], s_n=s_conv[1], s_u=s_conv[2], sup_e=s_up[0], sup_n=s_up[1], sup_u=s_up[2],
                         sstr_e=s_str[0], sstr_n=s_str[1], sstr_u=s_str[2]))
    g = pd.DataFrame(rows)
    g.to_csv(OUT / 'bostock_fault_geometry.csv', index=False)
    return g


def ocean_constituents(geo):
    grids = {z: np.load(SETUP / f'ocean_cells_EOT20_gauge_{z}.npz') for z in ('inner', 'outer')}
    rows = []
    for _, f in geo.iterrows():
        zf = f.dep_cat_km * 1e3
        n = np.array([f.n_e, f.n_n, f.n_u])
        dirs = {'tau': np.array([f.s_e, f.s_n, f.s_u]), 'tau_updip': np.array([f.sup_e, f.sup_n, f.sup_u]),
                'tau_strike': np.array([f.sstr_e, f.sstr_n, f.sstr_u]), 'sigma_n': n}
        tot = {k: 0 for k in dirs}
        for z, G in grids.items():
            LA, LO, ar, F = G['lat'], G['lon'], G['area'], G['F']
            dkm, az = gc_dist_az(f.lat, f.lon, LA, LO)
            dx = -dkm * 1e3 * np.sin(np.radians(az)); dy = -dkm * 1e3 * np.cos(np.radians(az))
            K = boussinesq_enu(dx, dy, zf)
            w = RHO_W * G0 * ar
            S = {cmp: (K[cmp] * w) @ F.T for cmp in COMPS}
            for nm, sv in dirs.items():
                tot[nm] = tot[nm] + traction(S, n, sv)
        for nm, val in tot.items():
            for ic, c in enumerate(MODEL_CONST):
                rows.append(dict(fid=int(f.fid), comp=nm, const=c, re=val[ic].real, im=val[ic].imag))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / 'bostock_ocean_constituents.csv', index=False)
    return T


def hourly(geo, T):
    t_all = np.arange(T0, T1, np.timedelta64(1, 'h'), dtype='datetime64[ns]')
    fids = geo.fid.values; nf = len(fids)
    # ---- ocean
    oc = {}
    for comp in ('tau', 'tau_updip', 'tau_strike', 'sigma_n'):
        q = T[T.comp == comp]
        hc = {c: q[q.const == c].set_index('fid').reindex(fids).eval('re + 1j*im').values for c in MODEL_CONST}
        parts = []
        for y in range(2003, 2015):
            tt = t_all[(t_all >= np.datetime64(f'{y}-01-01')) & (t_all < np.datetime64(f'{y+1}-01-01'))]
            if len(tt):
                parts.append(predict(hc, tt))
        oc[comp] = np.concatenate(parts, axis=1)
    # ---- body
    lat = geo.lat.values; lon = geo.lon.values; D = B.D
    slat = np.concatenate([lat, lat + D, lat - D, lat, lat]); slon = np.concatenate([lon, lon, lon, lon + D, lon - D])
    bd = {k: np.zeros((nf, len(t_all)), np.float32) for k in ('tau', 'sigma_n', 'tau_updip', 'tau_strike')}
    for y in range(2003, 2015):
        sel = np.where((t_all >= np.datetime64(f'{y}-01-01')) & (t_all < np.datetime64(f'{y+1}-01-01')))[0]
        if len(sel) == 0:
            continue
        t0 = _time.time()
        ds = C.SET_displacements(slon, slat, t_all[sel], type='time series', standard='datetime', variable=['N', 'E', 'R'])
        uN = ds['N'].values.reshape(5, nf, -1); uE = ds['E'].values.reshape(5, nf, -1); uR = ds['R'].values.reshape(5, nf, -1)
        for i in range(nf):
            ee, nn, en = B.surface_strain(uN[:, i], uE[:, i], uR[:, i], lat[i], D, D)
            See, Snn, Sen = B.stress_from_strain(ee, nn, en)
            r = geo.iloc[i]; n = np.array([r.n_e, r.n_n, r.n_u])
            bd['sigma_n'][i, sel] = B.project(See, Snn, Sen, n, n)
            bd['tau'][i, sel] = B.project(See, Snn, Sen, n, [r.s_e, r.s_n, r.s_u])
            bd['tau_updip'][i, sel] = B.project(See, Snn, Sen, n, [r.sup_e, r.sup_n, r.sup_u])
            bd['tau_strike'][i, sel] = B.project(See, Snn, Sen, n, [r.sstr_e, r.sstr_n, r.sstr_u])
        print(y, f'{_time.time()-t0:.1f}s', flush=True)
    df = pd.DataFrame({'fid': np.repeat(fids, len(t_all)).astype(np.int16), 'time': np.tile(t_all, nf)})
    for comp in ('tau', 'sigma_n', 'tau_updip', 'tau_strike'):
        df[f'body_{comp}_Pa'] = bd[comp].ravel()
        df[f'ocean_{comp}_Pa'] = oc[comp].ravel().astype(np.float32)
        df[f'{comp}_Pa'] = (df[f'body_{comp}_Pa'] + df[f'ocean_{comp}_Pa']).astype(np.float32)
    df['coulomb_0p4'] = (df.tau_Pa + 0.4 * df.sigma_n_Pa).astype(np.float32)
    df.to_parquet(OUT / 'tidal_stress_bostock_hourly.parquet', index=False)
    return df


def main():
    geo = geometry()
    T = ocean_constituents(geo)
    acc = {}
    # (a) ocean constituents vs setup
    S = pd.read_csv(SETUP / 'ocean_constituent_stress.csv')
    S = S[(S.model == 'EOT20') & (S.fill == 'gauge') & (S.depth == 'cat')].groupby(['fid', 'comp', 'const'])[['re', 'im']].sum()
    M = T.set_index(['fid', 'comp', 'const']).join(S, rsuffix='_setup', how='inner')
    rel = np.hypot(M.re - M.re_setup, M.im - M.im_setup) / np.hypot(M.re_setup, M.im_setup).clip(lower=1e-3)
    big = M[np.hypot(M.re_setup, M.im_setup) > 10]
    relb = np.hypot(big.re - big.re_setup, big.im - big.im_setup) / np.hypot(big.re_setup, big.im_setup)
    acc['ocean_const_n_compared'] = int(len(M)); acc['ocean_const_max_rel_diff_amp_gt_10Pa'] = float(relb.max())
    acc['geometry_vs_setup_max_dip_diff'] = float((geo.set_index('fid').dip - pd.read_csv(SETUP / 'family_fault_geometry.csv').set_index('fid').dip).abs().max())
    df = hourly(geo, T)
    # (b) hourly 2014-01 vs setup
    s = pd.read_parquet(SETUP / 'tidal_stress_svi_hourly.parquet', columns=['fid', 'time', 'tau_Pa', 'sigma_n_Pa', 'body_tau_Pa', 'ocean_tau_Pa'],
                        filters=[('time', '<', pd.Timestamp('2014-02-01'))])
    j = df[df.time >= pd.Timestamp('2014-01-01')].merge(s, on=['fid', 'time'], suffixes=('', '_setup'))
    for c in ('tau_Pa', 'sigma_n_Pa', 'body_tau_Pa', 'ocean_tau_Pa'):
        acc[f'hourly_2014_01_max_abs_diff_{c}'] = float((j[c] - j[f'{c}_setup']).abs().max())
    acc['hourly_rows_compared'] = int(len(j)); acc['n_families'] = int(df.fid.nunique()); acc['hours'] = int(df.time.nunique())
    acc['tau_std_median_Pa'] = float(df.groupby('fid').tau_Pa.std().median())
    acc['pass'] = bool(acc['ocean_const_max_rel_diff_amp_gt_10Pa'] < 1e-3 and acc['hourly_2014_01_max_abs_diff_tau_Pa'] < 1.0
                       and acc['hourly_2014_01_max_abs_diff_sigma_n_Pa'] < 1.0)
    json.dump(acc, open(OUT / 'step1_acceptance.json', 'w'), indent=1)
    print(json.dumps(acc, indent=1))


if __name__ == '__main__':
    main()
