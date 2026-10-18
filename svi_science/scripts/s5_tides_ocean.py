"""S5 setup, step 3: ocean-tide LOADING stress on the plate interface at each SVI family.

Method
  * Ocean tide elevations (pure ocean tide, relative to sea floor) for the 10 constituents
    q1 o1 p1 k1 n2 m2 s2 k2 2n2 m4 from EOT20 (primary; SEANOE, free) or GOT5.5 (check; NASA, free).
  * Integration grid: 0.02 deg cells within 330 km of the family-cloud centre (48.6N, 123.6W)
    ("inner"), 0.125 deg cells from 330 to 700 km ("outer", far-field estimate). Ocean cells =
    not Natural-Earth-10m land and not NE-10m lakes.
  * Elevation on each ocean cell: normalised bilinear interpolation of the valid model nodes.
    Ocean cells the model does not resolve (inland Salish Sea channels, Puget Sound, inlets) are
    filled by inverse-distance (k=8, power 2) from the model-valid cells PLUS the tide-gauge
    harmonic constants (NOAA + CHS; s5_tides_gauges.py) -> fill='gauge' (primary), or from the
    nearest model-valid cell only -> fill='nearest' (sensitivity variant).
  * Each cell is a vertical point force P = rho_w g eta A_cell on an elastic, homogeneous,
    isotropic HALF-SPACE (Boussinesq; nu = 0.25; stresses independent of shear modulus).
    Stress tensor at the family hypocentre (catalogue depth; Slab2 depth as a variant).
  * Constituent-domain: complex stress per constituent = sum over cells; projected onto the
    Slab2 fault plane -> sigma_n (tension +) and tau (positive = promotes hanging-wall slip in
    the 236 deg convergence direction; also up-dip and strike components).
  * Hourly series 2014-01-01..2026-09-30 by harmonic prediction (pyTMD nodal corrections, FES
    convention, minor constituents inferred from the majors by admittance).
Validation in this script: (1) Boussinesq checks (force balance, circular-load closed form,
equilibrium by finite differences, thrust sign test); (2) model vs tide-gauge constants.
Outputs (results/s5_setup/): boussinesq_validation.csv, model_vs_gauges.csv,
ocean_constituent_stress.csv, ocean_stress_hourly.parquet, ocean_cells_<model>_<fill>.npz
"""
import sys, numpy as np, pandas as pd, xarray as xr
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
import shapely
import cartopy.io.shapereader as shpreader
from s5_tides_common import (OUT, TIDES, CENTER, MODEL_CONST, RHO_W, G0, predict, gc_dist_az)

NU = 0.25
DG_IN, R_IN, DG_OUT, R_OUT = 0.02, 330.0, 0.125, 700.0
BOX = (40.0, 57.0, -138.0, -108.0)  # lat0, lat1, lon0, lon1 for model subsetting


# ----------------------------------------------------------------------------- Boussinesq
def boussinesq_enu(dx, dy, z, nu=NU):
    """Stress (Pa, tension +) per unit downward point force (N) at the surface origin, at an
    observation point offset (dx east, dy north) [m] and depth z>0 [m]. ENU, U up.
    Johnson (1985, Contact Mechanics, eq. 3.21) cylindrical solution rotated to ENU."""
    r = np.hypot(dx, dy)
    rho = np.sqrt(r**2 + z**2)
    with np.errstate(invalid='ignore', divide='ignore'):
        c = np.where(r > 0, dx / r, 1.0); s = np.where(r > 0, dy / r, 0.0)
    k = 1 / (2 * np.pi)
    srr = k * ((1 - 2 * nu) / (rho * (rho + z)) - 3 * z * r**2 / rho**5)
    stt = -k * (1 - 2 * nu) * (1 / (rho * (rho + z)) - z / rho**3)
    szz = -3 * k * z**3 / rho**5
    srz = -3 * k * r * z**2 / rho**5          # (r, z-down) component
    return dict(EE=srr * c**2 + stt * s**2, NN=srr * s**2 + stt * c**2, UU=szz,
                EN=(srr - stt) * c * s, EU=-srz * c, NU=-srz * s)


COMPS = ['EE', 'NN', 'UU', 'EN', 'EU', 'NU']


def traction(S, n, s):
    """S: dict comp -> array (...); n, s 3-vectors (E,N,U). returns s.S.n"""
    T = {('E', 'E'): 'EE', ('N', 'N'): 'NN', ('U', 'U'): 'UU', ('E', 'N'): 'EN', ('N', 'E'): 'EN',
         ('E', 'U'): 'EU', ('U', 'E'): 'EU', ('N', 'U'): 'NU', ('U', 'N'): 'NU'}
    ax = 'ENU'
    return sum(s[i] * S[T[(ax[i], ax[j])]] * n[j] for i in range(3) for j in range(3))


def validate_boussinesq():
    rows = []
    z = 30e3
    # (1) force balance: integral of sigma_zz over a horizontal plane = -P
    h = 500.0; x = np.arange(-600e3, 600e3, h) + h / 2
    X, Y = np.meshgrid(x, x)
    S = boussinesq_enu(X, Y, z)
    rows.append(dict(test='force_balance_int_szz_dA', value=S['UU'].sum() * h * h, expected=-1.0))
    rows.append(dict(test='int_sEU_dA (should be 0 by symmetry)', value=S['EU'].sum() * h * h, expected=0.0))
    # (2) uniform circular load p=1 Pa radius a, on axis at depth z (Johnson/Poulos & Davis)
    for a in (10e3, 30e3, 100e3):
        m = np.hypot(X, Y) <= a
        Sd = {k: (v.sum() * h * h) for k, v in boussinesq_enu(-X[m], -Y[m], z).items()}
        q = z / np.sqrt(a**2 + z**2)
        rows.append(dict(test=f'disk a={a/1e3:.0f}km sigma_z', value=Sd['UU'], expected=-(1 - q**3)))
        rows.append(dict(test=f'disk a={a/1e3:.0f}km sigma_r', value=Sd['EE'],
                         expected=-0.5 * ((1 + 2 * NU) - 2 * (1 + NU) * q + q**3)))
    # (3) equilibrium: div(S) = 0 by central differences at a generic point
    p0 = np.array([7e3, -4e3, 25e3]); e = 5.0
    def S_at(dx, dy, dz):
        return boussinesq_enu(np.array(p0[0] + dx), np.array(p0[1] + dy), np.array(p0[2] - dz))  # dz = +up
    d = lambda comp, ax: ((S_at(*(e * np.eye(3)[ax]))[comp] - S_at(*(-e * np.eye(3)[ax]))[comp]) / (2 * e))
    scale = abs(S_at(0, 0, 0)['UU']) / p0[2]
    for i, (a1, a2, a3) in enumerate([('EE', 'EN', 'EU'), ('EN', 'NN', 'NU'), ('EU', 'NU', 'UU')]):
        rows.append(dict(test=f'equilibrium row {"ENU"[i]} (relative)',
                         value=float((d(a1, 0) + d(a2, 1) + d(a3, 2)) / scale), expected=0.0))
    # (4) fault sign test: thrust dipping 10 deg toward +E (hanging wall above), horizontal
    #     E-W compression -> tau(up-dip) > 0 and sigma_n < 0
    dip = np.radians(10)
    n = np.array([np.sin(dip), 0, np.cos(dip)]); s_up = np.array([-np.cos(dip), 0, np.sin(dip)])
    Sc = {k: 0.0 for k in COMPS}; Sc['EE'] = -1.0
    rows.append(dict(test='thrust sign: E-W compression tau_updip', value=traction(Sc, n, s_up), expected=np.sin(dip) * np.cos(dip)))
    rows.append(dict(test='thrust sign: E-W compression sigma_n', value=traction(Sc, n, n), expected=-np.sin(dip)**2))
    v = pd.DataFrame(rows)
    v['ok'] = np.isclose(v.value, v.expected, rtol=2e-2, atol=2e-3)
    v.to_csv(OUT / 'boussinesq_validation.csv', index=False)
    print(v.to_string()); return v


# ----------------------------------------------------------------------------- models
def read_model(name):
    H = {}
    for c in MODEL_CONST:
        if name == 'EOT20':
            ds = xr.open_dataset(TIDES / f'EOT20/ocean_tides/{c.upper()}_ocean_eot20.nc')
            lat, lon = ds.lat.values, ds.lon.values
        else:
            ds = xr.open_dataset(TIDES / f'GOT/GOT5.5/ocean_tides/{c}.nc')
            lat, lon = ds.latitude.values.astype(float), ds.longitude.values.astype(float)
        lon = np.where(lon > 180, lon - 360, lon)
        iy = (lat >= BOX[0]) & (lat <= BOX[1]); ix = (lon >= BOX[2]) & (lon <= BOX[3])
        amp = ds.amplitude.values[np.ix_(iy, ix)].astype(float) / 100.0
        ph = ds.phase.values[np.ix_(iy, ix)].astype(float)
        hc = amp * np.exp(-1j * np.radians(ph))
        hc[~np.isfinite(hc) | (amp <= 0) | (amp > 50)] = np.nan
        o = np.argsort(lon[ix])
        H[c] = hc[:, o]; LAT = lat[iy]; LON = lon[ix][o]
    return LAT, LON, H


def sample(LAT, LON, H2, plat, plon):
    """normalised bilinear interpolation over valid nodes; NaN if no valid node."""
    v = np.isfinite(H2)
    f = lambda arr: RegularGridInterpolator((LAT, LON), arr, bounds_error=False, fill_value=0.0)((plat, plon))
    wsum = f(v.astype(float))
    num = f(np.where(v, H2.real, 0)) + 1j * f(np.where(v, H2.imag, 0))
    out = num / np.where(wsum > 1e-6, wsum, np.nan)
    return out


# ----------------------------------------------------------------------------- grid
def ocean_mask(lat, lon):
    bbox = (lon.min() - 0.5, lat.min() - 0.5, lon.max() + 0.5, lat.max() + 0.5)
    geoms = []
    for nm in ('land', 'lakes'):
        rdr = shpreader.Reader(shpreader.natural_earth('10m', 'physical', nm))
        for g in rdr.geometries():
            g = shapely.clip_by_rect(g, *bbox)
            if not g.is_empty:
                geoms.append(g)
    U = shapely.union_all(geoms)
    shapely.prepare(U)
    return ~shapely.contains_xy(U, lon, lat)


def make_grid(dg, rmin, rmax):
    dlat = rmax / 111.0 + dg; dlon = rmax / (111.0 * np.cos(np.radians(CENTER[0] + dlat))) + dg
    la = np.arange(CENTER[0] - dlat, CENTER[0] + dlat, dg) + dg / 2
    lo = np.arange(CENTER[1] - dlon, CENTER[1] + dlon, dg) + dg / 2
    LA, LO = np.meshgrid(la, lo, indexing='ij'); LA = LA.ravel(); LO = LO.ravel()
    d, _ = gc_dist_az(*CENTER, LA, LO)
    m = (d > rmin) & (d <= rmax)
    LA, LO = LA[m], LO[m]
    oc = ocean_mask(LA, LO)
    LA, LO = LA[oc], LO[oc]
    area = (dg * 111195.0) ** 2 * np.cos(np.radians(LA))
    return LA, LO, area


# Inland waters (Juan de Fuca Strait east of the Cape Flattery-Carmanah Pt mouth, Haro Strait,
# San Juan Is., Puget Sound, Hood Canal, Strait of Georgia, Discovery Is.): both global models
# fail there (M2 at Victoria 0.97 m modelled vs 0.37 m observed; phase errors 60-130 deg), so
# inside this polygon the elevation field is interpolated from the tide gauges only. The polygon
# edges run over land (Vancouver Island spine, Olympic Peninsula) except at the strait mouth
# and the Johnstone Strait end.
INLAND = shapely.Polygon([(-124.75, 48.61), (-124.40, 48.75), (-124.30, 49.05), (-124.60, 49.30),
                          (-125.10, 49.75), (-125.60, 50.10), (-125.90, 50.60), (-121.0, 50.60),
                          (-121.0, 46.80), (-123.5, 46.80), (-123.5, 47.95), (-124.2, 48.05),
                          (-124.68, 48.30), (-124.74, 48.40)])


def inland(lat, lon):
    return shapely.contains_xy(INLAND, lon, lat)


def _xy(lat, lon):
    return np.c_[(lon - CENTER[1]) * 111.195 * np.cos(np.radians(CENTER[0])), (lat - CENTER[0]) * 111.195]


def gauge_idw(g, lat, lon, k=6):
    """complex IDW (power 2) from gauge constants g (one constituent) to points."""
    src = g.amp_m.values * np.exp(-1j * np.radians(g.phase_deg.values))
    tree = cKDTree(_xy(g.lat.values, g.lon.values))
    dd, ii = tree.query(_xy(lat, lon), k=min(k, len(g)))
    w = 1.0 / np.maximum(dd, 0.5) ** 2
    return (w * src[ii]).sum(1) / w.sum(1)


def build_field(model, LA, LO, gauges, fill):
    """fill='gauge': inland polygon from gauges only; elsewhere model, gaps by IDW from model
    cells + gauges.  fill='nearest': pure model, gaps by nearest model cell (no gauges)."""
    LAT, LON, H = read_model(model)
    xy = _xy(LA, LO)
    inl = inland(LA, LO)
    F = np.zeros((len(MODEL_CONST), len(LA)), complex); src_flag = None
    for ic, c in enumerate(MODEL_CONST):
        v = sample(LAT, LON, H[c], LA, LO)
        ok = np.isfinite(v)
        if fill == 'gauge':
            g = gauges[gauges.const == c]
            gin = g[inland(g.lat.values, g.lon.values)]
            if inl.any():
                v[inl] = gauge_idw(gin, LA[inl], LO[inl])
            ok_model = ok & ~inl
            gap = ~ok & ~inl
            if gap.any():
                gxy = _xy(g.lat.values, g.lon.values)
                sx = np.r_[xy[ok_model], gxy]
                sv = np.r_[v[ok_model], g.amp_m.values * np.exp(-1j * np.radians(g.phase_deg.values))]
                dd, ii = cKDTree(sx).query(xy[gap], k=8)
                w = 1.0 / np.maximum(dd, 0.5) ** 2
                v[gap] = (w * sv[ii]).sum(1) / w.sum(1)
            flag = np.where(inl, 2, np.where(ok, 0, 1))
        else:
            if (~ok).any():
                dd, ii = cKDTree(xy[ok]).query(xy[~ok], k=1)
                v[~ok] = v[ok][ii]
            flag = np.where(ok, 0, 1)
        if src_flag is None:
            src_flag = flag
        F[ic] = v
    return F, src_flag   # flag 0 = model, 1 = gap-filled, 2 = inland gauge field


def gauge_loo(gauges):
    """leave-one-out check of the inland gauge field."""
    rows = []
    for c in MODEL_CONST:
        g = gauges[(gauges.const == c)]
        g = g[inland(g.lat.values, g.lon.values)].reset_index(drop=True)
        for j in range(len(g)):
            pred = gauge_idw(g.drop(j), g.lat.values[j:j + 1], g.lon.values[j:j + 1])[0]
            obs = g.amp_m[j] * np.exp(-1j * np.radians(g.phase_deg[j]))
            rows.append(dict(const=c, id=g.id[j], name=g['name'][j], lat=g.lat[j], lon=g.lon[j], amp_obs=abs(obs),
                             amp_pred=abs(pred), vec_misfit_m=abs(obs - pred), rel_misfit=abs(obs - pred) / abs(obs)))
    t = pd.DataFrame(rows); t.to_csv(OUT / 'gauge_field_loo.csv', index=False)
    return t


def model_vs_gauges(gauges):
    rows = []
    for model in ('EOT20', 'GOT5.5'):
        LAT, LON, H = read_model(model)
        gs = gauges.drop_duplicates('id')
        for c in MODEL_CONST:
            vv = sample(LAT, LON, H[c], gs.lat.values, gs.lon.values)
            # nearest valid node if bilinear has none
            valid = np.argwhere(np.isfinite(H[c]))
            tree = cKDTree(np.c_[LAT[valid[:, 0]], LON[valid[:, 1]] * np.cos(np.radians(48.6))])
            dd, ii = tree.query(np.c_[gs.lat.values, gs.lon.values * np.cos(np.radians(48.6))])
            near = H[c][valid[ii, 0], valid[ii, 1]]
            for j, (_, g) in enumerate(gs.iterrows()):
                gg = gauges[(gauges.id == g.id) & (gauges.const == c)]
                if len(gg) == 0:
                    continue
                obs = gg.amp_m.values[0] * np.exp(-1j * np.radians(gg.phase_deg.values[0]))
                mod = vv[j] if np.isfinite(vv[j]) else near[j]
                rows.append(dict(model=model, id=g.id, name=g['name'], source=g.source, lat=g.lat, lon=g.lon, inland=bool(inland(np.array([g.lat]), np.array([g.lon]))[0]),
                                 const=c, model_node='bilinear' if np.isfinite(vv[j]) else f'nearest {dd[j]*111:.0f}km',
                                 amp_obs=abs(obs), amp_mod=abs(mod),
                                 dphase_deg=((np.degrees(np.angle(obs)) - np.degrees(np.angle(mod)) + 180) % 360) - 180,
                                 vec_misfit_m=abs(obs - mod)))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / 'model_vs_gauges.csv', index=False)
    return t


# ----------------------------------------------------------------------------- main
def main():
    validate_boussinesq()
    geo = pd.read_csv(OUT / 'family_fault_geometry.csv')
    gauges = pd.read_csv(OUT / 'gauge_constants.csv', dtype={'id': str})
    # drop gauges in restricted basins / rivers (behind tidal rapids or river-dominated):
    # Victoria Gorge + Portage Inlet, Sechelt Inlet (Porpoise Bay), lower Fraser River
    bad = gauges.name.str.contains('Gorge|Portage Inlet|Porpoise Bay|Sechelt|New Westminster|Woodward|Mission|Deas', case=False)
    print('gauges dropped (restricted basins/rivers):', sorted(gauges[bad].name.unique()))
    gauges = gauges[~bad]
    mv = model_vs_gauges(gauges)
    print(mv[mv.const.isin(['m2', 'k1'])].groupby(['model', 'inland', 'const'])[['vec_misfit_m', 'amp_obs']].median().round(3))
    loo = gauge_loo(gauges)
    print('inland gauge-field leave-one-out:')
    print(loo.groupby('const')[['vec_misfit_m', 'rel_misfit', 'amp_obs']].median().round(3), flush=True)

    grids = {'inner': make_grid(DG_IN, -1, R_IN), 'outer': make_grid(DG_OUT, R_IN, R_OUT)}
    for z, (LA, LO, ar) in grids.items():
        print(f'{z}: {len(LA)} ocean cells', flush=True)
    variants = [('EOT20', 'gauge'), ('EOT20', 'nearest'), ('GOT5.5', 'gauge')]
    fields = {}
    for model, fill in variants:
        for z, (LA, LO, ar) in grids.items():
            F, filled = build_field(model, LA, LO, gauges, fill)
            fields[(model, fill, z)] = F
            np.savez_compressed(OUT / f'ocean_cells_{model}_{fill}_{z}.npz', lat=LA, lon=LO, area=ar, F=F,
                                filled=filled, consts=np.array(MODEL_CONST))
            print(model, fill, z, 'cell source fractions (0 model,1 gap-fill,2 inland gauges):', np.bincount(filled, minlength=3) / len(filled), flush=True)

    rows = []
    for depth_opt in ('cat', 'slab2'):
        for i, f in geo.iterrows():
            zf = (f.dep_cat_km if depth_opt == 'cat' else f.slab2_depth_km) * 1e3
            n = np.array([f.n_e, f.n_n, f.n_u])
            dirs = {'tau': np.array([f.s_e, f.s_n, f.s_u]), 'tau_updip': np.array([f.sup_e, f.sup_n, f.sup_u]),
                    'tau_strike': np.array([f.sstr_e, f.sstr_n, f.sstr_u]), 'sigma_n': n}
            for z, (LA, LO, ar) in grids.items():
                dkm, az = gc_dist_az(f.lat, f.lon, LA, LO)
                dx = -dkm * 1e3 * np.sin(np.radians(az)); dy = -dkm * 1e3 * np.cos(np.radians(az))
                K = boussinesq_enu(dx, dy, zf)
                w = RHO_W * G0 * ar
                for model, fill in variants:
                    if depth_opt == 'slab2' and (model, fill) != ('EOT20', 'gauge'):
                        continue
                    F = fields[(model, fill, z)]
                    S = {cmp: (K[cmp] * w) @ F.T for cmp in COMPS}   # (nconst,) complex, Pa per constituent
                    for nm, sv in dirs.items():
                        val = traction(S, n, sv)
                        for ic, c in enumerate(MODEL_CONST):
                            rows.append(dict(fid=int(f.fid), model=model, fill=fill, depth=depth_opt, zone=z,
                                             comp=nm, const=c, re=val[ic].real, im=val[ic].imag))
            if i % 20 == 0:
                print(depth_opt, i, flush=True)
    T = pd.DataFrame(rows)
    T['amp_Pa'] = np.hypot(T.re, T.im); T['phase_deg'] = (-np.degrees(np.arctan2(T.im, T.re))) % 360
    T.to_csv(OUT / 'ocean_constituent_stress.csv', index=False)

    # ---- hourly series: primary variant EOT20 + gauge fill, catalogue depth, inner+outer
    P_ = T[(T.model == 'EOT20') & (T.fill == 'gauge') & (T.depth == 'cat')]
    P_ = P_.groupby(['fid', 'comp', 'const'])[['re', 'im']].sum().reset_index()
    fids = np.sort(P_.fid.unique())
    t_all = np.arange('2014-01-01', '2026-10-01', np.timedelta64(1, 'h'), dtype='datetime64[ns]')
    out = {}
    for comp in ('tau', 'tau_updip', 'tau_strike', 'sigma_n'):
        q = P_[P_.comp == comp]
        hc = {c: (q[q.const == c].set_index('fid').reindex(fids).eval('re + 1j*im')).values for c in MODEL_CONST}
        out[comp] = np.concatenate([predict(hc, t_all[(t_all >= np.datetime64(f'{y}-01-01')) & (t_all < np.datetime64(f'{y+1}-01-01'))])
                                    for y in range(2014, 2027)], axis=1)
    df = pd.DataFrame({'fid': np.repeat(fids, len(t_all)).astype(np.int16), 'time': np.tile(t_all, len(fids))})
    for comp, arr in out.items():
        df[f'ocean_{comp}_Pa'] = arr.ravel().astype(np.float32)
    df.to_parquet(OUT / 'ocean_stress_hourly.parquet', index=False)
    print(df.groupby('fid')[['ocean_tau_Pa', 'ocean_sigma_n_Pa']].std().describe().round(1))


if __name__ == '__main__':
    main()
