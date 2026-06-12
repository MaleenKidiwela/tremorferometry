#!/usr/bin/env python
"""NULL TEST for the multi-window joint inversion.

Question: does the southern "hot zone" temporal RMS of 0.1–0.16% (and up to 0.30% in the 40–42 band,
1–2 stations) reflect real fault velocity change, or is it noise amplified through a sparse design matrix?

Method
------
Plant ZERO everywhere (zero true fault field).  Forward-model every (cell, station, window, month) datum
that exists in the REAL dataset using the same kernels and coverage pattern, but add synthetic noise
drawn from the REAL per-(station, window) residual statistics (fit from the real inversion's residuals;
see below).  Invert with the SAME machinery (joint LSQR + Laplacian + site terms, same lam_f/lam_s).
Repeat for N_REAL ≥ 20 noise realisations.

Per-(station, window) noise std is estimated from the real inversion residuals: for each station|window
combination, pool residuals across all months where that pair appears, then take the standard deviation.
This captures the realistic noise level without assuming it equals any global constant.

Outputs
-------
  fault_tomography/inversion/null_test_cache.npz  — assembled G, months, cell metadata (skip CSV reload)
  fault_tomography/inversion/null_test_results.npz — per-realisation recovered MF, summary statistics
  figures/null_test_southern.png                  — null RMS distribution vs observed, by lat band + coverage
  fault_tomography/PHASE_A_RESULTS.md             — appended results section (append-only)
"""
import os, sys, time
import numpy as np, pandas as pd
from scipy import sparse
from scipy.sparse.linalg import lsqr
from scipy.spatial import cKDTree
from scipy.interpolate import griddata
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

sys.path.insert(0, 'fault_tomography/kernels')
from kernel import kernel_singlescatter

# ── Constants (must match invert_multiwin.py exactly) ──────────────────────────────────────────
GRID = 0.10; BETA = 3.5; LSTAR = 40.0; CCMIN = 0.7
_SFX = os.environ.get('SFX', '_cal_des')   # SFX=_calT_des for the true-scale 35-station run
WINS = [('13', 1.0, 3.0, '1to3'+_SFX), ('24', 2.0, 4.0, '2to4'+_SFX), ('35', 3.0, 5.0, '3to5'+_SFX)]
STA = {
 'B927':(49.2188,-124.8113),'NLLB':(49.2271,-123.9882),'B928':(48.834,-125.134),'PGC':(48.6498,-123.4521),
 'B011':(48.65,-123.448),'B004':(48.202,-124.427),'B013':(47.813,-122.9108),'HDW':(47.6490,-123.0530),
 'GNW':(47.5641,-122.8250),'B014':(47.5133,-123.8125),'B941':(46.9868,-122.219),'B018':(46.9795,-123.0203),
 'B020':(46.3827,-123.8445),'B201':(46.3033,-122.2648),'B204':(46.136,-122.169),'B023':(46.1112,-123.0787),
 'B022':(45.9546,-123.931),'B026':(45.3094,-123.8231),'COLT':(45.17044,-122.438152),'COR':(44.5855,-123.3046),
 'B028':(44.4937,-122.9638),'B030':(43.9713,-122.7717),'B032':(43.668,-123.3923),'B033':(43.2917,-123.1245),
 'B036':(42.5058,-123.3817),'B040':(41.8308,-122.4205),'B039':(41.4667,-122.4847),'B935':(40.4787,-123.5732),
 'B017':(46.9960,-123.5575),'B001':(48.0431,-123.1314),'B005':(48.0596,-123.5034),'B003':(48.0623,-124.1416),
 'B045':(40.4360,-124.0008),'B932':(40.2825,-124.2245),'B049':(40.2403,-123.8225),  # +7 new 2026-06-10
}
lam_f, lam_s = 0.4, 0.05     # same regularisation as invert_multiwin.py
N_REAL = 25                   # number of noise realisations
_TAG   = os.environ.get('TAG', '')   # e.g. TAG=_calT35 -> separate cache/result/real-inversion
CACHE  = f'fault_tomography/inversion/null_test_cache{_TAG}.npz'
RESULT = f'fault_tomography/inversion/null_test_results{_TAG}.npz'
REAL_NPZ = os.environ.get('REAL_NPZ', 'fault_tomography/inversion/fault_4d_multiwin.npz')

# ── Step 1: Load or build cache (geometry + G matrices + monthly data lists) ───────────────────
if os.path.exists(CACHE):
    print(f'[null_test] loading cache from {CACHE}')
    c = np.load(CACHE, allow_pickle=True)
    Gf          = c['Gf']           # (n_cat, M)
    Gs_data     = c['Gs_data']      # COO data for Gs
    Gs_row      = c['Gs_row']
    Gs_col      = c['Gs_col']
    Gs_shape    = tuple(c['Gs_shape'])
    cells_lat   = c['cells_lat']
    cells_lon   = c['cells_lon']
    cells_nsta  = c['cells_nsta']
    months      = list(c['months'])
    key_cell    = c['key_cell']     # (n_cat,) -> cell index
    key_sta     = c['key_sta']      # (n_cat,) -> sta|win label
    sw_labels   = list(c['sw_labels'])
    M           = int(c['M'])
    B           = int(c['B'])
    # monthly index lists: for each month j, which cat-rows are active?
    month_rows  = c['month_rows'].tolist()   # list of arrays (object array)
    month_dd    = c['month_dd'].tolist()     # list of dvv arrays in %
    print(f'  {M} cells, {B} site terms, {len(months)} months loaded from cache')
else:
    t0 = time.time()
    print('[null_test] building geometry + G (first run -- will cache) ...')
    # ── 1a. load data tensor (same as invert_multiwin.py) ──
    parts = []
    for s, (sla, slo) in STA.items():
        for wk, w1, w2, suf in WINS:
            f = f'data/daily_dvv_{s}_{suf}.csv'
            if not os.path.exists(f): continue
            d = pd.read_csv(f); d = d[d.cc_max > CCMIN]
            if not len(d): continue
            pre  = d.patch.astype(str).str.split('__').str[0]
            plat = pre.str.split('_').str[0].astype(float)
            plon = pre.str.split('_', n=1).str[1].astype(float)
            d = d.assign(clat=(np.round(plat/GRID)*GRID).round(3),
                         clon=(np.round(plon/GRID)*GRID).round(3))
            d['cell'] = d.clat.astype(str) + '_' + d.clon.astype(str)
            d['ym']   = pd.to_datetime(d.date).dt.to_period('M').astype(str)
            d['station'] = s; d['sta_lat'] = sla; d['sta_lon'] = slo; d['win'] = wk
            g = d.groupby(['cell', 'station', 'win', 'ym']).agg(
                dvv=('dvv', 'mean'), clat=('clat', 'first'), clon=('clon', 'first'),
                sta_lat=('sta_lat', 'first'), sta_lon=('sta_lon', 'first')).reset_index()
            parts.append(g)
    ten = pd.concat(parts, ignore_index=True)
    print(f'  tensor: {len(ten):,} rows | {ten.cell.nunique()} cells | {ten.station.nunique()} sta '
          f'| {time.time()-t0:.0f}s')

    # ── 1b. cells + depth ──
    cells = ten[['cell','clat','clon']].drop_duplicates('cell').reset_index(drop=True)
    slab  = pd.read_csv('data/cas_slab2_input_04-18.csv', low_memory=False)[['lat','lon','depth','etype']].dropna(subset=['depth'])
    slab  = slab[(~slab.etype.isin(['BA','TO'])) & (slab.depth.between(10,70))]
    q = cells[['clat','clon']].values
    dl = griddata(slab[['lat','lon']].values, slab.depth.values, q, 'linear')
    dn = griddata(slab[['lat','lon']].values, slab.depth.values, q, 'nearest')
    cells['depth_km'] = np.clip(np.where(np.isnan(dl), dn, dl), 12, 55)
    nst = ten.groupby('cell')['station'].nunique().rename('n_stations')
    cells = cells.merge(nst, on='cell')
    lat0, lon0 = cells.clat.mean(), cells.clon.mean()
    def proj(lat, lon):
        return (np.asarray(lon)-lon0)*111.0*np.cos(np.radians(lat0)), (np.asarray(lat)-lat0)*111.0
    cells['x'], cells['y'] = proj(cells.clat.values, cells.clon.values)
    M = len(cells); cidx = {c: i for i, c in enumerate(cells.cell)}
    Xk = cells[['x','y']].values; Zk = cells.depth_km.values

    # ── 1c. data pairs + Gf ──
    cat = ten[['cell','station','win','clat','clon','sta_lat','sta_lon']].drop_duplicates(
            ['cell','station','win']).reset_index(drop=True)
    cat['a'] = cat.cell.map(cidx); cat['depth'] = cat.cell.map(cells.set_index('cell').depth_km)
    cat['ax'], cat['ay'] = proj(cat.clat.values, cat.clon.values)
    cat['bx'], cat['by'] = proj(cat.sta_lat.values, cat.sta_lon.values)
    WB = {wk: (w1, w2) for wk, w1, w2, _ in WINS}
    Gf = np.zeros((len(cat), M))
    for i, r in cat.iterrows():
        if i % 200 == 0: print(f'  kernel {i}/{len(cat)}...', flush=True)
        w1, w2 = WB[r.win]
        K = kernel_singlescatter(Xk, Zk, np.array([r.ax, r.ay]), r.depth,
                                 np.array([r.bx, r.by]), BETA, w1, w2, ell=LSTAR)
        Gf[i] = K
    Gf = -Gf
    sw    = (cat.station + '|' + cat.win)
    swu   = sorted(sw.unique()); swidx = {x: i for i, x in enumerate(swu)}; B = len(swu)
    Gs_col_arr = sw.map(swidx).values
    Gs_row_arr = np.arange(len(cat))
    Gs = sparse.csr_matrix((np.ones(len(cat)), (Gs_row_arr, Gs_col_arr)), shape=(len(cat), B))
    key = cat[['cell','station','win']].copy(); key['row'] = np.arange(len(cat))
    key_cell = cat['a'].values
    key_sta  = sw.map(swidx).values

    # ── 1d. monthly index lists ──
    months = sorted(ten.ym.unique())
    month_rows, month_dd = [], []
    for E in months:
        te = ten[ten.ym == E][['cell','station','win','dvv']]
        mm = key.merge(te, on=['cell','station','win'])
        if len(mm) < 60 or mm.station.nunique() < 10:
            month_rows.append(np.array([], dtype=int))
            month_dd.append(np.array([], dtype=float))
        else:
            month_rows.append(mm.row.values.astype(int))
            month_dd.append(mm.dvv.values * 100.0)

    # ── 1e. save cache ──
    print(f'[null_test] saving cache to {CACHE}...')
    np.savez(CACHE,
             Gf=Gf, Gs_data=np.ones(len(cat)), Gs_row=Gs_row_arr, Gs_col=Gs_col_arr,
             Gs_shape=np.array([len(cat), B]),
             cells_lat=cells.clat.values, cells_lon=cells.clon.values,
             cells_nsta=cells.n_stations.values,
             months=np.array(months), key_cell=key_cell, key_sta=key_sta,
             sw_labels=np.array(swu),
             M=np.int32(M), B=np.int32(B),
             month_rows=np.array(month_rows, dtype=object),
             month_dd=np.array(month_dd, dtype=object))
    print(f'  saved.  total setup: {time.time()-t0:.0f}s')
    # Alias variables to the names used in the cache-load branch
    Gs_data  = np.ones(len(cat)); Gs_row = Gs_row_arr; Gs_col = Gs_col_arr
    Gs_shape = (len(cat), B)
    cells_lat  = cells.clat.values
    cells_lon  = cells.clon.values
    cells_nsta = cells.n_stations.values
    sw_labels  = swu

# ── Reconstruct Gs as csr (works whether we loaded from cache or just built it) ──────────────
Gs = sparse.csr_matrix((Gs_data, (Gs_row, Gs_col)), shape=Gs_shape)

# ── Laplacian ────────────────────────────────────────────────────────────────────────────────
def graph_laplacian(P, k=6):
    tree = cKDTree(P); m = len(P); rr, cc, vv = [], [], []
    for i in range(m):
        _, nn = tree.query(P[i], k+1)
        for j in nn[1:]:
            rr += [i, i]; cc += [i, j]; vv += [1.0, -1.0]
    return sparse.coo_matrix((vv, (rr, cc)), shape=(m, m)).tocsr()

# Build cell coordinates for Laplacian (need x,y not lat/lon)
lat0_g = np.mean(cells_lat); lon0_g = np.mean(cells_lon)
cx = (cells_lon - lon0_g) * 111.0 * np.cos(np.radians(lat0_g))
cy = (cells_lat - lat0_g) * 111.0
L = graph_laplacian(np.column_stack([cx, cy]))
regf = sparse.hstack([lam_f * L, sparse.csr_matrix((M, B))])
regs = sparse.hstack([sparse.csr_matrix((B, M)), lam_s * sparse.identity(B)])

# ── Step 2: Estimate per-(station, window) noise std from real inversion residuals ──────────────
print('[null_test] estimating per-(station,window) noise from real inversion residuals ...')
z_real = np.load(REAL_NPZ, allow_pickle=True)
MF_real  = z_real['MF']   # (M, n_months_total)
SITE_real = z_real['SITE'] # (B, n_months_total)
months_real = list(z_real['months'])
# Build a month index map: months_real position -> our month index
real_month_idx = {m: i for i, m in enumerate(months_real)}

# For each ok month in the real inversion, compute residuals = dd - (Gf @ mf + Gs @ site)
sw_labels_str = [str(x) for x in sw_labels]
residuals_by_sw = {sw: [] for sw in sw_labels_str}
n_months_used = 0
for j, E in enumerate(months):
    rows = month_rows[j]; dd = month_dd[j]
    if len(rows) == 0: continue
    ri = real_month_idx.get(E, -1)
    if ri < 0: continue
    mf   = MF_real[:, ri]
    site = SITE_real[:, ri]
    if np.any(np.isnan(mf)) or np.any(np.isnan(site)): continue
    pred = Gf[rows] @ mf + Gs[rows] @ site
    resid = dd - pred
    # assign residuals to station|window labels
    sw_idx_rows = key_sta[rows]   # (len(rows),) -> site term index
    for sw_i, r in zip(sw_idx_rows, resid):
        residuals_by_sw[sw_labels_str[sw_i]].append(r)
    n_months_used += 1

print(f'  residuals from {n_months_used} ok months')
# Per-(station,window) noise std
noise_std = np.zeros(B)
for i, sw in enumerate(sw_labels_str):
    rv = np.array(residuals_by_sw[sw])
    if len(rv) > 5:
        noise_std[i] = np.std(rv)
    else:
        noise_std[i] = np.nanmedian([np.std(np.array(residuals_by_sw[s])) for s in sw_labels_str if len(residuals_by_sw[s]) > 5])

print(f'  noise_std per (sta,win): min={noise_std.min():.3f} median={np.median(noise_std):.3f} max={noise_std.max():.3f} %')

# ── Step 3: Null test — N_REAL realisations ──────────────────────────────────────────────────
print(f'[null_test] running {N_REAL} null realisations (zero truth + real coverage + real noise) ...')
m_true = np.zeros(M)   # ZERO truth everywhere

all_MF = []   # each element: (M, n_ok_months)
ok_month_list = [j for j, rows in enumerate(month_rows) if len(rows) > 0]
n_ok = len(ok_month_list)

t_start = time.time()
for real_i in range(N_REAL):
    rng = np.random.RandomState(real_i * 17 + 3)
    MF_null = np.full((M, n_ok), np.nan)
    for jj, j in enumerate(ok_month_list):
        rows = month_rows[j]; dd_true = np.zeros(len(rows))  # zero fault field -> Gf @ 0 = 0
        # noise per datum: use the station|window noise_std
        sw_idx_rows = key_sta[rows]
        noise = rng.randn(len(rows)) * noise_std[sw_idx_rows]
        dd_synth = dd_true + noise
        top  = sparse.hstack([sparse.csr_matrix(Gf[rows]), Gs[rows]])
        Aaug = sparse.vstack([top, regf, regs]).tocsr()
        baug = np.concatenate([dd_synth, np.zeros(M), np.zeros(B)])
        sol  = lsqr(Aaug, baug, atol=1e-7, btol=1e-7, iter_lim=2000)[0]
        MF_null[:, jj] = sol[:M]
    all_MF.append(MF_null)
    elapsed = time.time() - t_start
    eta = elapsed / (real_i + 1) * (N_REAL - real_i - 1)
    print(f'  realisation {real_i+1}/{N_REAL}  elapsed {elapsed/60:.1f}m  ETA {eta/60:.1f}m', flush=True)

# ── Step 4: Compute null temporal RMS per cell per realisation ───────────────────────────────
print('[null_test] computing null RMS statistics ...')
null_cell_rms = np.array([np.nanstd(mf, axis=1) for mf in all_MF])  # (N_REAL, M)
# Real observed temporal RMS (recompute over same ok months to be consistent)
real_ok_months_idx = [j for j, rows in enumerate(month_rows) if len(rows) > 0]
# Use the real inversion's MF at these same months
real_mf_at_ok = np.full((M, n_ok), np.nan)
for jj, j in enumerate(real_ok_months_idx):
    E = months[j]
    ri = real_month_idx.get(E, -1)
    if ri >= 0:
        real_mf_at_ok[:, jj] = MF_real[:, ri]   # may contain NaN; nanstd handles it
obs_cell_rms = np.nanstd(real_mf_at_ok, axis=1)

# Save results
print(f'[null_test] saving results to {RESULT} ...')
np.savez(RESULT,
         null_cell_rms=null_cell_rms,
         obs_cell_rms=obs_cell_rms,
         cells_lat=cells_lat, cells_lon=cells_lon, cells_nsta=cells_nsta,
         noise_std=noise_std, sw_labels=np.array(sw_labels_str),
         N_REAL=np.int32(N_REAL), n_ok_months=np.int32(n_ok))

# ── Step 5: Summary statistics by lat band × coverage class ─────────────────────────────────
print('\n[null_test] === RESULTS ===')
lat_bands  = [(40,42,'40-42'), (42,44,'42-44'), (44,46,'44-46'), (46,48,'46-48'), (48,50,'48-50')]
cov_classes = [(1,1,'1'), (2,2,'2'), (3,3,'3'), (4,99,'>=4')]

summary_rows = []
for lo, hi, bname in lat_bands:
    for n1, n2, cname in cov_classes:
        mask = (cells_lat >= lo) & (cells_lat < hi) & (cells_nsta >= n1) & (cells_nsta <= n2)
        if mask.sum() == 0: continue
        obs  = obs_cell_rms[mask]
        null = null_cell_rms[:, mask]          # (N_REAL, n_cells_in_class)
        null_flat = null.ravel()               # pool all realisations + cells in class
        obs_med  = np.median(obs)
        obs_90th = np.percentile(obs, 90)
        null_med = np.median(null_flat)
        null_95th = np.percentile(null_flat, 95)
        null_99th = np.percentile(null_flat, 99)
        null_max  = null_flat.max()
        # Per-cell percentile: where does observed sit in its null distribution?
        per_cell_pctile = np.array([
            100 * np.mean(null[:, k] <= obs[k]) for k in range(mask.sum())])
        pctile_med = np.median(per_cell_pctile)
        row = dict(band=bname, nsta=cname, n_cells=mask.sum(),
                   obs_median=obs_med, obs_90th=obs_90th,
                   null_median=null_med, null_95th=null_95th, null_99th=null_99th, null_max=null_max,
                   obs_pctile_in_null_med=pctile_med)
        summary_rows.append(row)
        print(f'  {bname} nsta={cname:3s}: n={mask.sum():3d}  '
              f'obs_med={obs_med:.3f}%  null_95={null_95th:.3f}%  null_max={null_max:.3f}%  '
              f'obs_pctile_in_null={pctile_med:.0f}th')

# ── Step 6: Figure ────────────────────────────────────────────────────────────────────────────
print('[null_test] making figure ...')
fig = plt.figure(figsize=(16, 12))
gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

band_axes = {
    '40-42': fig.add_subplot(gs[0, 0]),
    '42-44': fig.add_subplot(gs[0, 1]),
    '44-46': fig.add_subplot(gs[0, 2]),
    '46-48': fig.add_subplot(gs[1, 0]),
    '48-50': fig.add_subplot(gs[1, 1]),
}
cov_colors = {'1': '#d62728', '2': '#ff7f0e', '3': '#2ca02c', '>=4': '#1f77b4'}
cov_labels = {'1': '1 sta', '2': '2 sta', '3': '3 sta', '>=4': '≥4 sta'}

for lo, hi, bname in lat_bands:
    ax = band_axes[bname]
    for n1, n2, cname in cov_classes:
        mask = (cells_lat >= lo) & (cells_lat < hi) & (cells_nsta >= n1) & (cells_nsta <= n2)
        if mask.sum() == 0: continue
        obs  = obs_cell_rms[mask]
        null = null_cell_rms[:, mask].ravel()
        col  = cov_colors[cname]
        # KDE-like histogram of null
        bins = np.linspace(0, max(null.max(), obs.max()) * 1.1 + 0.01, 60)
        ax.hist(null, bins=bins, density=True, alpha=0.3, color=col, label=f'null {cov_labels[cname]}')
        # Observed values as vertical ticks
        for ov in obs:
            ax.axvline(ov, color=col, lw=1.2, alpha=0.7)
        # Mark median observed
        ax.axvline(np.median(obs), color=col, lw=2.5, ls='--')
    ax.set_title(f'Lat {bname}°N', fontsize=11, fontweight='bold')
    ax.set_xlabel('temporal RMS (%)', fontsize=9)
    ax.set_ylabel('density', fontsize=9)
    ax.legend(fontsize=7, loc='upper right')
    ax.tick_params(labelsize=8)

# Summary text panel
ax_txt = fig.add_subplot(gs[1, 2])
ax_txt.axis('off')
lines = ['NULL TEST VERDICT', '---------------------------------', '']
for r in summary_rows:
    flag = '(!!)' if r['obs_pctile_in_null_med'] > 95 else '(ok)'
    lines.append(f"{flag} {r['band']} n={r['nsta']:3s}: obs_med={r['obs_median']:.3f}% "
                 f"null_95={r['null_95th']:.3f}%  pctile={r['obs_pctile_in_null_med']:.0f}th")
ax_txt.text(0.02, 0.98, '\n'.join(lines), transform=ax_txt.transAxes,
            fontsize=7, va='top', family='monospace')

fig.suptitle('NULL TEST: zero-truth inversion with real coverage + real noise levels\n'
             'Histograms = null RMS distribution; dashed lines = observed median; '
             'ticks = per-cell observed values', fontsize=11)
outfig = 'figures/null_test_southern.png'
plt.savefig(outfig, dpi=130, bbox_inches='tight')
print(f'saved {outfig}')

# ── Step 7: Append to PHASE_A_RESULTS.md ──────────────────────────────────────────────────────
results_block = f"""
## NULL TEST: Is the southern "hot zone" temporal RMS real or a coverage artifact? (2026-06-10)

**Question:** The real multi-window inversion shows elevated temporal RMS in the sparse south
(42–44°N: median 0.09–0.15%, 1–2 stations; 40–42°N: median 0.29%, 1–2 stations).  Does this
reflect genuine fault velocity variability, or is it noise amplified through a sparse design matrix?

**Method:**  `fault_tomography/inversion/null_test.py`
- Planted fault field: ZERO everywhere, all months.
- Forward-modelled every real (cell, station, window, month) datum using the real kernels and
  real coverage pattern.
- Added synthetic noise drawn per datum from the REAL per-(station, window) residual std, estimated
  by pooling residuals from the real inversion across all solved months
  (noise_std range: {noise_std.min():.3f}–{noise_std.max():.3f}%, median {np.median(noise_std):.3f}%).
- Inverted with identical joint LSQR + Laplacian + site-term machinery (λ_f={lam_f}, λ_s={lam_s}).
- Repeated {N_REAL} realisations (different RNG seeds).

**Noise std used (per station|window pair):**
- Estimated from real inversion residuals (dd − Ĝf mf − Ĝs site), pooled across all ok months.
- Min: {noise_std.min():.3f}%, Median: {np.median(noise_std):.3f}%, Max: {noise_std.max():.3f}%.

**Results — null temporal RMS vs observed, by latitude band and station coverage:**

| Lat band | nsta | n cells | obs median (%) | null 95th pct (%) | null max (%) | obs pctile in null |
|---|---|---|---|---|---|---|
"""

for r in summary_rows:
    results_block += (f"| {r['band']} | {r['nsta']} | {r['n_cells']} | "
                      f"{r['obs_median']:.3f} | {r['null_95th']:.3f} | {r['null_max']:.3f} | "
                      f"{r['obs_pctile_in_null_med']:.0f}th |\n")

results_block += f"""
**Verdict:**

"""
# Generate verdict based on results
verdict_lines = []
for r in summary_rows:
    lo_str, hi_str = r['band'].split('-')
    lo_f, hi_f = float(lo_str), float(hi_str)
    is_sparse = (r['nsta'] in ['1','2']) and (hi_f <= 44)
    if not is_sparse: continue
    if r['obs_pctile_in_null_med'] > 95:
        verdict_lines.append(
            f"- **{r['band']} {r['nsta']}-sta (n={r['n_cells']}):** obs median {r['obs_median']:.3f}% "
            f"sits at {r['obs_pctile_in_null_med']:.0f}th pctile of the null — "
            f"SIGNIFICANT (exceeds null 95th = {r['null_95th']:.3f}%). Likely real signal or unmodelled systematic.")
    else:
        verdict_lines.append(
            f"- **{r['band']} {r['nsta']}-sta (n={r['n_cells']}):** obs median {r['obs_median']:.3f}% "
            f"sits at {r['obs_pctile_in_null_med']:.0f}th pctile of the null "
            f"(null 95th = {r['null_95th']:.3f}%, null max = {r['null_max']:.3f}%). "
            f"CANNOT BE DISTINGUISHED from noise amplification — likely a COVERAGE ARTIFACT.")

results_block += '\n'.join(verdict_lines)
results_block += f"""

**Overall conclusion:**
The null test directly answers the open question from the synthetic resolution test (§4, 2026-06-09).
Cells where the observed temporal RMS lies within the null distribution (obs_pctile ≤ 95th) are
consistent with pure noise amplification through the sparse design matrix and should NOT be
interpreted as real fault velocity change. Cells/classes where the observed RMS significantly
exceeds the null distribution warrant further investigation.

Figure: `figures/null_test_southern.png`
Result archive: `{RESULT}`
"""

with open('fault_tomography/PHASE_A_RESULTS.md', 'a') as fh:
    fh.write(results_block)
print('[null_test] appended to fault_tomography/PHASE_A_RESULTS.md')
print('[null_test] DONE.')
