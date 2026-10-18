"""S5 shared machinery: tidal phases/stresses at event times, whole-day-shift surrogates, grouped statistics.
Rules: results/s5/PLAN_S5.md.
"""
import numpy as np, pandas as pd
import pyTMD.constituents as PC
from s5_tides_common import design

NS_H = np.int64(3600 * 10**9)
NS_D = np.int64(86400 * 10**9)
FIT_CONSTS = ['q1', 'o1', 'p1', 'k1', 'n2', 'm2', 's2', 'k2', '2n2', 'm4', 'j1', 'oo1', 'mu2', 'nu2', 'l2', 't2',
              '2q1', 'rho1', 'm1', 'mf', 'mm']
PH_CONSTS = ['m2', 'o1', 'k1', 's2', 'n2']
SPEED_DEG_H = {'m2': 28.9841042, 'o1': 13.9430356, 'k1': 15.0410686, 's2': 30.0, 'n2': 28.4397295}
_EPOCH = np.datetime64('1992-01-01T00:00:00', 'ns')


def mjd(times):
    return (np.asarray(times, dtype='datetime64[ns]') - _EPOCH) / np.timedelta64(1, 'D') + 48622.0


def astro_th(times, consts=PH_CONSTS):
    """theta_astro = V + u (radians) per constituent, shape (n, nc)."""
    pu, pf, G = PC.arguments(mjd(times), consts, corrections='FES')
    return np.radians(G) + pu


class Stress:
    """Hourly stress arrays per family + harmonic constants for phases."""

    def __init__(self, df, comps):
        df = df.sort_values(['fid', 'time'])
        self.fids = np.sort(df.fid.unique())
        self.fidx = {f: i for i, f in enumerate(self.fids)}
        t = df.time[df.fid == self.fids[0]].values.astype('datetime64[ns]')
        self.t0 = t[0]; self.nh = len(t)
        assert np.all(np.diff(t) == np.timedelta64(1, 'h'))
        self.A = {c: df[c].values.reshape(len(self.fids), self.nh).astype(np.float32) for c in comps}
        self.std = {c: self.A[c].std(1) for c in comps}
        # harmonic fit of each component per family -> arg(hc) per phase constituent
        X = np.c_[design(t, FIT_CONSTS), np.ones(len(t))]
        pinv = np.linalg.pinv(X)
        n = len(FIT_CONSTS)
        self.hc = {}
        for c in comps:
            cf = self.A[c].astype(float) @ pinv.T                    # (nf, 2n+1)
            self.hc[c] = {k: cf[:, i] + 1j * cf[:, n + i] for i, k in enumerate(FIT_CONSTS)}
        # daily u (nodal angle) grid for exact shifted phases
        days = np.arange(t[0].astype('datetime64[D]') - 4000, t[-1].astype('datetime64[D]') + 4000, dtype='datetime64[D]')
        self.day0 = days[0].astype('datetime64[ns]')
        pu, pf, G = PC.arguments(mjd(days.astype('datetime64[ns]')), PH_CONSTS, corrections='FES')
        self.u_day = pu                                               # (ndays, nc) radians

    def value(self, comp, fi, tns):
        """linear interpolation of hourly comp at int64 ns times for family indices fi."""
        h = (tns - self.t0.astype(np.int64)) / NS_H
        h0 = np.floor(h).astype(np.int64); w = (h - h0).astype(np.float32)
        ok = (h0 >= 0) & (h0 + 1 < self.nh)
        h0c = np.clip(h0, 0, self.nh - 2)
        A = self.A[comp]
        v = A[fi, h0c] * (1 - w) + A[fi, h0c + 1] * w
        v[~ok] = np.nan
        return v


class EventSet:
    """events: DataFrame with fid, t (datetime64[ns]); optional label columns (int codes, -1 = excluded)."""

    def __init__(self, ev, S, comps, phase_comp='tau_Pa'):
        ev = ev[ev.fid.isin(S.fids)].reset_index(drop=True)
        self.ev = ev; self.S = S; self.comps = comps
        self.fi = ev.fid.map(S.fidx).values.astype(np.int64)
        self.tns = ev.t.values.astype('datetime64[ns]').astype(np.int64)
        self.day = (self.tns // NS_D).astype(np.int64)
        th = astro_th(ev.t.values)
        self.th = {}
        for j, c in enumerate(PH_CONSTS):
            self.th[c] = th[:, j]
        self.phase_comp = phase_comp
        # family arg(hc) for each phase constituent (theta = th_astro + arg(hc))
        self.arg = {c: np.angle(S.hc[phase_comp][c])[self.fi] for c in PH_CONSTS}
        self.arg_sn = {c: np.angle(S.hc['sigma_n_Pa'][c])[self.fi] for c in ('m2', 'o1')} if 'sigma_n_Pa' in S.hc else {}
        # (family, day) blocks
        key = self.fi * 10**6 + (self.day - self.day.min())
        self.pair_key, self.pair_idx = np.unique(key, return_inverse=True)
        self.pair_f = self.pair_key // 10**6
        self.pair_day = self.pair_key % 10**6 + self.day.min()
        nday = np.bincount(self.pair_f, minlength=len(S.fids))
        self.shiftable = nday[self.pair_f] >= 3

    def quantities(self, delta_days):
        """per-event quantity arrays at times shifted by delta_days (int array, per event)."""
        q = {}
        tns = self.tns + delta_days.astype(np.int64) * NS_D
        di = (self.day + delta_days - self.S.day0.astype(np.int64) // NS_D).astype(np.int64)
        d0 = (self.day - self.S.day0.astype(np.int64) // NS_D).astype(np.int64)
        for j, c in enumerate(PH_CONSTS):
            w = np.radians(SPEED_DEG_H[c]) * 24.0
            th = self.th[c] + w * delta_days + (self.S.u_day[di, j] - self.S.u_day[d0, j])
            ph = th + self.arg[c]
            q[f'cos_{c}'] = np.cos(ph); q[f'sin_{c}'] = np.sin(ph)
            if c in self.arg_sn:
                ps = th + self.arg_sn[c]
                q[f'cos_sn_{c}'] = np.cos(ps); q[f'sin_sn_{c}'] = np.sin(ps)
        for c in self.comps:
            q[c] = self.S.value(c, self.fi, tns)
            q[c + '_z'] = q[c] / self.S.std[c][self.fi]
        return q

    def surrogate_deltas(self, rng):
        """whole-day shift: permute each family's active days; events move with their day."""
        npair = len(self.pair_key)
        r = rng.random(npair)
        order = np.lexsort((r, self.pair_f))           # within-family random order
        base = np.lexsort((self.pair_day, self.pair_f))  # within-family day order (pairs already sorted by key)
        newday = np.empty(npair, np.int64)
        newday[base] = self.pair_day[order]
        d = newday - self.pair_day
        d[~self.shiftable] = 0
        return d[self.pair_idx]


def grouped(q, labels, nfull=None):
    """sums and counts of each quantity per label code; labels: dict name -> (codes int array, ncodes)."""
    out = {}
    for ln, (codes, nc) in labels.items():
        ok = codes >= 0
        cnt = np.bincount(codes[ok], minlength=nc).astype(float)
        d = {'n': cnt}
        for k, v in q.items():
            vv = v[ok]; good = np.isfinite(vv)
            d[k] = np.bincount(codes[ok][good], weights=vv[good], minlength=nc)
            if not (k.startswith('cos_') or k.startswith('sin_')):
                d[k + '_n'] = np.bincount(codes[ok][good], minlength=nc).astype(float)
                d[k + '_sq'] = np.bincount(codes[ok][good], weights=vv[good] ** 2, minlength=nc)
        out[ln] = d
    return out


def _hists(q, labels, specs):
    out = {}
    for nm, (fn, edges, ln) in specs.items():
        codes, nc = labels[ln]
        v = fn(q); ok = (codes >= 0) & np.isfinite(v)
        b = np.clip(np.digitize(v[ok], edges) - 1, 0, len(edges) - 2)
        out[nm] = np.bincount(codes[ok] * (len(edges) - 1) + b, minlength=nc * (len(edges) - 1)).reshape(nc, -1)
    return out


def run(es, labels, nsurr=1000, seed=5, hist_specs=None):
    """observed + surrogate grouped sums. returns obs, surr (label -> qty -> (nsurr, nc)), hist obs, hist surr mean."""
    rng = np.random.default_rng(seed)
    hist_specs = hist_specs or {}
    q0 = es.quantities(np.zeros(len(es.ev), np.int64))
    obs = grouped(q0, labels)
    hobs = _hists(q0, labels, hist_specs)
    hsum = {k: np.zeros(v.shape) for k, v in hobs.items()}
    surr = {ln: {k: np.zeros((nsurr, len(v))) for k, v in d.items()} for ln, d in obs.items()}
    for s in range(nsurr):
        q = es.quantities(es.surrogate_deltas(rng))
        g = grouped(q, labels)
        for ln, d in g.items():
            for k, v in d.items():
                surr[ln][k][s] = v
        for k, v in _hists(q, labels, hist_specs).items():
            hsum[k] += v
    return obs, surr, hobs, {k: v / nsurr for k, v in hsum.items()}


def summarize(obs, surr, ln, names=None, const='m2'):
    """per label code: n, C (obs mean cos), C_surr mean/sd, excess, a(%) and theta0, p one-sided."""
    o, s = obs[ln], surr[ln]
    rows = []
    nc = len(o['n'])
    for i in range(nc):
        n = o['n'][i]
        if n == 0:
            continue
        r = dict(label=ln, code=i, name=names[i] if names is not None else i, n=int(n))
        for c in ([const] if isinstance(const, str) else const):
            Co, So = o[f'cos_{c}'][i] / n, o[f'sin_{c}'][i] / n
            Cs, Ss = s[f'cos_{c}'][:, i] / n, s[f'sin_{c}'][:, i] / n
            ex, ey = Co - Cs.mean(), So - Ss.mean()
            r.update({f'C_{c}': Co, f'S_{c}': So, f'Csurr_{c}': Cs.mean(), f'Csurr_sd_{c}': Cs.std(),
                      f'Cex_{c}': ex, f'a_{c}_pct': 200 * np.hypot(ex, ey), f'theta0_{c}_deg': np.degrees(np.arctan2(ey, ex)),
                      f'p_C_{c}': (1 + np.sum(Cs >= Co)) / (1 + len(Cs)),
                      f'p2_C_{c}': (1 + np.sum(np.abs(Cs - Cs.mean()) >= abs(Co - Cs.mean()))) / (1 + len(Cs))})
        for comp in ('tau_Pa', 'sigma_n_Pa', 'tau_updip_Pa', 'coulomb_0p4', 'body_tau_Pa', 'ocean_tau_Pa'):
            if comp not in o:
                continue
            no = o[comp + '_n'][i]
            mo = o[comp][i] / no
            ms = s[comp][:, i] / np.maximum(s[comp + '_n'][:, i], 1)
            var_s = (s[comp + '_sq'][:, i] / np.maximum(s[comp + '_n'][:, i], 1) - ms ** 2).mean()
            dt = mo - ms.mean()
            short = comp.replace('_Pa', '')
            r.update({f'dmean_{short}_Pa': dt, f'dmean_{short}_sd': ms.std(),
                      f'p_{short}': (1 + np.sum(ms >= mo)) / (1 + len(ms)),
                      f'asig_{short}_kPa': (var_s / dt / 1e3) if dt != 0 else np.nan})
        rows.append(r)
    return pd.DataFrame(rows)


def diff_test(obs, surr, ln, i, j, const='m2'):
    """difference of C between codes i and j: obs, surrogate sd, two-sided p, MDD (2.8 sd)."""
    o, s = obs[ln], surr[ln]
    ci = o[f'cos_{const}'][i] / o['n'][i]; cj = o[f'cos_{const}'][j] / o['n'][j]
    si = s[f'cos_{const}'][:, i] / o['n'][i]; sj = s[f'cos_{const}'][:, j] / o['n'][j]
    ei = ci - si.mean(); ej = cj - sj.mean()
    d = ei - ej
    ds = (si - si.mean()) - (sj - sj.mean())
    p = (1 + np.sum(np.abs(ds) >= abs(d))) / (1 + len(ds))
    return dict(Cex_i=ei, Cex_j=ej, diff=d, sd=ds.std(), p_two=p, MDD_C=2.8 * ds.std(), MDD_a_pct=200 * 2.8 * ds.std(),
                ratio=ei / ej if ej != 0 else np.nan)


def rise_perm(es, dt_days, fe_codes, q_cos, first=(-1e9, 1.0), late=(2.0, 6.0), nperm=1000, seed=7):
    """P2: C(late) - C(first day) with whole-day permutation of day positions within family-episode.
    dt_days: per-event days since onset; fe_codes: family-episode code (-1 = excluded); q_cos: per-event cos(theta)."""
    rng = np.random.default_rng(seed)
    ok = fe_codes >= 0
    fe = fe_codes[ok]; dt = dt_days[ok]; cs = q_cos[ok]; day = es.day[ok]
    key = fe.astype(np.int64) * 10**6 + (day - day.min())
    pk, pidx = np.unique(key, return_inverse=True)
    pfe = pk // 10**6; pday = pk % 10**6
    # keep only family-episodes with >= 2 days
    nd = np.bincount(pfe)
    use = nd[pfe[pidx]] >= 2
    def stat(shift_days):
        d2 = dt + shift_days
        f = (d2 >= first[0]) & (d2 < first[1]) & use; l = (d2 >= late[0]) & (d2 < late[1]) & use
        if f.sum() == 0 or l.sum() == 0:
            return np.nan, 0, 0
        return cs[l].mean() - cs[f].mean(), int(f.sum()), int(l.sum())
    obs, nf, nl = stat(np.zeros(len(dt)))
    base = np.lexsort((pday, pfe))
    null = []
    for _ in range(nperm):
        order = np.lexsort((rng.random(len(pk)), pfe))
        newday = np.empty(len(pk), np.int64); newday[base] = pday[order]
        null.append(stat((newday - pday)[pidx])[0])
    null = np.array(null)
    return dict(diff=obs, n_first=nf, n_late=nl, perm_mean=np.nanmean(null), perm_sd=np.nanstd(null),
                p_one=(1 + np.sum(null >= obs)) / (1 + np.isfinite(null).sum()))


def rayleigh_p(cos, sin):
    n = len(cos); R = np.hypot(cos.mean(), sin.mean()); Z = n * R ** 2
    return float(np.exp(-Z) * (1 + (2 * Z - Z ** 2) / (4 * n) - (24 * Z - 132 * Z ** 2 + 76 * Z ** 3 - 9 * Z ** 4) / (288 * n ** 2))), float(R)
