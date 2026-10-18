"""Q8 RTR debug iteration D2: back-zone surge detector.

Per episode: onset map (4-km bins, 10th-percentile time). At each 30-min bin tau the SLIPPED zone = bins whose onset is
earlier than tau - LAG_H; the FRONT = median position of bins with onset in [tau - 24 h, tau]. Activity in the slipped
zone is counted per 30-min bin and compared with its own +-24 h mean (same slipped set). Surge bins (Poisson p < P_BIN
and count >= 3) are merged (gap <= 1 bin) into runs; each run's slipped-zone events give speed (Theil-Sen), span
(10-90 %), direction. RTR = run lasting >= 1.5 h, span >= 8 km, |v| >= 1 km/h, moving AWAY from the front.
Score S = -log10 Poisson P(run count | run background).
"""
import numpy as np
import pandas as pd
from scipy import stats

ONSET_BIN = 4.0
ONSET_MIN_EV = 5
DTB = 0.5           # h
LAG_H = 12.0
HALF_BG = 24.0
P_BIN = 0.01
MIN_DUR = 1.5
MIN_SPAN = 8.0
MIN_V = 1.0


def onset_map(t, x):
    b = np.floor(x / ONSET_BIN).astype(int)
    df = pd.DataFrame(dict(b=b, t=t))
    g = df.groupby("b").t
    on = g.quantile(0.10)[g.size() >= ONSET_MIN_EV]
    return on


def detect(t, x, T_end, unit=None, lag_h=LAG_H, main_dir=None, onset_override=None):
    on = onset_map(t, x) if onset_override is None else onset_override
    bins = on.index.values; onv = on.values
    eb = np.floor(x / ONSET_BIN).astype(int)
    nb = int(np.ceil(T_end / DTB))
    it = np.clip((t / DTB).astype(int), 0, nb - 1)
    # count matrix: time bin x onset bin
    bidx = {b: i for i, b in enumerate(bins)}
    ib = np.array([bidx.get(b, -1) for b in eb])
    ok = ib >= 0
    Cm = np.zeros((nb, len(bins))); np.add.at(Cm, (it[ok], ib[ok]), 1)
    taus = (np.arange(nb) + 0.5) * DTB
    slipped = onv[None, :] < (taus[:, None] - lag_h)              # (nb, nbins)
    c = (Cm * slipped).sum(1)
    hw = int(HALF_BG / DTB)
    cs = np.vstack([np.zeros(len(bins)), np.cumsum(Cm, 0)])
    bg = np.zeros(nb)
    for i in range(nb):
        a, b = max(0, i - hw), min(nb, i + hw + 1)
        bg[i] = ((cs[b] - cs[a]) * slipped[i]).sum() / (b - a)
    p = stats.poisson.sf(c - 1, np.maximum(bg, 0.05))
    surge = (c >= 3) & (p < P_BIN)
    runs = []
    i = 0
    while i < nb:
        if not surge[i]:
            i += 1; continue
        j = i
        while j + 1 < nb and (surge[j + 1] or (j + 2 < nb and surge[j + 2])):
            j += 1 if surge[j + 1] else 2
        runs.append((i, j)); i = j + 1
    out = []
    for i, j in runs:
        t0, t1 = i * DTB, (j + 1) * DTB
        if t1 - t0 < MIN_DUR:
            continue
        tau_mid = (i + j + 1) / 2 * DTB
        sl = set(bins[onv < (t0 + DTB / 2 - lag_h)])
        m = (t >= t0) & (t < t1) & np.isin(eb, list(sl))
        if m.sum() < 5:
            continue
        te, xe = t[m], x[m]
        v = stats.theilslopes(xe, te)[0]
        span = np.percentile(xe, 90) - np.percentile(xe, 10)
        fb = bins[(onv >= t0 - 24) & (onv <= t0)]
        front = (np.median(fb) + 0.5) * ONSET_BIN if len(fb) else np.nan
        away = np.sign(v) == np.sign(np.median(xe) - front) if np.isfinite(front) else False
        n_run = c[i:j + 1].sum(); b_run = bg[i:j + 1].sum()
        S = -stats.poisson.logsf(n_run - 1, max(b_run, 0.05)) / np.log(10)
        nf = len(np.unique(unit[m])) if unit is not None else np.nan
        rho = stats.spearmanr(te, xe)[0]
        k = max(len(te) // 4, 2)
        kin = span >= MIN_SPAN and abs(v) >= MIN_V
        opp = (main_dir is not None) and (np.sign(v) == -np.sign(main_dir))
        cls = "RTR" if (kin and opp) else ("forward" if kin else "surge_only")
        # slipped-zone positions in the run (reference for "starts at the front")
        ref_d = np.abs(xe - front) if np.isfinite(front) else np.array([np.nan])
        out.append(dict(t0=t0, t1=t1, t_first=te.min(), t_last=te.max(), N_ev=int(m.sum()), n_run=n_run, bg_run=b_run, S=min(S, 999.0),
                        v_ev=v, span_ev=span, rho=rho, x_start=float(np.median(xe[np.argsort(te)][:k])),
                        x_end=float(np.median(xe[np.argsort(te)][-k:])), front=front, n_fam=nf, cls=cls, away=bool(away), ref_front_dist_med=float(np.median(ref_d)),
                        idx=np.nonzero(m)[0]))
    return out


def main_direction(t, x):
    """Theil-Sen of daily median position (km/day) over days with >= 1 % of the episode's events."""
    d = pd.DataFrame(dict(day=np.floor(t / 24.0), x=x))
    g = d.groupby("day").x; n = g.size(); m = g.median()[n >= 0.01 * len(x)]
    if len(m) < 3:
        return np.nan
    return float(stats.theilslopes(m.values, m.index.values)[0])
