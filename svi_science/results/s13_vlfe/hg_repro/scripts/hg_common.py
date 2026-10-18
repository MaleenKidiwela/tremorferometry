"""hg_repro shared constants and loaders (rules: ../PLAN.md)."""
import os
import numpy as np
import pandas as pd

ROOT = '/home/jovyan/tremorferometry'
S13 = f'{ROOT}/svi_science/results/s13_vlfe'
HG = f'{S13}/hg_repro'
TB, FIG, LOGD, MF = f'{HG}/tables', f'{HG}/figures', f'{HG}/logs', f'{HG}/data/mf'
DISP_HG, DISP_V3 = f'{HG}/data/disp1sps', f'{S13}/data/disp1sps_v3'
STAS = ['CN.NLLB', 'CN.VGZ', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'TA.A04D', 'TA.B05D', 'TA.D03D']
COMPS = ['N', 'E', 'Z']
CHANS = [f'{s}.{c}' for s in STAS for c in COMPS]          # 27 channels, index = 3*station + comp
T0 = pd.Timestamp('2011-06-01')                              # sample 0 = 2011-06-01 00:00:00 UTC, 1 sps
T1 = pd.Timestamp('2013-12-01')
NDAY = (T1 - T0).days                                        # 914
N = NDAY * 86400
L = 90                                                       # template length (s) = CMT window
M = N - L + 1                                                # number of CC lags
BAND = (0.02, 0.05)
SEED = 20261015
# paper periods (UTC dates, [start, end))
PER = {'ETS2011': ('2011-07-23', '2011-09-06'), 'interETS': ('2011-09-06', '2012-08-30'),
       'interETS_b': ('2011-10-01', '2012-08-30'), 'ETS2012': ('2012-08-30', '2012-10-11'),
       'pre2011': ('2011-07-01', '2011-07-23'), 'post2012': ('2012-10-11', '2012-12-01'),
       'all18': ('2011-07-01', '2012-12-01'), 'win2011': ('2011-07-01', '2011-10-01')}


def tidx(t):
    """Sample index of a UTC time."""
    return int(round((pd.Timestamp(t) - T0).total_seconds()))


def itime(i):
    return T0 + pd.to_timedelta(np.asarray(i), unit='s')


def per_idx(name):
    a, b = PER[name]
    return tidx(a), tidx(b)


def templates():
    T = pd.read_csv(f'{S13}/tables/hg2018_table1.csv', comment='#')
    T['time'] = pd.to_datetime(T.time)
    return T[T.template <= 5].reset_index(drop=True)


def day_file(sta, d):
    rel = f'{sta}/{d.year}/{d.dayofyear:03d}.npz'
    for root in (DISP_V3, DISP_HG):
        p = f'{root}/{rel}'
        if os.path.isfile(p):
            return p
    return None


def tz_pt(t):
    return pd.Timestamp(t).tz_localize('UTC').tz_convert('America/Los_Angeles')
