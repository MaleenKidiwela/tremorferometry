#!/usr/bin/env python
"""Amendment Q1: station-era table (days on disk, broadband days, >= 180 broadband days, data on the 14 template-event
days). -> tables/era_station_year.csv, tables/era_event_days.csv"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import s13_common as C

STAS = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'CN.LZB', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2',
        'TA.A04D', 'TA.B05D', 'TA.D03D', 'C8.SHVB', 'C8.KHVB', 'C8.MGCB', 'C8.TWBB', 'C8.SHDB', 'UW.SQM', 'UW.HDW', 'UW.OTR',
        'UW.WISH', 'XU.BS11', 'XU.W040', 'XU.W020', 'XU.PL11', 'XU.W060', 'US.NLWA', 'CC.OBSR', 'TA.D04D', 'UW.ELW', 'UW.LEBA', 'UW.RATT']


def bb(chans):
    return any(c[:2] in ('BH', 'HH') for c in str(chans).split(','))


def main():
    d = pd.read_csv(f'{K.S13}/tables/inventory_station_year.csv')
    d = d[d.year.str.match(r'^\d{4}$') & d.sta.isin(STAS)].copy(); d['year'] = d.year.astype(int)
    d['bb'] = d.chans.apply(bb)
    rows = []
    for (s, y), g in d.groupby(['sta', 'year']):
        rows.append(dict(sta=s, year=y, days_any=int(g.n_files.max()), days_bb=int(g[g.bb].n_files.max()) if g.bb.any() else 0,
                         chans=';'.join(sorted(set(g.chans.astype(str))))))
    T = pd.DataFrame(rows)
    T['bb_ge180'] = T.days_bb >= 180
    T.to_csv(f'{K.TB}/era_station_year.csv', index=False)
    W = T.pivot_table(index='sta', columns='year', values='days_bb', fill_value=0).reindex(STAS).fillna(0).astype(int)
    W.to_csv(f'{K.TB}/era_station_year_bbdays_wide.csv')
    # event days
    from obspy import read
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    ev = []
    for r in V.itertuples():
        day = str(pd.Timestamp(r.epoch, unit='s').date())
        o = dict(cand_id=r.cand_id, role=r.role, date=day)
        for s in STAS:
            fs = C.day_files(s, day)
            band = ''
            for f in fs:
                try:
                    st = read(f, headonly=True)
                    bands = {tr.stats.channel[:2] for tr in st}
                    if bands & {'BH', 'HH'}:
                        band = 'BB'; break
                    band = band or '/'.join(sorted(bands))
                except Exception:
                    continue
            o[s] = band or '-'
        ev.append(o)
    E = pd.DataFrame(ev)
    E.to_csv(f'{K.TB}/era_event_days.csv', index=False)
    vet = E[E.role == 'vetted']
    n_bb = {s: int((vet[s] == 'BB').sum()) for s in STAS}
    print(W.to_string())
    print('broadband on the 8 vetted-event days:', n_bb)


if __name__ == '__main__':
    main()
