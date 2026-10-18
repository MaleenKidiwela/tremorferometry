"""S5 setup: tide-gauge harmonic constants around SVI (validation + inland-water fill).

* NOAA CO-OPS harmonic constants (public API, Greenwich phase lag, metres).
* Canadian Hydrographic Service (DFO IWLS public API): hourly tide PREDICTIONS (wlp) for 2023,
  harmonically analysed here (same nodal convention as the model prediction, pyTMD 'FES').
* Convention check: the same analysis applied to NOAA hourly predictions at Port Angeles 2023
  must return NOAA's published constants.
Output: results/s5_setup/gauge_constants.csv, gauge_convention_check.csv
"""
import json, time, urllib.request
import numpy as np, pandas as pd
from s5_tides_common import OUT, CENTER, MAJOR, harmonic_fit, gc_dist_az

CONSTS = MAJOR + ['m4']
FIT = CONSTS + ['j1', 'oo1', 'mu2', 'nu2', 'l2', 't2', 'm3', 'ms4', 'm6', '2q1', 'rho1', 'm1', 'eps2', 'lambda2', 'mn4', 'mf', 'mm', 'sa', 'ssa']
RMAX = 330.0


def get(url, tries=3):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            err = e; time.sleep(2)
    raise err


def noaa():
    st = get('https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations.json?type=harcon')['stations']
    rows = []
    for s in st:
        d, _ = gc_dist_az(*CENTER, s['lat'], s['lng'])
        if d > RMAX or (s['lat'] < 46.35 and s['lng'] > -124.05):   # Columbia River estuary excluded
            continue
        try:
            hc = get(f"https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/{s['id']}/harcon.json?units=metric")
        except Exception as e:
            print('NOAA fail', s['id'], e); continue
        for c in hc.get('HarmonicConstituents', []):
            name = c['name'].lower().replace('2n2', '2n2')
            if name in CONSTS and c['amplitude'] is not None:
                rows.append(dict(source='NOAA', id=s['id'], name=s['name'], lat=s['lat'], lon=s['lng'], dist_km=d,
                                 const=name, amp_m=c['amplitude'], phase_deg=c['phase_GMT'], fit_resid_m=np.nan))
    return pd.DataFrame(rows)


def fetch_dfo_wlp(sid, year=2023):
    vals = []
    for m in range(1, 13):
        a = f'{year}-{m:02d}-01T00:00:00Z'
        b = f'{year + (m == 12)}-{(m % 12) + 1:02d}-01T00:00:00Z'
        d = get(f'https://api-iwls.dfo-mpo.gc.ca/api/v1/stations/{sid}/data?time-series-code=wlp&from={a}&to={b}&resolution=SIXTY_MINUTES')
        vals += [(x['eventDate'], x['value']) for x in d]
    df = pd.DataFrame(vals, columns=['t', 'v']).drop_duplicates('t')
    t = pd.to_datetime(df.t).dt.tz_localize(None).values.astype('datetime64[ns]')
    return t, df.v.values.astype(float)


def dfo():
    st = get('https://api-iwls.dfo-mpo.gc.ca/api/v1/stations')
    rows = []
    for s in st:
        d, _ = gc_dist_az(*CENTER, s['latitude'], s['longitude'])
        codes = [x['code'] for x in s.get('timeSeries', [])]
        if d > RMAX or 'wlp' not in codes:
            continue
        try:
            t, v = fetch_dfo_wlp(s['id'])
        except Exception as e:
            print('DFO fail', s['code'], s['officialName'], e); continue
        if len(t) < 8000:
            print('DFO short', s['code'], s['officialName'], len(t)); continue
        hc, res = harmonic_fit(t, v, FIT)
        for c in CONSTS:
            rows.append(dict(source='CHS', id=s['code'], name=s['officialName'], lat=s['latitude'], lon=s['longitude'],
                             dist_km=d, const=c, amp_m=abs(hc[c]), phase_deg=(-np.degrees(np.angle(hc[c]))) % 360,
                             fit_resid_m=res))
        print('DFO', s['code'], s['officialName'], f"M2 {abs(hc['m2']):.3f} K1 {abs(hc['k1']):.3f} res {res:.3f}", flush=True)
    return pd.DataFrame(rows)


def convention_check():
    u = ('https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?product=predictions&station=9444090'
         '&begin_date={}&end_date={}&datum=MSL&interval=h&units=metric&time_zone=gmt&format=json')
    vals = []
    for m in range(1, 13):
        b = f'2023{m:02d}01'; e = (pd.Timestamp(b) + pd.offsets.MonthEnd(0)).strftime('%Y%m%d')
        vals += get(u.format(b, e))['predictions']
    df = pd.DataFrame(vals)
    t = pd.to_datetime(df.t).values.astype('datetime64[ns]'); v = df.v.astype(float).values
    hc, res = harmonic_fit(t, v, FIT)
    pub = get('https://api.tidesandcurrents.noaa.gov/mdapi/prod/webapi/stations/9444090/harcon.json?units=metric')
    pub = {c['name'].lower(): c for c in pub['HarmonicConstituents']}
    rows = []
    for c in CONSTS:
        if c in pub:
            dph = ((-np.degrees(np.angle(hc[c])) - pub[c]['phase_GMT'] + 180) % 360) - 180
            rows.append(dict(const=c, amp_fit=abs(hc[c]), amp_noaa=pub[c]['amplitude'], dphase_deg=dph))
    out = pd.DataFrame(rows); out.to_csv(OUT / 'gauge_convention_check.csv', index=False)
    print('convention check (Port Angeles 2023 NOAA predictions re-analysed):'); print(out.round(3))


if __name__ == '__main__':
    convention_check()
    g = pd.concat([noaa(), dfo()], ignore_index=True)
    g.to_csv(OUT / 'gauge_constants.csv', index=False)
    print(g[g.const == 'm2'].groupby('source').size())
