#!/usr/bin/env python
"""rf03_rf: P-SV-SH receiver functions per event, rfpy (Audet) conventions, Gosselin et al. 2020 settings.

Per event: merge; remove sensitivity; rotate ->ZNE with the StationXML azimuths; detrend/taper; low-pass 1 Hz
(2-pole zero-phase); resample to 5 Hz on an exact grid tied to the iasp91 P time; SNR on Z (rfpy calc_snr:
0.05-1 Hz, RMS 30 s after vs 30 s before P, dB); rotate Z-R-T -> P-SV-SH (rfpy 'PVH', surface vp 6.0, vs 3.5);
Wiener deconvolution of P from SV and SH with the pre-event P-mode noise in the denominator (rfpy 'wiener',
dts = 120 s: signal P-5..P+110 s, noise P-120..P-5 s); RF band-passed 0.05-0.5 Hz (main) and 0.05-0.7 Hz (variant).
Instrument response is not deconvolved: the RF is a spectral ratio of components of the same sensor, so it cancels
(sensitivity + orientation are applied).  All events are kept with their SNR (rank, do not discard); the >5 dB rule
is applied downstream.
Output: data/rf_teleseis/rf/<NET>.<STA>_rf.npz  (t, V05, H05, V07, H07, P05 + meta columns)
Usage: python -P rf03_rf.py <NET.STA> [...]  (or no args = all stations with a DONE download marker)
"""
import os, sys, glob, warnings
import numpy as np, pandas as pd
from multiprocessing import Pool
from obspy import read, read_inventory, UTCDateTime

warnings.filterwarnings("ignore")
ROOT = "/home/jovyan/tremorferometry"
D = f"{ROOT}/data/rf_teleseis"
OUTD = f"{D}/rf"
SR, DTS = 5.0, 120.0
VP, VS = 6.0, 3.5
TKEEP = (-10.0, 40.0)   # RF lag window stored (s)
INV = {}


def bandpass(x, f1, f2, sr=SR):
    from obspy.signal.filter import bandpass as bp
    return bp(x, f1, f2, sr, corners=2, zerophase=True)


def one(args):
    key, r = args
    fn = f"{D}/raw/{key}/{r['evid']}.mseed"
    try:
        st = read(fn)
    except Exception:
        return None
    try:
        st.merge(method=1, fill_value=0)
        tP = UTCDateTime(str(r["tP"]))
        inv = INV[key]
        st.remove_sensitivity(inv)
        st.rotate("->ZNE", inventory=inv)
        st = st.select(component="[ZNE]")
        if len(st) != 3:
            return None
        g0, g1 = tP - DTS - 10, tP + DTS + 10
        for tr in st:
            if tr.stats.starttime > g0 or tr.stats.endtime < g1:
                return None
        st.trim(g0 - 5, g1 + 5)
        st.detrend("demean"); st.detrend("linear"); st.taper(0.05, max_length=5)
        st.filter("lowpass", freq=1.0, corners=2, zerophase=True)
        tg = np.arange(-DTS, DTS, 1 / SR)  # exact grid around P
        Z = N = E = None
        for tr in st:
            x = tr.times() + (tr.stats.starttime - tP)
            y = np.interp(tg, x, tr.data.astype(float))
            c = tr.stats.channel[-1]
            if c == "Z": Z = y
            elif c == "N": N = y
            elif c == "E": E = y
        # zero-filled gaps make flat segments; reject if > 2 % of samples are exactly zero on any component
        if any(np.mean(c == 0) > 0.02 for c in (Z, N, E)):
            return None
        # SNR on Z (rfpy calc_snr)
        def snr(y):
            yy = bandpass(y - y.mean(), 0.05, 1.0)
            s = yy[(tg >= 0) & (tg <= 30)]; n = yy[(tg >= -30) & (tg <= 0)]
            return 10 * np.log10(np.mean(s ** 2) / np.mean(n ** 2))
        baz = np.deg2rad(r["baz_true"])
        R = -N * np.cos(baz) - E * np.sin(baz)
        T = E * np.cos(baz) - N * np.sin(baz)  # obspy NE->RT convention
        snrz, snrr = snr(Z), snr(R)
        p = r["slow"]
        qp = np.sqrt(1 / VP ** 2 - p ** 2); qs = np.sqrt(1 / VS ** 2 - p ** 2)
        m11 = p * VS * VS / VP; m12 = -(1 - 2 * VS * VS * p * p) / (2 * VP * qp)
        m21 = (1 - 2 * VS * VS * p * p) / (2 * VS * qs); m22 = p * VS
        Pm = -m11 * R + m12 * Z
        Vm = -m21 * R + m22 * Z
        Hm = -T / 2.0
        # rfpy Wiener deconvolution
        nn = int(round((DTS - 5.0) * SR)) + 1
        sig = (tg >= -5.0) & (tg < DTS - 10.0)
        noi = (tg >= -DTS) & (tg < -5.0)

        def win(y, m):
            w = y[m].copy()
            w = w - w.mean()
            nt = len(w); k = min(int(2.0 * SR), int(0.05 * nt))  # taper 5 % max 2 s
            if k > 1:
                h = np.hanning(2 * k)
                w[:k] *= h[:k]; w[-k:] *= h[k:]
            out = np.zeros(nn); out[:min(nn, nt)] = w[:nn]
            return out
        Ps, Vs_, Hs, Pn = win(Pm, sig), win(Vm, sig), win(Hm, sig), win(Pm, noi)
        Fp, Fv, Fh, Fn = (np.fft.fft(x, n=nn) for x in (Ps, Vs_, Hs, Pn))
        Spp = np.real(Fp * np.conj(Fp)); Snn = np.real(Fn * np.conj(Fn))
        den = Spp + Snn
        rfP = np.fft.ifftshift(np.real(np.fft.ifft(Spp / den)))
        rfV = np.fft.ifftshift(np.real(np.fft.ifft(Fv * np.conj(Fp) / den)))
        rfH = np.fft.ifftshift(np.real(np.fft.ifft(Fh * np.conj(Fp) / den)))
        lag = (np.arange(nn) - nn // 2) / SR
        keep = (lag >= TKEEP[0]) & (lag <= TKEEP[1])
        out = dict(evid=r["evid"], snr=snrz, snrr=snrr)
        for tag, f2 in (("05", 0.5), ("07", 0.7)):
            out["V" + tag] = bandpass(rfV, 0.05, f2)[keep]
            out["H" + tag] = bandpass(rfH, 0.05, f2)[keep]
            out["P" + tag] = bandpass(rfP, 0.05, f2)[keep]
        out["lag"] = lag[keep]
        return out
    except Exception as e:  # noqa
        return None


def station(key):
    os.makedirs(OUTD, exist_ok=True)
    net, sta = key.split(".")
    A = pd.read_csv(f"{D}/arrivals_{key}.csv", keep_default_na=False)
    have = {os.path.basename(f)[:-6] for f in glob.glob(f"{D}/raw/{key}/*.mseed")}
    A = A[A.evid.isin(have)].reset_index(drop=True)
    INV[key] = read_inventory(f"{D}/raw/{key}.xml")
    # back-azimuth recomputed here: rf00 (first run) stored the event->station azimuth in 'baz' by mistake
    from obspy.geodetics import gps2dist_azimuth
    S = pd.read_csv(f"{D}/stations.csv").set_index(["net", "sta"]).loc[(net, sta)]
    A["baz_true"] = [gps2dist_azimuth(S.lat, S.lon, a, b)[1] for a, b in zip(A.latitude, A.longitude)]
    with Pool(8) as P:
        res = P.map(one, [(key, r) for _, r in A.iterrows()], chunksize=20)
    ok = [x for x in res if x is not None]
    M = A.set_index("evid").loc[[x["evid"] for x in ok]].reset_index()
    M["baz_rf00"] = M["baz"]; M["baz"] = M["baz_true"]
    M["snr"] = [x["snr"] for x in ok]; M["snrr"] = [x["snrr"] for x in ok]
    arr = {k: np.array([x[k] for x in ok]) for k in ("V05", "H05", "P05", "V07", "H07", "P07")}
    np.savez_compressed(f"{OUTD}/{key}_rf.npz", lag=ok[0]["lag"], **arr)
    M.to_csv(f"{OUTD}/{key}_meta.csv", index=False)
    print(key, "events with data", len(A), "RFs", len(ok), "SNR>5dB", int((M.snr > 5).sum()), flush=True)


if __name__ == "__main__":
    keys = sys.argv[1:]
    if not keys:
        keys = [os.path.basename(f)[:-5] for f in sorted(glob.glob(f"{D}/raw/*.DONE"))]
    for k in keys:
        station(k)
