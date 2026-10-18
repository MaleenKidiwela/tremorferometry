#!/usr/bin/env python
"""E4 step 0 verification — which slot of the SeisBench stacks' "Z12" is North?

Test (Merlin's): for family stacks with a clear, rectilinear P arrival, the horizontal particle
motion of P is along the source-receiver great circle, so the azimuth of the horizontal principal
axis must agree with `path_backazimuth_deg` modulo 180 deg.  Compute that azimuth twice -- once
assuming slot1=N/slot2=E, once assuming slot1=E/slot2=N -- and report the median |difference|.
Expectation from the audit: ~9 deg for 1=N, ~40 deg for 1=E.

Writes results/e4_polarization_check.csv
"""
from __future__ import annotations
import numpy as np, pandas as pd, h5py
from scipy.signal import butter, sosfiltfilt

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
FS = 20.0
P0 = 800            # every trace is aligned on its own P at sample 800


def bp(x, fs, f1=1.0, f2=8.0):
    sos = butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def main():
    m = pd.read_csv(f"{SD}/metadata.csv")
    with h5py.File(f"{SD}/waveforms.hdf5", "r") as f:
        print("[hdf5] data_format:", {k: f["data_format"][k][()] for k in f["data_format"]})
        rows = []
        for r in m.itertuples():
            if not np.isfinite(r.trace_P_snr_db) or r.trace_P_snr_db <= 6:
                continue
            b, idx = r.trace_name.split("$")
            p = idx.split(",")
            w = np.asarray(f["data"][b][int(p[0]), :3, :int(p[2].lstrip(":"))], float)
            if w.shape[1] < P0 + 40:
                continue
            w = bp(w, FS)
            seg = w[:, P0 - 5:P0 + 25]           # -0.25 .. +1.25 s around P
            C = seg @ seg.T / seg.shape[1]
            ev, V = np.linalg.eigh(C)
            ev = ev[::-1]; V = V[:, ::-1]
            rect = 1.0 - (ev[1] + ev[2]) / (2 * ev[0]) if ev[0] > 0 else 0.0
            if rect <= 0.7:
                continue
            u = V[:, 0]
            # slots 1,2 of the eigenvector = the two horizontals
            az_1N = np.degrees(np.arctan2(u[2], u[1])) % 180.0   # 1=N,2=E: az=atan2(E,N)
            az_1E = np.degrees(np.arctan2(u[1], u[2])) % 180.0   # 1=E,2=N
            baz = float(r.path_backazimuth_deg) % 180.0
            d1N = min(abs(az_1N - baz), 180 - abs(az_1N - baz))
            d1E = min(abs(az_1E - baz), 180 - abs(az_1E - baz))
            rows.append(dict(sta=r.station_code, fam=r.source_family_id,
                             psnr=r.trace_P_snr_db, rect=rect, baz=baz,
                             az_1N=az_1N, az_1E=az_1E, d_1N=d1N, d_1E=d1E))
    D = pd.DataFrame(rows)
    D.to_csv(f"{RES}/e4_polarization_check.csv", index=False)
    print(f"[pol] {len(D)} stacks with P snr > 6 dB and rectilinearity > 0.7")
    for lab, c in (("slot1 = N (Z,N,E)", "d_1N"), ("slot1 = E (Z,E,N)", "d_1E")):
        v = D[c].to_numpy()
        print(f"  {lab}: median |az - baz| = {np.median(v):6.1f} deg   "
              f"frac<25deg = {np.mean(v < 25):.2f}   frac<45deg = {np.mean(v < 45):.2f}")
    print(f"[write] {RES}/e4_polarization_check.csv")


if __name__ == "__main__":
    main()
