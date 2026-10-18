#!/usr/bin/env python
"""S12c step 1: predicted zero-offset reflection lags of the LVL top/bottom under PGC, VGZ, SNB.

Method. The RF study picked LVL top / bottom DEPTHS (k=0 harmonic) with its own 1-D model (rf04: LAY_VP, Vp/Vs 1.75)
and teleseismic slowness. Undo that conversion to recover the observed Ps delay times t_Ps (p = 0.06 s/km), then
predict zero-offset two-way times in several velocity models that all honour the same t_Ps:
  PP = 2 int dz/Vp ; SS = 2 int dz/Vs ; PS = int dz/Vp + int dz/Vs (P down, S up, = RF PpPs) .
Models: (a) RF model (Vp/Vs 1.75 everywhere); (b) USGS CVM v1.7 south_VI profile, depth rescaled so that its Ps time to the
pick equals the observed t_Ps; (c) as (a) but the LVL (between top and bottom) has Vp/Vs 2.0 / 2.3 (high-pore-pressure
LVL), depth of bottom re-solved to keep the observed bottom t_Ps. Range over models = prediction band.
Also: delay change of everything below the LVL top for +2.5 % Vs (and Vp change 0 or half) in a 3.5 km layer at the LVL top.
"""
import json, numpy as np, pandas as pd
ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s12_lvl/acf_pgc/tables"
DZ = 0.01
LAY_Z = np.array([0, 5, 15, 25, 35, 45]); LAY_VP = np.array([5.8, 6.3, 6.7, 6.9, 7.2, 7.8]); K0 = 1.75
cvm = pd.read_csv(f"{ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv")
P = 0.06

def rf_model(z):
    i = np.searchsorted(LAY_Z, z, side="right") - 1
    a = LAY_VP[i]; return a, a / K0

def cvm_model(z):
    i = np.searchsorted(cvm.ztop_km.values, z, side="right") - 1
    return cvm.vp.values[i], cvm.vs.values[i]

def times(model, zmax, p=P, lvl=None):
    z = np.arange(0, zmax, DZ) + DZ / 2
    a, b = model(z); a = a.astype(float).copy(); b = b.astype(float).copy()
    if lvl is not None:  # (ztop, kappa) : inside [ztop, zmax] set Vs = Vp/kappa
        m = z >= lvl[0]; b[m] = a[m] / lvl[1]
    tps = np.sum((np.sqrt(1 / b**2 - p**2) - np.sqrt(1 / a**2 - p**2)) * DZ)
    tpp = 2 * np.sum(DZ / a); tss = 2 * np.sum(DZ / b); return tps, tpp, tss, (tpp + tss) / 2

def solve_depth(model, tps_target, lvl=None, zlo=5, zhi=90):
    for _ in range(60):
        zm = (zlo + zhi) / 2
        if times(model, zm, lvl=lvl)[0] < tps_target: zlo = zm
        else: zhi = zm
    return zm

rows = []; summ = {}
for s in ["CN.PGC", "CN.VGZ", "CN.SNB"]:
    r = json.load(open(f"{ROOT}/svi_science/results/s12_lvl/rf_gosselin/tables/full/{s}_result.json"))
    zt, zb = r["z_top"], r["z_bot"]
    tps_t = times(rf_model, zt)[0]; tps_b = times(rf_model, zb)[0]
    for nm, model, lv in [("rf_k1.75", rf_model, None), ("cvm_v1.7", cvm_model, None),
                          ("rf_lvl_k2.0", rf_model, 2.0), ("rf_lvl_k2.3", rf_model, 2.3)]:
        z1 = solve_depth(model, tps_t)
        z2 = solve_depth(model, tps_b, lvl=None if lv is None else (z1, lv))
        for lab, zz, l in [("top", z1, None), ("bottom", z2, None if lv is None else (z1, lv))]:
            tps, tpp, tss, tps0 = times(model, zz, lvl=l)
            rows.append(dict(station=s, model=nm, interface=lab, depth_km=round(zz, 2), t_Ps=round(tps, 2),
                             PP=round(tpp, 2), PS=round(tps0, 2), SS=round(tss, 2)))
    # forward target: +2.5 % Vs in a 3.5 km layer at the LVL top (Vs there ~ Vp/kappa)
    vp_top = rf_model(zt + 1)[0]
    tgt = {}
    for k in [1.75, 2.0, 2.3]:
        vs = vp_top / k; dt_ss = -2 * 3.5 / vs * 0.025; dt_pp_half = -2 * 3.5 / vp_top * 0.0125
        tgt[f"kappa{k}"] = dict(vs_lvl=round(vs, 2), dSS_ms=round(1e3 * dt_ss, 1), dPP_ms_if_dVp_half=round(1e3 * dt_pp_half, 1),
                                dPS_ms=round(1e3 * (dt_ss + dt_pp_half) / 2, 1))
    summ[s] = dict(z_top_rf=zt, z_bot_rf=zb, tPs_top=round(tps_t, 2), tPs_bot=round(tps_b, 2), target_delays=tgt)
T = pd.DataFrame(rows); T.to_csv(f"{OUT}/predicted_lags.csv", index=False)
band = T.groupby(["station", "interface"]).agg(PP_lo=("PP", "min"), PP_hi=("PP", "max"), PS_lo=("PS", "min"),
                                                PS_hi=("PS", "max"), SS_lo=("SS", "min"), SS_hi=("SS", "max")).reset_index()
band.to_csv(f"{OUT}/predicted_lag_bands.csv", index=False)
json.dump(summ, open(f"{OUT}/predicted_targets.json", "w"), indent=1)
pd.set_option("display.width", 200); print(T.to_string(index=False)); print(band.to_string(index=False)); print(json.dumps(summ, indent=1))
