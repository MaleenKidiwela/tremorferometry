"""S6 Deviation 1: matched-filter amplitude ratio per detection x station against the family x station reference
stack (3-C, 2-8 Hz velocity, 20 Hz, S window -0.5..+3.5 s, best lag within +-0.25 s). Noise = the same projection on
two pre-arrival windows. Output results/s6/mf_station.parquet. Rules: results/s6/PLAN_S6.md (Deviation 1)."""
import glob
import numpy as np, pandas as pd

ROOT = "/home/jovyan/tremorferometry"; OUT = f"{ROOT}/svi_science/results/s6"
STAS = ["CN.NLLB", "CN.PFB", "CN.PGC", "CN.SNB", "CN.VGZ", "PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007",
        "PB.B009", "PB.B010", "PB.B011", "PB.B926", "UW.LRIV"]
S0, L, LAG = 280, 80, 5          # S window starts at sample 280 (-0.5 s), 80 samples (4 s), lag +-5 samples (0.25 s)
NWIN = [5, 115]                  # noise windows start (-14.25 s, -8.75 s), same length and lag search


def main():
    T = pd.read_parquet(f"{OUT}/targets.parquet")
    T["head"] = (T.P_same >= 0.75) & T.keep_xfam8 & (T.kind == "event")
    fam = dict(zip(T.det_id, T.fam)); head = set(T.det_id[T["head"]]); ev05 = set(T.det_id[(T.kind == "event")])
    files = sorted(glob.glob(f"{OUT}/env/wave_*.npz"))
    # pass 1: stacks (headline; and P>=0.5 events as fallback)
    acc = {}; accb = {}
    for f in files:
        z = np.load(f); W = z["W"][:, :, S0 - LAG:S0 + L + LAG].astype(np.float64)
        for i, (d, s) in enumerate(zip(z["det_id"], z["sta"])):
            if d not in ev05: continue
            k = (fam[d], int(s)); x = W[i]
            ok = np.isfinite(x).all(axis=1)
            if not ok.any(): continue
            for D_, use in ((accb, True), (acc, d in head)):
                if not use: continue
                a = D_.setdefault(k, [np.zeros_like(x), np.zeros(3)])
                a[0][ok] += x[ok]; a[1][ok] += 1
    stk = {}
    for k, (sm, n) in acc.items():
        if n.max() >= 20: stk[k] = (sm[:, LAG:LAG + L] / np.maximum(n, 1)[:, None], n.max(), "headline")
    for k, (sm, n) in accb.items():
        if k not in stk and n.max() >= 5: stk[k] = (sm[:, LAG:LAG + L] / np.maximum(n, 1)[:, None], n.max(), "P050_fallback")
    pd.DataFrame([dict(fam=k[0], sta=STAS[k[1]], stack_rms=float(np.sqrt(np.mean(v[0] ** 2))), n_stack=int(v[1]), stack_src=v[2])
                  for k, v in stk.items()]).to_csv(f"{OUT}/family_station_stacks.csv", index=False)
    # pass 2: projections
    rows = []
    for f in files:
        z = np.load(f); W = z["W"].astype(np.float64)
        for i, (d, s) in enumerate(zip(z["det_id"], z["sta"])):
            k = (fam[d], int(s))
            if k not in stk: continue
            Tk, nst, src = stk[k]; x = W[i]
            ok = np.isfinite(x).all(axis=1) & (np.abs(Tk).sum(axis=1) > 0)
            if ok.sum() < 2: continue
            Tc = Tk[ok]; xc = x[ok]; tt = (Tc * Tc).sum()

            SW = np.lib.stride_tricks.sliding_window_view(xc, L, axis=1)

            def proj(st):
                seg = SW[:, st - LAG:st + LAG + 1, :]
                num = np.einsum("clk,ck->l", seg, Tc); j = int(np.argmax(num))
                return num[j] / tt, num[j] / np.sqrt((seg[:, j] ** 2).sum() * tt + 1e-60)
            a, cc = proj(S0)
            n1, _ = proj(NWIN[0]); n2, _ = proj(NWIN[1])
            rows.append((d, STAS[int(s)], a, cc, n1, n2, int(ok.sum()), nst, src))
    M = pd.DataFrame(rows, columns=["det_id", "sta", "a_mf", "cc_mf", "an1", "an2", "ncomp_mf", "n_stack", "stack_src"])
    M["sig_noise"] = np.sqrt((M.an1 ** 2 + M.an2 ** 2) / 2)
    M.to_parquet(f"{OUT}/mf_station.parquet")
    print(len(M), M.describe().T)


if __name__ == "__main__":
    main()
