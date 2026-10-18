"""S3 post-hoc (NOT pre-registered; labelled as such in RESULTS_S3.md): is the episode-to-episode patch pattern explained by
the episode's direction (start end)? Mean leave-pair-out residual correlation for same-direction vs opposite-direction
pairs, with a direction-label shuffle control. Also NW-ward vs SE-ward share of the NW half of each region."""
import json
import numpy as np
import pandas as pd
import importlib.util

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s3"
spec = importlib.util.spec_from_file_location("p", f"{ROOT}/svi_science/scripts/s3_patches.py"); P = importlib.util.module_from_spec(spec); spec.loader.exec_module(P)
rng = np.random.default_rng(3)


def main():
    K = pd.read_csv(f"{OUT}/episodes_kinematics.csv")
    res = {}
    for reg in ["SVI", "NW"]:
        Cdf, mid = P.counts(reg)
        C = Cdf.to_numpy().astype(float)
        k = K[K.region == reg].set_index("episode").loc[Cdf.columns]
        d = np.where(k.migrates, np.sign(k.lfe_speed), 0).astype(int)
        n = C.shape[1]; cor = {}
        for i in range(n):
            for j in range(i + 1, n):
                if abs((mid[j] - mid[i]) / np.timedelta64(1, "D")) < 90:
                    continue
                ri, rj = P.resid(C, i, j); cor[(i, j)] = np.corrcoef(ri, rj)[0, 1]

        def stat(dd):
            same = [c for (i, j), c in cor.items() if dd[i] != 0 and dd[i] == dd[j]]
            opp = [c for (i, j), c in cor.items() if dd[i] != 0 and dd[j] != 0 and dd[i] != dd[j]]
            return np.mean(same) - np.mean(opp), np.mean(same), np.mean(opp), len(same), len(opp)
        obs = stat(d)
        perm = np.array([stat(rng.permutation(d))[0] for _ in range(5000)])
        # consecutive-pair direction alternation
        dm = d[d != 0]
        alt = int(np.sum(dm[1:] != dm[:-1]))
        palt = np.array([np.sum((x := rng.permutation(dm))[1:] != x[:-1]) for _ in range(5000)])
        res[reg] = dict(same_minus_opposite=float(obs[0]), mean_corr_same_direction=float(obs[1]), mean_corr_opposite_direction=float(obs[2]),
                        n_same=obs[3], n_opp=obs[4], p_direction_shuffle=float((np.sum(perm >= obs[0]) + 1) / 5001),
                        direction_sequence="".join({1: "N", -1: "S", 0: "."}[x] for x in d),
                        n_direction_changes_consecutive=alt, n_consecutive=len(dm) - 1,
                        p_more_alternation_than_shuffle=float((np.sum(palt >= alt) + 1) / 5001))
    json.dump(res, open(f"{OUT}/post_direction_modes.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
