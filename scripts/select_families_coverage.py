"""Coverage-balanced per-station family selection for the densify budget (hard cap TARGET).
mode 'label' (BOREHOLE, in-domain picker): LFE-labeled only, rank by P(LFE).
mode 'snr'   (FLEET broadband): rank by family SNR (validated predictive of causality survival, AUC~0.68 at
             CLRS). The fleet picker's family-stack P(LFE) is ANTI-predictive (AUC 0.40) -> NOT used; causality
             is the sole verdict downstream. No picker CSV needed.
Dedup near-duplicate templates within 0.05deg bin (cc>=0.9, keep higher-rank). Writes <disc>_sel300.summary.csv.
Usage: python select_families_coverage.py <disc_prefix> <picker_csv_or_-> <target> [label|snr]"""
import sys, numpy as np, pandas as pd

disc, pk_csv, TARGET = sys.argv[1], sys.argv[2], int(sys.argv[3])
mode = sys.argv[4] if len(sys.argv) > 4 else "label"
s = pd.read_csv(disc + ".summary.csv")
if mode == "label":
    pk = pd.read_csv(pk_csv)[["fam", "pred", "p_lfe"]]
    s = s.merge(pk, left_on="family_id", right_on="fam")
    s = s[s.pred == "LFE"].copy(); RANK = "p_lfe"
else:                                    # 'snr' (fleet): rank by SNR from the cluster summary
    s = s.copy(); RANK = "snr"
d = np.load(disc + ".npz", allow_pickle=True)
s["bl"] = (s.lat / 0.05).round().astype(int); s["bo"] = (s.lon / 0.05).round().astype(int)


def maxcc(a, b, mx=10):
    a = a - a.mean(); b = b - b.mean()
    a /= (np.linalg.norm(a) + 1e-12); b /= (np.linalg.norm(b) + 1e-12)
    n = min(len(a), len(b)); best = 0.0
    for sh in range(-mx, mx + 1):
        v = np.dot(a[sh:n], b[:n - sh]) if sh >= 0 else np.dot(a[:n + sh], b[-sh:n])
        best = max(best, abs(v))
    return best


# dedup within bin: cc>=0.9 -> drop the lower-rank duplicate
drop = set()
for _, g in s.sort_values(RANK, ascending=False).groupby(["bl", "bo"]):
    fams = g.family_id.tolist()
    tmpl = {f: np.asarray(d[f], float) for f in fams if f in d.files}
    for i in range(len(fams)):
        if fams[i] in drop or fams[i] not in tmpl:
            continue
        for j in range(i + 1, len(fams)):
            if fams[j] in drop or fams[j] not in tmpl:
                continue
            if maxcc(tmpl[fams[i]], tmpl[fams[j]]) >= 0.9:
                drop.add(fams[j])
s = s[~s.family_id.isin(drop)].copy(); n_dedup = len(drop)

# coverage-balanced: top-rank family per occupied bin, then fill by rank, HARD CAP at TARGET
s = s.sort_values(RANK, ascending=False)
nbins = s.groupby(["bl", "bo"]).ngroups
sel = s.groupby(["bl", "bo"]).head(1)                     # coverage: best-rank per bin
if len(sel) > TARGET:
    sel = sel.sort_values(RANK, ascending=False).head(TARGET)
elif len(sel) < TARGET:
    fill = s[~s.family_id.isin(sel.family_id)].head(TARGET - len(sel))
    sel = pd.concat([sel, fill]).drop_duplicates("family_id").head(TARGET)
out = disc + "_sel300.summary.csv"
drop_cols = [c for c in ["fam", "pred", "p_lfe", "bl", "bo"] if c in sel.columns]
sel.drop(columns=drop_cols).to_csv(out, index=False)
print(f"{disc.split('/')[-1]} [{mode}]: {len(s)+n_dedup} fam -> dedup -{n_dedup} -> {nbins} bins -> "
      f"SELECTED {len(sel)} by {RANK} -> {out}")
