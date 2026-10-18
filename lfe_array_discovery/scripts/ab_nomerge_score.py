"""Score no-merge templates with ab_famscore, two site groupings: NW-original (B013/W030/B943 = one site)
and all stations < 3 km apart (ab_mute.sites, as used for SW/NVI).  usage: python ab_nomerge_score.py stacks.npz out.csv"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_common as A, ab_mute as M, ab_famscore as F
npz, out = sys.argv[1], sys.argv[2]
sub = sys.argv[3] if len(sys.argv) > 3 else "NW"
d0 = F.score(npz).rename(columns={"n_imp_sites": "n_sites_orig"})
st = A.POOL[sub]; sid = M.sites(st)
F.SITE = {s: f"G{k}" for s, k in zip(st, sid)}
d3 = F.score(npz)[["grp", "n_imp_sites"]].rename(columns={"n_imp_sites": "n_sites_3km"})
d = d0.merge(d3, on="grp"); d.to_csv(out, index=False)
for c in ["n_sites_orig", "n_sites_3km"]:
    print(c, ">=3:", int((d[c] >= 3).sum()), " ==2:", int((d[c] == 2).sum()), " of", len(d))
