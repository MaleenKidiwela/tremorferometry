"""Catalogue v1.3 = v1.2 + NW event-sharing flag columns (2026-10-07).

Reads data/catalogue/families_v1_2.csv (unchanged) and lfe_extension_wa/results/nw_sharing/ outputs; writes
data/catalogue/families_v1_3.csv. No row removed, no v1.2 column changed (asserted at the end).

New columns (Royer_NW_WA rows in the sharing study; empty for every other row):
  nw_sharing_tested       True for the 132 NW families in the sharing study; False for the 9 NW rows left out of it
                          (cross-region dup non-representatives); empty for other regions
  nw_source_group         pre-registered E3 source group (source_groups_E3.csv, G###; singletons have their own id)
  nw_array_group          array-confirmed group (source_groups_E3_array_confirmed.csv, C###)
  nw_copy_of              family_uid of the array group's primary, for non-primary members of multi-family array groups
  nw_link_type            same_source_array | inseparable_on_permanent | none (precedence in that order)
  nw_sharing_max          max E3 sharing fraction s (variant D shared repeats / smaller family's repeats) over all partners
  nw_sharing_max_partner  family_uid of that partner
  nw_sharing_max_n_shared shared repeats behind that s
Primary rule (nw_copy_of): best location_quality (A > B > C), then most array-era discovery detections
(n_detections), then smallest erh, then family_uid.
Usage: python lfe_extension_wa/scripts/catalogue_v1_3.py [--force]
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
SD = f"{ROOT}/lfe_extension_wa/results/nw_sharing"
V12 = f"{ROOT}/data/catalogue/families_v1_2.csv"
V13 = f"{ROOT}/data/catalogue/families_v1_3.csv"
NEW = ["nw_sharing_tested", "nw_source_group", "nw_array_group", "nw_copy_of", "nw_link_type", "nw_sharing_max",
       "nw_sharing_max_partner", "nw_sharing_max_n_shared"]

if os.path.exists(V13) and "--force" not in sys.argv:
    sys.exit(f"{V13} exists; refusing to overwrite (pass --force)")

cat = pd.read_csv(V12, low_memory=False)
meta = pd.read_csv(f"{SD}/nw_families_meta.csv")
uid = meta.set_index("fid").family_uid
G = pd.read_csv(f"{SD}/source_groups_E3.csv")
C = pd.read_csv(f"{SD}/source_groups_E3_array_confirmed.csv")
L = pd.read_csv(f"{SD}/links_E3_classified.csv")
P = pd.read_csv(f"{SD}/pair_table.csv", usecols=["era", "fid_i", "fid_j", "n_D_small", "s"], low_memory=False)
assert meta.family_uid.isin(cat.family_uid).all() and len(meta) == 132
assert len(L) == 58 and int(L.same_source.sum()) == 38

F = meta[["fid", "family_uid"]].copy()
F["nw_sharing_tested"] = True
F["nw_source_group"] = F.fid.map(G.set_index("fid").group)
F["nw_array_group"] = F.fid.map(C.set_index("fid").group)

# link type per family
ss = set(L[L.same_source].fid_i) | set(L[L.same_source].fid_j)
anyl = set(L.fid_i) | set(L.fid_j)
F["nw_link_type"] = np.where(F.fid.isin(ss), "same_source_array",
                             np.where(F.fid.isin(anyl), "inseparable_on_permanent", "none"))

# primary of each multi-family array group
q = cat.set_index("family_uid")
F["_qual"] = F.family_uid.map(q.location_quality).map({"A": 0, "B": 1, "C": 2}).fillna(3)
F["_ndet"] = F.family_uid.map(q.n_detections).fillna(0)
F["_erh"] = F.family_uid.map(q.erh).fillna(np.inf)
csize = F.groupby("nw_array_group").fid.transform("size")
prim = (F[csize > 1].sort_values(["nw_array_group", "_qual", "_ndet", "_erh", "family_uid"],
                                 ascending=[True, True, False, True, True])
        .groupby("nw_array_group").family_uid.first())
F["nw_copy_of"] = np.where(csize > 1, F.nw_array_group.map(prim), None)
F.loc[F.nw_copy_of == F.family_uid, "nw_copy_of"] = None

# max sharing fraction (E3, all partners)
P = P[P.era == "E3"]
B = pd.concat([P.rename(columns={"fid_i": "fid", "fid_j": "partner"}),
               P.rename(columns={"fid_j": "fid", "fid_i": "partner"})]).dropna(subset=["s"])
B = B.sort_values(["fid", "s", "n_D_small"], ascending=[True, False, False]).drop_duplicates("fid").set_index("fid")
F["nw_sharing_max"] = F.fid.map(B.s)
F["nw_sharing_max_partner"] = F.fid.map(B.partner).map(uid)
F["nw_sharing_max_n_shared"] = F.fid.map(B.n_D_small)

out = cat.merge(F[["family_uid"] + NEW], on="family_uid", how="left", validate="one_to_one")
nw = out.region == "Royer_NW_WA"
out.loc[nw & out.nw_sharing_tested.isna(), "nw_sharing_tested"] = False
out.to_csv(V13, index=False)

# ---- checks
a = pd.read_csv(V12, low_memory=False)
b = pd.read_csv(V13, low_memory=False)
assert len(a) == len(b) == 610
assert list(b.columns) == list(a.columns) + NEW
pd.testing.assert_frame_equal(a, b[a.columns])
assert b.loc[b.region != "Royer_NW_WA", NEW].isna().all().all()
t = b[b.nw_sharing_tested == True]
assert len(t) == 132 and (b.region == "Royer_NW_WA").sum() == 141
assert t.nw_source_group.nunique() == 93 and t.nw_array_group.nunique() == 100
assert t.nw_copy_of.notna().sum() == 32 and t.nw_copy_of.dropna().isin(t.family_uid).all()
assert set(t.nw_link_type) <= {"same_source_array", "inseparable_on_permanent", "none"}

# ---- counts for README
act = set(pd.read_csv(f"{SD}/headline_N_by_era.csv").query("era == 'E3' and N > 0").fid)
ta = t[t.family_uid.isin(uid[uid.index.isin(act)])]
rm = meta.set_index("family_uid").royer_match.astype(bool)
t = t.assign(royer=t.family_uid.map(rm))
Ln = L[~L.same_source]
both = ss & (set(Ln.fid_i) | set(Ln.fid_j))
print("link type:", t.nw_link_type.value_counts().to_dict(), "| same_source families that also have a permanent-only link:", len(both))
print("multi array groups:", int((t.groupby("nw_array_group").size() > 1).sum()),
      "families in them:", int(t.nw_array_group.map(t.nw_array_group.value_counts()).gt(1).sum()), "copies:", int(t.nw_copy_of.notna().sum()))
print("distinct: permanent", t.nw_source_group.nunique(), "(E3-active", ta.nw_source_group.nunique(), "of", len(ta), ")",
      "| array", t.nw_array_group.nunique(), "(E3-active", ta.nw_array_group.nunique(), ")")
print("royer match families", int(t.royer.sum()), "distinct array groups", t[t.royer].nw_array_group.nunique(),
      "distinct permanent groups", t[t.royer].nw_source_group.nunique())
print("sharing_max >= 0.3:", int((t.nw_sharing_max >= 0.3).sum()), "NaN:", int(t.nw_sharing_max.isna().sum()))
print("primaries:", t[t.nw_copy_of.isna() & t.nw_array_group.map(t.nw_array_group.value_counts()).gt(1)].location_quality.value_counts().to_dict())
print("OK", V13)
