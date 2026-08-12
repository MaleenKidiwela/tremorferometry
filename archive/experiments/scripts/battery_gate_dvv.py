#!/usr/bin/env python
"""Battery gate AFTER causality + dv/v (the user-mandated second gate, 2026-07-11).

Rationale: the Tier-1 trust battery false-negatived when run on the FULL cap-off family set (that set is
~2/3 matched-filter ringing, so real ~ fake). Run instead on the CAUSALITY-CERTIFIED subset, where the
input is no longer noise-dominated -> the battery is a properly-powered, independent second gate. Keep
families that are BOTH causality-reliable AND not battery-FAIL; recompute the Z dv/v over that doubly-
gated set. NEVER destroys the causality dv/v -- writes *_GATED outputs + a battery_gate summary so the
gate's impact is auditable before it is trusted fleet-wide.

Usage: battery_gate_dvv.py <STA> <TAG>
Needs: data/<sta>_causality_cert.csv, data/family_trust_tier1_<STA>.csv, data/daily_dvv_<TAG>_Z_2to4.csv
"""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

STA, TAG = sys.argv[1], sys.argv[2]; sta = STA.lower()

# --- causality-certified families ---
C = pd.read_csv(f"data/{sta}_causality_cert.csv")
cert = set(C[C.reliable].fam.astype(str))

# --- battery verdicts (real families only) ---
B = pd.read_csv(f"data/family_trust_tier1_{STA}.csv")
B = B[B.kind == "F"].copy(); B["fam"] = B.fam.astype(str)
verd = dict(zip(B.fam, B.verdict))

cert_scored = [f for f in cert if f in verd]
n_trusted = sum(verd[f] == "TRUSTED" for f in cert_scored)
n_undet   = sum(verd[f] == "UNDETERMINED" for f in cert_scored)
n_fail    = sum(verd[f] == "FAIL" for f in cert_scored)
n_noscore = len(cert) - len(cert_scored)          # certified but battery couldn't score (too few dets)

# doubly-gated: certified AND not battery-FAIL (a not-scored certified fam is KEPT, conservative)
gated = set(f for f in cert if verd.get(f, "UNDETERMINED") != "FAIL")
trusted_only = set(f for f in cert if verd.get(f) == "TRUSTED")

def zseries(fams):
    df = pd.read_csv(f"data/daily_dvv_{TAG}_Z_2to4.csv")
    df["patch"] = df.patch.astype(str)
    df = df[df.patch.isin(fams) & (df.cc_max >= 0.6)]; df["date"] = pd.to_datetime(df.date)
    med = df.groupby("date").dvv.median(); n = df.groupby("date").patch.nunique()
    return med[n >= 3].rolling(15, center=True, min_periods=5).median()

z_caus, z_gate = zseries(cert), zseries(gated)
z_gate.mul(100).rename("dvv_pct").to_csv(f"data/daily_dvv_{TAG}_Z_2to4_GATED.csv")

fig, ax = plt.subplots(figsize=(13, 4.3))
ax.plot(z_caus.index, z_caus.values*100, lw=1.0, color="0.6",  label=f"causality only ({len(cert)} fam)")
ax.plot(z_gate.index, z_gate.values*100, lw=1.6, color="#b21f2d", label=f"causality ∩ battery ({len(gated)} fam)")
ax.axhline(0, color="0.6", lw=0.6); ax.set_ylim(-1, 1); ax.grid(alpha=0.2)
ax.set_ylabel("dv/v (%)"); ax.legend(loc="upper right", frameon=False)
ax.set_title(f"{STA} certified dv/v — causality vs battery-gated (Z, 2-4 s, origin-anchored)")
fig.tight_layout(); fig.savefig(f"figures/{sta}_battery_gated_dvv.png", dpi=130)

summ = dict(station=STA, tag=TAG, causality_cert=len(cert),
            battery_scored=len(cert_scored), battery_trusted=n_trusted,
            battery_undet=n_undet, battery_fail=n_fail, cert_not_scored=n_noscore,
            gated_n=len(gated), gated_std_pct=round(float(z_gate.std()*100), 3),
            causality_std_pct=round(float(z_caus.std()*100), 3),
            trusted_only_n=len(trusted_only))
json.dump(summ, open(f"data/{sta}_battery_gate_summary.json", "w"), indent=2)
print(f"[{STA}] BATTERY GATE: cert={len(cert)} | within cert -> TRUSTED {n_trusted} / UNDET {n_undet} "
      f"/ FAIL {n_fail} / not-scored {n_noscore}")
print(f"[{STA}] doubly-gated (cert & !FAIL) n={len(gated)} std={summ['gated_std_pct']}% "
      f"vs causality-only std={summ['causality_std_pct']}%")
print(f"[{STA}] BATTERY_GATE_DONE {json.dumps(summ)}")
