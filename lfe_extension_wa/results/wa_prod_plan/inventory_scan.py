"""Read-only inventory of Washington station-days on disk (raw miniSEED + 20/40 Hz caches)."""
import os, re, sys, collections
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
roster = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/step4b/setup/stations.csv")
trial = set(pd.read_csv(f"{ROOT}/lfe_extension_wa/results/setup/stations.csv").netsta)
R = set(roster.netsta)
# wider WA list: permanent stations in the WA station tables
TEMP = {"XU", "XG", "ZK", "Z3", "7D", "ZU", "YW", "XD", "ZZ", "TA", "XA", "YB", "3F", "X4", "PO", "C8"}
wide = set()
for f in ["NW_WA", "SW_WA", "OLY2010"]:
    s = pd.read_csv(f"{ROOT}/data/extension_stations/{f}_stations.csv")
    wide |= set(s.net + "." + s.sta)

ROOTS = {
    "raw": ["lfe_array_discovery/traces", "lfe_array_discovery/traces_sealed", "lfe_array_discovery/traces_eq",
            "lfe_extension/traces", "lfe_catalog_repair/traces", "lfe_extension_wa/traces", "lfe_picker/traces",
            "lfe_south/traces", "data/waveforms"],
    "c20": ["lfe_extension_wa/cache20hz", "lfe_extension/cache20hz", "lfe_extension/cache20hz_r1",
            "lfe_picker/cache20hz", "lfe_catalog_repair/cache20hz", "lfe_array_discovery/cache20hz",
            "lfe_south/cache20hz"],
    "c40": ["lfe_array_discovery/cache40hz"],
}
pat = re.compile(r"^[A-Z0-9]{1,2}\.[A-Z0-9]{2,5}$")
rows = []
for kind, rl in ROOTS.items():
    for r in rl:
        base = f"{ROOT}/{r}"
        if not os.path.isdir(base):
            continue
        for ns in os.listdir(base):
            if not pat.match(ns) or ns not in (R | wide):
                continue
            p = f"{base}/{ns}"
            if not os.path.isdir(p):
                continue
            for y in os.listdir(p):
                py = f"{p}/{y}"
                if not (y.isdigit() and os.path.isdir(py)):
                    continue
                for e in os.scandir(py):
                    if e.name.startswith("."):
                        continue
                    try:
                        if e.is_symlink():
                            rp = os.path.realpath(e.path)
                            sz = os.path.getsize(rp) if os.path.exists(rp) else -1
                            link = 1
                        else:
                            rp = e.path; sz = e.stat().st_size; link = 0
                    except OSError:
                        continue
                    doy = e.name.split(".")[0]
                    rows.append((kind, r, ns, int(y), doy, sz, link, rp))
df = pd.DataFrame(rows, columns=["kind", "root", "netsta", "year", "doy", "bytes", "link", "real"])
df["in_roster"] = df.netsta.isin(R)
df["in_trial42"] = df.netsta.isin(trial)
df.to_parquet(sys.argv[1])
print(len(df))
