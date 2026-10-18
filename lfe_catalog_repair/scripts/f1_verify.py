#!/usr/bin/env python
"""F1 — consistency verification of the seven f1_*.csv tables + the report tables.

Read-only.  Prints (i) PASS/FAIL for every cross-table invariant and (ii) the
markdown tables that go into results/f1_fleet_inventory.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

R = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
YEARS = list(range(2009, 2027))
REGIONS = ["R1_N_VancIsl", "R2_S_VancIsl_JdF", "R3_Puget_Olympic", "R4_N_Oregon",
           "R5_C_Oregon", "R6_S_Oregon", "R7_NorCal_MTJ"]

sta = pd.read_csv(f"{R}/f1_stations.csv")
sy = pd.read_csv(f"{R}/f1_station_years.csv")
dens = pd.read_csv(f"{R}/f1_cell_density.csv")
regy = pd.read_csv(f"{R}/f1_regions.csv")
rsta = pd.read_csv(f"{R}/f1_region_stations.csv")
cost = pd.read_csv(f"{R}/f1_cost_estimate.csv")
gt = pd.read_csv(f"{R}/f1_groundtruth.csv")

ok = []


def chk(name, cond, extra=""):
    ok.append(bool(cond))
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}{(' — ' + extra) if extra else ''}")


print("=" * 96)
print("CONSISTENCY CHECKS")
print("=" * 96)
for n, d in [("f1_stations", sta), ("f1_station_years", sy), ("f1_cell_density", dens),
             ("f1_regions", regy), ("f1_region_stations", rsta),
             ("f1_cost_estimate", cost), ("f1_groundtruth", gt)]:
    chk(f"{n} non-empty", len(d) > 0, f"{len(d)} rows, {d.shape[1]} cols")

skey = set(zip(sta.net, sta.sta))
chk("station keys unique in f1_stations", len(skey) == len(sta))
chk("every f1_region_stations station is in f1_stations",
    set(zip(rsta.net, rsta.sta)) <= skey,
    f"{len(set(zip(rsta.net, rsta.sta)) - skey)} orphans")
chk("every f1_station_years station is in f1_stations",
    set(zip(sy.net, sy.sta)) <= skey)
chk("station_years = stations x years",
    len(sy) == len(sta) * len(YEARS), f"{len(sy)} vs {len(sta)}x{len(YEARS)}")
chk("years 2009-2026 complete in f1_station_years", sorted(sy.year.unique()) == YEARS)
chk("years 2009-2026 complete in f1_cell_density", sorted(dens.year.unique()) == YEARS)
chk("years 2009-2026 complete in f1_regions", sorted(regy.year.unique()) == YEARS)
ncell = dens.cell.nunique()
chk("cell_density = cells x years x 2 rosters",
    len(dens) == ncell * len(YEARS) * 2, f"{ncell} cells")
chk("regions table = 7 regions x years x 2 rosters",
    len(regy) == 7 * len(YEARS) * 2)
chk("region names identical in all three region tables",
    sorted(regy.region.unique()) == REGIONS == sorted(rsta.region.unique())
    == sorted(dens.region.unique()))
chk("n60 <= n100 <= n150 in every cell-year", bool(((dens.n60 <= dens.n100)
                                                    & (dens.n100 <= dens.n150)).all()))
chk("n60_bb <= n60 and n100_bb <= n100", bool(((dens.n60_bb <= dens.n60)
                                               & (dens.n100_bb <= dens.n100)).all()))
f = dens[dens.roster == "fleet"].set_index(["cell", "year"]).sort_index()
a = dens[dens.roster == "all"].set_index(["cell", "year"]).sort_index()
chk("fleet counts <= all-FDSN counts in every cell-year",
    bool((f.n100 <= a.n100).all() and (f.n60 <= a.n60).all()))
chk("az_gap in [0,360]", bool(dens.az_gap_deg.between(0, 360).all()
                              and dens.az_gap_60_deg.between(0, 360).all()))
chk("cells assigned to a region carry a finite slab depth 30-45 km",
    bool(dens.slab_km.between(30, 45).all()))

# region membership vs the per-year member count in f1_regions
bad = []
for reg in REGIONS:
    mem = set(zip(rsta[rsta.region == reg].net, rsta[rsta.region == reg].sta))
    for tag in ("all", "fleet"):
        for y in YEARS:
            live = sy[(sy.year == y) & sy.has_3c]
            if tag == "fleet":
                live = live[live.roster == "project"]
            n = len(set(zip(live.net, live.sta)) & mem)
            r = regy[(regy.region == reg) & (regy.year == y) & (regy.roster == tag)]
            if int(r.n_member_stations.iloc[0]) != n:
                bad.append((reg, tag, y, int(r.n_member_stations.iloc[0]), n))
chk("f1_regions.n_member_stations reproducible from f1_region_stations x f1_station_years",
    not bad, f"{len(bad)} mismatches {bad[:3]}")

# region x year fractions reproducible from the cell table
bad2 = []
for reg in REGIONS:
    for tag in ("all", "fleet"):
        for y in YEARS:
            d = dens[(dens.region == reg) & (dens.year == y) & (dens.roster == tag)]
            r = regy[(regy.region == reg) & (regy.year == y) & (regy.roster == tag)].iloc[0]
            for col, val in (("frac_ge4_60", (d.n60 >= 4).mean()),
                             ("frac_ge6_60", (d.n60 >= 6).mean()),
                             ("frac_ge3", (d.n100 >= 3).mean()),
                             ("n_cells", len(d))):
                if not np.isclose(float(r[col]), float(val)):
                    bad2.append((reg, tag, y, col, float(r[col]), float(val)))
chk("f1_regions fractions reproducible from f1_cell_density", not bad2,
    f"{len(bad2)} mismatches {bad2[:3]}")

# cost table vs station_years
for scope, sub in (("project_fleet", sy[sy.roster == "project"]), ("project_plus_extra", sy)):
    tot = cost[(cost.scope == scope) & (cost.year == "TOTAL")].iloc[0]
    styr = float(sub[sub.has_3c].days.clip(0, 366).sum() / 365.0)
    chk(f"cost TOTAL station-years matches f1_station_years [{scope}]",
        abs(float(tot.station_years_effective) - styr) < 0.2,
        f"{tot.station_years_effective} vs {styr:.1f}")
    chk(f"cost TOTAL n_stations_3c matches f1_station_years [{scope}]",
        int(tot.n_stations_3c) == sub[sub.has_3c].groupby(["net", "sta"]).ngroups)
    yr = cost[(cost.scope == scope) & (cost.year != "TOTAL")]
    chk(f"cost per-year rows sum to TOTAL [{scope}]",
        abs(yr.raw_download_GB_native.sum() - float(tot.raw_download_GB_native)) < 1.0)

chk("ground-truth rows cover all 7 regions x every catalog",
    gt.groupby("catalog").region.nunique().eq(7).all())
print(f"\n{sum(ok)}/{len(ok)} checks PASS")

print("\n" + "=" * 96)
print("REPORT TABLES")
print("=" * 96)

print("\n### (i) station table — project roster by network and class")
pj = sta[sta.roster == "project"]
print(pd.crosstab(pj.net, pj.klass, margins=True).to_string())
print("\n  project stations by role")
print(pj.role.value_counts().to_string())
print("\n  project stations that already produced a dv/v product / causality cert / inversion pair")
print(pj[["has_dvv", "has_cert", "in_inversion"]].sum().to_string())
print(f"\n  non-project ('extra') FDSN 3-C stations kept: {int((sta.roster == 'extra').sum())}")
print(pd.crosstab(sta[sta.roster == "extra"].net,
                  sta[sta.roster == "extra"].klass, margins=True).to_string())

print("\n### (i) 3-C >=20 Hz station-years per year, by roster and class")
sy3 = sy[sy.has_3c].copy()
print(sy3.pivot_table(index="year", columns=["roster", "klass"], values="sta",
                      aggfunc="count").fillna(0).astype(int).to_string())
print("\n  effective station-years (days/365) and orientation mix per year [project]")
p = sy3[sy3.roster == "project"]
t = p.groupby("year").agg(n_stations=("sta", "count"),
                          eff_station_years=("days", lambda s: round(s.clip(0, 366).sum() / 365, 1)),
                          med_days=("days", "median"))
t["n_ZNE"] = p[p.orient == "ZNE"].groupby("year").size()
t["n_Z12"] = p[p.orient == "Z12"].groupby("year").size().reindex(t.index).fillna(0).astype(int)
t["n_borehole"] = p[p.is_borehole].groupby("year").size().reindex(t.index).fillna(0).astype(int)
print(t.to_string())

print("\n### (ii) source-cell station density per year [project fleet]")
# NB pivot_table defaults to MEAN; these aggregate over the cells of a region, so the
# statistic must be stated explicitly (f1_regions.med_n* are medians).
for col in ("n60", "n100", "n150"):
    for how in ("median", "mean"):
        print(f"\n  {how} {col} per cell")
        print(dens[dens.roster == "fleet"].pivot_table(index="year", columns="region",
                                                       values=col, aggfunc=how)
              .round(1).to_string())

print("\n### (iii) region x year locatable fractions [project fleet]")
for col in ("frac_ge3_60", "frac_ge4_60", "frac_ge6_60", "frac_ge4_60_bb",
            "frac_ge4_60_gap180"):
    print(f"\n  {col}")
    print(regy[regy.roster == "fleet"].pivot_table(index="year", columns="region",
                                                   values=col).round(2).to_string())

print("\n### (iii) region summary, median over 2019-2026 (the dense era)")
for tag in ("fleet", "all"):
    print(f"\n  [{tag}]")
    s = regy[(regy.roster == tag) & regy.year.between(2019, 2026)].groupby("region").agg(
        n_cells=("n_cells", "median"), n_stations=("n_member_stations", "median"),
        med_n60=("med_n60", "median"), med_n60_bb=("med_n60_bb", "median"),
        med_n100=("med_n100", "median"), med_nearest_km=("med_nearest_km", "median"),
        med_az_gap_60=("med_az_gap_60", "median"), med_az_gap_100=("med_az_gap", "median"),
        bore_frac_60=("med_borehole_frac_60", "median"),
        frac_ge4_60=("frac_ge4_60", "median"), frac_ge6_60=("frac_ge6_60", "median"),
        frac_ge4_60_bb=("frac_ge4_60_bb", "median"))
    print(s.round(2).to_string())

print("\n### (iv) cost")
with pd.option_context("display.width", 250, "display.max_colwidth", 400):
    print(cost[cost.year != "TOTAL"].drop(columns=["note"]).to_string(index=False))
    print()
    print(cost[cost.year.isin(["TOTAL", "ASSOC", "LOCATE"])].to_string(index=False))

print("\n### (v) ground truth per region")
print(gt.pivot_table(index="region", columns="catalog", values="n_events_or_families",
                     aggfunc="sum").fillna(0).astype(int).to_string())
print("\n  catalog totals and notes")
print(gt.drop_duplicates("catalog")[["catalog", "source", "n_total", "note"]].to_string(index=False))

# ---------------------------------------------------------------- ground-truth detail
print("\n### (v) ground-truth detail: quality cuts and where the catalogs actually sit")
ext = "/home/jovyan/tremorferometry/lfe_catalog_repair/results/external"
lin = pd.read_csv(f"{ext}/lin_catalog/EQloc_001_0.1_3_S.csv",
                  usecols=["lat", "lon", "depth", "residual", "N"])
print(f"  Lin2024: {len(lin)} events; lat {lin.lat.min():.2f}-{lin.lat.max():.2f}, "
      f"lon {lin.lon.min():.2f}-{lin.lon.max():.2f}")
q = lin[(lin.residual < 0.5) & (lin.N >= 4)]
print(f"  Lin2024 after the paper's own quality cut (residual<0.5 & N>=4): {len(q)} "
      f"({100*len(q)/len(lin):.1f}%)")
bins = [(n, a, b) for n, a, b in
        [("R1_N_VancIsl", 48.70, 50.0), ("R2_S_VancIsl_JdF", 47.90, 48.70),
         ("R3_Puget_Olympic", 46.60, 47.90), ("R4_N_Oregon", 45.20, 46.60),
         ("R5_C_Oregon", 43.60, 45.20), ("R6_S_Oregon", 42.00, 43.60),
         ("R7_NorCal_MTJ", 40.00, 42.00)]]
for n, a, b in bins:
    m = (lin.lat >= a) & (lin.lat < b)
    mq = (q.lat >= a) & (q.lat < b)
    print(f"    {n:18s} all {int(m.sum()):8d}   quality-cut {int(mq.sum()):8d}")
print("  Lin lat histogram (0.25 deg, top 8 bins):")
h = (lin.lat // 0.25 * 0.25).value_counts().head(8).sort_index()
print("    " + ", ".join(f"{k:.2f}N:{v}" for k, v in h.items()))

bos = pd.read_csv(f"{ext}/bostock/bostock2015_sloweq_FULL_2003-2013.csv", usecols=["lat", "lon"])
print(f"  Bostock2015: {len(bos)} events, {len(bos.drop_duplicates(['lat','lon']))} distinct "
      f"(lat,lon) [published family count 130 -> some families share rounded hypocentres]; "
      f"lat {bos.lat.min():.2f}-{bos.lat.max():.2f}")
duc = pd.read_csv("/home/jovyan/tremorferometry/catalogs/ducellier_catalog/data/Plourde_2015/"
                  "templates_list.txt", sep=r"\s+", skiprows=1,
                  names=["temp", "fam", "lat", "lon", "dep", "eh", "ez", "n"]
                  ).drop_duplicates("fam")
print(f"  Ducellier/Plourde located families: {len(duc)}; lat {duc.lat.min():.2f}-"
      f"{duc.lat.max():.2f} -> {int(((duc.lat >= 40) & (duc.lat < 42)).sum())} inside R7, "
      f"{int((duc.lat >= 42).sum())} north of 42N, {int((duc.lat < 40).sum())} south of 40N")
sb = pd.read_csv(f"{ext}/other/seisbench_lfestacks_cascadia/metadata.csv", low_memory=False)
print(f"  SeisBench LFE stacks: {len(sb)} traces, "
      f"{sb.source_family_id.nunique()} families, "
      f"{sb.station_code.nunique()} stations; split "
      f"{sb.split.value_counts().to_dict() if 'split' in sb else 'n/a'}")

# ------------------------------------------------- (vi) per-region readiness summary
print("\n### (vi) region readiness — first year each criterion is met [project fleet]")
rows = []
for reg in REGIONS:
    d = regy[(regy.roster == "fleet") & (regy.region == reg)].sort_values("year")
    da = regy[(regy.roster == "all") & (regy.region == reg)].sort_values("year")

    def first(dd, col, thr):
        s = dd[dd[col] >= thr]
        return int(s.year.iloc[0]) if len(s) else None

    m = rsta[rsta.region == reg]
    rows.append(dict(
        region=reg, n_cells=int(d.n_cells.iloc[0]),
        member_all=len(m), member_fleet=int((m.roster == "project").sum()),
        member_fleet_bb=int(((m.roster == "project")
                             & m.klass.isin(["broadband", "broadband_borehole"])).sum()),
        y_ge4_60_all_cells=first(d, "frac_ge4_60", 0.999),
        y_ge6_60_90pct=first(d, "frac_ge6_60", 0.90),
        y_ge4_60_bb_90pct=first(d, "frac_ge4_60_bb", 0.90),
        y_ge6_60_90pct_allFDSN=first(da, "frac_ge6_60", 0.90),
        frac_ge6_60_2012=float(d[d.year == 2012].frac_ge6_60.iloc[0]),
        frac_ge6_60_2024=float(d[d.year == 2024].frac_ge6_60.iloc[0])))
print(pd.DataFrame(rows).to_string(index=False))
print("\n  member-station class mix per region (project fleet only)")
mm = rsta[rsta.roster == "project"]
print(pd.crosstab(mm.region, mm.klass, margins=True).to_string())
