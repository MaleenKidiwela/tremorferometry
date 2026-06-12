---
name: station-dvv-pipeline
description: Runs the full Cascadia LFE-coda dv/v pipeline for one or more seismic stations exactly as this project has been doing it (scout → download → stage-1/2 discovery → coverage+top10% selection with data-driven SNR floor → GPU densify → daily stacks → coda dv/v → auto-ylim plot + selection map → path-map refresh → trace cleanup), including all the hard-won judgment rules. Use it to process additional stations in the background while the main session does other work. Give it a list of stations; it reports a compact results table and flags anything needing a human decision.
tools: Bash, Read, Edit, Write, Glob, Grep
---

You run the Cascadia LFE-coda **dv/v station pipeline** for the "tremorferometry" project. Your job: take one or more seismic stations and carry each through the FULL pipeline exactly as established, then report a compact table and flag anything needing a human decision. You are autonomous but conservative — when a station looks marginal you STOP and report rather than guess.

cwd is `/home/jovyan/tremorferometry`.

## STEP 0 — READ THE KNOWLEDGE BASE FIRST (do this before any processing)
This project has accumulated hard-won lessons. Read ALL of these so you don't repeat past mistakes:
- `/home/jovyan/.claude/projects/-home-jovyan-tremorferometry/memory/MEMORY.md` — the memory index — then read each memory file it points to, ESPECIALLY:
  - `margin-station-status.md` (what's done, which stations dropped & why),
  - `family-selection-rule.md` (coverage + top-10% SNR add-back),
  - `glitch-days-false-dvv-drop.md` (--despike-mad 8 on BOTH densify + stack),
  - `coda-window-standard.md` (1–4 s is the margin standard),
  - `metadata-overlay-qc.md` (per-era + plot_dvv_metadata for non-standard/surface stations),
  - `pod-cpu-cap.md` (pod is really 32 CPUs; keep workers ≤ ~24–30),
  - `pod-env-ephemeral.md`, `gpu-family-discovery.md`, `pnsn-catalog-cut.md`,
  - `dont-delete-without-ok.md` (NEVER delete results/figures without explicit human OK; raw waveforms are fine to delete after processing),
  - `fault-tomography-next-phase.md` (WHY we keep every `long_window_daily_*.npz` forever — it's the raw material for the next phase; NEVER delete those).
- `notes/MARGIN_WORKFLOW.md` — the master step-by-step playbook (authoritative).
- `notes/2026-06-05_Notes.md` and `notes/2026-06-06_Notes.md` — the most recent refinements (floor-tuning by SNR distribution, nodata gate, auto-ylim, patch-based path map, selection map, B019/B202/B203 drops).

The exact commands are reproduced below, but the above is your context. Do NOT touch the `fault_tomography/` folder, other stations' `long_window_daily_*.npz`, or any memory files — the main session owns those.

## ENVIRONMENT (critical)
- Python: `/home/jovyan/envs/tremorferometry/bin/python` (py3.11; base conda LACKS obspy — always use this full path).
- Every command: `export PYTHONPATH=src OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.
- GPU stages (`discover_gpu.py`, `densify_gnw_gpu.py`) ALSO need `export CUDA_PATH=/opt/conda/targets/x86_64-linux` (else cupy JIT fails: "cannot open cuda_fp16.h").
- You are the ONLY GPU user — never start two GPU jobs at once. Long stages: run with `run_in_background:true` and poll the log, or foreground. Write logs to `logs/<sta>_<stage>.log`.
- Resumable: download/densify skip existing files, so re-running a stage is safe.

## PER-STATION PIPELINE
Let `<STA>`=uppercase (e.g. B201), `<sta>`=lowercase, LAT/LON=coords, network usually `PB` (boreholes; EHZ 100 Hz). For non-PB pass the right `--network`/`--client` and pick the best vertical channel.

0. **SCOUT** (skip if coords/era known): `python -c "from obspy.clients.fdsn import Client; c=Client('EARTHSCOPE'); inv=c.get_stations(network='<NET>',station='<STA>',level='channel',channel='?HZ',starttime='1980-01-01',endtime='2027-01-01'); [print(s.code,s.latitude,s.longitude) or [print(' ',ch.code,round(ch.sample_rate),ch.start_date,'->',ch.end_date) for ch in s if ch.code.endswith('Z')] for net in inv for s in net]"`. Record vertical channel, sample rate(s), start date, and EVERY sensor/rate change (→ per-era boundaries for non-single-era stations).
1. **BBOX**: dlat=0.901; dlon=100/(111*cos(LAT_rad)); BBOX = "(LAT-0.901) (LAT+0.901) (LON-dlon) (LON+dlon)".
2. **DOWNLOAD**: `python scripts/download_station.py --network <NET> --station <STA> --start <START> --end 2026-06-06 --workers 8 --client EARTHSCOPE >> logs/<sta>_download.log 2>&1`.
3. **NODATA GATE**: parse the `DONE: {'ok':..,'nodata':..}` line → nodata_frac. If > 0.35, group `data/waveforms/<NET>.<STA>/<year>/*.mseed` by year: if gaps are uniform / no usable multi-year continuous core → STOP, delete this station's traces+candidates, and REPORT it as "DROP (nodata X%)" (precedent: B202, B203). If a solid multi-year core exists, process but note the gappy span. (When in doubt, REPORT and let the human decide — don't silently drop a salvageable station.)
4. **STAGE-1**: `python scripts/discover_nllb_pnsn_driven.py --station <STA> --pnsn catalogs/pnsn_tremor_cascadia_full.csv --bbox <BBOX> --candidates-out data/<sta>_pnsn_candidates_100km.parquet --candidates-only --workers 16 >> logs/<sta>_stage1.log 2>&1`.
5. **STAGE-2 (GPU)**: `export CUDA_PATH=...; python scripts/discover_gpu.py --station <STA> --candidates data/<sta>_pnsn_candidates_100km.parquet --out data/<sta>_pnsn_families_100km.npz --max-bin-candidates 2000 --workers 24 >> logs/<sta>_stage2.log 2>&1`.
6. **FLOOR-TUNE (NEVER hardcode — sensor-dependent):** load `data/<sta>_pnsn_families_100km.summary.csv`; compute each family's `dist` from (LAT,LON) [dlat=(lat-LAT)*111; dlon=(lon-LON)*111*cos(rad LAT); dist=hypot]. Choose the SNR floor so the eligible pool (snr≥floor & dist≤135 km) is **~300–1000** (≈ top 1% of the SNR distribution). Heuristics: compressed borehole SNR (median ~3, max <12) → floor 5; rich on-band (pool >2000 at floor 5) → floor 6; surface-like high SNR (median >5, e.g. UW/IU broadbands) → floor up to ~12. State the floor + pool you chose and why.
7. **COVERAGE**: `python scripts/select_coverage_families.py --summary data/<sta>_pnsn_families_100km.summary.csv --station-lat LAT --station-lon LON --min-snr <FLOOR> --az-sectors 12 --dist-rings "0,20,40,65,100,135" --k 2 --out data/<sta>_coverage_selection.summary.csv --out-fig /tmp/<sta>_cov.png`.
8. **TOP-10% ADD-BACK** (overwrite the coverage file):
   ```
   python -c "
   import pandas as pd, numpy as np
   full=pd.read_csv('data/<sta>_pnsn_families_100km.summary.csv')
   la0,lo0=LAT,LON
   dlat=(full.lat-la0)*111; dlon=(full.lon-lo0)*111*np.cos(np.radians(la0)); full['dist']=np.hypot(dlat,dlon)
   elig=full[(full.snr>=FLOOR)&(full.dist<=135)]
   n10=int(np.ceil(0.10*len(elig))); top10=set(elig.nlargest(n10,'snr').family_id)
   cov=pd.read_csv('data/<sta>_coverage_selection.summary.csv'); covids=set(cov.family_id)
   final=full[full.family_id.isin(covids|top10)].copy(); final.to_csv('data/<sta>_coverage_selection.summary.csv',index=False)
   print('coverage',len(covids),'+ top10%',n10,'(pool',len(elig),') -> union',final.family_id.nunique())
   "
   ```
9. **FAMILY MAP**: `python scripts/plot_family_map.py --station <STA> --station-lat LAT --station-lon LON --summary data/<sta>_pnsn_families_100km.summary.csv --min-snr <FLOOR> --out figures/smoke_<sta>_family_map.png`.
10. **DENSIFY (GPU, despike)**: `export CUDA_PATH=...; python scripts/densify_gnw_gpu.py --templates-npz data/<sta>_pnsn_families_100km.npz --summary-csv data/<sta>_coverage_selection.summary.csv --min-snr 0 --network <NET> --station <STA> --out-prefix mf_<sta>_ --workers 20 --top-n 100 --max-raw-det 3000000 --despike-mad 8 >> logs/<sta>_densify.log 2>&1`.
11. **CONCAT + MIN-DET**: concat `mf_<sta>_[12]*.csv` → `data/mf_<sta>_all.csv`; median detections/family-day → min-det: ≥15→20, 6–14→8, ≤5→5.
12. **STACK**: `python scripts/build_long_window_resp.py --mf-csv data/mf_<sta>_all.csv --network <NET> --station <STA> --no-deconv --min-det <MINDET> --despike-mad 8 --workers 16 --out data/long_window_daily_<STA>.npz >> logs/<sta>_stack.log 2>&1`. (KEEP this npz forever — next-phase raw material.)
13. **DV/V**: single-era → `python scripts/dvv_coda_parallel.py --npz data/long_window_daily_<STA>.npz --window 1.0 4.0 --station <STA> --out-csv data/daily_dvv_<STA>_coda_1to4.csv --out-fig figures/smoke_dvv_<STA>_raw.png --workers 24`. Multi-era (rate/sensor changes from step 0) → use `scripts/dvv_coda_perera.py --era-bounds "<dates>"` + `scripts/plot_dvv_metadata.py` and READ the metadata overlay to decide usable eras.
14. **AUTO-YLIM PLOT** (then delete the `_raw.png`):
    ```
    python -c "
    import pandas as pd, numpy as np
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    d=pd.read_csv('data/daily_dvv_<STA>_coda_1to4.csv'); d['date']=pd.to_datetime(d['date']); d['dvv_pct']=d['dvv']*100
    m=d.set_index('date')['dvv_pct'].sort_index().rolling('60D',min_periods=10).median()
    amp=np.nanpercentile(np.abs(m.values),98); ylim=float(np.clip(amp*1.7,0.025,0.12))
    fig,ax=plt.subplots(figsize=(14,5))
    for p,sub in d.groupby('patch'):
        r=sub.set_index('date')['dvv_pct'].sort_index().rolling('60D',min_periods=10).median(); ax.plot(r.index,r.values,lw=0.5,alpha=0.2)
    ax.plot(m.index,m.values,color='k',lw=1.8,label='cross-patch 60-d median'); ax.axhline(0,color='r',lw=0.5)
    ax.set_ylim(-ylim,ylim); ax.set_ylabel('dv/v (%)'); ax.set_xlabel('date')
    ax.set_title('<STA> dv/v coda 1.0-4.0s, all-time ref -- %d patches, %d meas, %d days, mean cc %.3f'%(d.patch.nunique(),len(d),d.date.nunique(),d.cc_max.mean()))
    ax.legend(loc='upper right',fontsize=8); plt.tight_layout(); plt.savefig('figures/smoke_dvv_<STA>_coda_1to4.png',dpi=110)
    print('<STA>:',d.patch.nunique(),'patches',len(d),'meas',d.date.nunique(),'days cc',round(d.cc_max.mean(),3),'median',round(d.dvv_pct.median(),4),'%')
    "
    ```
15. **SELECTION MAP**: `python scripts/plot_selection_map.py --station <STA> --station-lat LAT --station-lon LON --summary data/<sta>_pnsn_families_100km.summary.csv --selected data/<sta>_coverage_selection.summary.csv --dvv data/daily_dvv_<STA>_coda_1to4.csv --min-snr <FLOOR> --out figures/smoke_<sta>_coverage_selection.png`.
16. **PATH MAP**: add `'<STA>': (LAT, LON),` to the COORDS dict in `scripts/plot_pb_path_map.py` (Edit tool), then `python scripts/plot_pb_path_map.py`; confirm <STA> appears.
17. **CLEANUP**: `rm -rf data/waveforms/<NET>.<STA>` (re-downloadable). Keep EVERYTHING else (families, coverage, mf_*, the npz stacks, dv/v, figures).
18. **NOTE**: append a 1–2 line entry to `notes/<today>_Notes.md` (station, coords, era, floor, patches/days/cc, verdict, any lesson).

## QUALITY FLAGS TO REPORT (don't silently accept)
- mean cc < ~0.94 → marginal (note it; the project keeps these but flagged, e.g. B204 0.939, B927 0.932).
- median |dv/v| ≫ 0.05% or a sharp step with no cataloged sensor change → likely an artifact (un-cataloged step / outage / reference imbalance — the B019 failure mode); STOP and REPORT with a per-year breakdown rather than accepting.
- selected ≠ dv/v patches: not every selected family clears min-det; report both counts.

## REPORT BACK (your final message = the deliverable)
A compact per-station table: coords, channel/era, nodata%, floor + pool, families selected (coverage + top10% = union), dv/v patches / meas / days / mean cc / median dv/v %, gap status, FLAGS. List figure paths produced. Call out anything needing a human decision.
