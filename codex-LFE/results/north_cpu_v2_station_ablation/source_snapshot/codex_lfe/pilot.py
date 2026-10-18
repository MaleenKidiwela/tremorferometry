"""Reproducible CPU/GPU-capable pilot; no imports from the old LFE implementation.

Run from codex-LFE: python -m codex_lfe.pilot --config configs/cpu_pilot.yaml
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import time

import numpy as np
import pandas as pd
import yaml

from .evaluate import null_summary, score_catalog
from .matched_filter import (align_to_anchor, combine_sites, extract_candidates, normalized_cc,
                             resolve_candidates, shifted_sites)
from .provenance import environment, sha256, source_manifest, write_json
from .reference import choose_families, read_bostock
from .templates import build_station_template, cut
from .waveforms import load_day, runs

ROOT = Path(__file__).resolve().parents[1]
COLUMNS = ["time_s", "score_mad", "network_cc", "n_sites", "n_supporting_sites", "supporting_sites",
           "family_id", "time_basis", "realization", "threshold", "date", "split", "block"]


def validate_config(cfg):
    sets = [set(cfg[k]) for k in ["train_dates", "calibration_dates", "test_dates"]]
    if any(sets[i] & sets[j] for i in range(3) for j in range(i)):
        raise ValueError("training, calibration and test dates overlap")
    if any(not s for s in sets):
        raise ValueError("all three date splits required")
    if cfg['block_guard_s'] <= cfg['null_shift_max_s'] + cfg['refractory_s']:
        raise ValueError("block guard must exceed maximum shift plus refractory interval")
    if 2 * cfg['block_guard_s'] >= cfg['normalization_block_s']:
        raise ValueError("guard consumes entire normalization block")
    if len({s['site_id'] for s in cfg['stations'].values()}) < cfg['min_sites']:
        raise ValueError("insufficient independent sites")
    if cfg['anchor_station'] not in cfg['stations']:
        raise ValueError("anchor station is not in network")


def date_start(date):
    return pd.Timestamp(date, tz="UTC").timestamp()


def load_waves(cfg, date, manifests):
    waves = {}
    dt = pd.Timestamp(date)
    for station, meta in cfg['stations'].items():
        path = (ROOT / cfg['raw_root'] / station / str(dt.year) / f"{dt.dayofyear:03d}.mseed").resolve()
        if not path.exists():
            raise FileNotFoundError(f"Missing explicitly requested raw input: {path}")
        wave, manifest = load_day(path, ROOT / 'cache' / 'waveforms', meta['channels'], meta['site_id'],
                                  date_start(date), fs=cfg['fs'], alignment=cfg.get('fractional_alignment', 'linear'))
        waves[station] = wave
        manifests[f"{date}/{station}"] = manifest | {"valid_fraction": wave.valid.mean(axis=1).tolist()}
    return waves


def build_templates(cfg, reference, families, out, manifests):
    windows = {(family, station): [] for family in families for station in cfg['stations']}
    training_ids = []
    for date in cfg['train_dates']:
        print(f"[train] processing raw {date}", flush=True)
        waves = load_waves(cfg, date, manifests)
        d = reference[(reference.date == date) & reference.family_id.isin(families)]
        training_ids.extend(d.reference_id.tolist())
        for family, frame in d.groupby('family_id'):
            for station, wave in waves.items():
                for t in frame.time_s:
                    window = cut(wave, t - 10, int(70 * cfg['fs']),
                                 alignment=cfg.get('template_alignment', 'linear'))
                    if window is not None:
                        windows[family, station].append(window)
    templates, qc, anchors, archive = {}, [], {}, {}
    for family in families:
        temp = {}
        for station in cfg['stations']:
            template, row = build_station_template(windows[family, station], cfg['fs'],
                template_s=cfg['template_s'], min_events=cfg['min_training_events'],
                min_split_cc=cfg['min_training_split_cc'])
            row.update(family_id=family, station=station)
            qc.append(row)
            if template is not None:
                temp[station] = (template, row)
        anchor = temp.get(cfg['anchor_station'])
        if anchor is None or len({cfg['stations'][s]['site_id'] for s in temp}) < cfg['min_sites']:
            print(f"[train] {family}: insufficient coherent sites/anchor; excluded from detector", flush=True)
            continue
        anchors[family] = anchor[1]['peak_offset_from_published_s']
        templates[family] = {}
        for station, (template, row) in temp.items():
            offset = row['template_start_from_published_s'] - anchors[family]
            templates[family][station] = (template, offset)
            archive[f"{family}__{station}__template"] = template
            archive[f"{family}__{station}__offset_s"] = np.asarray(offset)
        print(f"[train] {family}: {len(temp)} stations; anchor offset {anchors[family]:.2f}s", flush=True)
    pd.DataFrame(qc).to_csv(out / 'template_qc.csv', index=False)
    np.savez(out / 'templates.npz', **archive)
    write_json(out / 'training.json', {"reference_ids": training_ids, "anchors_from_published_s": anchors,
                                      "families_selected": families, "families_usable": list(templates),
                                      "time_basis": f"arrival_anchor:{cfg['anchor_station']}:training_stack_peak",
                                      "origin_time_validated": False})
    return templates, anchors


def scan_date(cfg, date, split, templates, out, manifests, backend, device):
    print(f"[{split}] processing raw {date}", flush=True)
    waves = load_waves(cfg, date, manifests)
    fs = cfg['fs']
    n = int(86400 * fs)
    site_series = {}
    for family, members in templates.items():
        aligned = {}
        for station, (template, offset) in members.items():
            cc, valid = normalized_cc(waves[station].x, template, waves[station].valid, backend, device)
            aligned[station], _ = align_to_anchor(cc, valid, offset, fs, n)
        site_series[family] = combine_sites(aligned, {s: cfg['stations'][s]['site_id'] for s in aligned})
    del waves
    block_n = int(cfg['normalization_block_s'] * fs)
    guard = int(cfg['block_guard_s'] * fs)
    all_rows, exposures, intervals = [], [], []
    basis = f"arrival_anchor:{cfg['anchor_station']}:training_stack_peak"
    for block, first in enumerate(range(0, n, block_n)):
        last = min(n, first + block_n)
        core = np.zeros(last - first, bool)
        core[guard:len(core) - guard] = True
        # Shared offsets across ALL families; controls preserve intra-site coherence.
        all_sites = sorted({site for sites in site_series.values() for site in sites})
        rng = np.random.default_rng(np.random.SeedSequence([cfg['seed'], int(date.replace('-', '')), block]))
        shifts = [{s: 0 for s in all_sites}]
        for _ in range(cfg['null_replicates']):
            shifts.append({s: int(rng.choice([-1, 1]) * rng.integers(
                int(cfg['null_shift_min_s'] * fs), int(cfg['null_shift_max_s'] * fs) + 1)) for s in all_sites})
        block_rows = {t: {r: [] for r in range(len(shifts))} for t in cfg['thresholds']}
        for family, full_sites in site_series.items():
            local = {s: v[first:last] for s, v in full_sites.items()}
            for realization, shift in enumerate(shifts):
                sites = local if realization == 0 else shifted_sites(local, shift)
                # Statistics and ranked peak suppression are threshold-independent.
                # Compute once at the lowest threshold; preserve separate final
                # cross-family assignment at each threshold. Regression-tested
                # against independent calls at every operating point.
                base_rows, good = extract_candidates(sites, fs, date_start(date) + first / fs,
                    min(cfg['thresholds']), support_mad=cfg['site_support_mad'], min_sites=cfg['min_sites'],
                    refractory_s=cfg['refractory_s'], common_valid=core)
                for threshold in cfg['thresholds']:
                    rows = [row.copy() for row in base_rows if row['score_mad'] >= threshold]
                    for row in rows:
                        row.update(family_id=family, time_basis=basis, realization=realization,
                                   threshold=threshold, date=date, split=split, block=block)
                    block_rows[threshold][realization].extend(rows)
                    if threshold == cfg['thresholds'][0]:
                        exposures.append(dict(date=date, split=split, block=block, family_id=family,
                                              realization=realization, exposure_s=float(good.sum() / fs)))
                        if realization == 0:
                            for lo, hi in runs(good):
                                intervals.append(dict(family_id=family, start_s=date_start(date) + (first + lo) / fs,
                                                      end_s=date_start(date) + (first + hi) / fs, date=date, split=split))
        for threshold, realizations in block_rows.items():
            for realization, rows in realizations.items():
                resolved = resolve_candidates(pd.DataFrame(rows, columns=COLUMNS), cfg['cross_family_tolerance_s'],
                                               cfg['family_ambiguity_margin_mad'])
                all_rows.extend(resolved.to_dict('records'))
    candidates = pd.DataFrame(all_rows)
    if candidates.empty:
        candidates = pd.DataFrame(columns=COLUMNS)
    candidates.to_parquet(out / f'{split}_{date}_candidates.parquet', index=False)
    pd.DataFrame(exposures).to_parquet(out / f'{split}_{date}_exposure.parquet', index=False)
    pd.DataFrame(intervals, columns=['family_id', 'start_s', 'end_s', 'date', 'split']).to_parquet(
        out / f'{split}_{date}_eligible_intervals.parquet', index=False)
    print(f"[{split}] {date}: {len(candidates)} final candidate rows across thresholds + controls", flush=True)
    return candidates, pd.DataFrame(exposures), pd.DataFrame(intervals)


def summarize_split(cfg, split, candidates, exposure, intervals, reference, families, anchors):
    ref = reference[reference.date.isin(cfg[f'{split}_dates']) & reference.family_id.isin(families)].copy()
    ref['time_s'] += ref.family_id.map(anchors).fillna(0)
    ref['time_basis'] = f"arrival_anchor:{cfg['anchor_station']}:training_stack_peak"
    eligible = np.zeros(len(ref), bool)
    for family, ii in ref.groupby('family_id').indices.items():
        times = ref.iloc[ii].time_s.to_numpy()
        ints = intervals[intervals.family_id == family].sort_values('start_s') if len(intervals) else intervals
        if not len(ints):
            continue
        starts, ends = ints.start_s.to_numpy(), ints.end_s.to_numpy()
        j = np.searchsorted(starts, times, side='right') - 1
        eligible[ii] = (j >= 0) & (times < ends[np.maximum(j, 0)])
    ref = ref[eligible]
    exp = exposure.groupby('realization').exposure_s.sum().to_dict() if len(exposure) else {}
    results = {}
    for threshold in cfg['thresholds']:
        d = candidates[(candidates.threshold == threshold) & (candidates.realization == 0)]
        controls = [int(((candidates.threshold == threshold) & (candidates.realization == r)).sum())
                    for r in range(1, cfg['null_replicates'] + 1)]
        null = null_summary(len(d), exp.get(0, 0), controls,
                            [exp.get(r, 0) for r in range(1, cfg['null_replicates'] + 1)]) if exp.get(0, 0) > 0 and all(exp.get(r, 0) > 0 for r in range(1, cfg['null_replicates'] + 1)) else None
        results[f'{threshold:g}'] = {'reference': score_catalog(d, ref, cfg['reference_tolerance_s']),
                                     'null': null, 'n_ambiguous': int(d.family_id.isna().sum()),
                                     'reference_total_selected_groups': int(len(eligible)),
                                     'reference_eligible': int(eligible.sum()),
                                     'reference_excluded': int((~eligible).sum()),
                                     'exposure_unit': 'sum of family-seconds; merged-catalog null ratio is approximate if family exposure mix differs'}
    return results, ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', default='configs/cpu_pilot.yaml')
    ap.add_argument('--output', default=None)
    ap.add_argument('--backend', choices=['numpy', 'torch'], default='numpy')
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--workers', type=int, default=1, help='independent CPU station-day groups; test waits for frozen calibration')
    args = ap.parse_args()
    cfg_path = Path(args.config).resolve()
    cfg = yaml.safe_load(cfg_path.read_text())
    validate_config(cfg)
    if args.workers < 1 or args.workers > 8:
        raise ValueError('workers must be between 1 and 8')
    if args.device != 'cpu' and args.workers != 1:
        raise ValueError('use one worker per GPU until concurrent GPU memory use is profiled')
    # Resolve read-only source paths in the manifest so archived code can be rerun
    # from a different directory without silently addressing different inputs.
    cfg['raw_root'] = str((ROOT / cfg['raw_root']).resolve())
    cfg['reference_csv'] = str((ROOT / cfg['reference_csv']).resolve())
    if args.device != 'cpu':
        import torch
        if args.backend != 'torch' or not torch.cuda.is_available():
            raise RuntimeError('CUDA requested but unavailable, or backend is not torch')
    out = Path(args.output) if args.output else ROOT / 'results' / datetime.now(timezone.utc).strftime('cpu_pilot_%Y%m%dT%H%M%SZ')
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'output must be empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    manifests = {}
    refpath = (ROOT / cfg['reference_csv']).resolve()
    manifest = {'config': cfg, 'config_sha256': sha256(cfg_path), 'environment': environment(),
                'source_code': source_manifest(ROOT / 'codex_lfe'), 'reference_path': str(refpath),
                'reference_sha256': sha256(refpath), 'backend': args.backend, 'device': args.device,
                'workers': args.workers,
                'created_utc': datetime.now(timezone.utc).isoformat(), 'status': 'running'}
    write_json(out / 'manifest.json', manifest)
    for name, expected_hash in manifest['source_code'].items():
        src = ROOT / 'codex_lfe' / name
        if sha256(src) != expected_hash:
            raise RuntimeError('Source changed while taking experiment snapshot')
        dst = out / 'source_snapshot' / 'codex_lfe' / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
    write_json(out / 'source_snapshot' / 'resolved_config.json', cfg)
    reference = read_bostock(refpath)
    families = choose_families(reference, cfg['train_dates'], cfg['family_count'], cfg['family_separation_km'])
    templates, anchors = build_templates(cfg, reference, families, out, manifests)
    write_json(out / 'input_manifest.json', manifests)
    if not templates:
        manifest.update(status='failed_template_qc', elapsed_s=time.monotonic() - start)
        write_json(out / 'manifest.json', manifest)
        raise RuntimeError('No families passed training-only QC; no valid pilot detector')
    selection = None
    report = {'catalog_status': 'experimental_candidates_only', 'splits': {}}
    for split in ['calibration', 'test']:
        def run_date(date):
            return scan_date(cfg, date, split, templates, out, manifests, args.backend, args.device)

        # Date-local random seeds and disjoint outputs preserve results regardless
        # of scheduling. All calibration futures finish and the selection file is
        # written before any held-out date is submitted.
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            pieces = list(pool.map(run_date, cfg[f'{split}_dates']))
        c, e, i = [pd.concat([p[k] for p in pieces], ignore_index=True) for k in range(3)]
        summary, eligible_reference = summarize_split(cfg, split, c, e, i, reference, families, anchors)
        report['splits'][split] = summary
        eligible_reference.to_parquet(out / f'{split}_eligible_reference.parquet', index=False)
        write_json(out / f'{split}_scores.json', summary)
        write_json(out / 'input_manifest.json', manifests)
        if split == 'calibration':
            passing = [t for t in cfg['thresholds'] if summary[f'{t:g}']['null'] is not None and
                       summary[f'{t:g}']['null']['observed_count'] >= cfg['calibration_min_detections'] and
                       summary[f'{t:g}']['null']['control_to_observed_range'][1] <= cfg['calibration_max_control_ratio']]
            selection = min(passing) if passing else None
            report['selection'] = {'threshold': selection, 'diagnostic_threshold': cfg['diagnostic_threshold'],
                                   'calibration_gate_passed': selection is not None,
                                   'frozen_before_test': True,
                                   'warning': 'Passing empirical controls alone does not certify an LFE catalog.'}
            write_json(out / 'selection_before_test.json', report['selection'])
        chosen = selection if selection is not None else cfg['diagnostic_threshold']
        d = c[(c.threshold == chosen) & (c.realization == 0)].copy()
        d['event_id'] = [f'{split}_{k:07d}' for k in range(len(d))]
        d['quality_status'] = 'experimental_unverified'
        d['time_utc'] = pd.to_datetime(d.time_s, unit='s', utc=True)
        d.to_parquet(out / f'{split}_experimental_candidates.parquet', index=False)
    report['elapsed_s'] = time.monotonic() - start
    write_json(out / 'report.json', report)
    manifest.update(status='complete_experiment_not_validated_catalog', elapsed_s=report['elapsed_s'])
    write_json(out / 'manifest.json', manifest)
    print(f'[done] {out}; selection={selection}; elapsed={report["elapsed_s"]:.1f}s', flush=True)


if __name__ == '__main__':
    main()
