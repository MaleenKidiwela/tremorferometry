"""Post-fit calibration diagnostics; never changes frozen catalog scores."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .pilot import load_waves, date_start
from .matched_filter import normalized_cc, align_to_anchor, combine_sites, robust_scale
from .provenance import write_json, sha256


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    run, out = Path(args.run), Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('output must be empty')
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((run / 'manifest.json').read_text())['config']
    refs = pd.read_parquet(run / 'calibration_eligible_reference.parquet')
    archive = np.load(run / 'templates.npz')
    families = json.loads((run / 'training.json').read_text())['families_usable']
    manifests, rows = {}, []
    fs = cfg['fs']; block_n = int(cfg['normalization_block_s'] * fs)
    guard = int(cfg['block_guard_s'] * fs)
    radius = int(round(cfg['reference_tolerance_s'] * fs))
    for date in cfg['calibration_dates']:
        print('[diagnostic]', date, flush=True)
        waves = load_waves(cfg, date, manifests)
        n = int(86400 * fs)
        for family in families:
            aligned = {}
            for station in cfg['stations']:
                key = f'{family}__{station}'
                if key + '__template' not in archive:
                    continue
                cc, valid = normalized_cc(waves[station].x, archive[key + '__template'], waves[station].valid)
                aligned[station], _ = align_to_anchor(cc, valid, float(archive[key + '__offset_s']), fs, n)
            sites = combine_sites(aligned, {s: cfg['stations'][s]['site_id'] for s in aligned})
            data = np.asarray(list(sites.values()))
            sub = refs[(refs.date == date) & (refs.family_id == family)]
            for block, group in sub.groupby(((sub.time_s - date_start(date)) * fs // block_n).astype(int)):
                first = block * block_n; x = data[:, first:first + block_n]
                good = np.isfinite(x).all(axis=0)
                good[:guard] = False; good[-guard:] = False
                stats = [robust_scale(v[good]) for v in x]
                scale = robust_scale(x[:, good].mean(axis=0))
                if scale is None or any(s is None for s in stats):
                    continue
                z = (x.mean(axis=0) - scale[0]) / scale[1]
                support = np.sum([(v - s[0]) / s[1] >= cfg['site_support_mad'] for v, s in zip(x, stats)], axis=0)
                for ref in group.itertuples():
                    idx = int(round((ref.time_s - date_start(date)) * fs)) - first
                    ix = np.arange(max(0, idx-radius), min(len(z), idx+radius+1))
                    ix = ix[good[ix]]
                    if not len(ix):
                        continue
                    peak = ix[np.argmax(z[ix])]
                    supported = ix[support[ix] >= cfg['min_sites']]
                    rows.append(dict(reference_id=ref.reference_id, family_id=family, date=date,
                        max_network_mad=float(z[peak]), support_at_network_max=int(support[peak]),
                        max_supported_mad=float(z[supported].max()) if len(supported) else None,
                        n_sites=len(sites), lag_at_network_max_s=float((peak-idx)/fs)))
    frame = pd.DataFrame(rows)
    frame.to_parquet(out / 'calibration_reference_scores.parquet', index=False)
    summary = []
    for family, g in frame.groupby('family_id'):
        summary.append(dict(family_id=family, n_reference=len(g),
            median_max_network_mad=float(g.max_network_mad.median()),
            no_supported_sample=int(g.max_supported_mad.isna().sum()),
            network_at_least_9=int((g.max_network_mad >= 9).sum()),
            supported_at_least_9=int((g.max_supported_mad >= 9).sum())))
    pd.DataFrame(summary).to_csv(out / 'family_summary.csv', index=False)
    write_json(out / 'manifest.json', dict(detector_manifest_sha256=sha256(run/'manifest.json'),
        source_sha256=sha256(__file__), inputs=manifests, split='calibration_only',
        interpretation='Reference-centered maxima within existing match tolerance. Diagnostic upper bounds before peak suppression and family resolution, not recovered detections. Threshold 9 is the existing frozen operating point; no retuning performed.'))

if __name__ == '__main__':
    main()
