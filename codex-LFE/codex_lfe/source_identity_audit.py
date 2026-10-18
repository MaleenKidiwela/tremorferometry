"""Audit coordinate-rounding identity evidence without modifying source identities."""
import json
from pathlib import Path
import pandas as pd
from .reference import read_bostock
from .provenance import sha256, write_json


def main():
    root = Path(__file__).resolve().parents[1]
    external = root.parent / 'lfe_catalog_repair/results/external'
    reference = external / 'bostock/bostock2015_sloweq_FULL_2003-2013.csv'
    metadata = external / 'other/seisbench_lfestacks_cascadia/metadata.csv'
    training = root / 'results/north_full_family_training_recovered_v1/training.json'
    out = root / 'results/north_source_identity_audit_v1'
    out.mkdir(exist_ok=True)
    ref = read_bostock(reference).groupby('family_id').agg(
        latitude=('latitude', 'first'), longitude=('longitude', 'first'),
        depth_km=('depth_km', 'first'), event_count=('reference_id', 'size')).reset_index()
    meta = pd.read_csv(metadata)
    fields = ['source_latitude_deg', 'source_longitude_deg', 'source_depth_km']
    sources = meta[['source_family_id'] + fields].drop_duplicates()
    if sources.source_family_id.duplicated().any():
        raise ValueError('inconsistent source coordinates')
    for target, field in zip(['latitude', 'longitude', 'depth_km'], fields):
        sources[target] = sources[field].round(2)
    joined = ref.merge(sources, on=['latitude', 'longitude', 'depth_km'], how='left')
    multiplicity = joined.groupby('family_id').source_family_id.transform('count')
    joined['coordinate_evidence'] = multiplicity.map(
        lambda n: 'unique_exact_rounding' if n == 1 else ('unmatched' if n == 0 else 'ambiguous'))
    selected = list(json.loads(training.read_text())['anchors_from_published_s'])
    joined['in_frozen_training_bank'] = joined.family_id.isin(selected)
    joined.to_csv(out / 'coordinate_identity_evidence.csv', index=False)
    matched = joined[joined.coordinate_evidence == 'unique_exact_rounding']
    # A printed coordinate at a rounding tie can differ across formatting rules.
    # Preserve the strict join and report interval compatibility separately.
    compatible = []
    for row in ref.itertuples():
        mask = ((sources.source_latitude_deg - row.latitude).abs() <= .00500001) & \
               ((sources.source_longitude_deg - row.longitude).abs() <= .00500001) & \
               ((sources.source_depth_km - row.depth_km).abs() <= .00500001)
        ids = sources.loc[mask, 'source_family_id'].tolist()
        compatible.append({'family_id': row.family_id, 'compatible_source_ids': ids,
                           'unique_interval_match': len(ids) == 1,
                           'in_frozen_training_bank': row.family_id in selected})
    write_json(out / 'rounding_interval_evidence.json', compatible)
    phase_rows = meta[meta.station_code.isin(['PGC', 'SNB', 'VGZ'])].merge(
        matched[['family_id', 'source_family_id', 'in_frozen_training_bank']], on='source_family_id')
    phase_rows.to_csv(out / 'overlapping_station_phase_metadata.csv', index=False)
    previous = []
    for path in sorted((root / 'results').glob('method_interleaved_published*/source_mapping.json')):
        rows = pd.DataFrame(json.loads(path.read_text()))
        compare = rows.merge(joined[['family_id','source_family_id','coordinate_evidence']], on='family_id')
        compare['same_source_id'] = compare.published_source_id == compare.source_family_id
        compare.to_csv(out / (path.parent.name + '_comparison.csv'), index=False)
        previous.append({'path': str(path.relative_to(root)), 'rows': len(compare),
                         'previously_accepted': int(compare.accepted_mapping.sum()),
                         'exact_rounding_agreements': int(compare.same_source_id.sum())})
    report = {'status': 'coordinate_identity_evidence_only', 'reference_groups': len(ref),
              'published_stack_families': len(sources),
              'unique_exact_rounding_matches': len(matched),
              'unmatched_groups': joined.loc[joined.coordinate_evidence == 'unmatched', 'family_id'].tolist(),
              'unique_rounding_interval_matches': sum(x['unique_interval_match'] for x in compatible),
              'training_interval_matches': sum(x['unique_interval_match'] and x['in_frozen_training_bank'] for x in compatible),
              'training_groups': len(selected),
              'training_unique_matches': int(matched.in_frozen_training_bank.sum()),
              'overlapping_station_rows': len(phase_rows),
              'prior_work_reused': previous,
              'production_validated': False,
              'limitations': ['Coordinates alone do not establish author-verified family identity.',
                  'No source IDs or frozen evaluation results changed.',
                  'Published stacks may share training events and are not independent test truth.',
                  'Aligned stacks do not preserve absolute interstation arrival clocks.'],
              'inputs': {str(p): sha256(p) for p in [reference, metadata, training, Path(__file__)]}}
    write_json(out / 'report.json', report)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
