"""Training-only family distinguishability audit; does not merge or retune families."""
import json
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from .evaluate import match_events
from .templates import vector_cc
from .reference import read_bostock
from .provenance import write_json, sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    baseline = ROOT/'results/north_full_family_bank_2012_v1'
    training = ROOT/'results/north_full_family_training_recovered_v1'
    out = ROOT/'results/north_family_distinguishability_v1'
    out.mkdir(exist_ok=False)
    cfg = json.loads((baseline/'manifest.json').read_text())['config']
    info = json.loads((training/'training.json').read_text())
    anchors = info['anchors_from_published_s']
    reference_path = ROOT/cfg['reference_csv']
    ref = read_bostock(reference_path)
    train = ref[ref.date.isin(cfg['train_dates']) & ref.family_id.isin(anchors)].copy()
    train['time_s'] += train.family_id.map(anchors)
    train['time_basis'] = 'arrival_anchor:CN.PGC:training_stack_peak'
    groups = {f: g.reset_index(drop=True) for f,g in train.groupby('family_id')}
    templates = np.load(training/'templates.npz', allow_pickle=False)
    rows, sites = [], []
    # Fixed 0.25 s tolerance assesses packet coincidence at the frozen training
    # anchors; no pair-specific time-shift optimization or evaluation data.
    for a,b in combinations(sorted(anchors),2):
        matched = match_events(groups[a], groups[b], .25, same_family=False)
        pair_sites = []
        for station, meta in cfg['stations'].items():
            ak,bk = f'{a}__{station}',f'{b}__{station}'
            if ak+'__template' not in templates or bk+'__template' not in templates:
                continue
            x,y = templates[ak+'__template'],templates[bk+'__template']
            ao,bo = float(templates[ak+'__offset_s']),float(templates[bk+'__offset_s'])
            left,right = max(ao,bo),min(ao+x.shape[1]/cfg['fs'],bo+y.shape[1]/cfg['fs'])
            n = int(round((right-left)*cfg['fs']))
            if n < 4*cfg['fs']:
                continue
            ai,bi = int(round((left-ao)*cfg['fs'])),int(round((left-bo)*cfg['fs']))
            cc = vector_cc(x[:,ai:ai+n],y[:,bi:bi+n])
            item = dict(family_a=a,family_b=b,station=station,site=meta['site_id'],
                        overlap_s=n/cfg['fs'],fixed_anchor_cc=cc)
            sites.append(item);pair_sites.append(item)
        site_cc = pd.DataFrame(pair_sites).groupby('site').fixed_anchor_cc.mean() if pair_sites else pd.Series(dtype=float)
        rows.append(dict(family_a=a,family_b=b,n_a=len(groups[a]),n_b=len(groups[b]),
                         coincident_training_events=len(matched),
                         coincident_fraction_a=len(matched)/len(groups[a]),
                         coincident_fraction_b=len(matched)/len(groups[b]),
                         anchor_difference_a_minus_b_s=anchors[a]-anchors[b],
                         shared_sites=len(site_cc),
                         site_mean_cc=float(site_cc.mean()) if len(site_cc) else None,
                         site_min_cc=float(site_cc.min()) if len(site_cc) else None,
                         sites_cc_ge_08=int((site_cc>=.8).sum())))
    pairs = pd.DataFrame(rows).sort_values(['coincident_training_events','site_mean_cc'],ascending=False)
    pairs.to_csv(out/'training_pair_evidence.csv',index=False)
    pd.DataFrame(sites).to_csv(out/'training_pair_station_evidence.csv',index=False)
    # Descriptive flag only, declared here before consulting test assignments.
    flags = pairs[(pairs.coincident_fraction_a>=.8)&(pairs.coincident_fraction_b>=.8)&(pairs.sites_cc_ge_08>=3)]
    flags.to_csv(out/'training_indistinguishable_pair_flags.csv',index=False)
    evaluation = ROOT/'results/north_full_family_frozen_evaluation_v1'
    candidates = pd.read_parquet(evaluation/'test_experimental_candidates.parquet')
    flagged_sets = [set([r.family_a,r.family_b]) for r in flags.itertuples()]
    candidates['contains_training_flagged_pair'] = candidates.family_alternatives.map(
        lambda s: any(f.issubset(set(str(s).split(','))) for f in flagged_sets))
    candidates.to_parquet(out/'evaluation_assignment_diagnostic.parquet',index=False)
    summary = dict(status='complete_training_distinguishability_diagnostic',
                   pair_count=len(pairs),training_flagged_pairs=len(flags),
                   evaluation_candidates=len(candidates),
                   evaluation_ambiguous=int(candidates.family_id.isna().sum()),
                   evaluation_ambiguous_with_training_flag=int((candidates.family_id.isna()&candidates.contains_training_flagged_pair).sum()),
                   family_merges_applied=0,thresholds_changed=False,production_validated=False,
                   limitations=['Coincident source records and similar templates do not prove identical sources.',
                                'Training anchors are packet clocks, not validated origin times.',
                                'Frozen evaluation failure is unchanged; these dates are not a new independent test.'],
                   inputs={str(p):sha256(p) for p in [reference_path,training/'training.json',training/'templates.npz',baseline/'manifest.json',evaluation/'test_experimental_candidates.parquet',Path(__file__)]})
    write_json(out/'report.json',summary)
    print(flags.to_string(index=False));print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
