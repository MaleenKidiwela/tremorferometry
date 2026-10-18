"""Training-derived clock alignment and within-published-family candidate grouping.

Cross-family duplicates remain unresolved. These are review candidates, not a
final event catalog or independently measured hypocenters.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from .provenance import write_json,sha256

def clock_offsets(mapping,comparisons):
    offsets={};diagnostics=[]
    for family,table in mapping.groupby('published_family_id'):
        anchor=sorted(table.template_id)[0];offsets[anchor]=0.
        for template in sorted(table.template_id)[1:]:
            pair=comparisons[(comparisons.template_a==anchor)&(comparisons.template_b==template)]
            network=pair[pair.station=='MEAN_NETWORK'].iloc[0]
            lag=float(network.b_index_minus_a_index_s)
            sites=pair[pair.station!='MEAN_NETWORK']
            support=sites[(sites.b_index_minus_a_index_s-lag).abs().le(.50001)&sites.peak_cc.ge(.2)]
            if len(support)<3 or network.peak_cc<.25:raise ValueError('Insufficient training-template clock support')
            offsets[template]=lag
            diagnostics.append(dict(published_family_id=family,anchor=anchor,template=template,add_to_template_time_s=lag,supporting_sites=len(support),network_cc=float(network.peak_cc)))
    return offsets,diagnostics

def group_candidates(frame,mapping,offsets,tolerance_s=1.):
    if frame.empty:return frame.copy()
    frame=frame.copy();lookup=mapping.set_index('template_id')
    frame['source_template_id']=frame.family_id
    frame['published_family_id']=frame.source_template_id.map(lookup.published_family_id)
    frame['aligned_time_s']=frame.time_s+frame.source_template_id.map(offsets)
    if frame.aligned_time_s.isna().any():raise ValueError('Unmapped template clock')
    result=[]
    for (family,threshold,control),table in frame.groupby(['published_family_id','threshold_mad','control']):
        table=table.sort_values(['score_mad','aligned_time_s','source_template_id'],ascending=[False,True,True]).reset_index(drop=True)
        remaining=np.ones(len(table),dtype=bool);times=table.aligned_time_s.to_numpy()
        for i,row in table.iterrows():
            if not remaining[i]:continue
            near=remaining&(np.abs(times-times[i])<=tolerance_s)
            record=row.to_dict();record['template_alternatives']=','.join(sorted(set(table.loc[near,'source_template_id'])))
            record['template_detection_count']=int(near.sum());record['event_status']='candidate_requires_review'
            record['time_basis']='canonical_published_template_start_'+family
            record['cross_family_duplicates_resolved']=False
            loc=lookup.loc[record['source_template_id']]
            record.update(latitude=float(loc.latitude),longitude=float(loc.longitude),depth_km=float(loc.depth_km),location_source='Published Plourde family location; not relocated individual event')
            result.append(record);remaining[near]=False
    return pd.DataFrame(result).sort_values(['control','threshold_mad','aligned_time_s']).reset_index(drop=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');args=ap.parse_args();root=Path(args.directory)
    source=Path('results/southern_benchmark_inventory_v1/comparison_bank_source_families.csv')
    clock=Path('results/ducellier_template_clock_diagnostic_v1/clock_comparison.csv')
    mapping=pd.read_csv(source);comparisons=pd.read_csv(clock);offsets,diagnostics=clock_offsets(mapping,comparisons)
    candidates=pd.read_csv(root/'candidates.csv');grouped=group_candidates(candidates,mapping,offsets)
    grouped.to_csv(root/'within_family_review_candidates.csv',index=False)
    write_json(root/'candidate_grouping_provenance.json',dict(source_sha256=sha256(__file__),input_sha256=sha256(root/'candidates.csv'),mapping_sha256=sha256(source),clock_comparison_sha256=sha256(clock),add_to_template_time_s=offsets,clock_support=diagnostics,tolerance_s=1.,physical_origins_estimated=False,cross_family_duplicates_resolved=False,final_catalog=False))
    counts=grouped.groupby(['threshold_mad','control','published_family_id']).size().rename('within_family_candidates').reset_index()
    counts.to_csv(root/'within_family_counts.csv',index=False)
    print(counts[counts.control.eq(-1)].to_string(index=False))

if __name__=='__main__':main()
