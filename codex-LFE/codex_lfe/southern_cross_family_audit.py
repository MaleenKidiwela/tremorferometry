"""Check published cross-family template similarity before resolving candidates.

Existing within-family clock curves are reused; no event times tune the clock or
similarity criteria, and weak template similarities never merge source families.
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from .ducellier_template_clocks import lag_cc
from .published_pickle import read_published_pickle
from .provenance import sha256,write_json


def main():
    root=Path(__file__).resolve().parents[1];out=root/'results/southern_cross_family_audit_v1';out.mkdir(exist_ok=False)
    source=root/'results/southern_benchmark_inventory_v1'
    mapping_path=source/'comparison_bank_source_families.csv';mapping=pd.read_csv(mapping_path)
    plan_path=source/'targeted_reference_packet_plan.csv';plan=pd.read_csv(plan_path).drop_duplicates('station')
    write_json(out/'design_before_audit.json',dict(source_sha256=sha256(__file__),production_validated=False,
        comparison='Distinct published families only; previous within-family clock comparisons not repeated.',
        lag_range_s=[-10,10],minimum_overlap_s=15,fs=20,
        flag='At least three physical sites with signed component-mean CC >=0.6 at the same training-derived network lag, with each site peak within 0.5 s.',
        candidate_tolerance_s=1.,action='Flag cross-family ambiguity only. Never rewrite original family identities or assert physical origin times.'))
    cached={};inputs={};rows=[];summaries=[];curves={}
    for template in mapping.template_id:
        for row in plan.itertuples():
            paths=[source/'source/src/templates_v1'/template/f'{row.station.split(".")[1]}_{c.split(".")[-1]}.pkl' for c in row.channels.split(',')]
            if not all(p.exists() for p in paths):continue
            values=[read_published_pickle(p) for p in paths]
            if any(tr.stats.sampling_rate!=20 for tr,_ in values):raise ValueError('Unexpected sample rate')
            cached[template,row.station]=(np.array([tr.data for tr,_ in values]),[angle for _,angle in values])
            inputs.update({str(p):sha256(p) for p in paths})
    lookup=mapping.set_index('template_id').published_family_id.to_dict()
    for a,b in itertools.combinations(sorted(mapping.template_id),2):
        if lookup[a]==lookup[b]:continue
        station_curves={}
        for station in sorted(plan.station):
            if (a,station) not in cached or (b,station) not in cached:continue
            (x,anglex),(y,angley)=cached[a,station],cached[b,station]
            if anglex!=angley:raise ValueError('Incompatible template coordinate conventions')
            cc=np.array([lag_cc(x,y,lag) for lag in range(-200,201)])
            station_curves[station]=cc;curves[a+'__'+b+'__'+station]=cc
        if not station_curves:continue
        network=np.nanmean(list(station_curves.values()),axis=0);best=int(np.nanargmax(network));lag=(best-200)/20
        support=[]
        for station,cc in station_curves.items():
            peak=int(np.nanargmax(cc));site_lag=(peak-200)/20
            accepted=bool(cc[best]>=.6 and abs(site_lag-lag)<=.50001)
            if accepted:support.append(station)
            rows.append(dict(template_a=a,template_b=b,station=station,site_peak_cc=cc[peak],site_peak_lag_s=site_lag,
                network_lag_s=lag,site_cc_at_network_lag=cc[best],supports_flag=accepted))
        summaries.append(dict(template_a=a,template_b=b,family_a=lookup[a],family_b=lookup[b],
            network_cc=network[best],add_to_b_detection_to_align_a_s=lag,common_sites=len(station_curves),
            supporting_sites=','.join(support),n_supporting_sites=len(support),flagged=len(support)>=3))
    comparisons=pd.DataFrame(summaries);comparisons.to_csv(out/'cross_family_template_comparison.csv',index=False)
    pd.DataFrame(rows).to_csv(out/'station_comparison.csv',index=False);np.savez_compressed(out/'training_correlation_curves.npz',**curves)
    # Only now inspect saved candidate timing, after the training-only flags exist.
    candidates_path=root/'results/ducellier_multifamily_november23_holdout_edges_v2/candidates.csv'
    candidates=pd.read_csv(candidates_path);candidates=candidates[(candidates.control==-1)&(candidates.threshold_mad==8)].copy()
    pairs=[]
    for comparison in comparisons.itertuples():
        if not comparison.flagged:continue
        aa=candidates[candidates.family_id==comparison.template_a];bb=candidates[candidates.family_id==comparison.template_b]
        for ia,a in aa.iterrows():
            for ib,b in bb.iterrows():
                delta=b.time_s+comparison.add_to_b_detection_to_align_a_s-a.time_s
                if abs(delta)<=1:
                    pairs.append(dict(source_row_a=int(ia),source_row_b=int(ib),template_a=comparison.template_a,
                        template_b=comparison.template_b,aligned_time_difference_s=delta,action='cross_family_ambiguity_flag_only'))
    pd.DataFrame(pairs,columns=['source_row_a','source_row_b','template_a','template_b','aligned_time_difference_s','action']).to_csv(out/'candidate_ambiguity_pairs.csv',index=False)
    write_json(out/'report.json',dict(status='complete_training_template_cross_family_audit',production_validated=False,
        cross_family_template_pairs=len(comparisons),flagged_template_pairs=int(comparisons.flagged.sum()),
        candidate_ambiguity_pairs=len(pairs),families_merged=0,raw_bytes_downloaded=0,
        limitation='Failure to flag an alias does not establish that all candidate rows are unique physical events.',
        sources={str(p):sha256(p) for p in [Path(__file__),root/'codex_lfe/ducellier_template_clocks.py',mapping_path,plan_path,candidates_path]},template_inputs=inputs))
    print(comparisons.to_string(index=False),flush=True)


if __name__=='__main__':main()
