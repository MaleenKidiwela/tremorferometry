"""Source-packet coincidence with training-derived peaks, NOT phase accuracy."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .evaluate import match_events
from .provenance import sha256,write_json


def main():
    root=Path('results/southern_picker_transfer_v1')
    out=root/'reference_coincidence';out.mkdir(exist_ok=True)
    design=json.loads((root/'design_before_inference.json').read_text())
    training_path=Path('results/southern_training_complete_v1/report.json')
    training=json.loads(training_path.read_text())
    rp=Path('results/southern_reference_audit_v1/published_reference.parquet')
    ref=pd.read_parquet(rp)
    ref=ref[ref.date.isin(design['dates'])&ref.family_id.isin(training['four_site_families'])].copy()
    offsets={(r['station'],f):q['peak_offset_from_published_s'] for r in training['results'] for f,q in r['qc_by_family'].items() if q.get('usable') and f in training['four_site_families']}
    shifts=[0,-600,-300,300,600]
    write_json(out/'design.json',dict(thresholds=design['threshold_grid'],tolerances_s=[2,4],shifts_s=shifts,
        references=len(ref),training_sha256=sha256(training_path),reference_sha256=sha256(rp),
        interpretation='S-pick coincidence around January-derived station stack peaks; these are not independently labeled S onsets. Shifted-reference comparison is a diagnostic, not event FDR. No time offset fitted on these February observations.'))
    allp=[];allc=[]
    for date in design['dates']:
        for station in design['stations']:
            folder=root/f'{date}_{station}'
            allp.append(pd.read_parquet(folder/'picks.parquet'))
            allc.append(pd.read_parquet(folder/'coverage.parquet'))
    picks=pd.concat(allp);coverage=pd.concat(allc)
    scores=[];matches=[];eligible_rows=[]
    for station in design['stations']:
        r=ref[ref.family_id.map(lambda f:(station,f) in offsets)].copy()
        r['time_s']=r.time_s+r.family_id.map(lambda f:offsets[(station,f)])
        cv=coverage[(coverage.station==station)&(coverage.phase=='S')]
        for shift in shifts:
            rr=r.copy();rr['time_s']+=shift
            for tol in [2,4]:
                eligible=np.zeros(len(rr),bool)
                for a in cv.itertuples():eligible|=(rr.time_s.to_numpy()-tol>=a.start_s)&(rr.time_s.to_numpy()+tol<a.end_s)
                selected=rr[eligible].copy().reset_index(drop=True)
                selected['family_id']='packet';selected['time_basis']='training_stack_peak_proxy'
                eligible_rows.extend(dict(station=station,shift_s=shift,tolerance_s=tol,reference_id=a.reference_id) for a in selected.itertuples())
                for threshold in design['threshold_grid']:
                    p=picks[(picks.station==station)&(picks.phase=='S')&(picks.model_score>=threshold)].copy()
                    # Only nearby predictions can match; retain all counts separately.
                    total=len(p)
                    times=np.sort(selected.time_s.to_numpy())
                    left=np.searchsorted(times,p.time_s.to_numpy()-tol,side='left');right=np.searchsorted(times,p.time_s.to_numpy()+tol,side='right')
                    p=p[right>left].reset_index(drop=True)
                    p['family_id']='packet';p['time_basis']='training_stack_peak_proxy'
                    m=match_events(p,selected,tolerance_s=tol,same_family=True)
                    scores.append(dict(station=station,threshold=threshold,tolerance_s=tol,shift_s=shift,eligible_reference_packets=len(selected),matched_packets=len(m),S_picks_in_panel=total,packet_coincidence_fraction=len(m)/len(selected) if len(selected) else None))
                    for a in m.itertuples():matches.append(dict(station=station,threshold=threshold,tolerance_s=tol,shift_s=shift,reference_id=selected.iloc[a.reference_index].reference_id,residual_from_stack_peak_s=a.residual_s))
    pd.DataFrame(scores).to_csv(out/'station_scores.csv',index=False)
    pd.DataFrame(matches).to_csv(out/'matched_reference_packets.csv',index=False)
    pd.DataFrame(eligible_rows).to_csv(out/'eligible_reference_packets.csv',index=False)
    write_json(out/'report.json',dict(status='complete_reference_packet_coincidence_only',reference_events=len(ref),phase_accuracy_validated=False,production_picker_validated=False))
    print(pd.DataFrame(scores).query('threshold == 0.3 and tolerance_s == 2').to_string(index=False))

if __name__=='__main__':main()
