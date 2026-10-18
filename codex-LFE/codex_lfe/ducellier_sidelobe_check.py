"""Flag possible template-correlation echoes without deleting close real events."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from .published_pickle import read_published_pickle
from .matched_filter import normalized_cc
from .provenance import sha256,write_json

def main():
    source=Path('results/southern_benchmark_inventory_v1');out=Path('results/ducellier_sidelobe_diagnostic_v1');out.mkdir(exist_ok=True)
    bank=json.loads((source/'comparison_bank_design.json').read_text());plan=pd.read_csv(source/'targeted_reference_packet_plan.csv').drop_duplicates('station')
    peaks=[];arrays={}
    for family in bank['families']:
        curves=[]
        for row in plan.itertuples():
            if row.station=='PB.B039':continue # mimic fixed next-day detection panel
            paths=[source/'source/src/templates_v1'/family/f'{row.station.split(".")[1]}_{c.split(".")[-1]}.pkl' for c in row.channels.split(',')]
            if not all(p.exists() for p in paths):continue
            template=np.asarray([read_published_pickle(p)[0].data for p in paths]);n=template.shape[1]
            cc,_=normalized_cc(np.pad(template,((0,0),(n,n))),template)
            curves.append(cc[n-200:n+201])
        mean=np.mean(curves,axis=0);arrays[family]=mean
        idx,_=find_peaks(mean,height=.05)
        for i in idx:
            lag=(i-200)/20
            if 2<=abs(lag)<=10:peaks.append(dict(template_id=family,lag_s=lag,relative_self_cc=float(mean[i])))
    peak=pd.DataFrame(peaks);peak.to_csv(out/'template_sidelobes.csv',index=False);np.savez_compressed(out/'self_correlation_profiles.npz',**arrays)
    root=Path('results/ducellier_multifamily_november23_holdout_edges_v2')
    candidates=pd.read_csv(root/'candidates.csv');candidates=candidates[candidates.control==-1]
    flags=[]
    for family,table in candidates.groupby('family_id'):
        table=table.sort_values('time_s').reset_index(drop=True)
        for i,a in table.iterrows():
            for j in range(i+1,len(table)):
                b=table.iloc[j];gap=b.time_s-a.time_s
                if gap>10:break
                lobes=peak[(peak.template_id==family)&((peak.lag_s-gap).abs()<=.15)]
                if len(lobes):flags.append(dict(template_id=family,earlier_time_s=float(a.time_s),later_time_s=float(b.time_s),gap_s=float(gap),earlier_score_mad=float(a.score_mad),later_score_mad=float(b.score_mad),nearest_self_cc=float(lobes.relative_self_cc.max()),status='possible_echo_or_real_close_pair_requires_review'))
    pd.DataFrame(flags).to_csv(out/'candidate_pair_flags.csv',index=False)
    write_json(out/'provenance.json',dict(source_sha256=sha256(__file__),candidate_sha256=sha256(root/'candidates.csv'),station_withheld='PB.B039',lag_tolerance_s=.15,minimum_self_cc=.05,events_deleted=0,warning='Noiseless template autocorrelation is an algorithmic echo diagnostic, not a physical sensitivity injection or proof a nearby event is false.'))
    print(pd.DataFrame(flags).to_string(index=False))

if __name__=='__main__':main()
