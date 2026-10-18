"""Training-template-only clock diagnostics within published family 026."""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from .published_pickle import read_published_pickle
from .provenance import sha256,write_json

def lag_cc(a,b,lag):
    # Compare A[i] with B[i+lag]. A matching B detection starts lag/fs earlier.
    lo=max(0,-lag);hi=min(a.shape[1],b.shape[1]-lag)
    if hi-lo<300:return np.nan
    x=a[:,lo:hi];y=b[:,lo+lag:hi+lag]
    x=x-x.mean(axis=1,keepdims=True);y=y-y.mean(axis=1,keepdims=True)
    return float(np.mean(np.sum(x*y,axis=1)/(np.linalg.norm(x,axis=1)*np.linalg.norm(y,axis=1))))

def main():
    source=Path('results/southern_benchmark_inventory_v1');out=Path('results/ducellier_template_clock_diagnostic_v1');out.mkdir(exist_ok=True)
    groups=pd.read_csv(source/'comparison_bank_source_families.csv')
    plan=pd.read_csv(source/'targeted_reference_packet_plan.csv').drop_duplicates('station')
    rows=[];curves={};inputs={}
    for group,table in groups.groupby('published_family_id'):
        for a,b in itertools.combinations(sorted(table.template_id),2):
            pair=[]
            for row in plan.itertuples():
                values=[];orientations=[]
                for family in [a,b]:
                    paths=[source/'source/src/templates_v1'/family/f'{row.station.split(".")[1]}_{c.split(".")[-1]}.pkl' for c in row.channels.split(',')]
                    if not all(p.exists() for p in paths):break
                    traces=[read_published_pickle(p) for p in paths]
                    values.append(np.asarray([tr.data for tr,angle in traces]));orientations.append([angle for tr,angle in traces])
                    inputs.update({str(p):sha256(p) for p in paths})
                if len(values)!=2:continue
                if orientations[0]!=orientations[1]:raise ValueError('Cross-template orientation differs; cannot mix normalized channels')
                cc=np.asarray([lag_cc(*values,lag) for lag in range(-200,201)])
                best=int(np.nanargmax(cc));lag=best-200
                rows.append(dict(published_family_id=group,template_a=a,template_b=b,station=row.station,b_index_minus_a_index_s=lag/20,b_detection_minus_a_detection_s=-lag/20,peak_cc=float(cc[best])))
                curves[a+'__'+b+'__'+row.station]=cc;pair.append(cc)
            network=np.nanmean(pair,axis=0);best=int(np.nanargmax(network));lag=best-200
            rows.append(dict(published_family_id=group,template_a=a,template_b=b,station='MEAN_NETWORK',b_index_minus_a_index_s=lag/20,b_detection_minus_a_detection_s=-lag/20,peak_cc=float(network[best])))
    pd.DataFrame(rows).to_csv(out/'clock_comparison.csv',index=False);np.savez_compressed(out/'correlation_curves.npz',**curves)
    write_json(out/'provenance.json',dict(source_sha256=sha256(__file__),template_sha256=inputs,minimum_overlap_s=15,lag_range_s=[-10,10],reference_event_times_used=False,clock_offsets_applied=False,warning='Training-template diagnostic only. Conflicting station lags or weak template similarity prevent automatic clock unification.'))
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
