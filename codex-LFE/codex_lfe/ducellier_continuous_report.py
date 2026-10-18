"""Reproducible one-to-one published-reference comparison of saved hourly scans."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from .published_pickle import read_published_pickle
from .matched_filter import shifted_sites
from .evaluate import match_events
from .provenance import write_json,sha256

def main():
    ap=argparse.ArgumentParser();ap.add_argument('directory');args=ap.parse_args()
    root=Path(args.directory)
    design=json.loads((root/'design.json').read_text());report=json.loads((root/'report.json').read_text())
    start=pd.Timestamp(design['start']).timestamp()
    with np.load(root/'station_correlations.npz') as arrays:sites={k:arrays[k] for k in arrays.files}
    good=np.isfinite(np.asarray(list(sites.values()))).all(axis=0)
    for shifts in report['control_shifts_samples']:
        good &= np.isfinite(np.asarray(list(shifted_sites(sites,shifts).values()))).all(axis=0)
    source=Path('results/southern_benchmark_inventory_v1/source/data/Ducellier/catalogs/080424.05.060/catalog_2004_2011.pkl')
    ref=read_published_pickle(source).copy()
    ref['time_s']=pd.to_datetime(ref[['year','month','day','hour','minute']],utc=True).map(lambda t: t.timestamp())+ref.second
    ref=ref[(ref.time_s>=start)&(ref.time_s<start+design['duration_s'])].copy()
    indices=np.rint((ref.time_s.to_numpy()-start)*20).astype(int)
    ref['eligible_common_exposure']=[bool(good[i]) if 0<=i<len(good) else False for i in indices]
    ref['family_id']='080424.05.060';ref['time_basis']='published_template_start'
    ref.to_csv(root/'reference_exposure.csv',index=False)
    eligible=ref.loc[ref.eligible_common_exposure.astype(bool)].reset_index(drop=True)
    detections=pd.read_csv(root/'candidates.csv')
    detections['family_id']='080424.05.060';detections['time_basis']='published_template_start'
    rows=[]
    for threshold in design['threshold_mad']:
        d=detections[(detections.control==-1)&(detections.threshold_mad==threshold)].reset_index(drop=True)
        for tol in [.5,1.,2.]:
            matches=match_events(d,eligible,tol)
            matches.to_csv(root/f'matches_{threshold}_{tol:g}s.csv',index=False)
            rows.append(dict(threshold_mad=threshold,tolerance_s=tol,candidates=len(d),matched=len(matches),eligible_references=len(eligible),all_hour_references=len(ref)))
    pd.DataFrame(rows).to_csv(root/'reference_comparison.csv',index=False)
    write_json(root/'comparison_provenance.json',dict(source_sha256=sha256(source),scan_report_sha256=sha256(root/'report.json'),script_sha256=sha256(__file__),reference_is_exhaustive=False,phase_labels=False,warning='Published-reference agreement does not measure precision. Development selection and incomplete array prevent unbiased recall claims.'))
    print(pd.DataFrame(rows).to_string(index=False))

if __name__=='__main__':main()
