"""Consolidate bounded central diagnostics without promoting candidates to LFEs."""
from pathlib import Path
import json,html
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .provenance import write_json,sha256


def main():
    out=Path('results/central_validation_summary_v1');out.mkdir(exist_ok=True)
    records=[];inputs=[]
    for label,active,quiet in [
        ('September, central stitching','central_activity_association_v1','central_activity_matched_quiet_association_v1'),
        ('September, library averaging','central_activity_picker_panel_v1_standard_association_v1','central_activity_matched_quiet_v1_standard_association_v1'),
        ('March, central stitching','central_strong_association_v2','central_strong_quiet_association_v2')]:
        for kind,directory in [('active',active),('quiet',quiet)]:
            p=Path('results')/directory/'summary.csv';inputs.append(p);frame=pd.read_csv(p)
            for (threshold,sites),group in frame.groupby(['threshold','minimum_S_sites']):
                observed=int(group.loc[group.realization==0,'candidates'].iloc[0]);null=group.loc[group.realization>0,'candidates'].to_numpy()
                records.append(dict(experiment=label,panel=kind,method='neural_association',threshold=threshold,minimum_sites=int(sites),observed=observed,control_mean=float(null.mean()),control_min=int(null.min()),control_max=int(null.max()),control_counts=json.dumps(null.tolist()),n_controls=len(null),source=str(p)))
    for label,directory in [('September, horizontal coherence','central_cross_station_v1'),('March, horizontal coherence','central_strong_cross_station_v2')]:
        p=Path('results')/directory/'summary.csv';inputs.append(p);frame=pd.read_csv(p)
        for (panel,threshold),group in frame.groupby(['panel','threshold']):
            observed=int(group.loc[group.realization==0,'candidates'].iloc[0]);null=group.loc[group.realization>0,'candidates'].to_numpy()
            records.append(dict(experiment=label,panel='quiet' if 'quiet' in panel else 'active',method='horizontal_coherence',threshold=threshold,minimum_sites=3,observed=observed,control_mean=float(null.mean()),control_min=int(null.min()),control_max=int(null.max()),control_counts=json.dumps(null.tolist()),n_controls=len(null),source=str(p)))
    frame=pd.DataFrame(records);frame.to_csv(out/'comparison.csv',index=False)
    shown=frame[((frame.method=='neural_association')&(frame.threshold==.1)&(frame.minimum_sites==3))|((frame.method=='horizontal_coherence')&(frame.threshold==.5))].reset_index(drop=True)
    fig,ax=plt.subplots(figsize=(11,7))
    for i,r in enumerate(shown.itertuples()):
        null=json.loads(r.control_counts);ax.scatter(null,np.full(len(null),i),marker='|',s=180,color='#377eb8',alpha=.65)
        ax.scatter(r.observed,i,color='black',s=38,zorder=3)
    ax.set_yticks(range(len(shown)),[f'{r.experiment}: {r.panel}' for r in shown.itertuples()]);ax.invert_yaxis();ax.set_xlabel('Candidate seed count in the bounded hour (not confirmed LFEs)');ax.grid(axis='x',alpha=.2)
    ax.set_title('Black: observed; blue: eight independently shifted controls\nMethods have different search statistics; compare within each row')
    fig.tight_layout();fig.savefig(out/'controls.png',dpi=160);plt.close(fig)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Central discovery validation</title><style>body{font:16px system-ui;max-width:1200px;margin:2em auto}img{width:100%}td,th{padding:.4em}table{border-collapse:collapse}</style><h1>Central Cascadia discovery diagnostics</h1><p>These are bounded development diagnostics, not a catalog release. Tremor context is not individual-event or phase truth. Proximity groups within 1 km count as one site. The March v1 cross-station attempt was aborted because it counted two nearby boreholes separately; only corrected v2 results appear here.</p><img src="controls.png" alt="Observed candidate counts and shifted controls"><p>Thresholds belong to different statistics and are not interchangeable. Control counts are not measured false discoveries. A one-hour result cannot validate full-history sensitivity.</p>'+frame.drop(columns=['control_counts','source']).to_html(index=False)+'<p><a href="comparison.csv">All counts and input paths</a></p>')
    write_json(out/'manifest.json',dict(status='development_diagnostic_summary',production_validated=False,inputs={str(p):sha256(p) for p in inputs},source_sha256=sha256(__file__),confirmed_LFE_labels_created=0))
    print(frame.to_string(index=False),flush=True)
if __name__=='__main__':main()
