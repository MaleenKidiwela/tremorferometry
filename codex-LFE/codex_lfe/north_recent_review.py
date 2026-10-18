"""Retain the bounded northern comparison and review both strict-threshold seeds."""
from pathlib import Path
import json,html
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .provenance import sha256,write_json


def main():
    root=Path('results');out=root/'north_recent_review_v1';out.mkdir(exist_ok=False)
    comparisons=[];pages=[];sources=[Path(__file__)];bytes_total=0
    for panel in ['active','quiet']:
        source=root/f'north_recent_{panel}_picker_v1';assoc=root/f'north_recent_{panel}_association_v1'
        summary=pd.read_csv(assoc/'summary.csv');sources.extend([assoc/'summary.csv',source/'report.json'])
        bytes_total+=json.loads((source/'report.json').read_text())['waveform_bytes']
        for (threshold,sites),g in summary.groupby(['threshold','minimum_S_sites']):
            observed=int(g[g.realization==0].candidates.iloc[0]);control=g[g.realization>0].candidates.to_numpy()
            comparisons.append(dict(panel=panel,threshold=threshold,minimum_S_sites=sites,observed=observed,control_min=int(control.min()),control_max=int(control.max()),control_mean=float(control.mean()),control_counts=json.dumps(control.tolist()),exposure_s=int(g.exposure_s.iloc[0])))
        events=pd.read_parquet(assoc/'t0.3_n3_r0_events.parquet');assignments=pd.read_parquet(assoc/'t0.3_n3_r0_assignments.parquet');stations=pd.read_csv(source/'stations.csv')
        sources.extend([assoc/'t0.3_n3_r0_events.parquet',assoc/'t0.3_n3_r0_assignments.parquet'])
        for event in events.itertuples():
            fig,axes=plt.subplots(len(stations),2,figsize=(16,17),sharex=True)
            for row,station in enumerate(stations.id):
                p=source/f'{station}_processed.npz';sources.append(p)
                with np.load(p) as a:
                    fs=float(a['fs']);t=float(a['start_s'])+np.arange(a['waveform'].shape[1])/fs-event.time;mask=(t>=-10)&(t<=70)
                    for component,x in enumerate(a['waveform']):
                        scale=max(float(np.percentile(np.abs(x[mask]),95)),1e-20)
                        axes[row,0].plot(t[mask],np.where(a['valid'][component,mask],x[mask]/scale+3*component,np.nan),lw=.4)
                    for phase,prob in zip(['P','S'],a['probabilities']):axes[row,1].plot(t[mask],prob[mask],label=phase,lw=.8)
                picks=assignments[(assignments.event_idx==event.idx)&(assignments.station==station)]
                for pick in picks.itertuples():
                    for ax in axes[row]:ax.axvline(pick.time-event.time,color='black',ls='--',lw=.7)
                axes[row,0].set_ylabel(station);axes[row,1].set_ylim(0,1);axes[row,1].legend(loc='upper right')
            axes[-1,0].set_xlabel('Seconds from diagnostic association time');axes[-1,1].set_xlabel('Seconds from diagnostic association time')
            fig.suptitle(f'{panel}: unvalidated candidate {event.time_utc}\nZ/N/E components scaled separately; dashed lines mark assigned model picks')
            fig.tight_layout(rect=(0,0,1,.96));name=f'{panel}_{event.idx:03d}.png';fig.savefig(out/name,dpi=130);plt.close(fig)
            pages.append(f'<h2>{html.escape(panel)}: {html.escape(str(event.time_utc))}</h2><p>Review status: pending; no LFE or phase truth assigned.</p><img src="{name}" alt="Candidate waveforms and picker probabilities">')
    frame=pd.DataFrame(comparisons);frame.to_csv(out/'comparison.csv',index=False)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Northern August 2026 diagnostic</title><style>body{font:16px system-ui;max-width:1500px;margin:2em auto}img{width:100%}td,th{padding:.4em;text-align:left}</style><h1>Northern August 2026 diagnostic</h1><p>Two one-hour panels, eight physical sites, identical picker and association settings. PNSN activity is context, not individual LFE truth. Neither hour has a four-site candidate at threshold 0.3; both have one three-site candidate. Lower-threshold activity counts remain comparable to shifted controls. No seeds promoted.</p>'+frame.to_html(index=False)+''.join(pages))
    write_json(out/'report.json',dict(status='complete_unvalidated_panel_review',raw_bytes=bytes_total,production_validated=False,phase_labels_created=0,individual_review_status='pending',sources={str(p):sha256(p) for p in sorted(set(sources))}))
    print(frame.to_string(index=False))

if __name__=='__main__':main()
