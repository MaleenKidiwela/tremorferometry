"""Summarize predeclared central controls and render the strict-threshold seed."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .provenance import write_json,sha256

def main():
    out=Path('results/central_activity_comparison_v1');out.mkdir(exist_ok=True)
    rows=[]
    for panel,directory in [('active','central_activity_association_v1'),('quiet','central_activity_matched_quiet_association_v1')]:
        source=Path('results')/directory
        frame=pd.read_csv(source/'summary.csv')
        for (threshold,sites),group in frame.groupby(['threshold','minimum_S_sites']):
            observed=int(group.loc[group.realization==0,'candidates'].iloc[0]);null=group.loc[group.realization>0,'candidates'].to_numpy()
            rows.append(dict(panel=panel,threshold=threshold,minimum_S_sites=sites,observed=observed,control_mean=float(null.mean()),control_min=int(null.min()),control_max=int(null.max()),control_counts=json.dumps(null.tolist()),rank_tail_fraction=float((1+(null>=observed).sum())/(1+len(null))),exposure_s=3120))
    pd.DataFrame(rows).to_csv(out/'comparison.csv',index=False)
    source=Path('results/central_activity_picker_panel_v1');assoc=Path('results/central_activity_association_v1')
    events=pd.read_parquet(assoc/'t0.3_n3_r0_events.parquet');assignments=pd.read_parquet(assoc/'t0.3_n3_r0_assignments.parquet')
    stations=pd.read_csv(source/'stations.csv')
    for event in events.itertuples():
        fig,axes=plt.subplots(len(stations),2,figsize=(16,17),sharex=True)
        for row,station in enumerate(stations.id):
            a=np.load(source/(station+'_processed.npz'));fs=float(a['fs']);start=float(a['start_s']);t=start+np.arange(a['waveform'].shape[1])/fs-event.time
            mask=(t>=-10)&(t<=70)
            for component,x in enumerate(a['waveform']):
                scale=np.percentile(np.abs(x[mask]),95)
                axes[row,0].plot(t[mask],x[mask]/max(scale,1e-20)+3*component,lw=.4)
            for phase,prob in zip(['P','S'],a['probabilities']):axes[row,1].plot(t[mask],prob[mask],label=phase,lw=.8)
            picks=assignments[(assignments.event_idx==event.idx)&(assignments.station==station)]
            for pick in picks.itertuples():
                for ax in axes[row]:ax.axvline(pick.time-event.time,color='black',ls='--',lw=.7)
            axes[row,0].set_ylabel(station);axes[row,1].set_ylim(0,1);axes[row,1].legend(loc='upper right')
        axes[-1,0].set_xlabel('Seconds from diagnostic association time');axes[-1,1].set_xlabel('Seconds from diagnostic association time')
        fig.suptitle(f'Unvalidated central seed {event.time_utc}\nComponent amplitudes scaled separately; dashed lines are assigned model picks')
        fig.tight_layout(rect=(0,0,1,.96));fig.savefig(out/f'seed_{event.idx:03d}.png',dpi=140);plt.close(fig)
    write_json(out/'report.json',dict(status='complete_failed_discrimination_diagnostic',production_validated=False,conclusion='Neither central panel exceeds its station-shift controls. Do not promote these seeds to a template bank.',phase_truth_available=False,source_hashes={str(p):sha256(p) for p in [Path(__file__),Path('results/central_activity_association_v1/summary.csv'),Path('results/central_activity_matched_quiet_association_v1/summary.csv')]}))
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
