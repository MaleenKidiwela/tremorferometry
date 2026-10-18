"""Calibration-only fixed-lag packet review; no inferred event labels."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .pilot import load_waves
from .templates import cut, vector_cc
from .provenance import write_json, sha256


def main():
    run = Path('results/method_interleaved_mf_v1')
    source = Path('results/method_interleaved_miss_diagnosis_v1/matcher_calibration_reference_scores.parquet')
    out = Path('results/method_interleaved_packet_review_v1')
    cfg = json.loads((run/'manifest.json').read_text())['config']
    refs = pd.read_parquet(source)
    assert set(refs.date) <= set(cfg['calibration_dates'])
    assert refs.reference_id.is_unique
    out.mkdir(exist_ok=False)
    chosen = []
    for _, g in refs.groupby(['family_id', 'status'], sort=True):
        g = g.sort_values(['max_network_mad', 'reference_id'])
        chosen.append(g.iloc[len(g)//2].reference_id)
    write_json(out/'design_before_review.json', dict(split='calibration_only', reference_count=len(refs),
        example_ids=chosen, selection='median network score per family/status', shifts_s=[-120,0,120],
        lag_fitting=False, review_status='pending'))
    bank = np.load(run/'templates.npz')
    records, pages, inputs, arrays = [], [], {}, {}
    for date, group in refs.groupby('date', sort=True):
        print('packet review', date, flush=True)
        waves = load_waves(cfg, date, inputs)
        for ref in group.itertuples():
            panel = ref.reference_id in chosen
            if panel:
                fig, axes = plt.subplots(len(cfg['stations']), 3, figsize=(14,10), squeeze=False)
            for row, (station, meta) in enumerate(cfg['stations'].items()):
                key = f'{ref.family_id}__{station}'
                if key+'__template' not in bank:
                    if panel:
                        for ax in axes[row]:
                            ax.text(.1,.5,'No frozen template'); ax.set_axis_off()
                    continue
                template = bank[key+'__template']
                start = ref.reference_anchor_s + float(bank[key+'__offset_s'])
                packet = None
                for shift in [-120,0,120]:
                    x = cut(waves[station], start+shift, template.shape[1], alignment=cfg.get('template_alignment','linear'))
                    records.append(dict(reference_id=ref.reference_id, family_id=ref.family_id, date=date,
                        status=ref.status, station=station, site_id=meta['site_id'], shift_s=shift,
                        valid=x is not None, fixed_lag_cc=vector_cc(x,template) if x is not None else np.nan))
                    if shift == 0:
                        packet = x
                if panel:
                    arrays[f'{ref.reference_id}__{station}__template'] = template
                    if packet is not None:
                        arrays[f'{ref.reference_id}__{station}__packet'] = packet
                    for c, ax in enumerate(axes[row]):
                        for x,color,label in [(template,'#cf6b16','template'),(packet,'#2268a5','reference packet')]:
                            if x is not None:
                                y = x-x.mean(axis=1,keepdims=True)
                                y = y/(np.sqrt(np.mean(y*y)) or 1.)
                                ax.plot(np.arange(y.shape[1])/cfg['fs'],y[c],color=color,lw=.8,label=label)
                        ax.set_title(meta['channels'][c],fontsize=9); ax.grid(alpha=.2)
                        if row == len(cfg['stations'])-1:
                            ax.set_xlabel('Seconds from frozen station window start')
                    axes[row,0].legend(fontsize=7)
            if panel:
                fig.suptitle(f'{ref.reference_id} | {ref.family_id}\n{ref.status}; network maximum {ref.max_network_mad:.2f} MAD; fixed reference lag',fontsize=11)
                fig.tight_layout(rect=(0,0,1,.95))
                name = ref.reference_id+'.png'
                fig.savefig(out/name,dpi=110); plt.close(fig)
                pages.append(f'<h2>{ref.reference_id}</h2><img loading="lazy" src="{name}">')
    frame = pd.DataFrame(records)
    frame.to_parquet(out/'station_packet_scores.parquet',index=False)
    sites = frame.groupby(['reference_id','family_id','status','shift_s','site_id'],as_index=False).agg(cc=('fixed_lag_cc','mean'),valid=('valid','all'))
    sites.loc[~sites.valid,'cc'] = np.nan
    sites.to_csv(out/'site_packet_scores.csv',index=False)
    sites.groupby(['status','shift_s','site_id']).agg(n=('reference_id','size'),n_valid=('cc','count'),median_cc=('cc','median')).reset_index().to_csv(out/'summary.csv',index=False)
    np.savez_compressed(out/'review_arrays.npz',**arrays)
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Calibration packet review</title><style>body{font:16px system-ui;max-width:1450px;margin:2em auto}img{width:100%}</style><h1>Calibration waveform/template review — labels pending</h1><p>One median-score reference per source/status. Orange: frozen training template; blue: calibration waveform. Fixed reference timing plus saved station offsets; no lag optimization. Each three-component packet is demeaned and divided by its own common RMS; amplitude differences between template and event are removed. Native components, including borehole 1/2. Shifted windows can contain other events and are not noise labels. Reference agreement is not event validation.</p>'+''.join(pages))
    write_json(out/'manifest.json',dict(sources={str(p.resolve()):sha256(p) for p in [Path(__file__),run/'manifest.json',run/'templates.npz',source]},waveform_inputs=inputs,split='calibration_only',review_status='pending',reference_count=len(refs),panel_count=len(pages),station_rows=len(frame),interpretation='Fixed-lag waveform shape diagnostic, not detector scores or precision.'))
    print('Saved',len(pages),'panels and',len(frame),'station windows',flush=True)

if __name__ == '__main__':
    main()
