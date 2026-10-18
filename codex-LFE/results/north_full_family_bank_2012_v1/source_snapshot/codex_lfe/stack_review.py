"""Save separate-date training stacks for phase and timing review; no labels inferred."""
import argparse
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yaml

from .pilot import load_waves
from .provenance import sha256, write_json
from .reference import read_bostock
from .templates import cut, stack_windows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',required=True)
    ap.add_argument('--matcher',required=True)
    ap.add_argument('--output',required=True)
    args = ap.parse_args()
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError('output must be empty')
    out.mkdir(parents=True,exist_ok=True)
    cfg = yaml.safe_load(Path(args.config).read_text())
    training = json.loads((Path(args.matcher)/'training.json').read_text())
    families = training['families_usable']
    dates = sorted(cfg['train_dates'])
    halves = {d:'early' if i<len(dates)//2 else 'late' for i,d in enumerate(dates)}
    windows = {(f,s,h):[] for f in families for s in cfg['stations'] for h in ['early','late']}
    refs = read_bostock(cfg['reference_csv'])
    manifest,used = {},[]
    for date in dates:
        print('training-stack raw date',date,flush=True)
        waves = load_waves(cfg,date,manifest)
        for r in refs[(refs.date==date)&refs.family_id.isin(families)].itertuples():
            for s,w in waves.items():
                snippet = cut(w,r.time_s-10,int(70*cfg['fs']),alignment=cfg.get('template_alignment','linear'))
                if snippet is not None:
                    windows[r.family_id,s,halves[date]].append(snippet)
                    used.append({'reference_id':r.reference_id,'station':s,'half':halves[date]})
    arrays,records,pages = {},[],[]
    for family in families:
        fig,axes = plt.subplots(len(cfg['stations']),3,figsize=(15,10),sharex=True)
        for row,(station,info) in enumerate(cfg['stations'].items()):
            stacks = {}
            for half in ['early','late']:
                w = windows[family,station,half]
                records.append({'family_id':family,'station':station,'half':half,'n_windows':len(w)})
                if len(w)>=2:
                    stacks[half] = stack_windows(w)
                    arrays[f'{family}__{station}__{half}'] = stacks[half]
            for channel,seed in enumerate(info['channels']):
                ax = axes[row,channel]
                scale = max((np.max(np.abs(x[channel])) for x in stacks.values()),default=1.) or 1.
                for half,color in [('early','#2469b3'),('late','#d47818')]:
                    if half in stacks:
                        x = stacks[half][channel]
                        ax.plot(np.arange(len(x))/cfg['fs']-10,x/scale,lw=.7,color=color,
                                label=f'{half}: n={len(windows[family,station,half])}')
                ax.set_title(seed,fontsize=9)
                ax.set_xlim(-5,50);ax.set_ylim(-1.1,1.1);ax.grid(alpha=.2)
                if channel==0: ax.legend(fontsize=7,loc='upper right')
                if row==len(cfg['stations'])-1: ax.set_xlabel('Seconds relative to published reference time')
        fig.suptitle(f'{family}: training stacks on separate dates; P/S labels pending',fontsize=12)
        fig.tight_layout(rect=(0,0,1,.96))
        name=f'{family}.png';fig.savefig(out/name,dpi=130);plt.close(fig)
        pages.append(f'<h2>{html.escape(family)}</h2><img src="{name}" alt="Training stacks for {html.escape(family)}">')
    np.savez_compressed(out/'stacks.npz',**arrays)
    write_json(out/'counts.json',records)
    write_json(out/'manifest.json',{'sources':{str(Path(p).resolve()):sha256(p) for p in
        [args.config,Path(args.matcher)/'training.json',cfg['reference_csv'],__file__]},
        'waveform_inputs':manifest,'train_dates':dates,'date_halves':halves,
        'used_reference_station_records':used,'fs':cfg['fs'],'start_relative_s':-10.,
        'coordinates':'native channels','normalization':'common event RMS before stacking; plots scaled per channel only',
        'review_status':'pending','phase_labels_created':0,'absolute_origin_time_validated':False})
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Training stack review</title>'
        '<style>body{font:16px system-ui;max-width:1500px;margin:2em auto}img{width:100%}</style>'
        '<h1>Training stack phase review</h1><p>Blue and orange use separate training dates. '
        'Native channels, 1–8 Hz. Each plot channel has its own display scale; use the stored arrays '
        'for amplitude comparisons. Times are relative to the unresolved published clock. '
        'No P/S labels or origin corrections have been inferred.</p>'+''.join(pages))
    print('Saved',len(families),'multistation stack panels',flush=True)


if __name__=='__main__':
    main()
