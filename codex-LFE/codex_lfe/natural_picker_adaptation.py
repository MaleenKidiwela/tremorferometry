"""Bounded date-separated real-waveform adaptation using incomplete phase proxies.

This is an experiment, not a reviewed phase-label dataset or production model.
Unknown published sources mask their possible arrival intervals out of the loss.
"""
from pathlib import Path
import argparse, json, os
import numpy as np
import pandas as pd
import yaml
from scipy.signal import find_peaks
from .pilot import load_waves
from .reference import read_bostock
from .published_stack_clock import response_correct
from .picker_adaptation import cardinality
from .provenance import sha256, write_json

ROOT = Path(__file__).resolve().parents[1]
FS = 20
N = 1200


def targets_and_mask(start_s, known, unknown_times):
    """Supervise the central 30 s, conservatively excluding unknown phases."""
    t = start_s + np.arange(N) / FS
    target = np.zeros((2, N), np.float32)
    mask = np.zeros((2, N), bool)
    mask[:, 300:900] = True
    for p in range(2):
        for origin in unknown_times[p]:
            mask[p, (t >= origin) & (t <= origin + 60)] = False
        for arrival in known[p]:
            delta = t - arrival
            target[p] = np.maximum(target[p], np.exp(-.5 * (delta / .25) ** 2))
            mask[p, np.abs(delta) <= 1.5] = True
    mask[:, :300] = False
    mask[:, 900:] = False
    return target, mask


def prepare(out, cfg, transfer, references):
    from obspy import read_inventory
    inv = read_inventory(ROOT / 'cache/metadata/cn_pilot_2012.xml')
    inputs = {}; records = []; samples = {k: [] for k in ['train', 'validation', 'test']}
    for split, key in [('train', 'train_dates'), ('validation', 'calibration_dates'), ('test', 'test_dates')]:
        for date in cfg[key]:
            refs = references[references.date == date]
            day = refs.merge(transfer[['family_id', 'station', 'phase', 'phase_arrival_after_published_s']], on='family_id')
            day['arrival_s'] = day.time_s + day.phase_arrival_after_published_s
            waves = load_waves(cfg, date, inputs)
            for station, wave in waves.items():
                d = day[day.station == station]
                arrivals = [d[d.phase == p].arrival_s.to_numpy() for p in ['P', 'S']]
                known_families = [set(transfer[(transfer.station == station) & (transfer.phase == p)].family_id) for p in ['P', 'S']]
                unknown = [refs[~refs.family_id.isin(families)].time_s.to_numpy() for families in known_families]
                starts = set(int(np.floor(((t - wave.start_s) * FS - 300) / 600)) * 600 for times in arrivals for t in times)
                positive_starts = starts.copy()
                # Reference-free contexts are background comparisons, not true negatives.
                possible = np.arange(600, wave.x.shape[1] - 1500, 600)
                rng = np.random.default_rng(int(date.replace('-', '')))
                rng.shuffle(possible)
                backgrounds = set()
                for first in possible:
                    lo = wave.start_s + (first - 300) / FS
                    hi = wave.start_s + (first + 1500) / FS
                    if ((refs.time_s >= lo - 60) & (refs.time_s <= hi)).any(): continue
                    if not wave.valid[:, first-300:first+1500].all(): continue
                    backgrounds.add(int(first))
                    if len(backgrounds) == 128: break
                starts |= backgrounds
                for first in sorted(starts):
                    lo, hi = first - 300, first + 1500
                    if lo < 0 or hi > wave.x.shape[1] or not wave.valid[:, lo:hi].all(): continue
                    start_s = wave.start_s + first / FS
                    near_known = [a[(a >= start_s - 2) & (a < start_s + 62)] for a in arrivals]
                    near_unknown = [a[(a >= start_s - 60) & (a < start_s + 60)] for a in unknown]
                    y, mask = targets_and_mask(start_s, near_known, near_unknown)
                    if not mask.any(): continue
                    x = response_correct(wave.x[:, lo:hi], cfg['stations'][station]['channels'], inv, [date], FS)[:, 300:1500].astype(np.float32)
                    x -= x.mean(axis=1, keepdims=True)
                    amplitude = np.max(np.abs(x))
                    if not np.isfinite(x).all() or amplitude <= 0: raise ValueError('Invalid input waveform')
                    x /= amplitude
                    truth = [[float((t-start_s)*FS) for t in a if 15 <= t-start_s < 45] for a in near_known]
                    samples[split].append((x, y, mask, truth))
                    records.append(dict(split=split, sample_index=len(samples[split])-1, date=date, station=station,
                        context_start_s=wave.start_s+lo/FS, context_end_s=wave.start_s+hi/FS,
                        input_start_s=start_s, input_end_s=start_s+60,
                        reference_selected=first in positive_starts, reference_free_context=first in backgrounds,
                        known_P=len(truth[0]), known_S=len(truth[1]), supervised_samples=int(mask.sum())))
                print('prepared', split, date, station, 'total', len(samples[split]), flush=True)
    pd.DataFrame(records).to_csv(out/'samples.csv', index=False)
    write_json(out/'waveform_inputs.json', inputs)
    result = {}
    for split, values in samples.items():
        if not values: raise ValueError('Empty split')
        x, y, mask, truth = zip(*values)
        np.savez_compressed(out/f'{split}_dataset.npz', waveform=np.asarray(x), target=np.asarray(y), mask=np.asarray(mask))
        write_json(out/f'{split}_proxy_samples.json', truth)
        result[split] = (np.asarray(x), np.asarray(y), np.asarray(mask), truth)
    return result


def evaluate(model, data, device):
    import torch
    x, _, mask, truth = data
    model.eval(); predictions = []
    with torch.inference_mode():
        for first in range(0, len(x), 128):
            predictions.append(model(torch.from_numpy(x[first:first+128]).to(device)).cpu().numpy())
    predictions = np.concatenate(predictions)
    rows = []
    for p, phase in enumerate(['P', 'S']):
        matched = proxies = unmatched = exposure = 0
        for k in range(len(x)):
            # Do not score peaks created at a masked boundary or reuse one pick
            # for several dense reference arrivals.
            support = np.convolve(mask[k,p].astype(int), np.ones(41, int), mode='same') == 41
            peaks = find_peaks(predictions[k,p], height=.3, distance=20)[0]
            peaks = peaks[support[peaks]]
            times = [t for t in truth[k][p] if support[int(round(t))]]
            hits = cardinality(peaks, times, tolerance=10)
            matched += hits; proxies += len(times); unmatched += len(peaks)-hits; exposure += support.sum()/FS
        rows.append(dict(phase=phase, matched_proxies=matched, eligible_proxies=proxies,
            unmatched_peaks_in_scored_intervals=unmatched, scored_exposure_s=exposure,
            proxy_precision=matched/max(1,matched+unmatched), proxy_recall=matched/max(1,proxies),
            proxy_f1=2*matched/max(1,proxies+matched+unmatched)))
    return rows, predictions


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', default='results/north_natural_picker_adaptation_v1')
    args = ap.parse_args(); out = ROOT / args.output; out.mkdir(exist_ok=False)
    cfgpath = ROOT/'configs/north_published_phase_transfer_2012.yaml'
    cfg = yaml.safe_load(cfgpath.read_text())
    transfers = ROOT/'results/north_date_stability_completion_v1/station_transfers.csv'
    transfer = pd.read_csv(transfers); transfer = transfer[transfer.accepted_waveform_transfer]
    refpath = ROOT/cfg['reference_csv']; references = read_bostock(refpath)
    original = ROOT.parent/'lfe_catalog_repair/results/external/other/seisbench_lfedetect/jms'
    design = dict(status='running_bounded_natural_adaptation', production_validated=False,
        dates={k: cfg[v] for k,v in [('train','train_dates'),('validation','calibration_dates'),('test','test_dates')]},
        split='Disjoint UTC dates, shared source families. No claim of source holdout or dates untouched by earlier matched-filter work.',
        labels='All 32 previously passing date-stable phase transfers; unreviewed individual-arrival proxies.',
        masking='Supervise central 30 s. Unknown-family phase intervals [published time, +60 s] are excluded; known proxy +/-1.5 s restores supervision.',
        background='Up to 128 reference-free 90 s contexts per station/day, deterministic random order; no proof of negative truth.',
        preprocessing='Same cached native filtering, response correction, ZNE and 90->60 s contexts as the completed natural proxy comparison.',
        training=dict(initial_weights='original jms', epochs=8, steps_per_epoch=32, batch_size=64, learning_rate=3e-5,
            weight_decay=1e-4, positive_BCE_weight=10, gaussian_sigma_s=.25, augmentations='polarity and horizontal rotation', seed=20260914),
        selection='Baseline epoch zero included; select best validation macro proxy F1 at threshold .3 and .5 s one-to-one tolerance; evaluate test dates once afterward.',
        limitations=['Phase proxies are incomplete and not independent manual truth.',
            'Unknown arrival mask is a conservative regional 0-60 s assumption, not a verified travel-time bound.',
            'Passing this experiment alone cannot validate production or full-margin generalization.'],
        sources={str(p):sha256(p) for p in [Path(__file__),cfgpath,transfers,refpath,original.with_suffix('.pt'),
            ROOT/'codex_lfe/published_stack_clock.py',ROOT/'codex_lfe/picker_adaptation.py',ROOT/'cache/metadata/cn_pilot_2012.xml']})
    write_json(out/'design_before_training.json', design)
    data = prepare(out, cfg, transfer, references)
    import torch
    os.environ.setdefault('SEISBENCH_CACHE_ROOT', str(ROOT/'cache/seisbench'))
    import seisbench.models as sbm
    torch.manual_seed(20260914); torch.set_num_threads(2)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = sbm.LFEDetect.load(str(original), map_location='cpu', weights_only=True).to(device)
    baseline, _ = evaluate(model, data['validation'], device)
    best = np.mean([r['proxy_f1'] for r in baseline]); best_epoch = 0
    state = {k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
    opt = torch.optim.AdamW(model.parameters(), lr=3e-5, weight_decay=1e-4)
    rng = np.random.default_rng(20260914); history = []
    x, y, mask, _ = data['train']
    for epoch in range(1,9):
        model.train(); loss_sum = 0
        for _ in range(32):
            ix = rng.integers(len(x), size=64); batch = x[ix].copy()
            angles = rng.uniform(-np.pi,np.pi,len(ix)); c,s = np.cos(angles)[:,None],np.sin(angles)[:,None]
            n,e = batch[:,1].copy(),batch[:,2].copy(); batch[:,1]=c*n-s*e; batch[:,2]=s*n+c*e
            batch *= rng.choice([-1,1],len(ix))[:,None,None]
            # Rotation changes the global maximum: renormalize to model contract.
            batch /= np.max(np.abs(batch),axis=(1,2),keepdims=True)
            tensor=torch.from_numpy(batch).to(device); target=torch.from_numpy(y[ix]).to(device); valid=torch.from_numpy(mask[ix]).to(device)
            opt.zero_grad(set_to_none=True); logits=model(tensor,logits=True)
            loss=torch.nn.functional.binary_cross_entropy_with_logits(logits,target,reduction='none',pos_weight=torch.tensor(10.,device=device))[valid].mean()
            loss.backward(); torch.nn.utils.clip_grad_norm_(model.parameters(),1.); opt.step(); loss_sum += float(loss.detach())
        scores,_=evaluate(model,data['validation'],device); metric=float(np.mean([r['proxy_f1'] for r in scores]))
        if metric>best:
            best=metric;best_epoch=epoch;state={k:v.detach().cpu().clone() for k,v in model.state_dict().items()}
        history.append(dict(epoch=epoch,loss=loss_sum/32,validation=scores))
        write_json(out/'training_progress.json',dict(epoch=epoch,total_epochs=8,best_epoch=best_epoch,history=history))
        print('epoch',epoch,'loss',loss_sum/32,'validation macro proxy F1',metric,flush=True)
    torch.save(state,out/'selected_state_dict.pt');model.load_state_dict(state)
    selected,prob=evaluate(model,data['test'],device);np.savez_compressed(out/'selected_test_predictions.npz',probabilities=prob)
    model=sbm.LFEDetect.load(str(original),map_location='cpu',weights_only=True).to(device)
    baseline_test,prob=evaluate(model,data['test'],device);np.savez_compressed(out/'baseline_test_predictions.npz',probabilities=prob)
    rows=[dict(model=name,**r) for name,values in [('baseline',baseline_test),('selected',selected)] for r in values]
    pd.DataFrame(rows).to_csv(out/'test_metrics.csv',index=False)
    write_json(out/'report.json',dict(status='complete_date_separated_natural_adaptation',production_validated=False,
        individual_phase_truth_available=False,best_epoch=best_epoch,baseline_validation=baseline,history=history,test=rows,
        raw_bytes_downloaded=0,weights_sha256=sha256(out/'selected_state_dict.pt')))
    print(pd.DataFrame(rows).to_string(index=False),flush=True)


if __name__=='__main__':main()
