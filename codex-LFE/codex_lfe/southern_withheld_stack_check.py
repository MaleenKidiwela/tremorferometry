"""Fixed-time withheld-station stacks of existing southern candidate cohorts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from obspy import read,read_inventory,UTCDateTime
from .waveforms import preprocess_stream
from .neural import orientation_matrix
from .published_pickle import read_published_pickle
from .ducellier_waveform_review import rms_scaled
from .provenance import sha256,write_json


def separated_rows(frame, duration):
    """Choose nonoverlapping windows using detector score, never withheld data."""
    selected=[]
    for row in frame.sort_values(['score_mad','aligned_time_s'],ascending=[False,True]).itertuples():
        if all(abs(row.aligned_time_s-t)>=duration for _,t in selected):
            selected.append((row.Index,row.aligned_time_s))
    return frame.loc[[i for i,_ in selected]].sort_values('aligned_time_s')


def similarity(a,b):
    x=rms_scaled(a);y=rms_scaled(b)
    return float(np.mean(x*y))


def main():
    root=Path(__file__).resolve().parents[1];out=root/'results/southern_withheld_stack_check_v1';out.mkdir(exist_ok=False)
    candidate_path=root/'results/ducellier_multifamily_november23_holdout_edges_v2/within_family_review_candidates.csv'
    candidates=pd.read_csv(candidate_path);candidates=candidates[(candidates.control==-1)&(candidates.threshold_mad==8)]
    mapping_path=root/'results/southern_benchmark_inventory_v1/comparison_bank_source_families.csv';mapping=pd.read_csv(mapping_path)
    acquisition_path=root/'results/ducellier_november23_acquisition_v1/acquisition.json'
    acquisition=json.loads(acquisition_path.read_text());record=next(r for r in acquisition['records'] if r['station']=='PB.B039' and r['status']=='received')
    config=json.loads((acquisition_path.parent/'design.json').read_text());start=float(UTCDateTime(config['start']));duration=config['duration_s']
    write_json(out/'design_before_waveform_check.json',dict(status='running',production_validated=False,
        candidates='All previously fixed 8-MAD observed next-day candidate rows, grouped within published families.',
        heldout_station='PB.B039, excluded from detection',alignment='Training-derived canonical-template starts; no adjustment using held-out waveforms.',
        cohorts=['all candidate rows','primary: nonoverlapping 35 s windows selected by descending detector score then time'],
        statistic='Signed equal-component correlation of the mean of component-RMS-normalized windows with the canonical published template; no maximization over lag or sign.',
        controls='100 distinct common time shifts per family from signed integer 60-600 s, retaining temporal spacing and candidate counts; exact same operation as observed.',
        interpretation='Cohort-level independent-station evidence, not individual phase truth or a false-discovery fraction.',
        sources={str(p):sha256(p) for p in [Path(__file__),candidate_path,mapping_path,acquisition_path]}))
    response_path=root/'results/ducellier_roster_packets_v1/responses/PB.xml';inv=read_inventory(response_path)
    raw=Path(record['path'])
    if sha256(raw)!=record['sha256']:raise ValueError('Raw checksum changed')
    stream=read(str(raw));channels=record['channels']
    for tr in stream:
        if str(inv.get_response(tr.id,tr.stats.starttime))!=str(inv.get_response(tr.id,tr.stats.endtime)):raise ValueError('Response epoch crossing')
        tr.remove_response(inventory=inv,output='VEL',pre_filt=(.2,.5,10,15),water_level=80)
    wave=preprocess_stream(stream,channels,'PB.B039',start,duration+60,fs=20,band=(1.5,9))
    current=orientation_matrix([inv.get_orientation(c,UTCDateTime(start)) for c in channels])
    zne=np.linalg.solve(current,wave.x);valid=wave.valid.all(axis=0)
    # The saved detector already excluded the response-removal taper. Apply the
    # same conservative end guards to all newly shifted stack windows.
    valid[:int((duration+60)*.05*20)+200]=False
    valid[-int((duration+60)*.05*20)-200:]=False
    rng=np.random.default_rng(20260914);summaries=[];control_rows=[];members=[];inputs={str(raw):sha256(raw),str(response_path):sha256(response_path)}
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for family,frame in candidates.groupby('published_family_id'):
        canonical=sorted(mapping[mapping.published_family_id==family].template_id)[0]
        paths=[root/'results/southern_benchmark_inventory_v1/source/src/templates_v1'/canonical/f'B039_{c.split(".")[-1]}.pkl' for c in channels]
        values=[read_published_pickle(p) for p in paths];template=np.asarray([tr.data for tr,_ in values]);n=template.shape[1]
        target=orientation_matrix([angle for _,angle in values]);x=target@zne
        inputs.update({str(p):sha256(p) for p in paths})
        selected=separated_rows(frame,n/20)
        shifts=np.r_[np.arange(-600,-59),np.arange(60,601)];rng.shuffle(shifts)
        indices=np.rint((frame.aligned_time_s.to_numpy()-start)*20).astype(int)
        def covered(ix):return ix>=0 and ix+n<=len(valid) and valid[ix:ix+n].all()
        if not all(covered(i) for i in indices):raise ValueError('Observed stack lacks complete untapered support')
        shifts=[int(s) for s in shifts if all(covered(i+s*20) for i in indices)][:100]
        if len(shifts)!=100:raise ValueError('Insufficient matched-exposure shifts')
        archive={'template':template};family_rows=[]
        for cohort,table in [('all',frame),('nonoverlapping',selected)]:
            base=np.rint((table.aligned_time_s.to_numpy()-start)*20).astype(int)
            members.extend(dict(family=family,cohort=cohort,source_index=int(i),time_s=float(row.aligned_time_s)) for i,row in table.iterrows())
            scores=[]
            for realization,shift in enumerate([0]+shifts):
                windows=np.asarray([rms_scaled(x[:,i+shift*20:i+shift*20+n]) for i in base])
                stack=windows.mean(axis=0);cc=similarity(stack,template);scores.append(cc)
                control_rows.append(dict(family=family,cohort=cohort,realization=realization,shift_s=shift,n_windows=len(base),stack_cc=cc))
                if realization==0:archive[cohort+'_stack']=stack;archive[cohort+'_windows']=windows
            null=np.asarray(scores[1:]);row=dict(family=family,cohort=cohort,n_windows=len(base),observed_cc=scores[0],
                shifted_median=float(np.median(null)),shifted_max=float(null.max()),shifted_min=float(null.min()),
                rank_fraction_including_observed=float((1+(null>=scores[0]).sum())/101))
            summaries.append(row);family_rows.append(row)
        np.savez_compressed(out/(family+'.npz'),**archive)
        fig,axes=plt.subplots(3,1,figsize=(11,6),sharex=True,layout='constrained')
        for c,ax in enumerate(axes):
            ax.plot(np.arange(n)/20,rms_scaled(archive['nonoverlapping_stack'])[c],color='.25',lw=.8,label='Withheld-station candidate stack')
            ax.plot(np.arange(n)/20,rms_scaled(template)[c],color='#c33732',lw=.8,alpha=.8,label='Published template')
            ax.set_ylabel(channels[c].split('.')[-1])
        axes[0].legend(fontsize=8);axes[-1].set_xlabel('Seconds relative to canonical template start')
        row=family_rows[-1];fig.suptitle(f'{family}: {len(selected)} nonoverlapping windows at withheld PB.B039\nCC {row["observed_cc"]:.3f}; shifted maximum {row["shifted_max"]:.3f}; fixed timing and polarity')
        fig.savefig(out/(family+'.png'),dpi=130);plt.close(fig)
        print(family,family_rows,flush=True)
    pd.DataFrame(summaries).to_csv(out/'summary.csv',index=False);pd.DataFrame(control_rows).to_csv(out/'control_scores.csv',index=False);pd.DataFrame(members).to_csv(out/'cohort_members.csv',index=False)
    write_json(out/'report.json',dict(status='complete_withheld_station_cohort_check',production_validated=False,
        families=len(candidates.published_family_id.unique()),individual_events_certified=0,raw_bytes_downloaded=0,inputs=inputs,
        limitations=['Shifts are correlated controls, not independent experiments; rank fractions are descriptive.',
            'Four family comparisons and two related cohorts are reported without selective promotion.',
            'Published-template phase labels and independent single-event truth are unavailable.']))


if __name__=='__main__':main()
