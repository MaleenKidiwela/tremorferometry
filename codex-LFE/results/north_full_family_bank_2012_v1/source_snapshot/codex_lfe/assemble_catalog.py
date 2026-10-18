"""Assemble the frozen northern candidate catalog with explicit validation strata."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .evaluate import match_events, score_catalog, null_summary
from .provenance import sha256, write_json
from .reference import read_bostock


def main():
    base=Path('results/method_interleaved_mf_v1')
    quiet=Path('results/north_catalog_block_20100906_10_v1')
    active=Path('results/north_catalog_active_holdout_v1')
    cfg=json.loads((base/'manifest.json').read_text())['config']
    groups=json.loads((base/'training.json').read_text())['families_usable']
    sources=[(base,'calibration',cfg['calibration_dates'],'development_calibration'),
             (base,'test',cfg['test_dates'],'development_evaluation'),
             (quiet,'test',json.loads((quiet/'selection_before_test.json').read_text())['config']['test_dates'],'heldout_quiet'),
             (active,'test',json.loads((active/'selection_before_test.json').read_text())['config']['test_dates'],'heldout_active')]
    observed=[]; controls=[]; masks=[]; exposures=[]; inputs=[]; refs=[]
    seen=set()
    for run,split,dates,label in sources:
        if seen & set(dates): raise ValueError('Overlapping catalog blocks')
        seen.update(dates)
        rp=run/(f'{split}_eligible_reference.parquet' if run==base else 'eligible_reference.parquet')
        r=pd.read_parquet(rp);r['validation_split']=label;refs.append(r);inputs.append(rp)
        for d in dates:
            paths=[run/f'{split}_{d}_{suffix}.parquet' for suffix in ['candidates','eligible_intervals','exposure']]
            inputs.extend(paths)
            c,m,e=[pd.read_parquet(p) for p in paths]
            c=c[c.threshold==9].copy();c['validation_split']=label
            observed.append(c[c.realization==0]);controls.append(c[c.realization!=0])
            m['validation_split']=label;e['validation_split']=label;masks.append(m);exposures.append(e)
    out=Path('results/north_candidate_catalog_v1');out.mkdir(exist_ok=False)
    cat=pd.concat(observed,ignore_index=True).sort_values('time_s').reset_index(drop=True)
    ref=pd.concat(refs,ignore_index=True)
    assert not ref.reference_id.duplicated().any()
    assert not cat.duplicated(['date','time_s']).any()
    cat.insert(0,'candidate_id',[f'codex_north_v1_{i:06d}' for i in range(1,len(cat)+1)])
    cat['arrival_anchor_utc']=pd.to_datetime(cat.time_s,unit='s',utc=True).astype(str)
    cat['validation_status']='experimental_unverified'
    published=read_bostock(cfg['reference_csv'])
    families=published[published.family_id.isin(groups)][['family_id','latitude','longitude','depth_km']].drop_duplicates()
    families=families.rename(columns={c:'source_group_'+c for c in ['latitude','longitude','depth_km']})
    cat=cat.merge(families,on='family_id',how='left',validate='many_to_one')
    matched=match_events(cat,ref,cfg['reference_tolerance_s'],True)
    cat['matched_reference_id']=pd.Series([None]*len(cat),dtype=object)
    cat['reference_residual_s']=np.nan
    for row in matched.itertuples():
        cat.loc[row.detection_index,'matched_reference_id']=ref.iloc[row.reference_index].reference_id
        cat.loc[row.detection_index,'reference_residual_s']=row.residual_s
    ctl=pd.concat(controls,ignore_index=True);mask=pd.concat(masks,ignore_index=True);exp=pd.concat(exposures,ignore_index=True)
    reports={}
    for _,_,_,label in sources:
        c=cat[cat.validation_split==label];r=ref[ref.validation_split==label]
        e=exp[exp.validation_split==label].groupby('realization').exposure_s.sum()
        n=ctl[ctl.validation_split==label].groupby('realization').size()
        reports[label]={'reference':score_catalog(c,r), 'controls':null_summary(len(c),e[0],[int(n.get(i,0)) for i in range(1,17)],[float(e[i]) for i in range(1,17)])}
    family_rows=[]
    for f in groups:
        c=cat[cat.family_id==f];r=ref[ref.family_id==f]
        family_rows.append({'family_id':f,**score_catalog(c,r)})
    cat.to_csv(out/'candidate_catalog.csv',index=False);cat.to_parquet(out/'candidate_catalog.parquet',index=False)
    families.to_csv(out/'source_groups.csv',index=False)
    pd.DataFrame(family_rows).to_json(out/'family_metrics.json',orient='records',indent=2)
    mask.to_parquet(out/'coverage_intervals.parquet',index=False);exp.to_parquet(out/'control_exposure.parquet',index=False)
    ctl.to_parquet(out/'shifted_control_candidates.parquet',index=False)
    ref.to_parquet(out/'eligible_reference.parquet',index=False)
    by_date=cat.groupby('date').size().reindex(sorted(seen),fill_value=0).rename('candidate_count').reset_index()
    by_date.to_csv(out/'daily_counts.csv',index=False)
    write_json(out/'validation.json',reports)
    write_json(out/'manifest.json',{'status':'experimental_unverified','dates':sorted(seen),'candidate_count':len(cat),'source_group_count':len(groups),
        'frozen_threshold_mad':9,'time_basis':'arrival_anchor:CN.PGC:training_stack_peak',
        'inputs':{str(p.resolve()):sha256(p) for p in inputs+[Path(__file__),base/'templates.npz',base/'manifest.json']},
        'limitations':['Northern known-source candidates only, not a margin-wide release','Source coordinates are inherited group coordinates, not individual hypocentres','No calibrated magnitude or origin time','Reference agreement is not measured precision','Development dates must not be pooled as independent validation','Zero candidates does not establish absence of LFEs']})
    (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Northern LFE candidate catalog</title><style>body{font:15px system-ui;margin:2em}table{border-collapse:collapse}td,th{padding:.4em;border:1px solid #ddd}th{background:#eee}</style><h1>Northern LFE candidate catalog v1</h1><p>Experimental, unverified candidates. Frozen eight-second matcher, 9 MAD. Times are PGC arrival anchors; source coordinates describe reference groups. No event precision, calibrated magnitude, or margin-wide completeness claim.</p><p><a href="candidate_catalog.csv">Candidate CSV</a> · <a href="validation.json">Validation by split</a> · <a href="source_groups.csv">Source groups</a> · <a href="manifest.json">Provenance and limits</a></p><h2>Daily counts</h2>'+by_date.to_html(index=False)+'<h2>Candidates</h2>'+cat[['candidate_id','arrival_anchor_utc','family_id','score_mad','validation_split','matched_reference_id']].to_html(index=False))
    print(json.dumps({'candidate_count':len(cat),'dates':len(seen),'validation':reports},indent=2),flush=True)

if __name__=='__main__':main()
