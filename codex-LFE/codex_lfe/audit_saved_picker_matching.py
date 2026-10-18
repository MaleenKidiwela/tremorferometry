"""Quantify legacy matching bias using read-only saved picks and clock offsets."""
from pathlib import Path
import numpy as np
import pandas as pd
from .evaluate import match_events
from .provenance import write_json, sha256


def epoch(series):
    return pd.to_datetime(series, utc=True).astype('datetime64[ns, UTC]').astype('int64').to_numpy()/1e9


def main():
    root=Path('../lfe_catalog_repair/results')
    out=Path('results/prior_picker_test_audit_v2')
    bo=pd.read_parquet(root/'e2_gt_bostock.parquet')
    cal=pd.read_csv(root/'e2_stack_calib.csv')
    clock=float(pd.read_csv(root/'e2_bostock_lin_dt.csv').dt.median())
    sp=cal[cal.station_code.eq('PGC')].set_index('source_family_id').sp_obs_s
    bo['pred']=epoch(bo.time)+clock+bo.ttP_PGC+bo.fam_id.map(sp)
    refs={'bostock_stack':bo.pred.dropna().to_numpy(),
          'bostock_near40':bo[bo.epi_PGC<40].pred.dropna().to_numpy()}
    legacy=pd.read_csv(root/'e2_validation_pgc2010.csv')
    selected=legacy[legacy.phase.eq('S')&legacy['gt'].isin(refs)].copy()
    rows=[]
    for model,group in selected.groupby('model'):
        picks=pd.read_parquet(root/f'e2_picks_pgc2010_{model}.parquet')
        picks=picks[picks.phase.eq('S')].copy();picks['time_s']=epoch(picks.time)
        for r in group.itertuples():
            times=np.sort(picks[(picks.thr_first<=r.threshold)&(picks.prob>=r.threshold)].time_s.to_numpy())
            ref=np.sort(refs[r.gt]+r.offset_s)
            assert len(times)==r.n_picks and len(ref)==r.n_gt
            ix=np.searchsorted(times,ref)
            nearest=np.minimum(abs(times[np.clip(ix,0,len(times)-1)]-ref),abs(times[np.clip(ix-1,0,len(times)-1)]-ref))
            old_recall=float((nearest<=r.tol_s).mean())
            if not np.isclose(old_recall,r.recall,atol=1e-12,rtol=0):raise ValueError('Cannot reproduce old recall')
            # Remove only isolated pick vertices; maximum matching is unchanged.
            left=np.searchsorted(ref,times-r.tol_s,side='left')
            right=np.searchsorted(ref,times+r.tol_s,side='right')
            dt=pd.DataFrame({'time_s':times[right>left],'family_id':'S','time_basis':'legacy_offset_adjusted_S_prediction'})
            rt=pd.DataFrame({'time_s':ref,'family_id':'S','time_basis':'legacy_offset_adjusted_S_prediction'})
            n=len(match_events(dt,rt,tolerance_s=r.tol_s,same_family=True))
            rows.append(dict(model=model,reference=r.gt,threshold=r.threshold,tolerance_s=r.tol_s,
                n_reference=len(ref),n_picks=len(times),legacy_recall=old_recall,one_to_one_matches=n,
                one_to_one_reference_recall=n/len(ref),recall_inflation_percentage_points=100*(old_recall-n/len(ref)),
                legacy_reference_agreement=r.precision,one_to_one_reference_agreement=n/len(times),
                original_clock_offset_s=r.offset_s))
        print('Audited',model,flush=True)
    frame=pd.DataFrame(rows);frame.to_csv(out/'saved_S_picker_matching_comparison.csv',index=False)
    inputs=[root/'e2_gt_bostock.parquet',root/'e2_stack_calib.csv',root/'e2_bostock_lin_dt.csv',root/'e2_validation_pgc2010.csv']+[root/f'e2_picks_pgc2010_{m}.parquet' for m in selected.model.unique()]
    write_json(out/'matching_audit.json',dict(status='complete_scoring_audit_only',rows=len(frame),all_legacy_recalls_reproduced=True,
        max_recall_inflation_percentage_points=float(frame.recall_inflation_percentage_points.max()),
        sources={str(p):sha256(p) for p in inputs+[Path(__file__)]},
        limitations=['Preserves legacy evaluation-fitted offsets solely to isolate matching effect.',
                    'Reference arrivals are predicted, not independent phase labels.',
                    'Legacy inference coverage is not reconstructed here; no operational rate claimed.',
                    'No new inference, waveform acquisition or external-project writes.',
                    'P scoring against S references is excluded, not legitimized.']))
    print(frame[['model','reference','threshold','tolerance_s','legacy_recall','one_to_one_reference_recall']].to_string(index=False),flush=True)

if __name__=='__main__':main()
