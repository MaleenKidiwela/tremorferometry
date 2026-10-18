import pandas as pd
import pytest
from codex_lfe.picker_evaluation import score_picks as evaluate_picks, WINDOW_COLUMNS, LABEL_COLUMNS, BASIS


def score_picks(p,w,l,threshold,tolerance,split,coverage=None):
    if coverage is None:
        coverage=pd.DataFrame([{'station':s,'phase':phase,'start_s':0.,'end_s':3600.} for s in w.station.unique() for phase in ['P','S']])
    return evaluate_picks(p,w,l,threshold,tolerance,split,coverage)

def data(status='complete', sampling='random'):
    w=pd.DataFrame([['w1','CN.PGC',0.,3600.,'development',sampling,status,'reviewer-A']],columns=WINDOW_COLUMNS)
    labels=pd.DataFrame([['l1','w1','P',100.,.1,'lfe'],['l2','w1','S',110.,.1,'lfe']],columns=LABEL_COLUMNS)
    p=pd.DataFrame([['CN.PGC','P',100.1,.8,BASIS],['CN.PGC','S',110.2,.7,BASIS]],columns=['station','phase','time_s','model_score','time_basis'])
    return p,w,labels


def test_duplicate_picks_do_not_inflate_recovery():
    p,w,l=data();p=pd.concat([p,p.iloc[:1]],ignore_index=True)
    report,_=score_picks(p,w,l,.5,.5,'development')
    phase=report['results'][0]
    assert phase['phase']=='P' and phase['n_matched']==1
    assert phase['review_panel_precision']==.5 and phase['unmatched_picks_per_reviewed_station_hour']==1


@pytest.mark.parametrize('status',['pending','ambiguous','unusable'])
def test_unreviewed_windows_cannot_be_used_as_negatives(status):
    p,w,l=data(status)
    with pytest.raises(ValueError,match='no completely reviewed'):score_picks(p,w,l,.5,.5,'development')


def test_wide_uncertainty_excludes_window_not_just_reference_label():
    p,w,l=data();l.loc[0,'uncertainty_s']=2
    with pytest.raises(ValueError,match='no completely reviewed'):score_picks(p,w,l,.5,.5,'development')


def test_wrong_station_and_phase_cannot_earn_credit():
    p,w,l=data();p.loc[0,'phase']='S';p.loc[1,'station']='CN.VGZ'
    report,_=score_picks(p,w,l,.5,.5,'development')
    assert sum(r['n_matched'] for r in report['results'])==0
    assert sum(r['unmatched_predictions'] for r in report['results'])==1


def test_empty_labels_in_completed_review_are_real_negative_exposure():
    p,w,l=data();l=l.iloc[:0]
    report,_=score_picks(p,w,l,.5,.5,'development')
    assert sum(r['unmatched_predictions'] for r in report['results'])==2
    assert all(r['recall'] is None for r in report['results'])


def test_ordinary_earthquake_is_not_lfe_positive():
    p,w,l=data();l['event_class']='ordinary_earthquake'
    report,_=score_picks(p,w,l,.5,.5,'development')
    assert sum(r['n_matched'] for r in report['results'])==0


def test_station_overlap_across_splits_is_rejected():
    p,w,l=data();other=w.copy();other.window_id='w2';other.split='test';w=pd.concat([w,other],ignore_index=True)
    with pytest.raises(ValueError,match='overlapping'):score_picks(p,w,l,.5,.5,'development')


def test_enriched_and_random_rates_are_not_pooled():
    p,w,l=data();other=w.copy();other.window_id='w2';other.station='CN.SNB';other.sampling='reference_enriched';w=pd.concat([w,other],ignore_index=True)
    report,_=score_picks(p,w,l,.5,.5,'development')
    assert len(report['results'])==4


def test_missing_reviewer_is_not_completed_review():
    p,w,l=data();w.reviewer=''
    with pytest.raises(ValueError,match='reviewer'):score_picks(p,w,l,.5,.5,'development')


def test_missing_prediction_coverage_is_not_silent_negative():
    p,w,l=data()
    coverage=pd.DataFrame([{'station':'CN.PGC','phase':phase,'start_s':0.,'end_s':1000.} for phase in ['P','S']])
    with pytest.raises(ValueError,match='no completely reviewed'):score_picks(p,w,l,.5,.5,'development',coverage)
