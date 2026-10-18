import pandas as pd
from codex_lfe.north_extended_validation import family_basis,score_family_rows


def test_overlapping_family_clocks_cannot_steal_each_others_reference():
    # Two families' physically different anchor phases happen to have the same
    # numeric timestamp. Only the matching source group may claim the reference.
    candidates=pd.DataFrame([dict(candidate_id='b1',family_id='b',realization=0,time_s=125.,time_basis=family_basis('b'))])
    reference=pd.DataFrame([dict(reference_id='a1',family_id='a',time_s=100.)])
    spans=pd.DataFrame([dict(family_id=f,start_s=0.,end_s=200.) for f in ['a','b']])
    scores,matches=score_family_rows(candidates,reference,spans,dict(a=25.,b=20.))
    assert scores.n_eligible_references.sum()==1
    assert scores.n_reference_matches.sum()==0
    assert matches.empty


def test_one_to_one_matching_and_half_open_reference_coverage():
    candidates=pd.DataFrame([dict(candidate_id='a1',family_id='a',realization=0,time_s=125.,time_basis=family_basis('a'))])
    reference=pd.DataFrame([dict(reference_id='r'+str(i),family_id='a',time_s=t) for i,t in enumerate([100.,100.1,175.])])
    spans=pd.DataFrame([dict(family_id='a',start_s=0.,end_s=200.)])
    scores,matches=score_family_rows(candidates,reference,spans,dict(a=25.))
    assert scores.n_eligible_references.sum()==2
    assert scores.n_reference_matches.sum()==1
    assert len(matches)==1
