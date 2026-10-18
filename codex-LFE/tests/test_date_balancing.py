from codex_lfe.north_date_stability_completion import balanced_dates

def test_whole_date_balance_is_disjoint_and_deterministic():
    counts={'2012-09-07':40,'2012-09-10':30,'2012-09-13':20,'2012-09-16':10}
    plan=balanced_dates(counts)
    assert plan['n_a']==plan['n_b']==50
    assert not set(plan['a'])&set(plan['b'])
    assert set(plan['a'])|set(plan['b'])==set(counts)
    assert plan==balanced_dates(dict(reversed(list(counts.items()))))

def test_dominant_single_date_cannot_create_independent_half():
    assert balanced_dates({'a':100,'b':5,'c':4}) is None
    assert balanced_dates({'a':100,'b':0}) is None
