import numpy as np
import pandas as pd
import pytest

from codex_lfe.evaluate import match_events, null_summary, score_catalog
from codex_lfe.matched_filter import resolve_candidates


def frame(times, families=None, basis='origin_utc'):
    return pd.DataFrame(dict(time_s=times, family_id=families or ['a'] * len(times), time_basis=basis))


def test_duplicates_cannot_inflate_recall_or_agreement():
    d, r = frame([10, 10.1, 10.2]), frame([10])
    s = score_catalog(d, r)
    assert s['n_matched_same_family'] == 1
    assert s['same_family_recall'] == 1
    assert s['reference_agreement_fraction'] == pytest.approx(1 / 3)


def test_augmenting_path_beats_nearest_greedy():
    # First detection has two possible references; second has only the first.
    d, r = frame([0.4, -0.5]), frame([0, 1])
    assert len(match_events(d, r, 0.6)) == 2


def test_family_identity_and_incompatible_clocks():
    s = score_catalog(frame([10], ['b']), frame([10], ['a']))
    assert s['same_family_recall'] == 0
    assert s['any_family_recall'] == 1
    with pytest.raises(ValueError, match='time bases'):
        score_catalog(frame([10]), frame([10], basis='arrival'))


def test_empty_nan_and_missing_families():
    assert score_catalog(frame([]), frame([0]))['same_family_recall'] == 0
    assert score_catalog(frame([]), frame([]))['same_family_recall'] is None
    assert len(match_events(frame([0], [None]), frame([0], [None]))) == 0
    with pytest.raises(ValueError, match='nonfinite'):
        match_events(frame([np.nan]), frame([0]))


def test_final_assignment_changes_family_recall():
    raw = frame([10, 10.1], ['a', 'b'])
    raw['score_mad'] = [8, 8.2]
    final = resolve_candidates(raw, ambiguity_margin=0.5)
    assert len(final) == 1
    assert final.assignment_status.iloc[0] == 'ambiguous'
    assert score_catalog(final, frame([10], ['a']))['same_family_recall'] == 0
    assert score_catalog(raw, frame([10], ['a']))['same_family_recall'] == 1


def test_exposure_adjusted_controls_and_zero_count():
    s = null_summary(100, 1000, [10, 20], [500, 1000])
    assert s['control_to_observed_mean'] == pytest.approx(0.2)
    assert null_summary(0, 100, [0], [100])['control_to_observed_mean'] is None
    with pytest.raises(ValueError):
        null_summary(1, 0, [0], [100])
