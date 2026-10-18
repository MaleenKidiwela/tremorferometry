import numpy as np
import pytest

from codex_lfe.matched_filter import (align_to_anchor, combine_sites, extract_candidates,
                                     normalized_cc, shifted_sites)


def test_signed_normalized_cc_and_reference_loop():
    rng = np.random.default_rng(4)
    tpl = rng.normal(size=(3, 40))
    x = rng.normal(size=(3, 500))
    x[:, 200:240] = tpl * 3 + 17
    x[:, 350:390] = -tpl
    cc, good = normalized_cc(x, tpl)
    assert cc[200] == pytest.approx(1)
    assert cc[350] == pytest.approx(-1)
    expected = np.mean([np.corrcoef(x[k, 20:60], tpl[k])[0, 1] for k in range(3)])
    assert cc[20] == pytest.approx(expected)
    assert good.all()


def test_zero_and_invalid_windows_rejected():
    rng = np.random.default_rng(2)
    tpl = rng.normal(size=(3, 20))
    x = np.zeros((3, 200))
    assert not normalized_cc(x, tpl)[1].any()
    x = rng.normal(size=(3, 200))
    valid = np.ones_like(x, bool)
    valid[1, 80:100] = False
    cc, good = normalized_cc(x, tpl, valid)
    assert not good[61:100].any()
    assert np.isnan(cc[61:100]).all()


def test_moveout_sign_and_no_wraparound():
    x = np.zeros(100)
    x[50] = 1
    out, valid = align_to_anchor(x, np.ones(100, bool), 2, 10, 100)
    assert np.nanargmax(out) == 30
    assert not valid[80:].any()
    with pytest.raises(ValueError, match='grid'):
        align_to_anchor(x, np.ones(100, bool), 0.025, 20, 100)


def test_colocated_duplicates_cannot_increase_site_weight():
    rng = np.random.default_rng(1)
    a, b = rng.normal(size=(2, 100))
    sites = combine_sites({'a': a, 'a2': a, 'b': b}, {'a': 'A', 'a2': 'A', 'b': 'B'})
    np.testing.assert_array_equal(sites['A'], a)
    np.testing.assert_allclose(np.mean(list(sites.values()), axis=0), (a + b) / 2)


def test_site_shift_preserves_missing_data_and_does_not_wrap():
    x = np.arange(20, dtype=float)
    x[10] = np.nan
    y = shifted_sites({'A': x}, {'A': 3})['A']
    assert np.isnan(y[:3]).all() and np.isnan(y[13])
    assert y[3] == 0 and y[-1] == 16


def test_coherent_event_detected_single_site_glitch_rejected():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(4, 4000))
    x[:, 1000] += 20
    x[0, 2000] += 200
    sites = {str(i): row for i, row in enumerate(x)}
    rows, good = extract_candidates(sites, 20, 0, 8, min_sites=3, support_mad=2)
    times = [r['time_s'] for r in rows]
    assert 50 in times
    assert 100 not in times
    assert good.all()


def test_torch_cpu_parity():
    pytest.importorskip('torch')
    rng = np.random.default_rng(6)
    x, t = rng.normal(size=(3, 400)), rng.normal(size=(3, 40))
    a, av = normalized_cc(x, t)
    b, bv = normalized_cc(x, t, backend='torch')
    np.testing.assert_allclose(a, b, atol=1e-10)
    np.testing.assert_array_equal(av, bv)


def test_threshold_filtering_equals_independent_peak_selection():
    rng = np.random.default_rng(94)
    sites = {str(s): rng.normal(size=30000) for s in range(4)}
    for row in sites.values():
        row[[300, 310, 330, 3000, 3020, 10000]] += [4, 5, 4, 5, 9, 12]
    low, low_valid = extract_candidates(sites, 20, 0, 4, min_sites=3)
    for threshold in [5, 6, 7, 8, 9, 10]:
        direct, valid = extract_candidates(sites, 20, 0, threshold, min_sites=3)
        assert [r for r in low if r['score_mad'] >= threshold] == direct
        np.testing.assert_array_equal(low_valid, valid)
