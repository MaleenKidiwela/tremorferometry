"""A known-answer multistation experiment with real failure modes mixed together."""
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt

from codex_lfe.evaluate import score_catalog
from codex_lfe.matched_filter import align_to_anchor, extract_candidates, normalized_cc, resolve_candidates


def test_network_injections_gaps_site_artifacts_and_final_scoring():
    rng = np.random.default_rng(22)
    fs, n, length = 20, 20000, 160
    sos = butter(3, [1, 7], fs=fs, output='sos', btype='bandpass')
    templates = {f: {s: sosfiltfilt(sos, rng.normal(size=(3, length)), axis=1)
                     for s in ['A', 'B', 'C', 'D']} for f in ['family1', 'family2']}
    waves = {s: sosfiltfilt(sos, rng.normal(size=(3, n)), axis=1) for s in ['A', 'B', 'C', 'D']}
    valid = {s: np.ones((3, n), bool) for s in waves}
    offsets = {'A': -1, 'B': 0, 'C': 1.5, 'D': 3}
    events = [(100, 'family1'), (230, 'family2'), (400, 'family1'), (550, 'family2')]
    for time, family in events:
        for station in waves:
            i = int((time + offsets[station]) * fs)
            waves[station][:, i:i + length] += templates[family][station] * 20
    # A huge one-site template-shaped artifact must not count as a network event.
    i = int((700 + offsets['A']) * fs)
    waves['A'][:, i:i + length] += templates['family1']['A'] * 1000
    # An all-station recording gap cannot create a normalized-CC detection.
    for s in waves:
        valid[s][:, 15000:16000] = False
        waves[s][:, 15000:16000] = 0
    rows = []
    for family in templates:
        sites = {}
        for s in waves:
            cc, good = normalized_cc(waves[s], templates[family][s], valid[s])
            sites[s], _ = align_to_anchor(cc, good, offsets[s], fs, n)
        found, exposure = extract_candidates(sites, fs, 0, 8, support_mad=3, min_sites=3)
        for row in found:
            row.update(family_id=family, time_basis='synthetic_origin_utc')
        rows.extend(found)
        assert not exposure[15100:15900].any()
    final = resolve_candidates(pd.DataFrame(rows))
    ref = pd.DataFrame([dict(time_s=t, family_id=f, time_basis='synthetic_origin_utc') for t, f in events])
    score = score_catalog(final, ref, tolerance_s=0.1)
    assert score['n_detections'] == 4
    assert score['n_matched_same_family'] == 4
    assert score['same_family_recall'] == 1
