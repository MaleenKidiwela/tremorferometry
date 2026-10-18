from codex_lfe.temporal_waveform_probes import covered_seconds


def test_union_does_not_fill_gap_or_double_count_overlaps():
    assert covered_seconds([(-5, 10), (8, 15), (20, 40)], 0, 30) == 25
    assert covered_seconds([], 0, 30) == 0
    assert covered_seconds([(0, 15), (15, 30)], 0, 30) == 30
