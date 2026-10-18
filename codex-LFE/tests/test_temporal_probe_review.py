from codex_lfe.temporal_probe_review import sample_grid_complete


def test_subsample_edge_is_not_an_internal_gap():
    assert sample_grid_complete([(0.019, 30.019, .025)], 0, 30)
    assert not sample_grid_complete([(.051, 30.026, .025)], 0, 30)
    assert not sample_grid_complete([(.019, 10.019, .025), (10.044, 30.019, .025)], 0, 30)
    assert not sample_grid_complete([(.019, 29.99, .025)], 0, 30)
    assert sample_grid_complete([(.019, 10.019, .025), (10.019, 30.019, .025)], 0, 30)
    assert not sample_grid_complete([], 0, 30)
