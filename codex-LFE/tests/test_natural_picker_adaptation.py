import numpy as np
from codex_lfe.natural_picker_adaptation import targets_and_mask


def test_unknown_sources_mask_only_their_phase_and_known_arrival_restores_local_label():
    y,mask=targets_and_mask(100.,[[125.],[]],[[110.],[]])
    assert not mask[0,300:470].any()
    assert mask[0,500]
    assert y[0,500]==1
    assert not mask[0,600:900].any()
    assert mask[1,300:900].all()
    assert not mask[:,:300].any()
    assert not mask[:,900:].any()


def test_multiple_proxies_merge_soft_targets_without_exceeding_one():
    y,mask=targets_and_mask(0.,[[25.,25.05],[]],[[],[]])
    assert y.max()==1
    assert y[0,500]==1 and y[0,501]==1
    np.testing.assert_array_equal(mask[:,300:900],True)
