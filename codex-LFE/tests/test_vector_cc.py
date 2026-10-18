import numpy as np
import pytest
from codex_lfe.matched_filter import normalized_cc
from codex_lfe.templates import vector_cc


def test_vector_coefficients_match_direct_demeaned_windows():
    rng=np.random.default_rng(382)
    template=rng.normal(size=(3,37))*np.array([1.,3.,.2])[:,None]
    data=rng.normal(size=(3,211))*np.array([2.,4.,.1])[:,None]
    actual,valid=normalized_cc(data,template,channel_combination='vector')
    expected=np.array([vector_cc(data[:,i:i+37],template) for i in range(175)])
    assert valid.all()
    np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=0)


def test_vector_gap_and_polarity_reversal():
    rng=np.random.default_rng(482)
    template=rng.normal(size=(3,20))
    data=np.concatenate([template,-template,template],axis=1)
    valid=np.ones_like(data,dtype=bool);valid[1,45]=False
    cc,mask=normalized_cc(data,template,valid,channel_combination='vector')
    assert cc[0]==pytest.approx(1)
    assert cc[20]==pytest.approx(-1)
    assert not mask[40] and np.isnan(cc[40])


def test_vector_component_energy_is_not_equal_channel_weight():
    rng=np.random.default_rng(983)
    template=rng.normal(size=(3,100));template[0]*=100
    data=template.copy();data[1:]*=-1
    equal,_=normalized_cc(data,template)
    vector,_=normalized_cc(data,template,channel_combination='vector')
    assert equal[0]==pytest.approx(-1/3)
    assert vector[0]>.99


def test_reject_unknown_channel_combination():
    with pytest.raises(ValueError,match='channel combination'):
        normalized_cc(np.ones((3,20)),np.ones((3,4)),channel_combination='invalid')
