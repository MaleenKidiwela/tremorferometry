import numpy as np
from codex_lfe.shared_reference_panel import in_intervals


def test_union_membership_with_nested_intervals():
    times=np.array([-1.,0.,3.,8.,10.,12.,15.])
    intervals=np.array([[2.,4.],[0.,10.],[12.,15.]])
    assert in_intervals(times,intervals).tolist()==[False,True,True,True,False,True,False]
    assert not in_intervals(times,[]).any()
