import numpy as np
from obspy import Stream,Trace,UTCDateTime
from codex_lfe.response_validity import response_window_mask

def test_complete_template_avoids_taper_and_separate_trace_edges():
    traces=[]
    for channel in ['BHE','BHN','BHZ']:
        traces.append(Trace(np.zeros(2001),dict(network='XX',station='TEST',channel=channel,starttime=UTCDateTime(0),sampling_rate=20)))
    stream=Stream(traces);channels=[tr.id for tr in traces]
    good=response_window_mask(stream,channels,0,2000,20,701)
    assert not good[0] and not good[200] and good[400]
    assert not good[1200]  # template at 60 s ends at 95 s, within end taper/guard
    broken=stream.copy();broken[-1].stats.starttime+=40
    masked=response_window_mask(broken,channels,0,2000,20,701)
    assert not masked[400]
    assert np.all(masked<=good)
