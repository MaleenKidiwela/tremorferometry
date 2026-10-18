import numpy as np
from obspy import Stream, Trace, UTCDateTime
import pytest

from codex_lfe.waveforms import preprocess_stream


def trace(data, component, start=0.0, rate=40):
    return Trace(np.asarray(data, dtype=float), header=dict(network='XX', station='TEST', location='',
                 channel='BH' + component, sampling_rate=rate, starttime=UTCDateTime(start)))


IDS = ['XX.TEST..BHZ', 'XX.TEST..BHN', 'XX.TEST..BHE']


def test_components_and_fractional_start_align_to_utc():
    rate, duration = 40, 60
    waves = []
    # All are samples of the same physical-time signal, with different phase grids.
    for component, phase, amp in [('E', 0.0125, 3), ('Z', 0.0, 1), ('N', 0.00625, 2)]:
        t = phase + np.arange(rate * duration) / rate
        waves.append(trace(amp * np.sin(2 * np.pi * 2 * t), component, phase))
    w = preprocess_stream(Stream(waves), IDS, 'SITE', 0, duration, edge_s=2)
    use = w.valid.all(axis=0)
    truth = np.sin(2 * np.pi * 2 * np.arange(duration * 20) / 20)
    for c, amp in enumerate([1, 2, 3]):
        assert np.corrcoef(w.x[c, use], truth[use])[0, 1] > 0.999
        assert np.std(w.x[c, use]) / np.std(truth[use]) == pytest.approx(amp, rel=0.06)


def test_gaps_and_zero_runs_never_become_valid_filtered_data():
    t = np.arange(60 * 40) / 40
    st = Stream()
    for c in ['Z', 'N', 'E']:
        st += trace(np.sin(2 * np.pi * 2 * t[:800]), c)
        st += trace(np.sin(2 * np.pi * 2 * t[1200:]), c, 30)
    w = preprocess_stream(st, IDS, 'SITE', 0, 60, edge_s=2)
    assert not w.valid[:, 18 * 20:32 * 20].any()
    assert w.valid[:, 5 * 20:15 * 20].all()
    z = Stream([trace(np.zeros(2400), c) for c in ['Z', 'N', 'E']])
    assert not preprocess_stream(z, IDS, 'SITE', 0, 60, edge_s=2).valid.any()


def test_missing_channel_stays_missing():
    t = np.arange(2400) / 40
    st = Stream([trace(np.sin(2 * np.pi * 2 * t), c) for c in ['Z', 'N']])
    w = preprocess_stream(st, IDS, 'SITE', 0, 60, edge_s=2)
    assert not w.valid[2].any()
    assert w.valid[:2].any()


def test_overlapping_conflict_masked():
    t = np.arange(2400) / 40
    st = Stream([trace(np.sin(2 * np.pi * 2 * t), c) for c in ['Z', 'N', 'E']])
    st += trace(np.ones(400) * 999, 'N', 20)
    w = preprocess_stream(st, IDS, 'SITE', 0, 60, edge_s=2)
    assert not w.valid[1, 18 * 20:32 * 20].any()


def test_native_phase_change_rejected():
    st = Stream([trace(np.ones(800), 'Z'), trace(np.ones(800), 'Z', 20.005)])
    with pytest.raises(ValueError, match='phase changes'):
        preprocess_stream(st, IDS, 'SITE', 0, 60)
