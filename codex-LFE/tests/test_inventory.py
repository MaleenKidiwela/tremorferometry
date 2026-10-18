from obspy import UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel

from codex_lfe.inventory import triplet_epochs


def test_triplet_requires_overlapping_epochs_and_independent_axes():
    def channel(code, azimuth, dip, start, end):
        return Channel(code, '', 44, -123, 0, 0, azimuth=azimuth, dip=dip, sample_rate=40,
                       start_date=UTCDateTime(start), end_date=UTCDateTime(end))
    channels = [channel('BHZ', 0, -90, 0, 100), channel('BH1', 30, 0, 10, 90), channel('BH2', 120, 0, 20, 80)]
    inv = Inventory([Network('XX', stations=[Station('A', 44, -123, 0, channels=channels)])], 'test')
    table = triplet_epochs(inv, 0, 100)
    assert len(table) == 1
    assert table.duration_days_metadata.iloc[0] == 60 / 86400
    channels[2].start_date = UTCDateTime(95)
    channels[2].end_date = UTCDateTime(100)
    assert triplet_epochs(inv, 0, 100).empty
