"""Deterministic physical-site representatives chosen without waveform scores."""
from .reuse_association import site_groups

def representative_sites(stations,distance_km=1.):
    if stations.id.duplicated().any():raise ValueError('Duplicate station metadata')
    groups=site_groups(stations,distance_km)
    # site_groups names each connected proximity group by its smallest SEED ID.
    selected=stations[stations.id.map(groups)==stations.id].sort_values('id').copy()
    return selected,groups
