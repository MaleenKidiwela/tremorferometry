"""Conservative support mask for ObsPy's default fractional response taper.

Exclude any template window touching the taper or an additional 10-second edge
of a separately response-corrected trace. This uses only metadata, not signal
amplitude or detection outcomes. It does not certify the response model itself.
"""
import numpy as np

def response_window_mask(stream,channels,start_s,n_windows,fs,template_samples,taper_fraction=.05,edge_s=10):
    if not 0<=taper_fraction<1 or edge_s<0:raise ValueError('Invalid edge policy')
    valid=np.ones(n_windows,dtype=bool)
    span=(template_samples-1)/fs
    for channel in channels:
        channel_valid=np.zeros(n_windows,dtype=bool)
        for tr in stream.select(id=channel):
            rate=float(tr.stats.sampling_rate)
            # Two additional native samples cover taper rounding conventions.
            guard=(np.ceil(tr.stats.npts*taper_fraction/2)+2)/rate+edge_s
            lo=float(tr.stats.starttime)+guard
            hi=float(tr.stats.endtime)-guard-span
            a=max(0,int(np.ceil((lo-start_s)*fs)))
            b=min(n_windows,int(np.floor((hi-start_s)*fs))+1)
            if b>a:channel_valid[a:b]=True
        valid &= channel_valid
    return valid
