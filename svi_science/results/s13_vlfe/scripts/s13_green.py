"""S13: whole-space displacement Green's function for a moment-tensor point source (Aki & Richards 4.29 generalised to a
moment tensor; near + intermediate + far terms), with a crude free-surface factor 2.  Coordinates: x north, y east,
z down (m).  Used only as a SHAPE / sign template and an order-of-magnitude amplitude for 20-50 s waves; the absolute
moment scale carries a ~x2 uncertainty (no layering, no true free surface) -> relative moments are the product.
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt

VP, VS, RHO = 6900.0, 3950.0, 3000.0       # lower crust / upper mantle wedge at 30-40 km (CVM-ish)


def mt_ned(strike, dip, rake, m0=1.0):
    """Moment tensor in (north, east, down) from Aki & Richards (Box 4.4)."""
    s, d, r = np.radians([strike, dip, rake])
    Mxx = -m0 * (np.sin(d) * np.cos(r) * np.sin(2 * s) + np.sin(2 * d) * np.sin(r) * np.sin(s) ** 2)
    Mxy = m0 * (np.sin(d) * np.cos(r) * np.cos(2 * s) + 0.5 * np.sin(2 * d) * np.sin(r) * np.sin(2 * s))
    Mxz = -m0 * (np.cos(d) * np.cos(r) * np.cos(s) + np.cos(2 * d) * np.sin(r) * np.sin(s))
    Myy = m0 * (np.sin(d) * np.cos(r) * np.sin(2 * s) - np.sin(2 * d) * np.sin(r) * np.cos(s) ** 2)
    Myz = -m0 * (np.cos(d) * np.cos(r) * np.sin(s) - np.cos(2 * d) * np.sin(r) * np.cos(s))
    Mzz = m0 * np.sin(2 * d) * np.sin(r)
    return np.array([[Mxx, Mxy, Mxz], [Mxy, Myy, Myz], [Mxz, Myz, Mzz]])


def synth(M, src_ned, rec_ned, t, tau=1.0, vp=VP, vs=VS, rho=RHO, free_surface=2.0):
    """Displacement (3, nt) in N,E,Z-UP order at times t (s, origin at 0) for moment tensor M (N m) with a smooth
    step moment function M(t) = M * 0.5*(1+erf(t/tau)).  Returns N, E, Up."""
    from scipy.special import erf
    d = np.asarray(rec_ned, float) - np.asarray(src_ned, float)
    r = np.linalg.norm(d); g = d / r
    dt = t[1] - t[0]
    mf = lambda tt: 0.5 * (1 + erf(tt / tau))                      # noqa: E731  moment function (unit)
    mr = lambda tt: np.exp(-(tt / tau) ** 2) / (tau * np.sqrt(np.pi))  # noqa: E731  moment rate (unit)
    a, b = vp, vs
    # near-field time integral  I(t) = int_{r/a}^{r/b} tau' M(t-tau') dtau'
    taus = np.arange(r / a, r / b + dt, dt / 4)
    I = np.array([np.trapezoid(taus * mf(tt - taus), taus) for tt in t])
    u = np.zeros((3, t.size))
    gMg = g @ M @ g
    Mg = M @ g
    trM = np.trace(M)
    for n in range(3):
        AN = 15 * g[n] * gMg - 3 * g[n] * trM - 6 * Mg[n]
        AIP = 6 * g[n] * gMg - g[n] * trM - 2 * Mg[n]
        AIS = -(6 * g[n] * gMg - g[n] * trM - 3 * Mg[n])
        AFP = g[n] * gMg
        AFS = -(g[n] * gMg - Mg[n])
        u[n] = (AN * I / r ** 4
                + AIP * mf(t - r / a) / (a ** 2 * r ** 2)
                + AIS * mf(t - r / b) / (b ** 2 * r ** 2)
                + AFP * mr(t - r / a) / (a ** 3 * r)
                + AFS * mr(t - r / b) / (b ** 3 * r)) / (4 * np.pi * rho)
    u = u * free_surface
    return np.vstack([u[0], u[1], -u[2]])                          # N, E, Up


def bandpass(x, fs, f1=0.02, f2=0.05, order=4):
    sos = butter(order, [f1, f2], btype='band', fs=fs, output='sos')
    return sosfiltfilt(sos, x, axis=-1)


def ned_offset(lat0, lon0, lat, lon):
    """Flat-earth north/east offsets (m) of (lat,lon) from (lat0,lon0)."""
    R = 6371e3
    return (np.radians(lat - lat0) * R, np.radians(lon - lon0) * R * np.cos(np.radians(0.5 * (lat + lat0))))


if __name__ == '__main__':
    # order-of-magnitude table: peak 20-50 s displacement per 1e12 N m (Mw ~1.9) thrust at 35 km depth
    t = np.arange(-200, 400, 0.5)
    M = mt_ned(315, 12, 90, 1e12)
    for epi in (0, 20e3, 40e3, 60e3, 100e3, 150e3):
        for az in (45, 225):
            rec = (epi * np.cos(np.radians(az)), epi * np.sin(np.radians(az)), 0.0)
            u = bandpass(synth(M, (0, 0, 35e3), rec, t), 2.0)
            print(f'epi {epi/1e3:5.0f} km az {az:3d}: peak N/E/Up nm = ' +
                  ' '.join(f'{np.abs(c).max()*1e9:.3f}' for c in u))
    # explosion sanity check: far-field P outward, positive radial
    u = synth(np.eye(3) * 1e12, (0, 0, 0), (50e3, 0, 0), t, free_surface=1.0)
    print('explosion test, N displacement final static sign (+ expected):', np.sign(u[0][-1]))
