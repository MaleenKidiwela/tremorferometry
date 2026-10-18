"""Georeference a map panel from its frame tick marks.

Tick marks sit where graticule lines cross the axes frame. A lon tick on a
horizontal edge constrains: the meridian lon=L passes through (x_tick, y_edge).
A lat tick on a vertical edge: the parallel lat=P passes through (x_edge, y_tick).
We fit   page = s * R(theta) * proj(lon,lat) + t   (y flipped)
for several candidate projections and report per-tick residuals (pt and km).
"""
import numpy as np
from pyproj import Proj
from scipy.optimize import least_squares, brentq


def make_proj(kind, lon0, lat0, lat1=None, lat2=None):
    if kind == "lcc":
        return Proj(proj="lcc", lon_0=lon0, lat_0=lat0, lat_1=lat1, lat_2=lat2, ellps="WGS84")
    if kind == "lcc_sphere":
        return Proj(proj="lcc", lon_0=lon0, lat_0=lat0, lat_1=lat1, lat_2=lat2, R=6371000)
    if kind == "tmerc":
        return Proj(proj="tmerc", lon_0=lon0, lat_0=lat0, ellps="WGS84")
    if kind == "eqdc":
        return Proj(proj="eqdc", lon_0=lon0, lat_0=lat0, lat_1=lat1, lat_2=lat2, ellps="WGS84")
    if kind == "aea":
        return Proj(proj="aea", lon_0=lon0, lat_0=lat0, lat_1=lat1, lat_2=lat2, ellps="WGS84")
    if kind == "stere":
        return Proj(proj="stere", lon_0=lon0, lat_0=lat0, ellps="WGS84")
    if kind == "aeqd":
        return Proj(proj="aeqd", lon_0=lon0, lat_0=lat0, ellps="WGS84")
    raise ValueError(kind)


class PanelGeoref:
    """kind: projection family. lonlim/latlim: map limits (m_map convention:
    lon0 = mid lon, standard parallels = lat limits) used as starting values;
    lon0/lat1/lat2 are then freed in the fit."""

    def __init__(self, kind, lonlim, latlim, ticks, free_std=True, free_rot=False):
        self.kind = kind
        self.lonlim, self.latlim = lonlim, latlim
        self.ticks = ticks  # list of (type, value, edge_coord, tick_coord)
        self.free_std = free_std
        self.free_rot = free_rot

    def _proj(self, p):
        lon0 = p[0]
        lat1, lat2 = p[1], p[2]
        return make_proj(self.kind, lon0, (lat1 + lat2) / 2, lat1, lat2)

    def fwd(self, p, lon, lat, proj=None):
        proj = proj or self._proj(p)
        X, Y = proj(np.asarray(lon, float), np.asarray(lat, float))
        s, tx, ty, th = p[3], p[4], p[5], p[6]
        c, sn = np.cos(th), np.sin(th)
        x = s * (c * X - sn * Y) / 1000.0 + tx
        y = -s * (sn * X + c * Y) / 1000.0 + ty
        return x, y

    def inv(self, p, x, y):
        proj = self._proj(p)
        s, tx, ty, th = p[3], p[4], p[5], p[6]
        c, sn = np.cos(th), np.sin(th)
        a = (np.asarray(x) - tx) * 1000.0 / s
        b = -(np.asarray(y) - ty) * 1000.0 / s
        X = c * a + sn * b
        Y = -sn * a + c * b
        lon, lat = proj(X, Y, inverse=True)
        return lon, lat

    def residuals(self, p, full=False):
        proj = self._proj(p)
        res = []
        for typ, val, edge, tc in self.ticks:
            if typ == "lon":  # horizontal edge at y=edge, find lat
                f = lambda la: self.fwd(p, val, la, proj)[1] - edge
                try:
                    la = brentq(f, 30, 70)
                except ValueError:
                    res.append(50.0)
                    continue
                res.append(self.fwd(p, val, la, proj)[0] - tc)
            else:  # lat tick on vertical edge x=edge, find lon
                f = lambda lo: self.fwd(p, lo, val, proj)[0] - edge
                try:
                    lo = brentq(f, -140, -110)
                except ValueError:
                    res.append(50.0)
                    continue
                res.append(self.fwd(p, lo, val, proj)[1] - tc)
        res = np.array(res)
        # priors that pin non-free params
        pri = []
        if not self.free_std:
            pri += [(p[1] - self.latlim[0]) * 100, (p[2] - self.latlim[1]) * 100]
        if not self.free_rot:
            pri += [p[6] * 1e4]
        return np.concatenate([res, pri]) if pri else res

    def fit(self, box):
        lon0 = np.mean(self.lonlim)
        lat1, lat2 = self.latlim
        # initial scale: box width / projected width at mid lat
        pr = make_proj(self.kind, lon0, np.mean(self.latlim), lat1, lat2)
        Xa, Ya = pr(self.lonlim[0], np.mean(self.latlim))
        Xb, Yb = pr(self.lonlim[1], np.mean(self.latlim))
        s = (box[2] - box[0]) / ((Xb - Xa) / 1000.0)
        cx, cy = (box[0] + box[2]) / 2, (box[1] + box[3]) / 2
        Xc, Yc = pr(lon0, np.mean(self.latlim))
        p0 = np.array([lon0, lat1, lat2, s, cx - s * Xc / 1000, cy + s * Yc / 1000, 0.0])
        lb = [lon0 - 3, 35, 35, s * 0.5, -1e4, -1e4, -0.5]
        ub = [lon0 + 3, 60, 60, s * 2.0, 1e4, 1e4, 0.5]
        sol = least_squares(self.residuals, p0, bounds=(lb, ub),
                            x_scale=[0.5, 0.5, 0.5, s * 0.01, 1, 1, 0.01])
        self.p = sol.x
        r = self.residuals(self.p)[: len(self.ticks)]
        self.res_pt = r
        # km per pt at panel centre
        lo1, la1 = self.inv(self.p, cx, cy)
        lo2, la2 = self.inv(self.p, cx, cy + 1.0)
        self.km_per_pt = abs(la2 - la1) * 111.19
        return self
