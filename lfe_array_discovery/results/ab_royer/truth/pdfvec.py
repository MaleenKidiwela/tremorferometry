"""Low-level helpers for reading MATLAB-exported vector figures with pymupdf.

MATLAB/Illustrator exported every stroke as a thin *filled* polygon, so frames,
ticks and symbol outlines are all type 'f' paths. We split each path into
closed sub-paths and classify them geometrically.
"""
import numpy as np
import pymupdf


def subpaths(g, tol=0.02):
    """Split a drawing's items into sub-paths (lists of vertices)."""
    out, cur, last = [], [], None
    for it in g["items"]:
        if it[0] == "re":
            r = it[1]
            if cur:
                out.append(cur)
            out.append([(r.x0, r.y0), (r.x1, r.y0), (r.x1, r.y1), (r.x0, r.y1)])
            cur, last = [], None
            continue
        pts = [(q.x, q.y) for q in it[1:] if hasattr(q, "x")]
        if not pts:
            continue
        if last is not None and (abs(pts[0][0] - last[0]) > tol or abs(pts[0][1] - last[1]) > tol):
            out.append(cur)
            cur = []
        if not cur:
            cur.append(pts[0])
        cur.extend(pts[1:] if it[0] == "l" else [pts[-1]])
        last = pts[-1]
    if cur:
        out.append(cur)
    return out


def rgb(g):
    f = g.get("fill")
    return tuple(round(x, 2) for x in f) if f else None


def panels(page, minsize=60):
    """Axes boxes = white 4-segment filled rectangles (the axes background)."""
    res = []
    for g in page.get_drawings():
        r = g["rect"]
        if len(g["items"]) == 4 and rgb(g) == (1.0, 1.0, 1.0) and r.width > minsize and r.height > minsize:
            res.append((r.x0, r.y0, r.x1, r.y1))
    return res


def ticks(page, box, tol=0.35, maxlen=3.0):
    """Tick marks = thin filled rectangles (~0.3 x ~1.8 pt) touching a frame edge.
    Returns dict edge -> sorted list of tick-centre coordinate along the edge."""
    x0, y0, x1, y1 = box
    out = {"top": [], "bottom": [], "left": [], "right": []}
    for g in page.get_drawings():
        r = g["rect"]
        if r.x1 < x0 - 1 or r.x0 > x1 + 1 or r.y1 < y0 - 1 or r.y0 > y1 + 1:
            continue
        if rgb(g) in (None, (1.0, 1.0, 1.0)):
            continue
        for sp in subpaths(g):
            a = np.array(sp)
            if len(a) < 4 or len(a) > 6:
                continue
            w = a[:, 0].max() - a[:, 0].min()
            h = a[:, 1].max() - a[:, 1].min()
            cx, cy = a[:, 0].mean(), a[:, 1].mean()
            if w < 0.45 and 1.0 < h < maxlen:  # vertical tick
                if abs(a[:, 1].min() - y0) < tol and x0 - tol <= cx <= x1 + tol:
                    out["top"].append(float((a[:, 0].min() + a[:, 0].max()) / 2))
                elif abs(a[:, 1].max() - y1) < tol and x0 - tol <= cx <= x1 + tol:
                    out["bottom"].append(float((a[:, 0].min() + a[:, 0].max()) / 2))
            elif h < 0.45 and 1.0 < w < maxlen:  # horizontal tick
                if abs(a[:, 0].min() - x0) < tol and y0 - tol <= cy <= y1 + tol:
                    out["left"].append(float((a[:, 1].min() + a[:, 1].max()) / 2))
                elif abs(a[:, 0].max() - x1) < tol and y0 - tol <= cy <= y1 + tol:
                    out["right"].append(float((a[:, 1].min() + a[:, 1].max()) / 2))
    for k in out:
        v = sorted(out[k])
        # de-duplicate
        dd = []
        for x in v:
            if not dd or abs(x - dd[-1]) > 0.3:
                dd.append(x)
        out[k] = dd
    return out


def filled_quads(page, colours, wmin=2.5, wmax=5.0):
    """Coloured filled 4-segment diamonds. Returns list of (colour, cx, cy, w, h, verts)."""
    res = []
    for g in page.get_drawings():
        c = rgb(g)
        if c not in colours:
            continue
        r = g["rect"]
        if not (wmin <= r.width <= wmax and wmin <= r.height <= wmax + 1.5):
            continue
        sp = subpaths(g)
        if len(sp) != 1:
            continue
        a = np.array(sp[0])
        if a.shape[0] not in (4, 5):
            continue
        res.append(dict(colour=c, cx=float((r.x0 + r.x1) / 2), cy=float((r.y0 + r.y1) / 2),
                        w=float(r.width), h=float(r.height), nv=a.shape[0],
                        vx=a[:, 0].tolist(), vy=a[:, 1].tolist()))
    return res
