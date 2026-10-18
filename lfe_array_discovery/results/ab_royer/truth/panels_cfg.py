"""Tick-value assignment for the four map panels (A overview, B NVI, C SVI, D NW+SW).

Values were read from the text tick labels of supplementary Fig. S1 (real PDF
text with positions). Main-paper Fig. 4 draws its labels as glyph outlines (no
text), but its panels have identical tick layouts (checked numerically in
digitize_royer2014.py: normalised tick positions agree to <0.3 % of panel size)
and identical labels in the raster.

Interior ticks only (corner ticks are drawn clipped and are not detected).
Longitudes negative (W). Lists are in increasing page coordinate order.
"""

M = 1 / 60.0

PANELS = {
    "A": dict(
        lonlim=(-130, -120), latlim=(46, 52),
        top=[-130, -128, -126, -124, -122, -120],
        bottom=[-128, -126, -124, -122],
        left=[51, 50, 49, 48, 47, 46],  # page y increasing = lat decreasing
        right=[51, 50, 49, 48, 47, 46],
    ),
    "B": dict(
        lonlim=(-127 - 20 * M - 0.0, -126 + 0.0), latlim=(49 + 40 * M, 50 + 35 * M),
        top=[-127 - 20 * M, -127, -126 - 40 * M, -126 - 20 * M, -126],
        bottom=[-127 - 20 * M, -127, -126 - 40 * M, -126 - 20 * M, -126],
        left=[50.5, 50 + 20 * M, 50 + 10 * M, 50.0, 49 + 50 * M],
        right=[50.5, 50 + 20 * M, 50 + 10 * M, 50.0, 49 + 50 * M],
    ),
    "C": dict(
        lonlim=(-126, -122), latlim=(47.9, 49.8),
        top=[-126, -125, -124, -123, -122],
        bottom=[-125, -124, -123],
        left=[49 + 40 * M, 49 + 20 * M, 49.0, 48 + 40 * M, 48 + 20 * M, 48.0],
        right=[49 + 40 * M, 49 + 20 * M, 49.0, 48 + 40 * M, 48 + 20 * M, 48.0],
    ),
    "D": dict(
        lonlim=(-125, -122), latlim=(46.3, 48.4),
        top=[-125, -124.5, -124, -123.5, -123, -122.5, -122],
        bottom=[-124.5, -124, -123.5, -123, -122.5],
        left=[48.0, 47.5, 47.0, 46.5],
        right=[48.0, 47.5, 47.0, 46.5],
    ),
}


def build_ticks(cfg, box, tk):
    """Combine detected tick coords `tk` (edge -> list) with values in cfg."""
    x0, y0, x1, y1 = box
    out = []
    for edge, ecoord, typ in (("top", y0, "lon"), ("bottom", y1, "lon"),
                              ("left", x0, "lat"), ("right", x1, "lat")):
        vals, pos = cfg[edge], tk[edge]
        if len(vals) != len(pos):
            raise ValueError(f"{edge}: {len(vals)} values vs {len(pos)} ticks {pos}")
        for v, c in zip(vals, pos):
            out.append((typ, v, ecoord, c))
    return out
