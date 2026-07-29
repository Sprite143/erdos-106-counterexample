"""Independent verification, sharing no code with geom/certify.

The certifier and the search both rest on the separating-axis test and on the
q = |cos|+|sin| containment bound.  If either is wrong they are wrong together,
and a false counterexample sails through both.  So this file re-checks a
configuration from scratch by a completely different route:

  * corners computed directly from the rotation, then containment checked
    corner by corner against [0,1]^2;
  * pairwise overlap measured as the AREA of the polygon intersection,
    computed by Sutherland-Hodgman clipping and the shoelace formula --
    no separating axes anywhere;
  * a Monte-Carlo / raster double-check that no point is covered twice.

Agreement between two unrelated methods is the minimum bar for believing a
counterexample.
"""

from __future__ import annotations

import numpy as np


def square_corners(cx, cy, th, s):
    """Corners in CCW order, straight from the rotation matrix."""
    hh = s / 2.0
    loc = np.array([[-hh, -hh], [hh, -hh], [hh, hh], [-hh, hh]])
    c, si = np.cos(th), np.sin(th)
    R = np.array([[c, -si], [si, c]])
    return loc @ R.T + np.array([cx, cy])


def shoelace(P):
    if len(P) < 3:
        return 0.0
    x, y = P[:, 0], P[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def clip(subject, clipper):
    """Sutherland-Hodgman: intersect convex `subject` with convex `clipper`."""
    out = list(subject)
    n = len(clipper)
    for i in range(n):
        a, b = clipper[i], clipper[(i + 1) % n]
        edge = b - a
        if not out:
            return np.zeros((0, 2))
        inp, out = out, []
        for j in range(len(inp)):
            p, q = inp[j], inp[(j + 1) % len(inp)]
            sp = edge[0] * (p[1] - a[1]) - edge[1] * (p[0] - a[0])
            sq = edge[0] * (q[1] - a[1]) - edge[1] * (q[0] - a[0])
            if sp >= 0:
                out.append(p)
            if (sp > 0) != (sq > 0):
                t = sp / (sp - sq)
                out.append(p + t * (q - p))
    return np.array(out) if out else np.zeros((0, 2))


def verify(x, y, t, side, target, raster=4000, seed=0):
    n = len(x)
    polys = [square_corners(x[i], y[i], t[i], side[i]) for i in range(n)]

    # containment, corner by corner
    worst_out = 0.0
    for P in polys:
        worst_out = max(worst_out, float(np.max(np.maximum(-P, P - 1.0))))

    # pairwise intersection AREA
    worst_area = 0.0
    worst_pair = None
    for i in range(n):
        for j in range(i + 1, n):
            a = shoelace(clip(polys[i], polys[j]))
            if a > worst_area:
                worst_area, worst_pair = a, (i, j)

    # independent raster check: no pixel covered twice
    rng = np.random.default_rng(seed)
    pts = rng.random((raster, 2))
    cover = np.zeros(raster, dtype=int)
    for i in range(n):
        c, s_ = np.cos(-t[i]), np.sin(-t[i])
        dx, dy = pts[:, 0] - x[i], pts[:, 1] - y[i]
        lx = c * dx - s_ * dy
        ly = s_ * dx + c * dy
        cover += ((np.abs(lx) <= side[i] / 2) & (np.abs(ly) <= side[i] / 2)).astype(int)
    double = int((cover > 1).sum())

    total = float(np.sum(side))
    areas = float(np.sum(side ** 2))
    return {
        "total": total,
        "sum_of_areas": areas,
        "worst_containment_excess": worst_out,
        "worst_pair_overlap_area": worst_area,
        "worst_pair": worst_pair,
        "raster_double_covered": double,
        "raster_pts": raster,
        "beats": total > target,
    }


def report(x, y, t, side, target, label=""):
    r = verify(x, y, t, side, target)
    print(f"[independent]{' ' + label if label else ''}")
    print(f"  sum of sides             {r['total']:.15f}   target {target}")
    print(f"  sum of areas             {r['sum_of_areas']:.15f}   (must be <= 1)")
    print(f"  worst containment excess {r['worst_containment_excess']:+.3e}   (must be <= 0)")
    print(f"  worst pair overlap AREA  {r['worst_pair_overlap_area']:.3e}  pair {r['worst_pair']}")
    print(f"  raster double-covered    {r['raster_double_covered']} / {r['raster_pts']}")
    ok = (r["worst_containment_excess"] <= 1e-12
          and r["worst_pair_overlap_area"] <= 1e-14
          and r["raster_double_covered"] == 0)
    print(f"  VALID PACKING?           {ok}")
    print(f"  BEATS TARGET?            {r['beats']}")
    return r, ok
