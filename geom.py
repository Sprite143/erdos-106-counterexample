"""Exact geometry for Erdos #106: squares packed in the unit square.

A packing is n squares with pairwise disjoint interiors, each contained in
[0,1]^2.  Square i is (x_i, y_i, t_i, h_i): centre, angle, HALF side-length.
The quantity being maximised is  sum of side-lengths = 2 * sum h_i.

Two facts make everything cheap and exact:

  * Containment.  A square of half-side h tilted by t has axis-aligned
    half-width  q*h  in BOTH directions, where q = |cos t| + |sin t|.
    So containment is just  q*h <= min(x, 1-x, y, 1-y).

  * Disjointness.  Separating-axis theorem.  For two squares the only
    candidate axes are u_i, v_i, u_j, v_j (edge normals; each square
    contributes two).  Along any axis d the half-extent of square i is
    h_i*(|u_i.d| + |v_i.d|), which for these four axes collapses to either
    h or g*h with  g = |cos(t_j - t_i)| + |sin(t_j - t_i)|.

Sign convention throughout: a constraint is SATISFIED when the value is <= 0.
`pen[i,j] <= 0` means squares i and j have disjoint interiors.
"""

from __future__ import annotations

import numpy as np

AXES = 4  # u_i, v_i, u_j, v_j


def frames(t):
    """u = direction of one pair of edges, v = the perpendicular pair."""
    c, s = np.cos(t), np.sin(t)
    return np.stack([c, s], -1), np.stack([-s, c], -1)


def qfac(t):
    """Half-width of the axis-aligned bounding box, per unit half-side."""
    return np.abs(np.cos(t)) + np.abs(np.sin(t))


def gfac(t):
    """g[i,j]: extent of square j along an edge normal of i, per unit h_j."""
    dt = t[None, :] - t[:, None]
    return np.abs(np.cos(dt)) + np.abs(np.sin(dt))


def containment(x, y, t, h):
    """(n,4) array; every entry <= 0 iff all squares sit inside the unit square."""
    w = qfac(t) * h
    return np.stack([w - x, w - (1.0 - x), w - y, w - (1.0 - y)], -1)


def overlaps(x, y, t, h):
    """(4,n,n) separating-axis overlap values.

    ov[a,i,j] <= 0 means axis a separates squares i and j.  The pair is
    disjoint iff the minimum over a is <= 0.
    """
    u, v = frames(t)
    c = np.stack([x, y], -1)
    d = c[None, :, :] - c[:, None, :]          # d[i,j] = c_j - c_i
    g = gfac(t)
    hi, hj = h[:, None], h[None, :]

    du_i = np.einsum("ijk,ik->ij", d, u)
    dv_i = np.einsum("ijk,ik->ij", d, v)
    du_j = np.einsum("ijk,jk->ij", d, u)
    dv_j = np.einsum("ijk,jk->ij", d, v)

    return np.stack([
        hi + hj * g - np.abs(du_i),
        hi + hj * g - np.abs(dv_i),
        hj + hi * g - np.abs(du_j),
        hj + hi * g - np.abs(dv_j),
    ], 0)


def penetration(x, y, t, h):
    """(n,n) pairwise penetration; <= 0 on the off-diagonal iff disjoint."""
    return overlaps(x, y, t, h).min(axis=0)


def assignment(x, y, t, h):
    """Pick, for each pair, the axis that separates it best and that axis' sign.

    Returns (axis, sign) each (n,n).  Feeding these back into the LP keeps the
    current configuration feasible, which is what makes the outer loop
    monotone: the incumbent is always in the feasible set of the next LP.
    """
    ov = overlaps(x, y, t, h)
    axis = ov.argmin(axis=0)

    u, v = frames(t)
    c = np.stack([x, y], -1)
    d = c[None, :, :] - c[:, None, :]
    proj = np.stack([
        np.einsum("ijk,ik->ij", d, u),
        np.einsum("ijk,ik->ij", d, v),
        np.einsum("ijk,jk->ij", d, u),
        np.einsum("ijk,jk->ij", d, v),
    ], 0)
    chosen = np.take_along_axis(proj, axis[None], axis=0)[0]
    sign = np.where(chosen >= 0.0, 1.0, -1.0)
    return axis, sign


def violation(x, y, t, h):
    """Worst constraint violation of the whole configuration (<= 0 is feasible)."""
    n = len(x)
    pen = penetration(x, y, t, h)
    iu = np.triu_indices(n, 1)
    worst_pair = pen[iu].max() if n > 1 else -np.inf
    worst_box = containment(x, y, t, h).max()
    return max(worst_pair, worst_box)


def total(h):
    """Sum of side-lengths -- the objective of Erdos #106."""
    return 2.0 * float(np.sum(h))


def corners(x, y, t, h):
    """(n,4,2) polygon vertices, for plotting and for exact re-checking."""
    u, v = frames(t)
    hh = h[:, None]
    c = np.stack([x, y], -1)
    return np.stack([
        c + hh * (u + v), c + hh * (u - v),
        c - hh * (u + v), c - hh * (u - v),
    ], 1)
