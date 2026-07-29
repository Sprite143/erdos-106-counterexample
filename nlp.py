"""Smooth local polish over ALL variables, angles included.

The LP freezes angles; this does not.  With t restricted to [0, pi/2] every
absolute value in the constraint system has a determined sign:

    q_i = |cos t_i| + |sin t_i|  =  cos t_i + sin t_i          (t_i in [0,pi/2])
    g_ij = |cos D| + |sin D|     =  cos D + eps*sin D,  D = t_j - t_i

with eps = sign(t_j - t_i) frozen from the incoming configuration.  Freeze the
separating axis and its sign too and the whole thing is a smooth NLP, so SLSQP
applies.  Every feasible point is still a genuine packing, for the same reason
as in the LP: the separation rows are sufficient conditions.

Used in alternation with lp.polish -- LP jumps between combinatorial basins,
this slides angles downhill inside one.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

import geom

HALFPI = np.pi / 2.0


def _pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def polish(x, y, t, h, maxiter=200):
    """Returns (x, y, t, h, value).  Falls back to the input if SLSQP fails."""
    n = len(x)
    t = np.mod(t, HALFPI)
    axis, sign = geom.assignment(x, y, t, h)
    pairs = _pairs(n)
    a = np.array([axis[i, j] for i, j in pairs])
    sg = np.array([sign[i, j] for i, j in pairs])
    eps = np.array([1.0 if t[j] >= t[i] else -1.0 for i, j in pairs])
    I = np.array([i for i, _ in pairs])
    J = np.array([j for _, j in pairs])

    z0 = np.concatenate([x, y, t, h])

    def split(z):
        return z[:n], z[n:2 * n], z[2 * n:3 * n], z[3 * n:]

    def cons(z):
        X, Y, T, H = split(z)
        ct, st = np.cos(T), np.sin(T)
        q = ct + st
        box = np.concatenate([q * H - X, q * H + X - 1.0,
                              q * H - Y, q * H + Y - 1.0])

        D = T[J] - T[I]
        g = np.cos(D) + eps * np.sin(D)
        dx = np.where(a <= 1, ct[I], ct[J])
        dy = np.where(a <= 1, st[I], st[J])
        # v = (-sin, cos) for the odd axes
        vx = np.where(a <= 1, -st[I], -st[J])
        vy = np.where(a <= 1, ct[I], ct[J])
        use_v = (a == 1) | (a == 3)
        DX = np.where(use_v, vx, dx)
        DY = np.where(use_v, vy, dy)

        ci = np.where(a <= 1, 1.0, g)
        cj = np.where(a <= 1, g, 1.0)
        proj = (X[J] - X[I]) * DX + (Y[J] - Y[I]) * DY
        sep = ci * H[I] + cj * H[J] - sg * proj
        return -np.concatenate([box, sep])          # SLSQP wants >= 0

    res = minimize(
        lambda z: -np.sum(z[3 * n:]),
        z0,
        jac=lambda z: np.concatenate([np.zeros(3 * n), -np.ones(n)]),
        bounds=[(0, 1)] * (2 * n) + [(0, HALFPI)] * n + [(0, 0.5)] * n,
        constraints=[{"type": "ineq", "fun": cons}],
        method="SLSQP",
        options={"maxiter": maxiter, "ftol": 1e-14},
    )

    X, Y, T, H = split(res.x)
    if geom.violation(X, Y, T, H) > 1e-9:
        return x, y, t, h, geom.total(h)
    val = geom.total(H)
    if val <= geom.total(h):
        return x, y, t, h, geom.total(h)
    return X, Y, T, H, val
