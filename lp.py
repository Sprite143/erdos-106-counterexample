"""The linear program at the heart of the search.

Freeze the angles t, and freeze which separating axis (and which sign) handles
each pair.  Then

    maximise   sum_i h_i
    s.t.       q_i h_i  <=  x_i,  1 - x_i,  y_i,  1 - y_i          (containment)
               h_A + g_ij h_B  <=  sigma * (c_j - c_i) . d          (separation)
               h_i >= 0

is a linear program in (x, y, h).  Every feasible point of it is a genuine
packing -- the separation rows are a sufficient condition, never a relaxation --
so the LP optimum is always a real, achievable configuration.

The outer loop is:  solve LP -> re-derive the assignment from the solution ->
solve again.  The re-derived assignment is always satisfied by the solution we
just got, so the incumbent stays feasible and the objective can only go up.
That makes `polish` monotone and convergent.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix

import geom


def _rows(t, axis, sign):
    """Build the sparse constraint matrix for one frozen combinatorial class."""
    n = len(t)
    q = geom.qfac(t)
    g = geom.gfac(t)
    u, v = geom.frames(t)

    data, rows, cols, rhs = [], [], [], []
    r = 0

    def add(col, val):
        rows.append(r)
        cols.append(col)
        data.append(val)

    X = lambda i: i
    Y = lambda i: n + i
    H = lambda i: 2 * n + i

    # containment: q h - x <= 0, q h + x <= 1, likewise in y
    for i in range(n):
        for col, s, b in ((X(i), -1.0, 0.0), (X(i), 1.0, 1.0),
                          (Y(i), -1.0, 0.0), (Y(i), 1.0, 1.0)):
            add(col, s)
            add(H(i), q[i])
            rhs.append(b)
            r += 1

    # separation, one row per unordered pair
    for i in range(n):
        for j in range(i + 1, n):
            a = int(axis[i, j])
            sg = float(sign[i, j])
            d = (u[i] if a == 0 else v[i] if a == 1 else
                 u[j] if a == 2 else v[j])
            # h coefficients: axes 0,1 measure i at h_i and j at g*h_j
            if a in (0, 1):
                ci, cj = 1.0, g[i, j]
            else:
                ci, cj = g[i, j], 1.0
            add(H(i), ci)
            add(H(j), cj)
            # -sigma * ((x_j-x_i) dx + (y_j-y_i) dy)
            add(X(i), sg * d[0])
            add(X(j), -sg * d[0])
            add(Y(i), sg * d[1])
            add(Y(j), -sg * d[1])
            rhs.append(0.0)
            r += 1

    A = csr_matrix((data, (rows, cols)), shape=(r, 3 * n))
    return A, np.asarray(rhs)


def solve(t, axis, sign, hmax=0.5):
    """One LP.  Returns (x, y, h) or None if infeasible."""
    n = len(t)
    A, b = _rows(t, axis, sign)
    c = np.concatenate([np.zeros(2 * n), -np.ones(n)])   # maximise sum h
    bounds = [(0.0, 1.0)] * (2 * n) + [(0.0, hmax)] * n
    res = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method="highs")
    if not res.success:
        return None
    z = res.x
    return z[:n], z[n:2 * n], z[2 * n:]


def polish(x, y, t, h, iters=60, tol=1e-12):
    """Alternate LP / reassignment until the objective stops improving."""
    best = (x.copy(), y.copy(), h.copy())
    best_val = geom.total(h) if geom.violation(x, y, t, h) <= 1e-9 else -np.inf

    for _ in range(iters):
        axis, sign = geom.assignment(x, y, t, h)
        out = solve(t, axis, sign)
        if out is None:
            break
        x, y, h = out
        val = geom.total(h)
        if val > best_val + tol:
            best_val, best = val, (x.copy(), y.copy(), h.copy())
        else:
            best_val = max(best_val, val)
            break
    x, y, h = best
    return x, y, h, best_val
