"""Is the axis-parallel optimum a strict local max once tilting is allowed?

BKU proved the axis-parallel maximum g(k^2+1)=k is exactly k, and the LP
reproduces it.  Define

    F(t) = max over packings with angle vector t of sum of side-lengths

which the LP computes exactly (given enough restarts to find the right
assignment basin).  F(0) = k.  The conjecture says max_t F(t) = k.

This scans F along tilt directions leaving t=0: single squares tilted, pairs,
and random directions, across the full range up to 45 degrees.  If F ever
exceeds k the conjecture is false; if F drops monotonically in every direction
that is real evidence for it.
"""

from __future__ import annotations

import sys

import numpy as np

import geom
import lp
import refine
import search


def F(t, n, rng, restarts=40, k=None):
    """Best value achievable at this exact angle vector."""
    best = -np.inf
    bcfg = None
    for r in range(restarts):
        if k is not None and r % 2 == 0:
            x, y, h = search.grid_config(n, rng, k)
        else:
            x, y, h = search.seed_config(n, rng)
        xx, yy, hh, v = lp.polish(x, y, t, h)
        if v > best:
            best, bcfg = v, (xx, yy, t.copy(), hh)
    return best, bcfg


def scan(n, k, rng, restarts=40):
    degs = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]
    target = float(k)
    print(f"n={n}  k={k}   axis-parallel optimum = {target}   "
          f"(Cauchy-Schwarz ceiling {np.sqrt(n):.6f})")

    print("\n-- tilt ONE square --")
    row = []
    for d in degs:
        t = np.zeros(n)
        t[0] = np.radians(d)
        v, _ = F(t, n, rng, restarts, k)
        row.append(v)
        print(f"   {d:5.1f} deg   F = {v:.9f}   {v - target:+.2e}")

    print("\n-- tilt TWO squares together --")
    for d in degs[1:]:
        t = np.zeros(n)
        t[0] = t[1] = np.radians(d)
        v, _ = F(t, n, rng, restarts, k)
        print(f"   {d:5.1f} deg   F = {v:.9f}   {v - target:+.2e}")

    print("\n-- tilt ALL squares together --")
    for d in degs[1:]:
        t = np.full(n, np.radians(d))
        v, _ = F(t, n, rng, restarts, k)
        print(f"   {d:5.1f} deg   F = {v:.9f}   {v - target:+.2e}")

    print("\n-- random tilt directions, small magnitude --")
    worst = -np.inf
    for trial in range(12):
        d = rng.choice([1.0, 2.0, 5.0, 10.0])
        mask = rng.random(n) < 0.4
        t = np.zeros(n)
        t[mask] = np.radians(d) * rng.random(mask.sum())
        v, _ = F(t, n, rng, restarts, k)
        worst = max(worst, v)
        print(f"   trial {trial:2d}  |tilted|={mask.sum()}  max {d:4.1f} deg   "
              f"F = {v:.9f}   {v - target:+.2e}")
    print(f"\n   best over random directions: {worst:.9f}  ({worst - target:+.2e})")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    k = int(np.floor(np.sqrt(n)))
    rng = np.random.default_rng(20260729)
    scan(n, k, rng, restarts=int(sys.argv[2]) if len(sys.argv) > 2 else 40)
