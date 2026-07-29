"""Outer search.

The LP already solves centres and sizes exactly once the angles are fixed, so
the only thing left to search over is the angle vector t (n numbers, period
pi/2) plus which basin of the assignment combinatorics we land in.  That is a
low-dimensional continuous search, which is why this problem is tractable at
all.

`evaluate` does the inner work for one angle vector: several random position
seeds, each polished by the monotone LP loop, best kept.
`climb` is a (1+lambda) hill-climber over t with an exploration level that
rises while stuck -- same idea as the semicircle engine.
"""

from __future__ import annotations

import numpy as np

import geom
import lp

HALFPI = np.pi / 2.0


def seed_config(n, rng, spread=1.0):
    x = rng.uniform(0.05, 0.95, n)
    y = rng.uniform(0.05, 0.95, n)
    h = np.full(n, 0.01 * spread)
    return x, y, h


def grid_config(n, rng, k):
    """Positions on a k x k lattice plus extras jittered in -- the natural basin."""
    pts = [((i + 0.5) / k, (j + 0.5) / k) for i in range(k) for j in range(k)]
    while len(pts) < n:
        pts.append((rng.random(), rng.random()))
    pts = np.array(pts[:n]) + rng.normal(0, 0.02, (n, 2))
    return np.clip(pts[:, 0], 0.02, 0.98), np.clip(pts[:, 1], 0.02, 0.98), np.full(n, 0.01)


def evaluate(t, n, rng, restarts=6, k=None):
    """Best packing found for this angle vector."""
    best = (-np.inf, None)
    for r in range(restarts):
        if k is not None and r < max(1, restarts // 3):
            x, y, h = grid_config(n, rng, k)
        else:
            x, y, h = seed_config(n, rng)
        x, y, h, val = lp.polish(x, y, t, h)
        if val > best[0]:
            best = (val, (x, y, t.copy(), h))
    return best


def climb(n, rng, steps=400, lam=6, restarts=4, k=None, t0=None, log=None):
    """(1+lambda) hill-climb on the angle vector."""
    if t0 is None:
        t0 = rng.uniform(0, HALFPI, n)
    t = t0.copy()
    val, cfg = evaluate(t, n, rng, restarts, k)
    best_val, best_cfg, best_t = val, cfg, t.copy()

    explore = 0.15
    for step in range(steps):
        improved = False
        for _ in range(lam):
            cand = t.copy()
            # mix of one-angle kicks and whole-vector jitter
            if rng.random() < 0.6:
                idx = rng.integers(0, n, size=max(1, int(rng.integers(1, 3))))
                cand[idx] += rng.normal(0, explore, len(idx))
            else:
                cand += rng.normal(0, explore * 0.35, n)
            if rng.random() < 0.12:                      # occasional snap to 0
                cand[rng.integers(0, n)] = 0.0
            cand = np.mod(cand, HALFPI)
            cval, ccfg = evaluate(cand, n, rng, restarts, k)
            if cval > val + 1e-12:
                t, val, cfg = cand, cval, ccfg
                improved = True
        if val > best_val + 1e-12:
            best_val, best_cfg, best_t = val, cfg, t.copy()
        explore = max(0.02, explore * 0.93) if improved else min(0.9, explore * 1.12)
        if log and step % max(1, steps // 20) == 0:
            log(f"  step {step:4d}  val {val:.9f}  best {best_val:.9f}  explore {explore:.3f}")
    return best_val, best_cfg, best_t
