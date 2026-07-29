"""Parallel global hunt for a tilted packing beating k.

One island per core.  Each island does cheap wide sampling of angle vectors
(LP-only, few restarts) and promotes anything that looks alive to the full
refine + angle hill-climb.  Best-ever per island is checkpointed to JSON so a
kill never loses a find.

Angles are canonicalised by sorting: the objective is invariant under
permuting the squares, so sorted angle vectors cover the whole space with the
10!-fold redundancy removed.

Anything that beats the target is written to a separate `HIT_*.json` and
re-verified from scratch before it is believed.
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

import geom
import lp
import refine
import search

HALFPI = np.pi / 2.0


def sample_angles(n, rng):
    """Mixed prior: mostly-axis-parallel with a few tilts, through to fully random."""
    mode = rng.random()
    t = np.zeros(n)
    if mode < 0.30:                              # a few squares tilted
        m = int(rng.integers(1, max(2, n // 2)))
        idx = rng.choice(n, m, replace=False)
        t[idx] = rng.uniform(0, HALFPI, m)
    elif mode < 0.50:                            # small tilts everywhere
        t = rng.uniform(0, np.radians(rng.choice([2, 5, 10, 20])), n)
    elif mode < 0.65:                            # a shared tilt angle, common in packings
        a = rng.uniform(0, HALFPI)
        m = int(rng.integers(1, n + 1))
        t[rng.choice(n, m, replace=False)] = a
    elif mode < 0.80:                            # two-angle families
        a, b = rng.uniform(0, HALFPI, 2)
        pick = rng.random(n) < 0.5
        t = np.where(pick, a, b)
    else:                                        # fully random
        t = rng.uniform(0, HALFPI, n)
    return np.sort(t)


def island(args):
    wid, n, k, budget, seed, outdir = args
    rng = np.random.default_rng(seed)
    target = float(k)
    best = -np.inf
    best_cfg = None
    evals = 0
    promoted = 0
    t_end = time.time() + budget
    path = os.path.join(outdir, f"island_{wid}.json")

    while time.time() < t_end:
        t = sample_angles(n, rng)
        v, cfg = search.evaluate(t, n, rng, restarts=6, k=k)
        evals += 1

        # promote anything close enough to be worth the expensive treatment
        if v > target - 0.005:
            promoted += 1
            x, y, tt, h = cfg
            try:
                x, y, tt, h, v = refine.run(x, y, tt, h)
            except AssertionError:
                continue
            bv, bcfg, _ = search.climb(n, rng, steps=12, lam=4, restarts=4, k=k, t0=tt)
            if bcfg is not None and bv > v:
                x, y, tt, h = bcfg
                try:
                    x, y, tt, h, v = refine.run(x, y, tt, h)
                except AssertionError:
                    pass
            cfg = (x, y, tt, h)

        if v > best:
            best, best_cfg = v, cfg
            x, y, tt, h = cfg
            rec = {
                "worker": wid, "n": n, "k": k, "value": float(v),
                "excess": float(v - target),
                "violation": float(geom.violation(x, y, tt, h)),
                "area": float(np.sum((2 * h) ** 2)),
                "x": x.tolist(), "y": y.tolist(),
                "theta": tt.tolist(), "side": (2 * h).tolist(),
                "evals": evals, "promoted": promoted,
            }
            with open(path, "w") as f:
                json.dump(rec, f, indent=1)
            # A raw value above the target means nothing -- the search
            # routinely returns 3.0000000004 with a 6e-10 overlap.  Only an
            # exact-arithmetic certificate counts as a hit.
            if v > target + 1e-12:
                import certify
                c = certify.certify(x, y, tt, h, target)
                rec["certified_total"] = float(c["certified_total"])
                rec["rigorous_violation"] = float(c["worst_violation"])
                if c["beats_target"]:
                    rec["CERTIFIED"] = True
                    with open(os.path.join(outdir, f"HIT_{wid}_{v:.12f}.json"), "w") as f:
                        json.dump(rec, f, indent=1)

    return wid, best, evals, promoted


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 600.0
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 11
    k = int(np.floor(np.sqrt(n)))
    outdir = f"out_n{n}"
    os.makedirs(outdir, exist_ok=True)

    print(f"hunting n={n}  target > {k}   ceiling {np.sqrt(n):.6f}   "
          f"{workers} islands x {budget:.0f}s", flush=True)

    import multiprocessing as mp
    args = [(w, n, k, budget, 90210 + 7919 * w, outdir) for w in range(workers)]
    with mp.Pool(workers) as pool:
        results = pool.map(island, args)

    best = max(r[1] for r in results)
    tot = sum(r[2] for r in results)
    prom = sum(r[3] for r in results)
    print(f"\nislands done: {tot} angle vectors evaluated, {prom} promoted")
    for wid, b, e, p in results:
        print(f"  island {wid:2d}  best {b:.10f}  ({b - k:+.2e})  evals {e}")
    print(f"\nBEST n={n}: {best:.12f}   excess {best - k:+.3e}   "
          f"{'*** COUNTEREXAMPLE ***' if best > k + 1e-9 else 'no counterexample'}")


if __name__ == "__main__":
    main()
