"""Alternate the two local moves until neither helps.

lp.polish   -- exact in (centres, sizes), can hop between assignment basins.
nlp.polish  -- slides angles downhill, stuck in one basin.

Alternating gets both.  Both are monotone in the objective and both only ever
return genuine packings, so `run` never reports a value it cannot realise.
"""

from __future__ import annotations

import numpy as np

import geom
import lp
import nlp


def run(x, y, t, h, rounds=12, tol=1e-13):
    best = geom.total(h)
    for _ in range(rounds):
        x, y, h, v1 = lp.polish(x, y, t, h)
        x, y, t, h, v2 = nlp.polish(x, y, t, h)
        if v2 <= best + tol:
            best = max(best, v2)
            break
        best = v2
    if geom.violation(x, y, t, h) > 1e-9:
        raise AssertionError("refine produced an infeasible packing")
    return x, y, t, h, geom.total(h)
