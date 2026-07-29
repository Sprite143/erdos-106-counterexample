"""Validate the inner solver against the Baek-Koizumi-Ueoro theorem.

BKU (2024) determine the axis-parallel optimum exactly:
    g(k^2)          = k
    g(k^2+2c+1)     = k + c/k      for -k < c < k
So for a whole range of n there is a PROVEN target to hit.  The solver is run
with t = 0 and must reproduce every one of them.  Matching ~25 independent
proven values is far better evidence for the inner solver than one MILP that
does not converge.
"""
import numpy as np, search, geom
from fractions import Fraction as Fr

exact = {}
for k in range(1, 7):
    exact.setdefault(k * k, Fr(k))
    for c in range(-k + 1, k):
        exact.setdefault(k * k + 2 * c + 1, Fr(k) + Fr(c, k))

rng = np.random.default_rng(31337)
bad = 0
print(f"{'n':>4} {'BKU exact':>12} {'solver':>14} {'diff':>11}   status")
for n in sorted(x for x in exact if x <= 26):
    tgt = exact[n]
    k = max(1, int(np.floor(np.sqrt(n))))
    v, cfg = search.evaluate(np.zeros(n), n, rng, restarts=120, k=k)
    d = v - float(tgt)
    if d > 1e-7:
        st = "*** EXCEEDS PROVEN VALUE -- BUG ***"; bad += 1
    elif d < -1e-7:
        st = "solver short (search, not correctness)"
    else:
        st = "match"
    print(f"{n:>4} {float(tgt):>12.6f} {v:>14.9f} {d:>+11.2e}   {st}")
print(f"\n{bad} values exceeded -- {'CORRECTNESS BUG' if bad else 'no upper-bound violations'}")
