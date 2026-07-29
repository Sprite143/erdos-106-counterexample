"""Rigorous verification of a candidate packing.

Floating point cannot settle this problem.  The target is to beat exactly 3,
and the search routinely returns things like 3.0000000004 whose worst
constraint violation is 6e-10 -- i.e. the squares OVERLAP, and the excess is
noise.  Anything reported as a counterexample has to survive exact arithmetic.

Everything here runs in interval arithmetic over exact Fractions.  cos and sin
are enclosed by their alternating Taylor series, whose truncation error is
bounded by the first omitted term, so every enclosure is provably correct.

The shrink lemma that makes certification easy:

    reducing EVERY half-side by delta reduces every constraint value by at
    least delta,

because containment has coefficient q = |cos|+|sin| >= 1 on h, and separation
has coefficient (1 + g) >= 2 on the pair.  So if the worst rigorous violation
is V > 0, shrinking all squares by V yields a provably valid packing whose
total is smaller by exactly 2*n*V.  A candidate is a genuine counterexample
only if it still beats the target AFTER that shrink.
"""

from __future__ import annotations

from fractions import Fraction as Fr

TERMS = 24


class Iv:
    """Interval with exact Fraction endpoints.  All ops are outward-rounded."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Fr(lo)
        self.hi = Fr(lo if hi is None else hi)
        if self.lo > self.hi:
            self.lo, self.hi = self.hi, self.lo

    def __add__(self, o):
        o = as_iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = as_iv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def __mul__(self, o):
        o = as_iv(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Iv(min(p), max(p))

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return as_iv(o) - self

    def abs(self):
        if self.lo >= 0:
            return Iv(self.lo, self.hi)
        if self.hi <= 0:
            return Iv(-self.hi, -self.lo)
        return Iv(Fr(0), max(-self.lo, self.hi))

    def __repr__(self):
        return f"[{float(self.lo):.17g},{float(self.hi):.17g}]"


def as_iv(v):
    return v if isinstance(v, Iv) else Iv(v)


def _series(x: Iv, even: bool) -> Iv:
    """cos (even=True) or sin (even=False) by alternating Taylor series.

    For |x| <= pi/2 the terms decrease from the start, so truncating after
    TERMS terms leaves an error bounded by the magnitude of the first omitted
    term -- which is what gets added as the enclosure width.
    """
    lo = Iv(Fr(0))
    term_pow = Iv(Fr(1)) if even else x
    fact = Fr(1)
    total = Iv(Fr(0))
    start = 0 if even else 1
    k = start
    sign = 1
    i = 0
    p = Iv(Fr(1))
    f = Fr(1)
    # build sum_{i<TERMS} (-1)^i x^(start+2i) / (start+2i)!
    for i in range(TERMS):
        e = start + 2 * i
        if i == 0:
            p = x if start == 1 else Iv(Fr(1))
            f = Fr(1)
        else:
            p = p * x * x
            f = f * Fr((e - 1) * e)
        term = p * Iv(Fr(1, 1) / f)
        total = total + (term if sign > 0 else -term)
        sign = -sign
    # remainder bound: first omitted term
    e = start + 2 * TERMS
    xm = max(abs(x.lo), abs(x.hi))
    fe = Fr(1)
    for m in range(2, e + 1):
        fe *= m
    rem = Fr(xm) ** e / fe
    return Iv(total.lo - rem, total.hi + rem)


def cos_iv(x):
    return _series(as_iv(x), True)


def sin_iv(x):
    return _series(as_iv(x), False)


def certify(x, y, t, h, target):
    """Rigorously verify a packing and return its certified total.

    Returns a dict with the exact certified sum of side-lengths after shrinking
    the configuration by whatever the rigorous worst violation turns out to be.
    """
    n = len(x)
    X = [Iv(Fr(float(v))) for v in x]
    Y = [Iv(Fr(float(v))) for v in y]
    T = [Iv(Fr(float(v))) for v in t]
    H = [Iv(Fr(float(v))) for v in h]

    C = [cos_iv(v) for v in T]
    S = [sin_iv(v) for v in T]
    Q = [C[i].abs() + S[i].abs() for i in range(n)]

    worst = Fr(-10)

    # A square of zero side has empty interior, so it is disjoint from
    # everything by definition and contributes nothing.  Drop such squares
    # rather than letting the separating-axis test flag them as overlaps.
    live = [i for i in range(n) if Fr(float(h[i])) > 0]

    # containment: q*h - dist <= 0 for the four walls
    for i in live:
        qh = Q[i] * H[i]
        for d in (X[i], Iv(Fr(1)) - X[i], Y[i], Iv(Fr(1)) - Y[i]):
            worst = max(worst, (qh - d).hi)

    # separation: need SOME axis with a provably negative overlap
    for ii, i in enumerate(live):
        for j in live[ii + 1:]:
            cd = C[j] * C[i] + S[j] * S[i]          # cos(t_j - t_i)
            sd = S[j] * C[i] - C[j] * S[i]          # sin(t_j - t_i)
            g = cd.abs() + sd.abs()
            dx, dy = X[j] - X[i], Y[j] - Y[i]
            proj = [
                dx * C[i] + dy * S[i],
                dx * (-S[i]) + dy * C[i],
                dx * C[j] + dy * S[j],
                dx * (-S[j]) + dy * C[j],
            ]
            ext = [H[i] + g * H[j], H[i] + g * H[j],
                   H[j] + g * H[i], H[j] + g * H[i]]
            best = min((ext[a] - proj[a].abs()).hi for a in range(4))
            worst = max(worst, best)

    delta = max(worst, Fr(0))
    shrunk = [max(Fr(float(hi)) - delta, Fr(0)) for hi in h]
    certified = 2 * sum(shrunk)

    return {
        "worst_violation": worst,
        "delta": delta,
        "raw_total": 2 * sum(Fr(float(v)) for v in h),
        "certified_total": certified,
        "target": Fr(target),
        "beats_target": certified > Fr(target),
    }


def report(x, y, t, h, target, label=""):
    r = certify(x, y, t, h, target)
    print(f"[certify]{' ' + label if label else ''}")
    print(f"  raw total        {float(r['raw_total']):.15f}")
    print(f"  worst violation  {float(r['worst_violation']):+.3e}  (rigorous)")
    print(f"  shrink delta     {float(r['delta']):.3e}")
    print(f"  CERTIFIED total  {float(r['certified_total']):.15f}")
    print(f"  target           {float(r['target'])}")
    print(f"  VERDICT          {'COUNTEREXAMPLE' if r['beats_target'] else 'does not beat target'}")
    return r
