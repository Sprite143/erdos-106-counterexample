#!/usr/bin/env python3
"""
INDEPENDENT VERIFIER for the claimed counterexample to Erdos problem #106.

This file is deliberately written from the PROBLEM STATEMENT, not from the
search code that produced the configuration.  It shares no code, and no method,
with the solver:

  * the solver decides overlap with the separating-axis theorem;
    this file computes the exact AREA of the polygon intersection instead.
  * the solver checks containment with the bounding-box identity
    (|cos|+|sin|)*h <= dist-to-wall;
    this file simply checks all four corners lie in the unit square.
  * the solver uses numpy floats plus interval arithmetic;
    this file uses only exact rational arithmetic (fractions.Fraction) and
    the Python standard library.  There is no floating point anywhere in the
    verification path.

Exactness is possible because the tilted square uses a RATIONAL rotation:
(255, 32, 257) is a Pythagorean triple, so

    cos t = 255/257 ,  sin t = 32/257 ,  cos^2 + sin^2 = 1  exactly

is an exact rotation of the plane (t = 7.1527 degrees).  Every corner
coordinate is therefore an exact rational, and every test below is decided
with no rounding of any kind.

WHAT IS BEING CHECKED
---------------------
Erdos #106 asks for f(n), the maximum of the sum of the side-lengths of n
squares drawn inside the unit square with no common interior point (i.e.
pairwise disjoint interiors -- squares are allowed to touch).  The conjecture
is f(k^2+1) = k.  For k = 4 that reads f(17) = 4.

A valid packing therefore requires exactly two things:
  (1) every square lies inside the unit square, and
  (2) no two squares share an interior point.
Both are checked below, exactly.  Then the sum of side-lengths is compared
with 4 as exact rationals.

Usage:  python verify_independent.py [CERTIFICATE_n17.json]
Exit code 0 means verified, 1 means the claim FAILED verification.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction as F

# ----------------------------------------------------------------- rotation

COS = F(255, 257)
SIN = F(32, 257)

assert COS * COS + SIN * SIN == 1, "rotation is not exact"
assert 255 ** 2 + 32 ** 2 == 257 ** 2, "(255,32,257) is not a Pythagorean triple"


# ------------------------------------------------------------- exact polygons

def corners(cx, cy, side, tilted):
    """The four corners of a square, in counter-clockwise order, exactly.

    Built straight from the definition: take the four corners of an
    axis-aligned square of the given side about the origin, rotate by t if the
    square is tilted, then translate to the centre.
    """
    h = F(side, 2)
    local = [(-h, -h), (h, -h), (h, h), (-h, h)]
    out = []
    for lx, ly in local:
        if tilted:
            rx = COS * lx - SIN * ly
            ry = SIN * lx + COS * ly
        else:
            rx, ry = lx, ly
        out.append((cx + rx, cy + ry))
    return out


def signed_area(poly):
    """Exact shoelace.  Positive for counter-clockwise."""
    if len(poly) < 3:
        return F(0)
    s = F(0)
    for i in range(len(poly)):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % len(poly)]
        s += x1 * y2 - x2 * y1
    return s / 2


def clip_halfplane(poly, a, b):
    """Sutherland-Hodgman: keep the part of `poly` left of the directed line a->b.

    All arithmetic is rational.  The intersection parameter is only ever formed
    when the two endpoints lie strictly on opposite sides, so the denominator
    can never be zero.
    """
    if not poly:
        return []
    (ax, ay), (bx, by) = a, b
    ex, ey = bx - ax, by - ay

    def side(p):
        return ex * (p[1] - ay) - ey * (p[0] - ax)

    out = []
    n = len(poly)
    for i in range(n):
        p, q = poly[i], poly[(i + 1) % n]
        sp, sq = side(p), side(q)
        if sp >= 0:
            out.append(p)
        if (sp > 0 and sq < 0) or (sp < 0 and sq > 0):
            t = sp / (sp - sq)
            out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return out


def intersection_area(P, Q):
    """Exact area of the intersection of two convex polygons (CCW)."""
    poly = P
    for i in range(len(Q)):
        poly = clip_halfplane(poly, Q[i], Q[(i + 1) % len(Q)])
        if not poly:
            return F(0)
    return abs(signed_area(poly))


# ------------------------------------------------------------------- checking

def check(squares, target, label="", verbose=True):
    """squares: list of dicts with exact Fraction cx, cy, side and bool tilted."""
    problems = []
    polys = []
    for s in squares:
        if s["side"] < 0:
            problems.append(f"square {s['i']}: negative side {s['side']}")
        P = corners(s["cx"], s["cy"], s["side"], s["tilted"])
        if signed_area(P) < 0:
            problems.append(f"square {s['i']}: corners are not counter-clockwise")
        polys.append(P)

    # (1) containment: every corner inside the unit square
    for s, P in zip(squares, polys):
        for (x, y) in P:
            if x < 0 or x > 1 or y < 0 or y > 1:
                problems.append(
                    f"square {s['i']}: corner ({x}, {y}) is outside the unit square")

    # (2) pairwise disjoint interiors: exact intersection area must be zero
    worst = F(0)
    for i in range(len(polys)):
        for j in range(i + 1, len(polys)):
            a = intersection_area(polys[i], polys[j])
            if a > worst:
                worst = a
            if a != 0:
                problems.append(
                    f"squares {squares[i]['i']} and {squares[j]['i']} overlap, "
                    f"exact area {a} = {float(a):.3e}")

    total = sum(s["side"] for s in squares)

    if verbose:
        tag = f" [{label}]" if label else ""
        print(f"  squares checked          : {len(squares)}{tag}")
        print(f"  corners outside unit sq. : "
              f"{sum(1 for p in problems if 'outside' in p)}")
        print(f"  overlapping pairs        : "
              f"{sum(1 for p in problems if 'overlap' in p)}")
        print(f"  largest overlap area     : {worst} (exactly)")
        print(f"  sum of side-lengths      : {total}")
        print(f"                           = {float(total):.18f}")
        print(f"  target                   : {target}")
        for p in problems[:8]:
            print(f"    ! {p}")
    return (not problems), total, problems


# ------------------------------------------------------------------ self-test

def selftest():
    """Prove the verifier itself works before trusting its verdict."""
    print("SELF-TEST (the verifier must pass these before its verdict means anything)")
    ok = True

    def mk(i, cx, cy, side, tilted=False):
        return {"i": i, "cx": F(cx), "cy": F(cy), "side": F(side), "tilted": tilted}

    # 1. the 4x4 grid: 16 squares of side 1/4, valid, sums to exactly 4
    grid = [mk(4 * a + b, F(2 * a + 1, 8), F(2 * b + 1, 8), F(1, 4))
            for a in range(4) for b in range(4)]
    good, total, _ = check(grid, 4, verbose=False)
    r1 = good and total == 4
    print(f"  [{'PASS' if r1 else 'FAIL'}] 4x4 grid accepted, sum == 4 exactly "
          f"(got {total}, valid={good})")
    ok &= r1

    # 2. two squares that merely TOUCH must be accepted (the problem allows it)
    touch = [mk(0, F(1, 4), F(1, 2), F(1, 2)), mk(1, F(3, 4), F(1, 2), F(1, 2))]
    good, _, _ = check(touch, 1, verbose=False)
    print(f"  [{'PASS' if good else 'FAIL'}] edge-to-edge touching squares accepted")
    ok &= good

    # 3. a genuine overlap must be caught
    over = [mk(0, F(3, 10), F(1, 2), F(3, 10)), mk(1, F(1, 2), F(1, 2), F(3, 10))]
    good, _, probs = check(over, 1, verbose=False)
    r3 = (not good) and any("overlap" in p for p in probs)
    print(f"  [{'PASS' if r3 else 'FAIL'}] overlapping squares rejected")
    ok &= r3

    # 4. a square poking outside must be caught
    out = [mk(0, F(19, 20), F(1, 2), F(2, 5))]
    good, _, probs = check(out, 1, verbose=False)
    r4 = (not good) and any("outside" in p for p in probs)
    print(f"  [{'PASS' if r4 else 'FAIL'}] out-of-bounds square rejected")
    ok &= r4

    # 5. a TILTED overlap must be caught (exercises the rotation path)
    tover = [mk(0, F(2, 5), F(1, 2), F(3, 10)),
             mk(1, F(3, 5), F(1, 2), F(3, 10), tilted=True)]
    good, _, probs = check(tover, 1, verbose=False)
    r5 = (not good) and any("overlap" in p for p in probs)
    print(f"  [{'PASS' if r5 else 'FAIL'}] overlapping TILTED squares rejected")
    ok &= r5

    print(f"  self-test: {'all passed' if ok else 'FAILED -- do not trust the verdict'}\n")
    return ok


# ----------------------------------------------------------------------- main

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "CERTIFICATE_n17.json"

    if not selftest():
        print("VERIFIER IS BROKEN")
        return 1

    with open(path) as fh:
        cert = json.load(fh)

    print(f"CERTIFICATE: {path}")
    print(f"  claim    : {cert.get('claim')}")
    print(f"  rotation : cos={cert['rotation']['cos']} sin={cert['rotation']['sin']}"
          f"  ({cert['rotation']['theta_deg']:.6f} deg)\n")

    # rebuild from centre/side/tilted -- stored corners are NOT trusted
    squares = [{"i": s["i"], "cx": F(s["cx"]), "cy": F(s["cy"]),
                "side": F(s["side"]), "tilted": bool(s["tilted"])}
               for s in cert["squares"]]

    n = len(squares)
    k = 4
    print(f"VERIFYING n={n} against the conjecture f(k^2+1)=k with k={k}:")
    good, total, problems = check(squares, k, label=f"n={n}")

    print()
    if not good:
        print("RESULT: NOT A VALID PACKING -- the claim is refuted.")
        return 1

    print(f"  valid packing            : YES (exact, no floating point)")
    print(f"  exceeds {k} by            : {total - k} = {float(total - k):.6e}")
    print()
    if total > k:
        print(f"RESULT: VERIFIED.  {n} squares pack the unit square with total")
        print(f"        side-length {total} > {k}.")
        print(f"        Therefore f({n}) > {k}, so f(k^2+1) = k is FALSE at k = {k}.")
        return 0
    print(f"RESULT: valid packing, but total {total} does not exceed {k}.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
