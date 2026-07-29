#!/usr/bin/env python3
"""
CLOSED-FORM CONSTRUCTION of the Erdos #106 counterexample at n = 17.

The configuration is not a pile of 17 arbitrary numbers -- it is built from two
rationals (a, u) that solve a 2x2 linear system.  Everything else follows.

THE ANGLE.  Tilt the odd square by t with

    cos t = 255/257,   sin t = 32/257        (255^2 + 32^2 = 257^2, so exact)

Write  q = cos t + sin t = 287/257.

THE LAYOUT.  With a = width of the left two columns and u = height of the
middle band, and b = a - u:

  1. left column   : 4 squares of side 1/4       in x in [0, 1/4]
  2. second column : 3 squares of side a - 1/4   in x in [1/4, a], stacked from y=0
  3. bottom-right  : 4 squares of side (1-a)/2   filling x in [a,1], y in [0, 1-a]
  4. top-left      : 1 square  of side 3/4 - 2b  in x in [1/4, 1-2b], y in [1/4+2b, 1]
  5. top-right     : 2 squares of side b         in y in [1-b, 1]
  6. middle band   : 2 squares of side u         in y in [1-a, 1-b], x in [1-2u, 1]
  7. the odd one   : 1 square of side s = u/q, tilted by t, inscribed in the
                     box x in [1-3u, 1-2u], y in [1-a, 1-b]

Two squares poke into the odd square's box and are what make tilting pay:
  * square #6 (top of column 2) pokes up from the lower left, corner (a, 3a-3/4)
  * the top-left square pokes down from the upper left, corner (1-2b, 1/4+2b)

THE TWO EQUATIONS.  a and u are fixed by requiring the tilted square's two
left-hand edges to pass exactly through those two corners.  Both conditions are
LINEAR in (a, u):

    a(cos t + 4 sin t) + u cos t (3 - sin t / q)      =  cos t + (7/4) sin t
    a(-2 sin t - 3 cos t) + u(5 sin t + 2 cos t + cos^2 t / q) = -(3/4) cos t

THE SCORE.  Summing all seventeen side-lengths collapses to

    total = 3 + a + 2u + s          with s = u/q

which for this angle is the exact rational 2190452873/547596200, beating 4 by
exactly 68073/547596200 = 1.2431e-4.
"""

from fractions import Fraction as F

COS, SIN = F(255, 257), F(32, 257)
Q = COS + SIN
assert COS * COS + SIN * SIN == 1


def solve():
    """The 2x2 linear system for (a, u)."""
    A11 = COS + 4 * SIN
    A12 = COS * (3 - SIN / Q)
    B1 = COS + F(7, 4) * SIN
    A21 = -2 * SIN - 3 * COS
    A22 = 5 * SIN + 2 * COS + COS * COS / Q
    B2 = -F(3, 4) * COS
    det = A11 * A22 - A12 * A21
    a = (B1 * A22 - A12 * B2) / det
    u = (A11 * B2 - B1 * A21) / det
    return a, u


def build():
    """All 17 squares as (centre_x, centre_y, side, tilted) in exact rationals."""
    a, u = solve()
    b = a - u
    s = u / Q
    sq = []

    def add(cx, cy, side, tilted=False):
        sq.append((cx, cy, side, tilted))

    # 1. left column: four squares of side 1/4
    for k in range(4):
        add(F(1, 8), F(1, 8) + F(k, 4), F(1, 4))

    # 2. second column: three squares of side a - 1/4
    w = a - F(1, 4)
    for k in range(3):
        add(F(1, 4) + w / 2, w / 2 + k * w, w)

    # 3. bottom-right block: four squares of side (1-a)/2
    d = (1 - a) / 2
    for i in range(2):
        for j in range(2):
            add(a + d / 2 + i * d, d / 2 + j * d, d)

    # 4. top-left square of side 3/4 - 2b
    c = F(3, 4) - 2 * b
    add(F(1, 4) + c / 2, 1 - c / 2, c)

    # 5. top-right: two squares of side b
    add(1 - 2 * b + b / 2, 1 - b / 2, b)
    add(1 - b / 2, 1 - b / 2, b)

    # 6. middle band: two squares of side u
    add(1 - 2 * u + u / 2, 1 - a + u / 2, u)
    add(1 - u / 2, 1 - a + u / 2, u)

    # 7. the tilted square, inscribed in the box x:[1-3u,1-2u], y:[1-a,1-b]
    add(1 - F(5, 2) * u, 1 - a + u / 2, s, True)

    return sq, (a, u, b, s)


if __name__ == "__main__":
    import sys
    sq, (a, u, b, s) = build()

    print("closed form")
    print(f"  a = {a}")
    print(f"  u = {u}")
    print(f"  b = {b}")
    print(f"  s = {s}")
    total = sum(x[2] for x in sq)
    print(f"\n  total = {total}")
    print(f"        = {float(total):.18f}")
    print(f"  formula 3 + a + 2u + s = {3 + a + 2*u + s}  (matches: {total == 3 + a + 2*u + s})")
    print(f"  beats 4 by {total - 4} = {float(total - 4):.6e}\n")

    # verify with the independent checker (exact rational polygon clipping)
    sys.path.insert(0, ".")
    import verify_independent as V

    squares = [{"i": i, "cx": c[0], "cy": c[1], "side": c[2], "tilted": c[3]}
               for i, c in enumerate(sq)]
    good, tot, problems = V.check(squares, 4, label="closed form")
    print()
    print(f"  valid packing : {good}")
    print(f"  beats 4       : {tot > 4}")
    if problems:
        for p in problems[:6]:
            print("   !", p)
