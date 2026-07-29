#!/usr/bin/env python3
"""
WITNESS certificate: disjointness proved without any clipping algorithm.

The standing objection to a tight packing is that some particular geometry
routine might miss a microscopic overlap.  This file removes the routine from
the trust chain entirely.

For every pair of squares it exhibits an explicit SEPARATING AXIS together with
an exact rational gap.  Checking one witness needs no algorithm at all -- just
four multiplications and a comparison:

    for a unit axis d,  square i projects to an interval of half-width E_i,
    square j to one of half-width E_j, and their centres are  D = c_j - c_i
    apart.  If

        |D . d|  >=  E_i + E_j

    then the projections do not overlap in their interiors, so the squares
    cannot share an interior point.  The gap is |D . d| - (E_i + E_j) >= 0.

Because every square here is either axis-parallel or turned by the exact
rational angle (cos t, sin t) = (255/257, 32/257), each d is a rational vector
after clearing denominators, and every quantity below is an exact rational.
Nothing is rounded, so "microscopically overlapping" is not a state this
configuration can be in: a gap is either negative or it is not.

A gap of exactly 0 means the two squares TOUCH, which the problem permits --
it forbids only common INTERIOR points.
"""

from fractions import Fraction as F

from construct import build, COS, SIN


def axes_for(t_i, t_j):
    """The four candidate separating directions, as exact rational vectors.

    For a square the edge normals are its own two axes.  Scaling a direction by
    a positive constant changes nothing in the test, so the vectors are kept
    unnormalised and the extents are scaled to match.
    """
    def frame(t):
        return ((F(1), F(0)), (F(0), F(1))) if not t else ((COS, SIN), (-SIN, COS))
    ui, vi = frame(t_i)
    uj, vj = frame(t_j)
    return [ui, vi, uj, vj]


def extent(t_sq, half, d):
    """Half-width of a square's projection onto direction d (|d| = 1)."""
    if not t_sq:
        return half * (abs(d[0]) + abs(d[1]))
    u = (COS * d[0] + SIN * d[1])
    v = (-SIN * d[0] + COS * d[1])
    return half * (abs(u) + abs(v))


def separation_witness(A, B):
    """Best separating axis for two squares, exactly.  Returns (axis, gap)."""
    (ax, ay, sa, ta), (bx, by, sb, tb) = A, B
    D = (bx - ax, by - ay)
    best = None
    for k, d in enumerate(axes_for(ta, tb)):
        proj = abs(D[0] * d[0] + D[1] * d[1])
        gap = proj - (extent(ta, sa / 2, d) + extent(tb, sb / 2, d))
        if best is None or gap > best[1]:
            best = (k, gap)
    return best


def main():
    sq, (a, u, b, s) = build()
    n = len(sq)
    names = ["u_i", "v_i", "u_j", "v_j"]

    strict = touching = overlapping = 0
    worst = None
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            k, gap = separation_witness(sq[i], sq[j])
            if gap < 0:
                overlapping += 1
            elif gap == 0:
                touching += 1
            else:
                strict += 1
            if worst is None or gap < worst[0]:
                worst = (gap, i, j, k)
            rows.append((i, j, names[k], gap))

    print("PAIRWISE DISJOINTNESS  (separating-axis witnesses, exact rationals)")
    print(f"  pairs examined          : {len(rows)}")
    print(f"  strictly separated      : {strict}")
    print(f"  touching (gap exactly 0): {touching}   <- legal, only interiors must be disjoint")
    print(f"  OVERLAPPING (gap < 0)   : {overlapping}")
    print(f"  smallest gap            : {worst[0]}  (pair {worst[1]},{worst[2]} on axis {names[worst[3]]})")

    # containment: exact clearance of every square to every wall
    worst_wall = None
    out = 0
    for idx, (cx, cy, side, t) in enumerate(sq):
        half = side / 2
        w = half * ((COS + SIN) if t else F(1))     # half-width of the bounding box
        for name, clear in (("x-", cx - w), ("x+", 1 - cx - w),
                            ("y-", cy - w), ("y+", 1 - cy - w)):
            if clear < 0:
                out += 1
            if worst_wall is None or clear < worst_wall[0]:
                worst_wall = (clear, idx, name)
    print("\nCONTAINMENT  (exact clearance to each wall)")
    print(f"  walls violated          : {out}")
    print(f"  smallest clearance      : {worst_wall[0]}  (square {worst_wall[1]}, wall {worst_wall[2]})")

    total = sum(x[2] for x in sq)
    print(f"\nTOTAL side-length         : {total} = {float(total):.18f}")
    print(f"  minus 4                 : {total - 4} = {float(total - 4):.6e}")

    ok = overlapping == 0 and out == 0 and total > 4
    print(f"\nVERDICT: {'valid packing beating 4' if ok else 'FAILED'}")

    print("\n--- the ten tightest pairs, for hand-checking ---")
    rows.sort(key=lambda r: r[3])
    for i, j, ax, gap in rows[:10]:
        print(f"  squares {i:2d},{j:2d}  axis {ax}  gap {gap}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
