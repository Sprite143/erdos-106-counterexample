"""Render the two-panel diagram: the whole packing, and a detail of the tilt."""
from fractions import Fraction as F
import construct

sq, (a, u, b, s) = construct.build()
C, S = float(construct.COS), float(construct.SIN)
A, U, B = float(a), float(u), float(b)
TOT = float(sum(q[2] for q in sq))
SA = (1 - 2 * U) - A                      # best axis-parallel square in the same hole


def corners(cx, cy, side, t):
    h = float(side) / 2
    cx, cy = float(cx), float(cy)
    out = []
    for lx, ly in ((-h, -h), (h, -h), (h, h), (-h, h)):
        rx, ry = (C * lx - S * ly, S * lx + C * ly) if t else (lx, ly)
        out.append((cx + rx, cy + ry))
    return out


polys = [(corners(*q), float(q[2]), q[3]) for q in sq]

W, H, PAD, P = 1580, 900, 54, 750
ox1, oy1 = PAD, PAD
ox2, oy2 = PAD + P + 130, PAD
zx0, zx1, zy0, zy1 = 0.395, 0.665, 0.535, 0.785
zs = P / max(zx1 - zx0, zy1 - zy0)
zw, zh = (zx1 - zx0) * zs, (zy1 - zy0) * zs

p1 = lambda x, y: (ox1 + x * P, oy1 + (1 - y) * P)
p2 = lambda x, y: (ox2 + (x - zx0) * zs, oy2 + (zy1 - y) * zs)
pts_str = lambda pts, f: " ".join("%.4f,%.4f" % f(x, y) for x, y in pts)

o = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
     'font-family="Georgia,serif"><rect width="100%%" height="100%%" fill="#fff"/>' % (W, H)]

# ---- panel A: the whole packing
for pts, side, t in polys:
    o.append('<polygon points="%s" fill="%s" stroke="%s" stroke-width="2"/>'
             % (pts_str(pts, p1), "#c0392b33" if t else "#5b7fa622",
                "#b8391f" if t else "#43607f"))
    cx = sum(q[0] for q in pts) / 4
    cy = sum(q[1] for q in pts) / 4
    tx, ty = p1(cx, cy)
    o.append('<text x="%.1f" y="%.1f" font-size="15" fill="#33475c" text-anchor="middle" '
             'font-family="Consolas,monospace">%.4f</text>' % (tx, ty + 5, side))
o.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#111" stroke-width="3"/>'
         % (ox1, oy1, P, P))
o.append('<text x="%d" y="%d" font-size="23" fill="#111">17 squares in the unit square, '
         'one tilted 7.1527&#176;</text>' % (ox1, oy1 + P + 36))
o.append('<text x="%d" y="%d" font-size="19" fill="#b8391f" font-family="Consolas,monospace">'
         'total = 2190452873/547596200 = %.15f &gt; 4</text>' % (ox1, oy1 + P + 64, TOT))

# ---- panel B: the detail
o.append('<clipPath id="cp"><rect x="%d" y="%d" width="%.1f" height="%.1f"/></clipPath>'
         % (ox2, oy2, zw, zh))
o.append('<g clip-path="url(#cp)">')
for pts, side, t in polys:
    o.append('<polygon points="%s" fill="%s" stroke="%s" stroke-width="2.5"/>'
             % (pts_str(pts, p2), "#c0392b33" if t else "#5b7fa618",
                "#b8391f" if t else "#8ba3bb"))
alt = [(A, 1 - A), (A + SA, 1 - A), (A + SA, 1 - A + SA), (A, 1 - A + SA)]
o.append('<polygon points="%s" fill="none" stroke="#2c6e49" stroke-width="2.5" '
         'stroke-dasharray="8 6"/>' % pts_str(alt, p2))
x1, y1 = p2(A, zy0)
x2, y2 = p2(A, zy1)
o.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#b8391f" stroke-width="2" '
         'stroke-dasharray="6 5"/>' % (x1, y1, x2, y2))
for px, py in ((A, 3 * A - 0.75), (1 - 2 * B, 0.25 + 2 * B)):
    cx, cy = p2(px, py)
    o.append('<circle cx="%.1f" cy="%.1f" r="7" fill="#b8391f"/>' % (cx, cy))
o.append('</g>')
o.append('<rect x="%d" y="%d" width="%.1f" height="%.1f" fill="none" stroke="#111" stroke-width="3"/>'
         % (ox2, oy2, zw, zh))
o.append('<text x="%d" y="%d" font-size="23" fill="#111">detail: why the tilt gains</text>'
         % (ox2, oy2 + zh + 36))
o.append('<text x="%d" y="%d" font-size="17" fill="#2c6e49" font-family="Consolas,monospace">'
         'dashed green: largest axis-parallel square here, side %.9f</text>'
         % (ox2, oy2 + zh + 64, SA))
o.append('<text x="%d" y="%d" font-size="17" fill="#b8391f" font-family="Consolas,monospace">'
         'red: tilted, side %.9f &#8212; reaches left past the line</text>'
         % (ox2, oy2 + zh + 88, float(s)))
o.append('<text x="%d" y="%d" font-size="15" fill="#555" font-family="Consolas,monospace">'
         'dots: the two neighbour corners its edges pass through</text>' % (ox2, oy2 + zh + 112))
o.append('</svg>')

open("../erdos106-repo/figure.svg", "w").write("\n".join(o))
print("total          %.15f" % TOT)
print("tilted side    %.9f" % float(s))
print("straight alt   %.9f" % SA)
print("gain           %.9f" % (float(s) - SA))
print("wrote erdos106-repo/figure.svg")
