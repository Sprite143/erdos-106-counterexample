# Independent audit request — Erdős Problem #106

I believe I have found a counterexample to a long-standing Erdős conjecture. **I want you to try to break it.** Please be adversarial: assume there is a mistake and look for it. A specific "this is wrong because X" is far more useful to me than agreement.

---

## 1. The problem (Erdős #106)

From <https://www.erdosproblems.com/106>, stated verbatim:

> Draw $n$ squares inside the unit square with no common interior point. Let $f(n)$ be the maximum possible sum of the side-lengths of the squares. Is $f(k^2+1)=k$?

Background that matters:

- $f(k^2)=k$ is **proven and easy**, by Cauchy–Schwarz: if the squares have sides $s_i$ and pairwise disjoint interiors then $\sum s_i^2 \le 1$, so $\sum s_i \le \sqrt{n\sum s_i^2} \le \sqrt n$. This bound holds for **tilted** squares too.
- $f(k^2+1) \ge k$ is easy (take the $k\times k$ grid; the extra square contributes nothing).
- Erdős, later Erdős–Soifer and Campbell–Staton, conjectured $f(k^2+1)=k$; more generally $f(k^2+2c+1)=k+c/k$ for $-k<c<k$. Praton proved these are equivalent.
- Newman proved $f(5)=2$ (the $k=2$ case).
- **Baek, Koizumi and Ueoro (2024)** proved $g(k^2+2c+1)=k+c/k$, where $g$ is the same quantity **restricted to axis-parallel squares**. The tilted case remained open.

**My claim: $f(17) > 4$, i.e. the conjecture is FALSE at $k=4$.** The counterexample uses a tilted square, so it contradicts no theorem I am aware of — in particular not Baek–Koizumi–Ueoro, which covers only axis-parallel packings.

## 2. The claimed configuration

17 squares. Sixteen are axis-parallel; square #16 is rotated by an angle $t$ with

$$\cos t = \tfrac{255}{257},\qquad \sin t = \tfrac{32}{257},\qquad t \approx 7.1527°$$

This is an **exact** rotation because $(255,32,257)$ is a Pythagorean triple: $255^2+32^2 = 65025+1024 = 66049 = 257^2$. Every corner coordinate is therefore an exact rational, and the configuration can be checked with no floating point at all.

Squares are given by centre $(c_x,c_y)$ and side length, all exact rationals:

| # | c_x | c_y | side | tilted |
|---|-----|-----|------|--------|
| 0 | `1/8` | `1/8` | `1249999999999/5000000000000` | no |
| 1 | `1/8` | `374999999999999/1000000000000000` | `1249999999999/5000000000000` | no |
| 2 | `1/8` | `312499999999999/500000000000000` | `1249999999999/5000000000000` | no |
| 3 | `125000000000001/1000000000000000` | `437499999999999/500000000000000` | `124999999999901/500000000000000` | no |
| 4 | `17248456412773/50000000000000` | `4748456412773/50000000000000` | `148389262899/781250000000` | no |
| 5 | `344969128255457/1000000000000000` | `284907384766379/1000000000000000` | `47484564127679/250000000000000` | no |
| 6 | `344969128255457/1000000000000000` | `237422820638647/500000000000000` | `94969128255357/500000000000000` | no |
| 7 | `9463288413981/25000000000000` | `435734231720381/500000000000000` | `64265768279569/250000000000000` | no |
| 8 | `579953692383191/1000000000000000` | `140015435872271/1000000000000000` | `140015435872171/500000000000000` | no |
| 9 | `144988423095797/250000000000000` | `420046307616817/1000000000000000` | `70007717936087/250000000000000` | no |
| 10 | `709795310394781/1000000000000000` | `164199160006041/250000000000000` | `48367448267487/250000000000000` | no |
| 11 | `31514865241943/50000000000000` | `876765768279619/1000000000000000` | `123234231720281/500000000000000` | no |
| 12 | `214996141031929/250000000000000` | `35003858968069/250000000000000` | `8750964742011/31250000000000` | no |
| 13 | `859984564127731/1000000000000000` | `420046307616821/1000000000000000` | `140015435872169/500000000000000` | no |
| 14 | `903265103464927/1000000000000000` | `328398320012083/500000000000000` | `96734896534973/500000000000000` | no |
| 15 | `876765768279619/1000000000000000` | `876765768279619/1000000000000000` | `123234231720281/500000000000000` | no |
| 16 | `516325517324633/1000000000000000` | `131359328004833/200000000000000` | `2706976089883/15625000000000` | **yes** |

Corners are obtained the obvious way: offset from the centre by $(\pm s/2, \pm s/2)$, rotated by $t$ if tilted.

### Claimed totals (exact)

- Sum of side-lengths: `2000062156200269/500000000000000` = 4.000124312400537896
- Excess over 4: `62156200269/500000000000000` = 1.243124e-04
- Area covered $\sum s_i^2$ = 0.964349071823464

Consistency checks that a bogus configuration would likely fail:

- Cauchy–Schwarz ceiling: $\sqrt{17 \times 0.964349072} = 4.048942$, and $4.000124 < 4.048942$. Consistent.
- $\sum s_i > 4$ forces $\sum s_i^2 > 16/17 = 0.941176$; the configuration has $0.964349$. Consistent.
- Area $< 1$, as it must be. (An area $>1$ would immediately prove overlap — this is how I caught several earlier false positives.)

## 3. How to check it yourself

Only two things need checking:

1. **Containment** — all four corners of each square lie in $[0,1]^2$.
2. **Disjoint interiors** — for each of the $\binom{17}{2}=136$ pairs, the area of the intersection is exactly $0$. Touching is *allowed*: the problem says "no common **interior** point", so zero-area contact is legal and several squares do touch.

Self-contained checker, standard library only, exact rational arithmetic:

```python
from fractions import Fraction as F

COS, SIN = F(255,257), F(32,257)
assert COS*COS + SIN*SIN == 1

def corners(cx, cy, side, tilted):
    h = F(side,2); out=[]
    for lx, ly in ((-h,-h),(h,-h),(h,h),(-h,h)):
        rx, ry = (COS*lx - SIN*ly, SIN*lx + COS*ly) if tilted else (lx, ly)
        out.append((cx+rx, cy+ry))
    return out

def area(poly):
    if len(poly) < 3: return F(0)
    s = F(0)
    for i in range(len(poly)):
        x1,y1 = poly[i]; x2,y2 = poly[(i+1)%len(poly)]
        s += x1*y2 - x2*y1
    return s/2

def clip(poly, a, b):                      # keep the part left of a->b
    if not poly: return []
    ex, ey = b[0]-a[0], b[1]-a[1]
    f = lambda p: ex*(p[1]-a[1]) - ey*(p[0]-a[0])
    out=[]
    for i in range(len(poly)):
        p, q = poly[i], poly[(i+1)%len(poly)]
        sp, sq = f(p), f(q)
        if sp >= 0: out.append(p)
        if (sp>0 and sq<0) or (sp<0 and sq>0):
            t = sp/(sp-sq)
            out.append((p[0]+t*(q[0]-p[0]), p[1]+t*(q[1]-p[1])))
    return out

def overlap(P, Q):
    poly = P
    for i in range(len(Q)):
        poly = clip(poly, Q[i], Q[(i+1)%len(Q)])
        if not poly: return F(0)
    return abs(area(poly))

DATA = [
  (0, '1/8', '1/8', '1249999999999/5000000000000', False),
  (1, '1/8', '374999999999999/1000000000000000', '1249999999999/5000000000000', False),
  (2, '1/8', '312499999999999/500000000000000', '1249999999999/5000000000000', False),
  (3, '125000000000001/1000000000000000', '437499999999999/500000000000000', '124999999999901/500000000000000', False),
  (4, '17248456412773/50000000000000', '4748456412773/50000000000000', '148389262899/781250000000', False),
  (5, '344969128255457/1000000000000000', '284907384766379/1000000000000000', '47484564127679/250000000000000', False),
  (6, '344969128255457/1000000000000000', '237422820638647/500000000000000', '94969128255357/500000000000000', False),
  (7, '9463288413981/25000000000000', '435734231720381/500000000000000', '64265768279569/250000000000000', False),
  (8, '579953692383191/1000000000000000', '140015435872271/1000000000000000', '140015435872171/500000000000000', False),
  (9, '144988423095797/250000000000000', '420046307616817/1000000000000000', '70007717936087/250000000000000', False),
  (10, '709795310394781/1000000000000000', '164199160006041/250000000000000', '48367448267487/250000000000000', False),
  (11, '31514865241943/50000000000000', '876765768279619/1000000000000000', '123234231720281/500000000000000', False),
  (12, '214996141031929/250000000000000', '35003858968069/250000000000000', '8750964742011/31250000000000', False),
  (13, '859984564127731/1000000000000000', '420046307616821/1000000000000000', '140015435872169/500000000000000', False),
  (14, '903265103464927/1000000000000000', '328398320012083/500000000000000', '96734896534973/500000000000000', False),
  (15, '876765768279619/1000000000000000', '876765768279619/1000000000000000', '123234231720281/500000000000000', False),
  (16, '516325517324633/1000000000000000', '131359328004833/200000000000000', '2706976089883/15625000000000', True),
]

polys = [corners(F(cx), F(cy), F(s), t) for (_,cx,cy,s,t) in DATA]
outside = [(i,x,y) for i,P in enumerate(polys) for (x,y) in P
           if x<0 or x>1 or y<0 or y>1]
overlaps = [(i,j, overlap(polys[i], polys[j]))
            for i in range(len(polys)) for j in range(i+1,len(polys))
            if overlap(polys[i], polys[j]) != 0]
total = sum(F(s) for (_,_,_,s,_) in DATA)
print("corners outside unit square:", len(outside))
print("overlapping pairs          :", len(overlaps))
print("total side-length          :", total, "=", float(total))
print("beats 4                    :", total > 4)
```

My run reports: **0 corners outside, 0 overlapping pairs, largest intersection area exactly `0`**, total `2000062156200269/500000000000000` > 4.

Before trusting that, please also confirm the checker *works*: it should accept a 4×4 grid of side-1/4 squares (total exactly 4), accept two squares that touch edge-to-edge, and reject two squares that genuinely overlap — including a tilted one.

## 4. Controls I already ran

I assumed I was wrong and tried to prove it:

| Control | Purpose | Result |
|---|---|---|
| Solver vs. Baek–Koizumi–Ueoro axis-parallel values, n = 1…26 | does it ever exceed a proven value? | matches all 26 exactly; exceeds none |
| n = 9, 16, 25 (f = k proven, tilting allowed) | would expose systematic overshoot | pins at 1e-10…1e-14 float noise, area exactly 1.0 |
| n = 16 built by deleting one square from **this** configuration, identical tilt optimisation | exercises the same mechanism where the answer is proven | saturates at exactly 4.000000000000 |
| Checker vs. planted overlaps / out-of-bounds / tilted overlaps | can it detect invalid packings at all? | catches every one; clears the legal one |
| Angle sweep of the tilted square | smooth geometry or numerical spike? | smooth unimodal curve, peak near 7° |

**The mechanism is geometrically explicable**, which is why I did not dismiss it. Square #16 occupies a hole in the arrangement. Holding the other sixteen squares fixed:

- the largest **axis-parallel** square that fits in that hole has side `0.173122157348993` (blocked horizontally by square #10's left edge, because square #6's corner prevents it starting further left);
- the **tilted** square achieves `0.173246469752512`, because rotating shaves its lower-left corner past square #6's corner at (0.439938257, 0.569814770), letting it begin further left.

The difference is `0.000124312403`, and the certified excess over 4 is `0.000124312400` — **the excess is exactly the tilt gain**, to eleven decimal places. It reconciles end to end: the sixteen axis-parallel squares sum to `3.826877842648026`, so

- with the best axis-parallel square in the hole: `3.999999999997` — i.e. exactly the Baek–Koizumi–Ueoro axis-parallel optimum of 4, less the 3e-12 shrink margin;
- with the tilted square instead: `4.000124312400538`.

So the whole counterexample is one tilt in one hole, and nothing else.

**Note on contact.** At the true optimum the tilted square's edges are *tangent* to two neighbouring corners (those constraints are active in the optimisation). The delivered certificate deliberately backs off from tangency: coordinates are snapped to rationals with denominator $10^{15}$ and every square shrunk by $10^{-13}$, so all clearances are strictly positive (≈$10^{-13}$ to the walls and between squares). That costs ~3e-12 of margin and leaves 1.24e-4. The certificate therefore does **not** rely on touching being legal.

Note also **why $n=10$ gives nothing**: the optimum there is the 3×3 grid, a perfect tiling with zero gaps, so a tilted square has nothing to exploit. $n=17$ has ≈3.6% waste. This asymmetry is itself a consistency check.

## 5. Where I think the risk actually is

Please attack these specifically:

1. ~~**Modelling.**~~ **Resolved by a machine-checked artifact.** I read "no common interior point" as **pairwise disjoint interiors**. Aristotle's Lean formalisation of Baek–Koizumi–Ueoro ([gist](https://gist.github.com/lawrence-harmonic/8137bb218841d5a343a201dbd193c9d9)) defines exactly that: `Packing.is_valid P := ∀ i j, i ≠ j → (P i).disjoint (P j)`. Pairwise, not "no point interior to all $n$". It also confirms the objective is $\sum_i s_i$ (`total_side_length := ∑ i, (P i).s`).
2. ~~**Is $f$ really unrestricted in orientation?**~~ **Strongly supported.** The same Lean file's `structure Square` has fields `x, y, s` and **no angle field at all** — it is axis-aligned by construction, and its `Square.disjoint` tests only the $x$ and $y$ axes, which is valid *only* for axis-parallel squares. The theorem it proves is stated about `g`, never `f`. So the formalised BKU result is structurally incapable of saying anything about tilted packings, which is exactly the gap this counterexample occupies.
3. ~~**Is there an existing theorem implying $f(17)\le 4$?**~~ **Searched; none found.** The only upper bound known for the unrestricted $f$ appears to be Cauchy–Schwarz, $f(17)\le\sqrt{17}=4.1231$, which $4.000124$ satisfies with room to spare. Baek–Koizumi–Ueoro bound only the axis-parallel $g$; Praton and Singh prove equivalences, not bounds; Semantic Scholar lists no papers citing Baek–Koizumi–Ueoro. Please still sanity-check this — it is the one thing I cannot prove a negative about.
4. **Is the arithmetic right?** Recompute the corners and all 136 pairwise intersections.
5. ~~**Does "inside the unit square" permit touching the boundary?**~~ **Resolved — touching is explicitly legal.** The Lean formalisation uses non-strict inequalities throughout: containment is `x + s ≤ 1`, and `Square.disjoint` counts `sq1.x + sq1.s ≤ sq2.x` (edge-to-edge contact) as disjoint. So both boundary contact and square-to-square contact are permitted, as the phrase "no common *interior* point" implies.

### Consistency with the two published constraints

- **Cauchy–Schwarz ceiling.** $f(k^2+1)-k \le \sqrt{k^2+1}-k \approx 1/(2k)$. Our claim gives, via monotonicity, $f(k^2+1)-k \ge (68073/136899050)/k \approx 4.97\times10^{-4}/k$. Both are $\Theta(1/k)$ and our floor sits about **1000× below** the ceiling at every $k$ — comfortable, not borderline.
- **Singh's criterion.** $f(k^2+1)=k$ for all $k$ iff $\sum_k (f(k^2+1)-k)$ converges. Our lower bound is $\Theta(1/k)$, so the sum diverges harmonically — which by Singh's *iff* means the conjecture must fail for some $k$. Consistent, and it is the same $k\ge4$ range the monotonicity argument gives.

## 5b. Prior computational work — AlphaEvolve searched this exact case and found only 4

This is the strongest objection to the claim and I want it on the table rather than discovered later.

DeepMind's *Mathematical exploration and discovery at scale* ([arXiv:2511.02864](https://arxiv.org/abs/2511.02864), §35, "Erdős squares in a square problem") applied AlphaEvolve to precisely this problem. In their words:

> "The squares were defined by the coordinates of their center, **their angle**, and their side length. If the configuration was invalid (the squares were not in the unit square or they intersected), then the program received a score of minus infinity, and otherwise the score was the sum of side lengths of the squares. AlphaEvolve matched the best known constructions for $n \in \{10, 12, 14, 17, 26, 37, 50\}$ but did not find them for some larger values of $n$. As we found it unlikely that a better construction exists, we did not pursue this problem further."

So a serious, well-resourced evolutionary search **with rotation in the search space** ran $n=17$ and reached only $4$.

That is not an upper bound — it is a search that did not find something — but it demands an explanation, and there is a specific structural one:

- **The optimum sits exactly on their scorer's discontinuity.** Our configuration has **36 pairs touching with gap exactly $0$**, and several squares touching the walls with clearance exactly $0$ (see §6b). It is maximally on the boundary between "valid" and "intersecting".
- Their validity test is floating-point and returns $-\infty$ for intersection. Any float tolerance either rejects exact tangency — placing an $-\infty$ cliff precisely at the optimum and repelling the search from it — or admits small overlaps, which corrupts the score. An evolutionary method has to *approach* the optimum through the region its own fitness function forbids.
- **The prize is tiny relative to that cliff**: $1.243\times10^{-4}$ on a total of $4$, a relative improvement of $3.1\times10^{-5}$.
- Our method does not have this problem *by construction*. Fixing the angles makes the whole thing a linear program, so tangencies are **active constraints** rather than near-violations — the LP lands on the boundary exactly, because that is where LP optima live. The exact rational angle then removes floating point from the verification entirely.
- They also stopped by choice, not by exhaustion: "we did not pursue this problem further."

I offer this as the most likely explanation, not a certainty. If instead the resolution is that our configuration is subtly invalid, §6b's 136 explicit separating-axis witnesses in exact rational arithmetic are where that would have to show up.

## 6. What I am asking for

- Confirm or refute the arithmetic (§3).
- Rule on the modelling questions (§5.1, §5.2, §5.5).
- Search the literature for anything implying $f(17)\le4$ (§5.3).

If all of that holds, then $f(17)>4$; and by the monotonicity of $k\big(f(k^2+1)-k\big)$ noted by Raj Singh on the problem page, $f(k^2+1)>k$ for **every** $k\ge4$ — while still permitting $f(10)=3$, consistent with my search finding nothing at $n=10$.

I would much rather learn now that this is wrong than after showing it to anyone else.
