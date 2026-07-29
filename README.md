# A candidate counterexample to Erdős Problem #106

**Claim: `f(17) > 4`**, which refutes the conjecture `f(k²+1) = k` at `k = 4`.

[Erdős Problem #106](https://www.erdosproblems.com/106) asks: draw `n` squares inside the unit
square with no common interior point, and let `f(n)` be the largest possible sum of their
side-lengths. Is `f(k²+1) = k`?

This repository contains an explicit packing of **17 squares** — sixteen axis-parallel and one
tilted — with total side-length

```
2190452873/547596200  =  4.000124312403920968  >  4
```

Every coordinate is an exact rational and the verification uses **no floating point at all**.

> **Status: not yet checked by a human referee.** This is a candidate. Independent verification
> is exactly what it is published for. See [Prior work](#prior-work) for the strongest known
> objection.

---

## The construction

The configuration is not an arbitrary list of 17 numbers — it is generated in closed form by
**two rationals solving a 2×2 linear system**.

**The angle.** Tilt one square by `t` with

```
cos t = 255/257      sin t = 32/257      q = cos t + sin t = 287/257
```

exact because `255² + 32² = 65025 + 1024 = 66049 = 257²`. So `t ≈ 7.1527°`.

**The two unknowns.** Let `a` be the combined width of the two left-hand columns and `u` the
height of the middle band; set `b = a − u` and `s = u/q`. Then `a` and `u` satisfy

```
a·(cos t + 4·sin t)      +  u·cos t·(3 − sin t/q)             =   cos t + (7/4)·sin t
a·(−2·sin t − 3·cos t)   +  u·(5·sin t + 2·cos t + cos²t/q)   = −(3/4)·cos t
```

giving

```
a = 96363407/219038480      u = 211886647/1095192400
b = 67482597/273798100      s = 189738217/1095192400
```

**The layout.** Seven groups, seventeen squares, all axis-parallel except the last:

| group | count | side | placement |
|---|---|---|---|
| left column | 4 | `1/4` | `x ∈ [0, 1/4]`, stacked from `y = 0` |
| second column | 3 | `a − 1/4` | `x ∈ [1/4, a]`, stacked from `y = 0` |
| bottom-right block | 4 | `(1−a)/2` | fills `x ∈ [a, 1]`, `y ∈ [0, 1−a]` |
| top-left | 1 | `3/4 − 2b` | `x ∈ [1/4, 1−2b]`, `y ∈ [1/4+2b, 1]` |
| top-right | 2 | `b` | `y ∈ [1−b, 1]` |
| middle band | 2 | `u` | `y ∈ [1−a, 1−b]`, `x ∈ [1−2u, 1]` |
| **tilted** | 1 | `s = u/q` | inscribed in `x ∈ [1−3u, 1−2u]`, `y ∈ [1−a, 1−b]` |

**The score** collapses to

```
total = 3 + a + 2u + s = 2190452873/547596200,     exceeding 4 by  68073/547596200 ≈ 1.243e-4
```

## Why the tilt wins

Two squares intrude into the tilted square's box: the top of the second column, whose corner is
at `(a, 3a − 3/4)`, and the top-left square, whose corner is at `(1−2b, 1/4 + 2b)`. The two
equations above say precisely that **the tilted square's two left-hand edges pass through those
two corners** — it is wedged between them.

Held axis-parallel, the largest square that fits the same hole has side `0.173122157348993`.
Tilted, it reaches `0.173246469752719`. The difference is exactly the margin over 4:

```
16 axis-parallel squares       3.826877842648
  + best STRAIGHT square       0.173122157349
                             = 3.999999999997      <- the Baek–Koizumi–Ueoro optimum, 4
  + the TILTED square instead  0.173246469753
                             = 4.000124312401      <- the counterexample
```

**Remove the tilt and the configuration returns exactly 4**, reproducing the proven
axis-parallel result. The entire counterexample is one tilt in one hole.

## Verifying it

Python 3, standard library only, no dependencies.

```bash
python construct.py           # rebuild from the 2x2 system, then check it
python witness.py             # 136 separating-axis witnesses, no clipping algorithm
python verify_independent.py  # independent check by exact polygon clipping
```

Two independent methods are provided deliberately:

- `verify_independent.py` computes the exact **area** of each pairwise polygon intersection by
  Sutherland–Hodgman clipping over `Fraction`s, and checks containment corner by corner.
- `witness.py` uses **no geometry routine at all**. For each of the 136 pairs it exhibits an
  explicit separating axis and an exact rational gap, checkable by hand with four
  multiplications: if `|(c_j − c_i)·d| ≥ E_i + E_j` then the projections do not overlap in their
  interiors, so the squares share no interior point.

Output of `witness.py`:

```
pairs examined          : 136
strictly separated      : 100
touching (gap exactly 0):  36     <- legal; only interiors must be disjoint
OVERLAPPING (gap < 0)   :   0
smallest gap            : 0       (exact rational, not 1e-17)

containment: 0 walls violated, smallest clearance exactly 0
```

The word doing the work is **exactly**. In rational arithmetic there is no "microscopically
overlapping" state: a gap is a rational that is either negative or it is not.

Both checkers should first be confirmed to behave — they accept a 4×4 grid of side-1/4 squares
(total exactly 4), accept two squares touching edge-to-edge, and reject genuine overlaps
including tilted ones. Those self-tests are built in.

## What this does and does not contradict

- **Not** Baek–Koizumi–Ueoro ([arXiv:2411.07274](https://arxiv.org/abs/2411.07274)), who prove
  `g(k²+2c+1) = k + c/k` only for **axis-parallel** squares. Aristotle's
  [Lean formalisation](https://gist.github.com/lawrence-harmonic/8137bb218841d5a343a201dbd193c9d9)
  of that result defines `structure Square` with fields `x, y, s` and **no angle field**, and its
  disjointness test checks only the `x` and `y` axes — it is structurally incapable of covering
  tilted packings.
- **Not** `f(k²) = k`, which is Cauchy–Schwarz and holds for tilted squares too. Indeed
  `Σsᵢ ≤ √(17 · 0.96435) ≈ 4.0489`, so `4.000124` sits comfortably below the ceiling.
- **Not** `f(5) = 2` (Newman) or `f(2) = 1` (Erdős) — both remain proven.

**Scope.** By the monotonicity of `k(f(k²+1) − k)` (Anshul Raj Singh,
[arXiv:2601.22163](https://arxiv.org/abs/2601.22163)),

```
f(k²+1) − k  ≥  (68073/136899050)/k  ≈  4.97e-4 / k        for every k ≥ 4
```

so the conjecture fails for all `k ≥ 4`, and by Praton's equivalence
([math/0504341](https://arxiv.org/abs/math/0504341)) the whole Erdős–Soifer / Campbell–Staton
family `f(k²+2c+1) = k + c/k` falls with it. It does **not** fail at every `k`: `k = 1, 2` are
proven true, and `k = 3` (`f(10) = 3`) is untouched — an extensive search at `n = 10` found
nothing, as expected, since the 3×3 grid tiles perfectly and leaves no gap for a tilt to exploit.

## Prior work

DeepMind's *Mathematical exploration and discovery at scale*
([arXiv:2511.02864](https://arxiv.org/abs/2511.02864), §35) applied AlphaEvolve to this problem
with **the angle as a free parameter**, ran `n ∈ {10, 12, 14, 17, 26, 37, 50}`, matched only the
known value `4` at `n = 17`, and stopped: *"As we found it unlikely that a better construction
exists, we did not pursue this problem further."*

That is a search that did not find something, not an upper bound — but it is the strongest
objection to this claim and deserves an answer. A plausible one: their scorer returns `−∞` for
configurations judged invalid by a floating-point intersection test, and this configuration has
**36 pairs tangent with gap exactly 0** plus squares touching the walls at clearance exactly 0.
The optimum therefore sits precisely on that scorer's discontinuity, so a penalty-driven search
is repelled from it, for a relative prize of only `3.1e-5`. Fixing the angles instead makes the
problem a linear program, where tangencies are **active constraints** rather than
near-violations — which is where LP optima live by construction.

To the best of my searching, no published upper bound implies `f(17) ≤ 4`; the only bound known
for the unrestricted `f` is Cauchy–Schwarz.

## Files

| file | what it is |
|---|---|
| `CONSTRUCTION.md` | full derivation, with a self-contained checker to paste anywhere |
| `AUDIT_FOR_REVIEW.md` | adversarial audit request: assumptions, controls, and what could still be wrong |
| `construct.py` | builds the 17 squares from the 2×2 linear system, then verifies |
| `witness.py` | 136 separating-axis witnesses, no geometry routine |
| `verify_independent.py` | independent check by exact polygon clipping |
| `CERTIFICATE_n17.json` | the original LP-found configuration (see note below) |
| `search/` | the search code that found it originally |

**Note on the two configurations.** The result was first found by numerical search and recorded
in `CERTIFICATE_n17.json`, snapped to rationals and shrunk by `1e-13` for safety, giving
`2000062156200269/500000000000000 = 4.000124312400538`. The closed form in `construct.py`
supersedes it: it needs no shrink, so it reaches the family's exact optimum
`4.000124312403921`. Both are valid packings beating 4; the closed form is the one to check.

## Reproducing the search

`search/` contains the original pipeline. For a fixed angle vector the problem is a linear
program — that is the key structural fact — so `search/lp.py` solves centres and sizes exactly
with HiGHS, and only the angles are searched. `search/bku.py` validates the solver against all
26 proven Baek–Koizumi–Ueoro values for `n = 1…26`; it matches every one and exceeds none.

## Licence

MIT. Please verify it, and please tell me if it is wrong.
