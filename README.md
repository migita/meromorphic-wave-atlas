# Meromorphic Wave Atlas

**[Open the public explorer](https://migita.github.io/meromorphic-wave-atlas/)**

An interactive atlas of exact travelling waves, by Alexander Migita. Select an equation family,
move its parameter slider, and select a lattice branch to see the corresponding solution immediately.
The explorer covers 13 admissible `(n,p)` pairs through ODE order six and 19 coefficient slices.

The solution panel offers a pole-free slice when the lattice has a bounded real oval, and a
meromorphic real-axis view with poles marked and the vertical display range clipped. Complex
solutions show separate real and imaginary parts. For real periodic profiles it also reports the
full real period, peak-to-trough amplitude, and period mean. The two Kawahara endpoint buttons
show the solitary pulse, the constant small-amplitude limit, and the nonconstant trigonometric
degeneration. Kawahara is displayed in physical normalization `v=−1680u`, matching the separate
physical-coordinate example; other profiles solve their displayed normalization.

This repository is self-contained: all assets are local, and the small exact torsion-divisor
constructor is included in `elliptic_divisor.py`. The static site is served from the root of `main`
by GitHub Pages, with `.nojekyll`; opening `index.html` directly also works offline.

## Build and check

```sh
python3 -m pip install -r requirements.txt
python3 build.py
python3 test_profiles.py       # Node.js needed for the browser math engine
node browser_check.mjs         # Node >=22 and a local Chromium browser
```

To rebuild only the interactive profile data and page, use `python3 profiles.py`.
The [profile validation record](data/profile_checks.json) checks the JavaScript output against
independent Taylor-jet substitutions into the full ODEs, as well as periods and analytic continuation
through zeros. The [exact algebra checks](data/verification.json) verify the underlying families.
Eleven sample jets too close to zeros for reliable fractional-power differentiation are explicitly
reported separately; their power identities and the relevant zero continuations are still checked.
Wave values and physical readouts in the interactive panel are numerical evaluations, not interval certificates.

## The atlas

Open the [offline interactive gallery](index.html), the [21-page PDF book](g2g3_atlas.pdf), or the
[overview image](figures/atlas_overview.png). Every individual plate is also supplied as PNG, SVG,
and PDF in [figures](figures/). The gallery lets you select a family, move its parameter, inspect
the selected lattice points, and read the corresponding equation coefficients.

For a concrete physical readout, see [Kawahara in wavelength–amplitude–mean coordinates](physical/README.md):
[plots and profiles](physical/kawahara_physical_coordinates.png), an [interactive 3-D curve](physical/kawahara_physical.html),
and the [solitary-pulse limit](physical/kawahara_pulse_limit.png).

The agreed range is **2 ≤ p ≤ 6, 2 ≤ n ≤ 7**. All **13** pairs satisfying the necessary pole
balance `q = p/(n−1) ∈ Z_{>0}` are represented. Six companion slices make **19 plates** in total.
The other 17 pairs are marked in the overview's coverage table. Here `p` is the order of the
travelling-wave ODE; for an evolution equation integrated once, the PDE has order `p+1`.

This is a standalone computational illustration of exact travelling-wave families. All algebraic
identities needed to reproduce the displayed curves are included in this repository.

## Pictures

| (n,p) | Main plate | Additional slices |
|---|---|---|
| (2,2) | [KdV: the free-energy line](figures/n2_p2_kdv.png) | A continuum for one conservative equation |
| (3,2) | [mKdV: a tilted energy line](figures/n3_p2_mkdv.png) | Character of order two |
| (2,3) | [Kuramoto–Sivashinsky](figures/n2_p3_ks.png) | The compatibility condition is `b² = 16a` |
| (4,3) | [Quartic flux: the threefold character](figures/n4_p3_quartic.png) | Dixon-type character space |
| (2,4) | [Kawahara: the returning loop](figures/n2_p4_kawahara.png) | [Kawahara–KS with an odd derivative](figures/n2_p4_mixed.png) |
| (3,4) | [Cubic fourth order: two lattice curves](figures/n3_p4_cubic_positive.png) | [Opposite constant-coefficient sign](figures/n3_p4_cubic_negative.png); [dispersive Swift–Hohenberg](figures/n3_p4_mixed.png) |
| (5,4) | [Quintic: a fourfold contact](figures/n5_p4_quintic.png) | [Even operator, character two](figures/n5_p4_even.png) |
| (2,5) | [Nikolaevskiy type: four algebraic branches](figures/n2_p5_quadratic_fifth.png) | Exact quartic condition and three recovered exceptional fibres |
| (6,5) | [Sextic: a fivefold character](figures/n6_p5_sextic.png) | The `a₄=1` chart; the scaling pole `t=3` is excluded |
| (2,6) | [Quadratic sixth order: a three-branch fold](figures/n2_p6_quadratic_sixth.png) | Exact cubic condition |
| (3,6) | [Cubic sixth order: twisted branches](figures/n3_p6_cubic_twisted.png) | [The ordinary branch at `a₄=0`](figures/n3_p6_ordinary.png) |
| (4,6) | [Quartic sixth order: a parabola](figures/n4_p6_quartic.png) | Ordinary wave `u=℘+1/168` |
| (7,6) | [Septic: a sixfold character](figures/n7_p6_septic.png) | [Even operator, character two](figures/n7_p6_even.png) |

The common physical reductions motivate the range. Quartic, sextic, and septic fluxes are included
as polynomial extensions so that the low-order table has no omitted admissible pair; the table
does not assert that every extension is a standard physical model. Equations with nonlinear derivative terms are outside the pure-power scope of this atlas.

## Reading the figures and counts

All planes use the standard invariants of

`X=℘(z;g₂,g₃), Y=℘′(z;g₂,g₃), Y²=4X³−g₂X−g₃`.

The grey cusp is `Δ=g₂³−27g₃²=0`; the pale region has `Δ>0`. Points off the cusp are genuine
elliptic lattices. Cusp points are degenerations, and the origin is the rational degeneration.
The invariants belong to the **full twisted period lattice Γ**, as in the structure paper.
A character of order `d` gives `d` poles in the ordinary period cell. In particular, the two cubic
fourth-order curves must not be combined by confusing an ordinary lattice with its index-two cover.

Colors encode the declared parameter, not discriminant sign. Diamonds mark real degeneration
parameters after checking the discriminant; resultant candidates are filtered, so other waves
at the same parameter are not automatically marked. Each figure has a separate enlarged view.
Finite view windows can crop large branches, especially in the order-five character normalization.
The CSVs retain all finite real points in the sampled parameter interval, including those outside
the displayed window. The gallery explicitly reports when a selected point is outside the view.

**A real point means real lattice invariants, not automatically a smooth, bounded real-valued
physical wave.** The profile, amplitude normalization, and choice of a real slice still matter.
The counts are of nonconstant meromorphic wave classes, modulo translation and phase, within
the declared slice. Complex conjugate lattices can leave this real plane without the corresponding
complex meromorphic waves ceasing to exist.

There are two distinct uses of the slider:

- **A fully fixed equation for each parameter.** Kawahara has two real lattice solutions for
  `|a|<13/6`, one degenerate solution at equality, and no real lattice solutions beyond it
  (two complex conjugate solutions remain). The even cubic fourth-order operator also has two
  character branches; in the `a₀=+1` slice they merge at `a=±5/2`. The other compatible families
  with a one-wave bound have one wave per displayed fixed equation.
- **An equation reconstructed along its compatibility locus.** For quadratic orders five and
  six, and the twisted cubic sixth-order slice, different algebraic roots at one slider value
  generally give different remaining coefficients. The small chart says **compatible coefficient
  choices**, not “four waves of one equation.” The full coefficients and integration constant
  are shown in the gallery and recorded in each CSV. Coalescing algebraic roots are counted once.

For `p=2` the two plates instead vary the first integral of **one conservative equation**. Each
such equation carries infinitely many waves. Mixed second-order exceptions can have meromorphic
solutions outside the elliptic/exponential/rational class; these are not claimed to be represented
by a lattice plane. Counts refer to the full set of coefficients in the stated construction, not just to one
projection of its coefficient locus.

The atlas covers **every pair in the agreed finite range**, with explicitly declared normalized
coefficient slices. It is not a global decomposition of every high-dimensional coefficient space,
nor a catalogue of all independent rational/exponential fronts off the elliptic compatibility locus.
For example, the separate KS fronts in the earlier picture are not extra points on its elliptic line.

## Equations and normalization

Except for the familiar physical KS formula, profiles have leading Laurent coefficient one:

`P(D)u − K uⁿ + C = 0,  K = (−1)^p (q)_p`.

The integration constant `C` is used only for quadratic nonlinearities, in the shifted form
with zero linear coefficient. To restore the paper's `vⁿ/n` convention, put `v=c u` with
`c^(n−1)=−nK`; its integration constant is `cC`. For odd nonlinearities this amplitude can be
complex. A dilation `uλ(z)=λ^q u(λz)` changes
`a_j → λ^(p−j)a_j`, `g₂ → λ⁴g₂`, `g₃ → λ⁶g₃`.

All exact formulas, operator coefficients, profiles/powers, intervals, and normalization choices
are in [models.json](data/models.json), generated by [models.py](models.py). The CSV operator column
is in increasing order `[a₀,a₁,…,a_{p−1},1]`; the next column is the additive constant. For the KS
plate it is the displayed physical equation `v‴+4v″+v′+v²/2+C=0`.

Three useful explicit examples:

1. **Kawahara:**
   `662158224 g₂² − 653016 a² g₂ + 1457 a⁴ − 28561 = 0`,
   `g₃ = a(31a²−42588g₂)/4745520`.
   Its exact coalescences are `(g₂,g₃)=(1/432, ±1/46656)`.
2. **Cubic fourth order:** write `b=a₂/60`, `γ=a₀/360−b²`. The ordinary branch is
   `(20γ,20b(b²−γ))`; the twisted branch is
   `(20b²−40γ/3,−8b³+80bγ/3)`.
3. **Quintic fourth order, `a₃=1`:**
   `g₂=(1600a²−240a−39)/120000`,
   `g₃=(40a−9)(1600a²−2160a+441)/216000000`,
   `Δ=(4a−1)(40a−11)^4/1600000000000`.
   The fourth-order contact is not a self-intersecting loop.

The new higher-order ordinary branches are obtained by exact substitution in the Weierstrass
coordinate ring. Characters of orders 3–6 use a torsion divisor `div(F)=mQ−mO`, with `u=F^(1/m)`.
The included [torsion constructor](elliptic_divisor.py) checks exact order
and constructs `F` by elliptic-curve addition. The Tate equations used before coefficient scaling are:

| Character order m | Weierstrass equation with Q=(0,0) |
|---|---|
| 3 | `y²+xy+t y=x³` |
| 4 | `y²+xy−t y=x³−t x²` |
| 5 | `y²+(1−t)xy−t y=x³−t x²` |
| 6 | `y²+(1−t)xy−t(1+t)y=x³−t(1+t)x²` |

Coordinates are converted to `Y²=4X³−g₂X−g₃` before plotting. Principal-part cancellation determines
the operator; a holomorphic function with the same nontrivial finite character must vanish, which
explains why those finite equations give a global wave. The formulas and exact-order checks are
in [derived_tate.json](data/derived_tate.json).

For quadratic order five, a linear readout denominator vanishes at `a=−1/7,−8/441,26/49`.
[boundaries.py](boundaries.py) returns to the undivided polynomial equations at these parameters
and recovers the finite fibres. They are included in the plots and verified separately. The
[boundary record](data/boundaries.json) also stores the exact elimination factors used to find
folds and cusp contacts of all three implicit families.

## Verification and reproduction

The [verification record](data/verification.json) contains exact checks of the full differential
equations at rational parameter values, including **all algebraic roots** of each implicit matching
polynomial, eight degenerate operators checked directly after rational-in-exponential substitution,
and a negative control that rejects a perturbed off-locus point. These are exact substitutions,
not floating-point residual fits. The generic formula derivations are retained separately.
The sampled curves and isolated high-degree boundary coordinates are numerical, not interval
certificates. The displayed counts are those of the stated algebraic constructions. Completeness beyond
these explicitly constructed families is a separate theoretical question.

Run everything from any directory:

```sh
python3 /path/to/g2g3_atlas_20261002/build.py
```

Or run the stages separately from this directory:

```sh
python3 derive.py ordinary
python3 derive.py twisted
python3 derive.py tate
python3 derive.py mixed
python3 boundaries.py
python3 verify.py
python3 render.py
```

Dependencies: Python, SymPy, NumPy, SciPy, and Matplotlib; see [requirements.txt](requirements.txt).
The gallery itself requires none of them and works directly from its HTML file. The optional
`node browser_check.mjs` checks all 19 selectors, 57 slider positions, the exact Kawahara
coalescences, links, and desktop/mobile layout. It uses a locally installed Chromium; set
`ATLAS_CHROME` if needed. See [browser_check.json](data/browser_check.json).

The mathematical conventions are the standard Weierstrass equation `Y²=4X³−g₂X−g₃` and the monic
pure-power wave equation displayed above. The construction files, exact identities, and
numerical profile checks document the specific solution families shown here. This explorer is
a computational illustration, not a stability analysis or a new theorem about all travelling waves.
