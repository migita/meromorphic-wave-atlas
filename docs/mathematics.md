# Mathematical scope and conventions

The atlas constructs exact meromorphic travelling-wave families in declared coefficient slices. The pure-power range is `2 ≤ p ≤ 6, 2 ≤ n ≤ 7`. Here `p` is the travelling-wave ODE order; an evolution equation integrated once generally has PDE order `p+1`.

The necessary pole balance is `q = p/(n−1) ∈ Z>0`. All thirteen admissible pairs are represented, with six companion slices making nineteen pure-power slices. The separate C4 constructor adds nonlinear derivative terms.

| Pair (n,p) | Main construction |
|---|---|
| (2,2) | KdV energy line |
| (3,2) | mKdV energy line, character two |
| (2,3) | Dispersive KS compatibility slice |
| (4,3) | Quartic flux, character three |
| (2,4) | Kawahara; companion with an odd derivative |
| (3,4) | Even cubic operator, ordinary and twisted branches; opposite linear sign and mixed companions |
| (5,4) | Quintic, character four; even companion |
| (2,5) | Nikolaevskiy-type quadratic family, with exceptional fibres |
| (6,5) | Sextic, character five; the scaling pole t=3 is excluded |
| (2,6) | Quadratic sixth-order compatibility branches |
| (3,6) | Twisted cubic sixth order; ordinary companion |
| (4,6) | Quartic sixth-order ordinary branch |
| (7,6) | Septic, character six; even companion |

Polynomial extensions complete the finite table; inclusion does not assert that every extension is a standard physical model.

## Lattice and profile

```text
X = ℘(z;g₂,g₃), Y = ℘′(z;g₂,g₃)
Y² = 4X³ − g₂X − g₃
Δ = g₂³ − 27g₃²
```

The invariants belong to the full twisted period lattice Γ. A character of order d gives d poles in an ordinary period cell. In particular, an ordinary lattice and its index-two cover must not be confused when comparing the cubic fourth-order branches.

Off the cusp, the lattice is elliptic. On `Δ=0`, it degenerates; the origin is the rational degeneration. Real invariants do not alone guarantee a smooth real-valued bounded profile. The slice, amplitude normalization and analytic continuation still matter. The engine explicitly continues finite characters through zeros rather than gluing principal roots.

For three real roots `e₁ ≥ e₂ ≥ e₃`, the bounded real oval is

```text
X(z) = e₃ + (e₂−e₃) sn²(√(e₁−e₃) z | m)
m = (e₂−e₃)/(e₁−e₃).
```

Here `m` is the Jacobi parameter, the square of the modulus. Its complementary parameter `mc=1−m` is stored separately in the precise Kawahara transition. The full profile period can be a multiple of the lattice’s real period. A collapsed real oval gives a constant slice, which has no primitive period. See [DLMF §23.6](https://dlmf.nist.gov/23.6#ii).

## Normalization

Except for the physical KS and Kawahara displays, pure-power waves have leading Laurent coefficient one:

```text
P(D)u − K uⁿ + C = 0, K = (−1)^p (q)_p.
```

The additive constant is used for quadratic nonlinearities in the shifted equation with zero linear coefficient, where declared. To restore the paper’s `vⁿ/n` convention, set `v=c u` with `c^(n−1)=−nK`; the constant becomes `cC`. For odd nonlinearities this amplitude can be complex.

The Kawahara display uses `v=−1680u`, giving `v⁗+a v″−v+v²/2=0`. KdV retains `u″−u−6u²=0`; its displayed solitary limit is `u=−¼ sech²(z/2)`. These normalizations are stated next to the equations.

A dilation `uλ(z)=λ^q u(λz)` changes `a_j → λ^(p−j)a_j`, `g₂ → λ⁴g₂`, `g₃ → λ⁶g₃`. The CSV operator column uses increasing order `[a₀,a₁,…,a_{p−1},1]`. Complete coefficient formulas and normalization choices are in [models.json](../data/models.json).

## What is counted

- Pure-power counts describe constructed nonconstant meromorphic wave classes modulo translation and phase in the declared slice.
- For KdV and mKdV, an energy coordinate varies within one conservative equation. Those equations carry continuous families of waves.
- Kawahara’s two real lattices exist for `|a|<13/6`, coalesce into one degenerate lattice at equality, and leave the real invariant plane beyond it. Complex conjugate lattices remain outside this plane.
- For quadratic orders five and six and twisted cubic order six, different algebraic roots at one slider value generally recover different remaining coefficients. The readout counts compatible coefficient choices, not multiple waves of one fixed equation.
- C4 counts displayed representatives of `v=α℘+β`, including separate amplitudes on a shared lattice. Free coordinates can parameterize a continuous family for one fixed equation. Multiple-pole elliptic and genus-two sectors are outside this constructor.

Curve selection snaps to solved frames. High-precision Kawahara limit states are separate from the original parameter grid. Their displayed profile retains the complementary Jacobi parameter even when a rounded coefficient cannot distinguish two nearby states.

## Algebraic records and limits

The exact Kawahara lattice conditions are

```text
662158224 g₂² − 653016 a²g₂ + 1457a⁴ − 28561 = 0
g₃ = a(31a² − 42588g₂)/4745520.
```

The two coalescences are `(g₂,g₃)=(1/432, ±1/46656)`. The smooth branch tends to `2−(35/12)sech⁴(z/(2√6))` at `a=−13/6` and to the constant `2` at `a=13/6`. At the constant end the limit of the preceding wavelengths is `2π√6`; that number is not assigned as the primitive period of a constant.

For the fixed Lax C4 equation with `a₂=0,a₀=−1,a₀₀=0`, the `α=−2` family has `g₂=½−15β²`, `g₃=35β³−β`. There is also an isolated `α=−6` representative at `(1/42,0)`. At `β=1/6` the regular profile is `½ sech²(z/2)`. The undivided coefficient equations and resonances are verified in [c4_formulas.json](../data/c4_formulas.json).

For quadratic order five, the vanishing readout denominators at `a=−1/7,−8/441,26/49` are handled with the original polynomial equations. The finite fibres and boundary elimination factors are in [boundaries.json](../data/boundaries.json).

The [exact verification record](../data/verification.json) checks full ODE substitutions, all algebraic roots of the tested implicit equations, degenerate operators and a perturbed negative control. [Profile checks](../data/profile_checks.json) use independent Taylor jets, periods and zero continuation; [C4 checks](../data/c4_checks.json) independently substitute the complete nonlinear ODE. The [Kawahara transition checks](../data/journey_checks.json) compare JavaScript against 110-digit references and independently verify the pulse and periodic averages.

The plotted values and individual readouts are numerical evaluations. A finite plot window is not a global census or an interval certificate. Completeness outside the constructed sectors and stability are separate theoretical questions.
