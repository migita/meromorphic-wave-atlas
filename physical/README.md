# Kawahara in wavelength–amplitude–mean coordinates

[Main figure](kawahara_physical_coordinates.png) · [Rotate the 3-D curve and select profiles](kawahara_physical.html)
· [Three-page PDF](kawahara_physical.pdf) · [Pulse comparison](kawahara_pulse_limit.png)
· [CSV](kawahara_physical.csv) · [Checks](checks.json)

This example maps the smooth real branch of

\[
v''''+a v''-v+\frac{v^2}{2}=0,\qquad -13/6<a<13/6,
\]

from the lattice plane into observable coordinates:

- **L:** the full spatial repeat distance of the profile;
- **A:** peak-to-trough amplitude, `max(v)−min(v)`;
- **M:** the average of `v` over one full period.

The three projections use the same curve and parameter colors. Profiles A–C are shown on a common
physical spatial axis; D marks a longer-period wave. The separate pulse figure shows periodic
depressions separating into a solitary pulse. Changing `a` changes the dispersion coefficient;
these are different equations in the displayed normalized family. The figures do not encode stability.

## The map

Take the larger real root `g₂(a)` of

\[
662158224g_2^2-653016a^2g_2+1457a^4-28561=0,
\qquad g_3=\frac{a(31a^2-42588g_2)}{4745520}.
\]

This branch has `Δ>0` in the open parameter interval. Let `e₁>e₂>e₃` be the roots of
`4X³−g₂X−g₃=0`, and put

\[
G=e_1-e_3,\qquad m=\frac{e_2-e_3}{G},\qquad
X(z)=e_3+(e_2-e_3)\operatorname{sn}^2(\sqrt G\,z\mid m).
\]

This is the bounded real oval of the Weierstrass curve. Here `m` is the **parameter** of the Jacobi
function, the square of its modulus. In these variables the real wave is

\[
v(z)=-1680\left[X(z)^2+\frac{a}{78}X(z)+b\right],
\qquad b=-\frac{g_2}{10}-\frac{31a^2+507}{851760}.
\]

Using complete elliptic integrals with parameter `m`,

\[
L=\frac{2K(m)}{\sqrt G},\quad
\langle X\rangle=e_1-G\frac{E(m)}{K(m)},\quad
\langle X^2\rangle=\frac{g_2}{12},
\]

\[
M=-1680\left[\frac{g_2}{12}+\frac{a}{78}\langle X\rangle+b\right].
\]

For the amplitude, evaluate the quadratic profile at `X=e₃,e₂` and also at its vertex
`X=−a/156` when that point lies between them. The largest minus smallest value is `A`.
Thus an interior extremum is retained; endpoint values alone would give an incorrect amplitude
on part of the branch.

The Jacobi/Weierstrass change of variables is the standard one in
[DLMF §23.6(ii)](https://dlmf.nist.gov/23.6#ii); the complete elliptic-integral conventions are in
[DLMF §19.2](https://dlmf.nist.gov/19.2). The Kawahara matching equations and profile are verified directly by the
[standalone calculation](kawahara.py).

## What the plots show

Wavelength is not monotone in `a`. Numerically its minimum on this branch is approximately
`L=12.0402258` at `a=−1.42207584`; the mean has a different minimum, about `0.53259997` at
`a=−1.88061096`. These are numerical optimizations of explicit formulas, not interval certificates.

At the small-amplitude end,

\[
a\to13/6:\quad (L,A,M)\to(2\pi\sqrt6,0,2).
\]

The limiting profile is constant, so the displayed finite wavelength is a **branch limit**, not
an assigned primitive period of a constant. Near this endpoint a full period contains two slightly
different oscillations: the leading small-amplitude harmonic has period `L/2`, while the full
nonconstant profile has period `L`. The arrows on the profile panels mark the full repeat distance.

At the other end,

\[
a\to-13/6:\quad L\to\infty,\qquad A\to35/12,\qquad M\to2,
\]

and the centered profiles converge to

\[
v_*(z)=2-\frac{35}{12}\operatorname{sech}^{4}\!\left(\frac{z}{2\sqrt6}\right).
\]

The integrated depression of this pulse is `∫(2−v*)dz=70√6/9`. Consequently

\[
M=2-\frac{70\sqrt6}{9L}+o(L^{-1}).
\]

The pulse is not placed at a finite wavelength in the plots. Its mean is quoted as the limit of
periodic averages. The wavelength axes are explicitly marked when logarithmic. All spatial and
amplitude units are those of the normalized equation; no mean subtraction or amplitude rescaling
is applied to the profiles.

## Checking and reproduction

Run `python3 kawahara.py` from this directory, or use its absolute path. Dependencies are those
of the atlas plus `plotly` for the offline interactive figure. No server is required to open the HTML.

The script checks the full ODE symbolically modulo the exact lattice equations and checks the
solitary pulse by a separate hyperbolic substitution. Seven independent numerical integrations
verify the mean formula and the averaged ODE identity `⟨v²⟩=2M`; dense profile sampling checks
the extrema. Endpoint checks verify the small-amplitude period and the pulse's shape and integrated
depression. Results are saved in `checks.json`.

Near `a=−13/6`, a finite but long period requires a coefficient exponentially close to the endpoint.
The calculation uses up to 135 decimal digits; the CSV retains both the distance to the endpoint
and high-precision decimal values of `a,g₂,g₃,m`. The convenient floating-point `a` column alone
cannot distinguish the longest-period states. The sampled wavelength range ends near 236; the
asymptotic endpoint remains at infinity.
