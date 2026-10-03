# Meromorphic Wave Atlas

**[Open the interactive atlas](https://migita.github.io/meromorphic-wave-atlas/)** · Alexander Migita

Exact travelling waves and the geometry of their elliptic lattices. Every graph in the explorer is interactive; no picture gallery, external font, CDN, or package installation is needed to view it. Open `index.html` directly to work offline.

## Two ways to explore

**Discover** starts with seven named examples: Kawahara, KdV, dispersive Kuramoto–Sivashinsky, Lax fifth-order KdV, Sawada–Kotera, Kaup–Kupershmidt, and a Swift–Hohenberg stationary reduction. Each declares its equation, normalization, constructed sector, and literature sources. All nineteen pure-power coefficient slices remain available in the full-atlas selector.

**Equation lab** lets you edit all six real C4 coefficients, including zero terms and fractions, or use the restricted Kawahara, dispersive KS, and KdV templates. The C4 equation is

```text
v⁗ + A v v″ + B (v′)² + C v³ + a₂ v″ + a₀ v + a₀₀ = 0.
```

The C4 constructor solves the undivided coefficient equations in the sector `v = α℘ + β`. It includes resonances with free `β`, `g₂`, or `g₃`, complex amplitudes on real lattices, and separate representatives sharing a lattice. Its counts describe this constructed sector. The other lab templates declare which coordinate varies and which coefficients are fixed; KdV’s energy varies within one fixed equation.

## Interactions

- Click or tap **anywhere along a lattice curve** to select a solved parameter and branch. Hover reveals its coordinates. Curve selection snaps to the actual solved frames.
- Drag the lattice to pan; use `+`, `−`, and **Fit** to zoom. With the plot focused, arrow keys select neighbouring parameters or branches; `0` fits the finite curves. Ctrl/Alt + wheel also zooms.
- Click the **wave** to pin a value, or type a coordinate into **Probe z**. Choose a common spatial window to compare profiles without rescaling their spatial coordinate.
- Use the parameter field for numbers or fractions. C4 evaluates the requested coordinate directly; other templates select the nearest precomputed frame and report a snap when needed.
- Choose the pole-free or meromorphic slice, compare the solitary limit, export wave or lattice CSV, and **Copy link** to preserve the equation, branch, slice, window, probe, zoom, and pulse-transition state. Custom coefficients are remembered when switching between Discover and the lab.

Finite-period windows are limited to twenty lattice periods so individual oscillations remain resolved. **Reset view** restores the plot controls without changing your equation or parameter.

### Kawahara: a periodic train becomes one pulse

The limit control follows the smooth branch of

```text
v⁗ + a v″ − v + v²/2 = 0
```

on a common spatial scale. The periodic depressions separate as `a → −13/6`, with the exact limit

```text
v*(z) = 2 − (35/12) sech⁴(z/(2√6)).
```

The 161 transition states are computed with 110 decimal digits before compilation. The complementary Jacobi parameter is stored explicitly, retaining finite periods when `a` and `m` would otherwise round to the endpoint in JavaScript. The last state is the exact infinite-period pulse. A localized pulse has no period mean; the mean of the preceding periodic family is a different quantity.

**Play transition** varies the equation parameter. It does not simulate time evolution. Wave values, readouts and sampled curves are numerical evaluations, not interval certificates or stability results.

### Dispersive KS

The named example is the dispersive slice

```text
v‴ + 4v″ + v′ + v²/2 + C = 0,
```

which satisfies `b² = 16μν` with `b = 4, μ = ν = 1`. Its smooth real periodic interval is `−18 < C < −8`. At `C = −18` the selected real oval tends to an asymmetric localized pulse; at `C = −8` it becomes constant. This identifies the displayed slice explicitly rather than suggesting generic elliptic waves for every KS equation. See [Eremenko’s meromorphic classification](https://arxiv.org/abs/nlin/0504053).

## Rebuild and verify

Quick rebuild of the page and precise Kawahara transition, using the existing verified atlas data:

```sh
python3 journeys.py
```

Recompute all interactive data and exact algebra:

```sh
python3 -m pip install -r requirements.txt
python3 build.py
```

The default build produces interactive HTML, data, CSVs and validation fixtures. For changes to profile construction alone, run `python3 profiles.py`.

```sh
python3 test_profiles.py    # independent Taylor jets, ODEs, characters and periods
python3 test_c4.py          # independent C4 lattice, resonance and full ODE checks
python3 test_journeys.py    # 110-digit references, pulse identity and window-independent readouts
node browser_check.mjs     # Chromium: real clicks, controls, links, animation and mobile layout
```

The browser check uses Node ≥22 and a local Chromium installation. It detects common Playwright cache layouts; otherwise set `ATLAS_CHROME` to the executable. Set `ATLAS_URL` to check the published page with the same suite.

Validation records: [exact algebra](data/verification.json), [profile checks](data/profile_checks.json), [C4 checks](data/c4_checks.json), [pulse-transition checks](data/journey_checks.json), [browser checks](data/browser_check.json).

Eleven original fractional-power jets near zeros are reported as ill conditioned for floating-point differentiation; their power identities and analytic continuation through zeros are still tested. Curve clicking never substitutes an interpolated off-locus point for a solved wave. Period means and amplitudes use a full real wave period, independently of the plot window. Singular real-axis slices have no physical period mean or finite peak-to-trough amplitude.

## Repository layout

| Files | Purpose |
|---|---|
| `gallery_template.html`, `atlas.css` | Accessible page structure and responsive styling |
| `atlas_ui.js`, `atlas_plot.js`, `atlas_cases.js` | Scene controls, interactive canvas plots and curated equations |
| `wave_math.js`, `c4_math.js`, `c4_ui.js` | Jacobi/Weierstrass profiles and the C4 constructor |
| `models.py`, `derive.py`, `boundaries.py`, `verify.py` | Exact families, exceptional fibres and algebraic checks |
| `render.py`, `profiles.py`, `journeys.py` | Interactive data, analytic profile metadata and page compilation |
| `data/` | Reproducible formulas, solved frames, CSVs and validation records |
| `physical/` | Independent Kawahara calculations and the [interactive observable-coordinate explorer](physical/kawahara_physical.html) |

The static publication artifacts in `figures/` and `g2g3_atlas.pdf` are retained as an archive. They are not loaded by the explorer. Rebuild them explicitly with `python3 render.py --figures` or `python3 c4.py --figures`.

GitHub Pages serves the repository root of `main`, with `.nojekyll`. The page embeds its mathematical data to support `file://` viewing and uses only local JavaScript and CSS assets.

For coefficient conventions, coverage, characters, and count semantics, see [Mathematical scope](docs/mathematics.md). For the physical Kawahara normalization and limiting observables, see [the independent calculation](physical/README.md).
