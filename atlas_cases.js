/* Curated scenes. Each names the actual reduction and its constructed sector. */
const CASES = [
  {
    id: "kawahara",
    slug: "n2_p4_kawahara",
    name: "Kawahara",
    caption: "A loop, then a solitary pulse",
    kicker: "01 / Fifth-order dispersion",
    parameter: -1.5,
    width: 60,
    description:
      "Two lattice branches trace a returning loop. Follow the smooth branch toward the cusp and a periodic train separates into one sech⁴ depression pulse.",
    context:
      "Travelling-wave ODE · v = −1680u · changing a changes the dispersion coefficient",
    sources: [
      ["Kawahara (1972)", "https://doi.org/10.1143/JPSJ.33.260"],
      [
        "Elliptic travelling waves · Mancas",
        "https://arxiv.org/abs/1612.07209",
      ],
      ["This exact pulse & normalization", "physical/README.md"],
    ],
    pulse: (z) => 2 - 35 / 12 / Math.cosh(z / (2 * Math.sqrt(6))) ** 4,
  },
  {
    id: "kdv",
    slug: "n2_p2_kdv",
    name: "Korteweg–de Vries",
    caption: "Cnoidal waves become a soliton",
    kicker: "02 / One equation, many waves",
    parameter: 0,
    width: 28,
    description:
      "Move along the energy line of one fixed KdV equation. As h approaches −1, the real period diverges and the cnoidal wave becomes a single sech² pulse.",
    context:
      "Travelling-wave ODE · the energy h varies; the equation stays fixed · displayed normalization u",
    sources: [
      [
        "Cnoidal waves as repeated solitary waves · Boyd",
        "https://doi.org/10.1137/0144066",
      ],
      ["Jacobi–Weierstrass relations · DLMF", "https://dlmf.nist.gov/23.6#ii"],
    ],
    pulse: (z) => -0.25 / Math.cosh(z / 2) ** 2,
  },
  {
    id: "ks",
    slug: "n2_p3_ks",
    name: "Dispersive KS",
    caption: "An asymmetric periodic wave",
    kicker: "03 / Kuramoto–Sivashinsky with dispersion",
    parameter: -13,
    width: 28,
    description:
      "An asymmetric smooth wave lives between C = −18 and −8. This dispersive KS slice satisfies the elliptic compatibility condition b² = 16μν. Its lower endpoint is an asymmetric localized pulse.",
    context:
      "Integrated travelling-wave ODE · ν = μ = 1, b = 4 · includes third-order PDE dispersion",
    sources: [
      [
        "Meromorphic KS travelling waves · Eremenko",
        "https://arxiv.org/abs/nlin/0504053",
      ],
    ],
    pulse: (z) => -6 + (15 / Math.cosh(z / 2) ** 2) * (1 - Math.tanh(z / 2)),
  },
  {
    id: "lax",
    slug: "c4",
    preset: "oneGap",
    name: "Lax fifth-order KdV",
    caption: "A free background, a sech² pulse",
    kicker: "04 / An integrable hierarchy",
    parameter: 0,
    width: 28,
    description:
      "A one-gap affine-Weierstrass family shares one fixed equation. The free background β moves its lattice along a cubic curve; β = 1/6 gives the exact pulse ½ sech²(z/2).",
    context: "C4 reduction · A = 10, B = 5, C = 10 · a₂ = 0, a₀ = −1, a₀₀ = 0",
    sources: [
      [
        "Integrals of nonlinear equations · Lax",
        "https://doi.org/10.1002/cpa.3160210503",
      ],
      ["C4 symbolic substitution", "data/c4_formulas.json"],
    ],
    pulse: (z) => 0.5 / Math.cosh(z / 2) ** 2,
  },
  {
    id: "sk",
    slug: "c4",
    preset: "sk",
    name: "Sawada–Kotera",
    caption: "Two affine elliptic branches",
    kicker: "05 / An integrable fifth-order equation",
    parameter: -2,
    width: 28,
    description:
      "The Sawada–Kotera core has two possible double-pole balances. Choose a branch to compare its lattice and profile, then take the coefficients into the equation lab.",
    context: "C4 reduction · A = 15, B = 0, C = 15 · affine-℘ sector",
    sources: [
      [
        "Cnoidal and travelling waves of SK & KK",
        "https://doi.org/10.1016/S0960-0779(02)00162-5",
      ],
      [
        "SK & KK integrable systems · Wang & Zhu",
        "https://arxiv.org/abs/2307.08196",
      ],
    ],
  },
  {
    id: "kk",
    slug: "c4",
    preset: "kk",
    name: "Kaup–Kupershmidt",
    caption: "Two balances, two wave shapes",
    kicker: "06 / Nonlinear derivative balance",
    parameter: -1,
    width: 36,
    description:
      "Different leading coefficients can produce very different profiles. Compare both affine branches on their real slices, or edit the nonlinear derivative terms yourself.",
    context: "C4 reduction · A = 10, B = 15/2, C = 20/3 · affine-℘ sector",
    sources: [
      [
        "Cnoidal and travelling waves of SK & KK",
        "https://doi.org/10.1016/S0960-0779(02)00162-5",
      ],
      [
        "SK & KK integrable systems · Wang & Zhu",
        "https://arxiv.org/abs/2307.08196",
      ],
    ],
  },
  {
    id: "swift",
    slug: "n3_p4_cubic_positive",
    name: "Swift–Hohenberg reduction",
    caption: "Ordinary and twisted geometry",
    kicker: "07 / Cubic fourth-order balance",
    parameter: 0,
    width: 36,
    description:
      "Two character branches belong to the same even cubic operator. Compare their real slices and follow their coalescence at a = ±5/2. The displayed stationary reduction uses the stated cubic sign.",
    context: "Stationary fourth-order reduction · u⁗ + a u″ + u − 120u³ = 0",
    sources: [
      [
        "Hydrodynamic fluctuations · Swift & Hohenberg",
        "https://doi.org/10.1103/PhysRevA.15.319",
      ],
    ],
    pulse: (z) => -0.125 / Math.cosh(z / (2 * Math.sqrt(2))) ** 2,
  },
];
