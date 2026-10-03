/* State, controls and linked views. Keep the mathematics in wave_math/C4. */
const $ = (id) => document.getElementById(id),
  family = $("family"),
  slider = $("position"),
  canvas = $("plot");
const fmt = AtlasPlot.fmt,
  superscript = (x) =>
    String(x).replace(/[0-9-]/g, (v) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[+v] || "⁻");
let current = 0,
  selected = 0,
  selectedBranch = null,
  preferredProfileMode = "regular",
  lastProfile = null,
  restoringScene = false;
let mode = "discover",
  activeCase = null,
  labMemory = "c4",
  journeyIndex = null,
  journeyPlaying = false,
  journeyTimer = null;
let latticeView = null,
  latticeFrame = null,
  waveFrame = null,
  pinnedZ = null,
  drag = null,
  profileCache = null;
let linkTimer = null;
let labC4Params = null,
  labParameter = null;
const model = () => ATLAS[current];
function effectiveFrame() {
  if (journeyIndex === null || model().slug !== JOURNEY.slug)
    return model().frames[+slider.value];
  const r = JOURNEY.records[journeyIndex],
    branch = r.distance === 0 ? 0 : 1;
  return {
    t: r.t,
    journey: r,
    points: [
      {
        g2: r.g2,
        g3: r.g3,
        delta: r.delta,
        branch,
        character: 1,
        operator: { 2: r.t, 0: -1 },
        constant: 0,
        profile: r.profile,
      },
    ],
  };
}
function sceneCase() {
  return (
    CASES.find((c) => c.id === activeCase && c.slug === model().slug) ||
    CASES.find((c) => c.slug === model().slug && c.slug !== "c4")
  );
}
function pulseReference() {
  if (model().kind === "c4") {
    const s = C4Panel.state();
    return C4.presetKey(s.coeff) === "lax" &&
      s.coeff.a2 === 0 &&
      s.coeff.a0 === -1 &&
      s.coeff.a00 === 0
      ? CASES.find((c) => c.id === "lax").pulse
      : null;
  }
  const c = sceneCase();
  if (!c?.pulse) return null;
  return c.pulse;
}
function setMode(value, keep = false) {
  if (mode === value && !keep) return;
  if (mode === "lab" && value !== "lab") {
    labMemory = model().slug;
    labParameter = effectiveFrame().t;
    if (model().kind === "c4") {
      const params = new URLSearchParams();
      C4Panel.toParams(params);
      labC4Params = params.toString();
    }
  }
  mode = value;
  $("discover-panel").hidden = mode !== "discover";
  $("lab-panel").hidden = mode !== "lab";
  for (const name of ["discover", "lab"]) {
    const tab = $(name + "-tab");
    tab.setAttribute("aria-selected", String(mode === name));
    tab.tabIndex = mode === name ? 0 : -1;
  }
  if (!keep) {
    if (mode === "lab") {
      if (labMemory === "c4" && labC4Params)
        C4Panel.restore(new URLSearchParams(labC4Params));
      choose(
        ATLAS.findIndex((m) => m.slug === labMemory),
        { keepMode: true },
      );
      if (labMemory !== "c4" && labParameter !== null)
        toParameter(labParameter);
    } else selectCase(activeCase || "kawahara");
  } else if (model().kind === "c4") {
    mode === "lab" ? C4Panel.show() : C4Panel.hide();
  }
}
function stopJourney() {
  journeyPlaying = false;
  if (journeyTimer) clearTimeout(journeyTimer);
  journeyTimer = null;
  $("journey-play").textContent = "Play transition";
  $("journey-play").setAttribute("aria-pressed", "false");
}
function exitJourney() {
  stopJourney();
  journeyIndex = null;
  $("journey-position").value = 0;
}
function equation(pt, m) {
  if (m.kind === "c4") return C4.equation(effectiveFrame().coefficients);
  if (!pt)
    return "No real lattice in the constructed sector at this parameter.";
  const variable = m.variable || "u",
    scale = m.display_scale || 1;
  const derivative = (j) =>
    variable +
    (j === 1
      ? "′"
      : j === 2
        ? "″"
        : j === 3
          ? "‴"
          : j === 4
            ? "⁗"
            : "⁽" + superscript(j) + "⁾");
  let out = derivative(m.p);
  for (let j = m.p - 1; j >= 0; j--) {
    const a = pt.operator[String(j)] || 0;
    if (Math.abs(a) < 1e-13) continue;
    out +=
      " " +
      (a < 0 ? "−" : "+") +
      " " +
      (Math.abs(a) === 1 ? "" : fmt(Math.abs(a)) + " ") +
      (j ? derivative(j) : variable);
  }
  const nonlinear = pt.physical ? 0.5 : -m.K / Math.pow(scale, m.n - 1),
    constant = (pt.constant || 0) * scale;
  out +=
    " " +
    (nonlinear < 0 ? "−" : "+") +
    " " +
    (Math.abs(nonlinear) === 0.5
      ? variable + "²/2"
      : (Math.abs(nonlinear) === 1 ? "" : fmt(Math.abs(nonlinear)) + " ") +
        variable +
        superscript(m.n));
  if (constant)
    out += " " + (constant < 0 ? "−" : "+") + " " + fmt(Math.abs(constant));
  return out + " = 0";
}
function draw() {
  if (!model()) return;
  latticeFrame = AtlasPlot.lattice(
    canvas,
    model(),
    effectiveFrame().points,
    selected,
    latticeView || model().limits,
  );
}
function drawWave() {
  const m = model(),
    frame = effectiveFrame(),
    point = frame.points[selected],
    wave = $("wave");
  if (!point) {
    lastProfile = null;
    waveFrame = null;
    profileCache = null;
    const c = wave.getContext("2d");
    wave.width = wave.clientWidth * (devicePixelRatio || 1);
    wave.height = wave.clientHeight * (devicePixelRatio || 1);
    c.scale(devicePixelRatio || 1, devicePixelRatio || 1);
    c.fillStyle = "#596c78";
    c.font = "12px system-ui";
    c.textAlign = "center";
    c.fillText(
      "No real lattice in this sector",
      wave.clientWidth / 2,
      wave.clientHeight / 2,
    );
    $("waveLegend").replaceChildren();
    $("profileStats").replaceChildren();
    $("profileNote").textContent =
      "These coefficients have no real lattice in the displayed construction. Complex lattices may exist outside this real plane.";
    $("wave-status").textContent = "No real lattice";
    $("wave-status").className = "wave-status singular";
    hideWaveHover();
    return;
  }
  const spec = point.profile,
    actual = WaveMath.actualMode(spec, preferredProfileMode);
  $("profileMode").options[0].disabled = spec.geometry.type !== "rect";
  $("profileMode").value = actual;
  $("wave-variable").textContent = spec.variable + "(z)";
  const fixed = $("window-mode").value === "fixed",
    span = C4.parseNumber($("window-span").value);
  let width = fixed && span > 0 ? span : null;
  const limited =
    width !== null && spec.geometry.period && width > 20 * spec.geometry.period;
  if (limited) {
    width = 20 * spec.geometry.period;
    $("window-span").value = fmt(width);
  }
  $("span-field").hidden = !fixed;
  const p =
    profileCache?.spec === spec &&
    profileCache.mode === actual &&
    profileCache.width === width
      ? profileCache.p
      : WaveMath.generate(spec, actual, 2, 1201, width);
  profileCache = { spec, mode: actual, width, p };
  if (frame.journey && actual === "regular") {
    p.period = frame.journey.period;
    p.amplitude = frame.journey.amplitude;
    p.mean = frame.journey.mean;
  }
  lastProfile = p;
  const reference = pulseReference();
  $("pulse-reference-field").hidden = !reference;
  waveFrame = AtlasPlot.wave(
    wave,
    p,
    spec,
    $("pulse-reference").checked && actual === "regular" ? reference : null,
    pinnedZ,
  );
  $("waveLegend").innerHTML = p.anyImag
    ? `<span><i class="swatch"></i>Re ${spec.variable}</span><span><i class="swatch imag"></i>Im ${spec.variable}</span>`
    : `<span><i class="swatch"></i>${spec.variable}(z)</span>`;
  if (reference && $("pulse-reference").checked && actual === "regular")
    $("waveLegend").insertAdjacentHTML(
      "beforeend",
      '<span><i class="swatch reference"></i>solitary limit</span>',
    );
  if (p.mean !== null)
    $("waveLegend").insertAdjacentHTML(
      "beforeend",
      '<span><i class="swatch mean"></i>mean</span>',
    );
  let note = p.empty
    ? "No finite samples in this window. Widen the window or choose the pole-free slice."
    : p.constant
      ? "Constant limit of the pole-free slice. A constant has no primitive period. The meromorphic real-axis slice shows the nonconstant degeneration."
      : p.mode === "axis"
        ? "Meromorphic real-axis profile. Vertical dashes mark poles; the vertical range is clipped. Period averages and physical amplitudes are undefined on this singular slice."
        : p.anyImag
          ? "Pole-free complex-valued solution: both real and imaginary parts are shown. Real lattice invariants do not imply a real-valued profile."
          : spec.geometry.period === null
            ? "Localized hyperbolic degeneration, displayed in a finite spatial window. The full real period is infinite."
            : "Smooth real periodic profile. The faint horizontal line marks its full-period average.";
  if (preferredProfileMode === "regular" && spec.geometry.type === "one")
    note = "This lattice has no bounded real oval. " + note;
  if (m.slug === "n2_p3_ks" && (frame.t < -18 - 1e-10 || frame.t > -8 + 1e-10))
    note += " The smooth periodic KS interval is −18 < C < −8.";
  if (limited)
    note +=
      " The window was limited to twenty lattice periods to keep individual waves resolved.";
  $("profileNote").textContent = note;
  $("wave-status").textContent = p.empty
    ? "No finite samples"
    : p.constant
      ? "Constant limit"
      : p.mode === "axis"
        ? "Meromorphic · poles"
        : p.anyImag
          ? "Smooth slice · complex"
          : p.period
            ? "Smooth periodic · real"
            : "Solitary limit · real";
  $("wave-status").className =
    "wave-status" + (p.mode === "axis" ? " singular" : "");
  const stat = (label, value) =>
    `<div><span>${label}</span><strong>${value}</strong></div>`;
  let stats = stat(
    "Full real period",
    p.constant
      ? "No primitive period"
      : p.period
        ? fmt(p.period)
        : spec.geometry.type === "rational"
          ? "Rational limit"
          : "∞",
  );
  if (p.amplitude !== null)
    stats += stat("Peak-to-trough amplitude", fmt(p.amplitude));
  if (p.mean !== null)
    stats += stat(p.constant ? "Constant value" : "Period mean", fmt(p.mean));
  if (spec.geometry.type === "rect" && !spec.geometry.node)
    stats += stat("Jacobi parameter m", fmt(spec.geometry.m));
  $("profileStats").innerHTML = stats;
  updateProbe();
}
function update() {
  const m = model(),
    frame = effectiveFrame();
  let found = frame.points.findIndex((p) => p.branch === selectedBranch);
  if (found < 0)
    found = Math.max(
      0,
      frame.points.findIndex((p) => p.delta > 0),
    );
  selected = found;
  if (frame.points[selected]) selectedBranch = frame.points[selected].branch;
  $("parameter-min").textContent = fmt(m.interval[0]);
  $("parameter-max").textContent = fmt(m.interval[1]);
  if (m.kind === "c4") C4Panel.onFrame(m, frame, frame.points[selected]);
  $("paramValue").value = fmt(frame.t);
  $("parameter-input").value = fmt(frame.t);
  $("parameter-input").title = "Selected parameter: " + String(frame.t);
  slider.setAttribute("aria-valuetext", fmt(frame.t));
  const n = frame.points.length;
  $("count").textContent = frame.journey
    ? "Smooth Kawahara branch shown · the other real lattice branch is available in the parameter sweep"
    : m.kind === "c4"
      ? C4Panel.count(frame)
      : m.fixed_equation
        ? "One fixed equation · a continuous family of energy levels"
        : m.count_kind === "equations"
          ? `${n} real compatible coefficient choice${n === 1 ? "" : "s"} at this parameter`
          : `${n} distinct real lattice solution${n === 1 ? "" : "s"} at this parameter`;
  const pointLabel = (p, i) =>
    m.kind === "c4"
      ? C4Panel.pointLabel(p)
      : m.slug === "n2_p4_kawahara"
        ? p.profile.geometry.node
          ? "Merged limit"
          : p.branch === 1
            ? "Smooth branch"
            : "Companion branch"
        : `Wave ${i + 1}${p.character > 1 ? " · χ" + superscript(p.character) : ""}`;
  $("branch-choices").innerHTML = frame.points
    .map(
      (p, i) =>
        `<button data-index="${i}" aria-pressed="${i === selected}">${pointLabel(p, i)}</button>`,
    )
    .join("");
  $("points").innerHTML = n
    ? "<table><thead><tr><th>Wave</th><th>g₂</th><th>g₃</th><th>Δ</th></tr></thead><tbody>" +
      frame.points
        .map((p, i) => {
          const scale = Math.max(
              Math.abs(p.g2 ** 3),
              Math.abs(27 * p.g3 * p.g3),
              1e-40,
            ),
            sign =
              frame.journey && frame.journey.distance !== 0
                ? "> 0"
                : Math.abs(p.delta) / scale < 1e-7
                  ? "0"
                  : p.delta > 0
                    ? "> 0"
                    : "< 0";
          return `<tr><td><button class="point" data-index="${i}">${i === selected ? "● " : ""}${pointLabel(p, i)}</button></td><td>${fmt(p.g2)}</td><td>${fmt(p.g3)}</td><td>${sign}</td></tr>`;
        })
        .join("") +
      "</tbody></table>"
    : "";
  for (const root of [$("points"), $("branch-choices")])
    root.querySelectorAll("button").forEach(
      (b) =>
        (b.onclick = () => {
          selectedBranch = frame.points[+b.dataset.index].branch;
          update();
        }),
    );
  if (m.kind === "c4") C4Panel.renderNotes(frame);
  $("equation").textContent = equation(frame.points[selected], m);
  $("shortcuts")
    .querySelectorAll("[data-parameter]")
    .forEach((b) =>
      b.setAttribute(
        "aria-pressed",
        String(
          Math.abs(frame.t - Number(b.dataset.parameter)) < 1e-10 &&
            journeyIndex === null,
        ),
      ),
    );
  if (frame.journey) {
    const r = frame.journey;
    $("journey-value").textContent = r.period
      ? "L = " + fmt(r.period)
      : "L → ∞";
    $("journey-note").textContent = r.distance
      ? "a = −13/6 + " +
        r.distance.toExponential(2) +
        ". Extra precision preserves the finite period; the orange curve is the exact solitary limit. This animates a parameter, not time evolution."
      : "a = −13/6 exactly. v(z) = 2 − (35/12) sech⁴(z/(2√6)). The period is infinite; the localized pulse has no period mean.";
  } else {
    $("journey-value").textContent = "";
    $("journey-note").textContent =
      "Use a common spatial scale to watch pulses separate. The orange curve is the exact solitary limit. This is a parameter transition.";
  }
  $("lab-coefficients").innerHTML =
    m.kind === "c4"
      ? ""
      : Object.entries(frame.points[selected]?.operator || {})
          .sort((a, b) => Number(b[0]) - Number(a[0]))
          .map(
            ([j, v]) =>
              `<span${m.slug === "n2_p4_kawahara" && j === "2" ? ' class="active"' : ""}>a${superscript(j)} = ${fmt(v)}</span>`,
          )
          .join("") +
        (m.slug === "n2_p3_ks"
          ? '<span class="active">C = ' + fmt(frame.t) + "</span>"
          : "");
  draw();
  drawWave();
  saveSceneLink();
}
function toParameter(value) {
  exitJourney();
  slider.value = AtlasPlot.closestFrame(model(), value);
  $("parameter-status").textContent = "";
  update();
}
function choose(index, options = {}) {
  if (index < 0 || index >= ATLAS.length) return;
  exitJourney();
  current = index;
  selected = 0;
  selectedBranch = null;
  preferredProfileMode = "regular";
  latticeView = null;
  pinnedZ = null;
  profileCache = null;
  activeCase = options.caseId || null;
  const m = model();
  if (!options.keepMode && m.kind === "c4" && !activeCase) setMode("lab", true);
  if (
    mode === "lab" &&
    !Array.from($("lab-family").options).some((o) => o.value === m.slug)
  )
    setMode("discover", true);
  family.value = String(index);
  const c = sceneCase();
  $("currentTitle").textContent =
    mode === "lab" && m.kind === "c4"
      ? "C4 coefficient explorer"
      : c?.name || m.title;
  $("case-kicker").textContent =
    mode === "lab"
      ? "Equation lab / Supported construction"
      : c?.kicker || "Full atlas / Pure-power coefficient slice";
  $("case-description").textContent =
    mode === "lab" && m.kind === "c4"
      ? "Choose coefficients, including zeros, and follow the affine elliptic representatives. A free coordinate changes the wave within a fixed equation."
      : c?.description ||
        "Explore this declared coefficient slice. Select a curve or branch to see its corresponding travelling profile.";
  $("badge").textContent =
    m.kind === "c4"
      ? "Affine-℘ sector · order 4"
      : `n = ${m.n} · p = ${m.p} · pole order ${m.q}`;
  $("equation-context").textContent =
    mode === "lab" && m.kind === "c4"
      ? "v = α℘ + β · real coefficients, real lattice invariants · solved amplitudes may be complex"
      : c?.context ||
        "Displayed normalization · exact polynomial compatibility conditions";
  $("paramLabel").textContent = m.parameter;
  slider.max = m.frames.length - 1;
  slider.value =
    m.kind === "c4" ? m.selectedIndex : Math.floor(m.frames.length * 0.5);
  $("parameter-min").textContent = fmt(m.interval[0]);
  $("parameter-max").textContent = fmt(m.interval[1]);
  $("parameter-status").textContent = "";
  $("window-mode").value = "auto";
  $("window-span").value = c?.width || 60;
  $("pulse-reference").checked = false;
  $("pulse-journey").hidden = m.slug !== JOURNEY.slug || mode !== "discover";
  $("lab-family").value = m.slug;
  $("lab-slice").hidden = m.kind === "c4";
  if (m.kind === "c4") {
    mode === "lab" ? C4Panel.show() : C4Panel.hide();
    C4Panel.shortcuts();
  } else {
    C4Panel.hide();
    $("notes").innerHTML = m.notes.map((n) => `<p>${n}</p>`).join("");
    const shortcuts =
      m.slug === "n2_p4_kawahara"
        ? [
            [-13 / 6, "Solitary pulse"],
            [-1.5, "Periodic train"],
            [0, "Zero dispersion"],
            [13 / 6, "Constant limit"],
          ]
        : m.slug === "n2_p3_ks"
          ? [
              [-18, "Asymmetric pulse"],
              [-13, "Smooth periodic"],
              [-8, "Constant limit"],
            ]
          : m.slug === "n2_p2_kdv"
            ? [
                [-1, "Solitary pulse"],
                [-0.95, "Near the pulse"],
                [0, "Cnoidal wave"],
                [1, "Constant limit"],
              ]
            : m.slug === "n3_p4_cubic_positive"
              ? [
                  [-2.5, "Solitary limit"],
                  [0, "Compare branches"],
                  [2.5, "Upper contact"],
                ]
              : [];
    $("shortcuts").innerHTML = shortcuts
      .map(([t, label]) => `<button data-parameter="${t}">${label}</button>`)
      .join("");
    $("shortcuts")
      .querySelectorAll("button")
      .forEach(
        (b) =>
          (b.onclick = () => {
            preferredProfileMode = "regular";
            toParameter(+b.dataset.parameter);
          }),
      );
    if (m.slug === JOURNEY.slug && mode === "discover") {
      const button = document.createElement("button");
      button.textContent = "Follow the pulse limit ↓";
      button.onclick = () => {
        selectJourney(0);
        $("pulse-journey").scrollIntoView({ block: "end", behavior: "smooth" });
        $("journey-position").focus({ preventScroll: true });
      };
      $("shortcuts").append(button);
    }
  }
  $("lab-constraint").textContent =
    m.slug === "n2_p4_kawahara"
      ? "Only the even dispersion coefficient a varies. Odd derivatives and the additive constant are zero; the linear coefficient is −1."
      : m.slug === "n2_p3_ks"
        ? "The elliptic slice requires b² = 16μν. Here b = 4 and μ = ν = 1. Only the integration constant C varies."
        : m.slug === "n2_p2_kdv"
          ? "One fixed conservative equation: u″ − u − 6u² = 0. The energy coordinate h chooses a wave; it does not change this equation."
          : "";
  $("case-sources").replaceChildren(
    ...(c?.sources || []).map(([name, url]) => {
      const a = document.createElement("a");
      a.href = url;
      a.textContent = name + " ↗";
      return a;
    }),
  );
  document
    .querySelectorAll(".case-card")
    .forEach((b) =>
      b.setAttribute(
        "aria-pressed",
        String(c?.id === b.dataset.case && mode === "discover"),
      ),
    );
  update();
  if (m.kind !== "c4" && c) toParameter(c.parameter);
}
function selectCase(id) {
  const c = CASES.find((c) => c.id === id) || CASES[0];
  setMode("discover", true);
  choose(
    ATLAS.findIndex((m) => m.slug === c.slug),
    { caseId: c.id, keepMode: true },
  );
  if (c.preset === "oneGap") C4Panel.oneGap();
  else if (c.preset) C4Panel.setPreset(c.preset);
  if (c.slug === "c4") {
    toParameter(c.parameter);
    // Curated scenes keep their own equation core; broader core changes belong
    // to the lab, where all coefficients and restrictions are visible.
    if (c.id === "lax") {
      $("c4-two-branches")?.remove();
    } else $("shortcuts").replaceChildren();
    const button = document.createElement("button");
    button.textContent = "Edit this equation →";
    button.onclick = () => {
      const t = effectiveFrame().t;
      setMode("lab", true);
      choose(current, { keepMode: true });
      toParameter(t);
    };
    $("shortcuts").append(button);
  }
}
function selectJourney(index) {
  stopJourney();
  journeyIndex = Math.max(
    0,
    Math.min(JOURNEY.records.length - 1, Math.round(index)),
  );
  $("journey-position").value = journeyIndex;
  slider.value = AtlasPlot.closestFrame(
    model(),
    JOURNEY.records[journeyIndex].t,
  );
  preferredProfileMode = "regular";
  $("window-mode").value = "fixed";
  $("window-span").value = JOURNEY.width;
  $("pulse-reference").checked = true;
  update();
}
function playJourney() {
  if (journeyPlaying) {
    stopJourney();
    return;
  }
  if (journeyIndex === null || journeyIndex === JOURNEY.records.length - 1)
    selectJourney(0);
  journeyPlaying = true;
  $("journey-play").textContent = "Pause";
  $("journey-play").setAttribute("aria-pressed", "true");
  const tick = () => {
    if (!journeyPlaying) return;
    const i = journeyIndex + 1;
    journeyIndex = Math.min(i, JOURNEY.records.length - 1);
    $("journey-position").value = journeyIndex;
    slider.value = AtlasPlot.closestFrame(
      model(),
      JOURNEY.records[journeyIndex].t,
    );
    update();
    if (i >= JOURNEY.records.length - 1) stopJourney();
    else journeyTimer = setTimeout(tick, 70);
  };
  journeyTimer = setTimeout(tick, 70);
}
function saveSceneLink(immediate = false) {
  if (restoringScene) return;
  const m = model(),
    frame = effectiveFrame(),
    p = new URLSearchParams({
      family: m.slug,
      parameter: String(frame.t),
      view: preferredProfileMode,
      mode,
    });
  if (activeCase) p.set("case", activeCase);
  if (frame.points[selected])
    p.set("branch", String(frame.points[selected].branch));
  if (m.kind === "c4") C4Panel.toParams(p);
  if (journeyIndex !== null) p.set("journey", String(journeyIndex));
  if ($("window-mode").value === "fixed") {
    p.set("window", "fixed");
    p.set("width", $("window-span").value);
  }
  if ($("pulse-reference").checked) p.set("reference", "1");
  if (pinnedZ !== null) p.set("probe", String(pinnedZ));
  if (latticeView) p.set("zoom", latticeView.flat().map(String).join(","));
  // Slider drags and the transition can produce hundreds of updates. Browsers
  // throttle History calls; write once after activity settles, and flush when
  // the user copies a link. Never let a delayed old scene overwrite a new hash.
  if (linkTimer) clearTimeout(linkTimer);
  const write = () => {
    linkTimer = null;
    const hash = "#" + p.toString();
    if (location.hash === hash) return;
    try {
      history.replaceState(null, "", hash);
    } catch (error) {
      /* Some offline viewers disallow history changes. */
    }
  };
  if (immediate) write();
  else linkTimer = setTimeout(write, 200);
}
function restoreSceneLink() {
  if (linkTimer) clearTimeout(linkTimer);
  linkTimer = null;
  const p = new URLSearchParams(location.hash.slice(1)),
    index = ATLAS.findIndex((m) => m.slug === p.get("family"));
  if (index < 0) return false;
  restoringScene = true;
  try {
    const c = CASES.find(
      (c) => c.id === p.get("case") && c.slug === ATLAS[index].slug,
    );
    setMode(
      p.get("mode") === "lab" || (ATLAS[index].kind === "c4" && !c)
        ? "lab"
        : "discover",
      true,
    );
    if (ATLAS[index].kind === "c4") C4Panel.restore(p);
    choose(index, { caseId: c?.id, keepMode: true });
    if (p.has("parameter") && Number.isFinite(Number(p.get("parameter"))))
      toParameter(Number(p.get("parameter")));
    if (model().slug === JOURNEY.slug && /^\d+$/.test(p.get("journey") || ""))
      selectJourney(Number(p.get("journey")));
    preferredProfileMode = p.get("view") === "axis" ? "axis" : "regular";
    if (p.has("branch") && Number.isFinite(Number(p.get("branch"))))
      selectedBranch = Number(p.get("branch"));
    if (p.get("window") === "fixed") {
      $("window-mode").value = "fixed";
      const span = C4.parseNumber(p.get("width"));
      if (span > 0) $("window-span").value = span;
    }
    $("pulse-reference").checked = p.get("reference") === "1";
    if (p.has("probe")) pinnedZ = C4.parseNumber(p.get("probe"));
    const zoom = (p.get("zoom") || "").split(",").map(Number);
    if (
      zoom.length === 4 &&
      zoom.every(Number.isFinite) &&
      zoom[0] < zoom[1] &&
      zoom[2] < zoom[3]
    )
      latticeView = [
        [zoom[0], zoom[1]],
        [zoom[2], zoom[3]],
      ];
    update();
  } finally {
    restoringScene = false;
  }
  saveSceneLink();
  return true;
}
function zoom(factor, center = null) {
  const limits = latticeView || model().limits;
  latticeView = limits.map((a, j) => {
    const mid = center ? center[j] : (a[0] + a[1]) / 2;
    return [mid + (a[0] - mid) * factor, mid + (a[1] - mid) * factor];
  });
  draw();
  saveSceneLink();
}
function fitLattice() {
  const m = model(),
    tracks =
      m.kind === "c4" && m.focus
        ? m.tracks.filter((t) => t.branch === selectedBranch)
        : m.tracks;
  const points = tracks
    .flatMap((t) => t.points.filter(Boolean))
    .concat(effectiveFrame().points.map((p) => [p.g2, p.g3]));
  if (!points.length) {
    latticeView = null;
  } else
    latticeView = [0, 1].map((axis) => {
      const vals = points.map((p) => p[axis]).filter(Number.isFinite),
        lo = Math.min(...vals),
        hi = Math.max(...vals),
        span =
          hi - lo || Math.max(Math.abs(hi) * 0.3, axis === 0 ? 0.01 : 0.001);
      return [lo - span * 0.08, hi + span * 0.08];
    });
  draw();
  saveSceneLink();
}
function position(event, element) {
  const r = element.getBoundingClientRect();
  return [event.clientX - r.left, event.clientY - r.top];
}
function tooltip(id, dotId, x, y, html, wrap) {
  const tip = $(id),
    dot = $(dotId);
  tip.innerHTML = html;
  tip.hidden = false;
  tip.style.left =
    Math.max(3, Math.min(wrap.clientWidth - tip.offsetWidth - 3, x + 13)) +
    "px";
  tip.style.top =
    Math.max(3, Math.min(wrap.clientHeight - tip.offsetHeight - 3, y - 13)) +
    "px";
  dot.hidden = false;
  dot.style.left = x + "px";
  dot.style.top = y + 4 + "px";
}
function hideLatticeHover() {
  $("lattice-tooltip").hidden = true;
  $("lattice-dot").hidden = true;
}
function hideWaveHover() {
  $("wave-tooltip").hidden = true;
  $("wave-dot").hidden = true;
  $("wave-crosshair").hidden = true;
}
canvas.onpointerdown = (e) => {
  const [x, y] = position(e, canvas);
  if (!latticeFrame?.inside(x, y)) return;
  drag = {
    x,
    y,
    limits: (latticeView || model().limits).map((a) => a.slice()),
    moved: false,
    pointerId: e.pointerId,
  };
  canvas.setPointerCapture(e.pointerId);
};
canvas.onpointermove = (e) => {
  const [x, y] = position(e, canvas);
  if (drag && e.pointerId === drag.pointerId) {
    const dx = x - drag.x,
      dy = y - drag.y;
    if (Math.hypot(dx, dy) > 4) drag.moved = true;
    if (drag.moved) {
      latticeView = drag.limits.map((a, j) => {
        const shift =
          (j === 0 ? -dx / latticeFrame.w : dy / latticeFrame.h) *
          (a[1] - a[0]);
        return [a[0] + shift, a[1] + shift];
      });
      draw();
      hideLatticeHover();
      canvas.style.cursor = "grabbing";
    }
    return;
  }
  const hit = AtlasPlot.pick(latticeFrame, x, y);
  canvas.style.cursor = hit ? "pointer" : "grab";
  if (!hit) {
    hideLatticeHover();
    return;
  }
  tooltip(
    "lattice-tooltip",
    "lattice-dot",
    hit.x,
    hit.y,
    `<span class="tooltip-label">Select this solved point</span>parameter = ${fmt(hit.t)}<br>g₂ = ${fmt(hit.point.g2)}<br>g₃ = ${fmt(hit.point.g3)}`,
    $("lattice-wrap"),
  );
};
canvas.onpointerup = (e) => {
  if (!drag) return;
  const moved = drag.moved;
  drag = null;
  if (canvas.hasPointerCapture(e.pointerId))
    canvas.releasePointerCapture(e.pointerId);
  if (moved) {
    saveSceneLink();
    return;
  }
  const [x, y] = position(e, canvas),
    hit = AtlasPlot.pick(latticeFrame, x, y, 22);
  if (hit) {
    exitJourney();
    slider.value = hit.frameIndex;
    selectedBranch = hit.branch;
    update();
  }
};
canvas.onpointercancel = () => {
  drag = null;
};
canvas.onpointerleave = hideLatticeHover;
canvas.ondblclick = () => {
  latticeView = null;
  draw();
  saveSceneLink();
};
canvas.addEventListener(
  "wheel",
  (e) => {
    if (!e.ctrlKey && !e.altKey) return;
    e.preventDefault();
    const [x, y] = position(e, canvas);
    if (latticeFrame?.inside(x, y))
      zoom(e.deltaY > 0 ? 1.15 : 1 / 1.15, [
        latticeFrame.xvalue(x),
        latticeFrame.yvalue(y),
      ]);
  },
  { passive: false },
);
canvas.onkeydown = (e) => {
  if (
    [
      "ArrowLeft",
      "ArrowRight",
      "ArrowUp",
      "ArrowDown",
      "+",
      "=",
      "-",
      "0",
    ].includes(e.key)
  )
    e.preventDefault();
  if (e.key === "+" || e.key === "=") zoom(0.75);
  else if (e.key === "-") zoom(1 / 0.75);
  else if (e.key === "0") {
    $("zoom-fit").click();
  } else if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
    exitJourney();
    slider.value = Math.max(
      0,
      Math.min(+slider.max, +slider.value + (e.key === "ArrowLeft" ? -1 : 1)),
    );
    update();
  } else if (e.key === "ArrowUp" || e.key === "ArrowDown") {
    const pts = effectiveFrame().points;
    if (pts.length) {
      selectedBranch =
        pts[
          (selected + (e.key === "ArrowUp" ? -1 : 1) + pts.length) % pts.length
        ].branch;
      update();
    }
  }
};
function waveHover(e) {
  if (!waveFrame || !lastProfile) return;
  const [x, y] = position(e, $("wave"));
  if (!waveFrame.inside(x, y)) {
    hideWaveHover();
    return;
  }
  const z = waveFrame.xvalue(x),
    pt = effectiveFrame().points[selected],
    v = WaveMath.at(pt.profile, z, lastProfile.mode);
  if (!v) {
    hideWaveHover();
    return;
  }
  const yy = waveFrame.py(v[0]);
  if (yy < waveFrame.P.t || yy > waveFrame.P.t + waveFrame.h) {
    hideWaveHover();
    return;
  }
  tooltip(
    "wave-tooltip",
    "wave-dot",
    x,
    yy,
    `<span class="tooltip-label">Click to pin this value</span>z = ${fmt(z)}<br>${pt.profile.variable} = ${fmt(v[0])}${lastProfile.anyImag ? " + " + fmt(v[1]) + "i" : ""}`,
    $("wave-wrap"),
  );
  const cross = $("wave-crosshair");
  cross.hidden = false;
  cross.style.left = x + "px";
  cross.style.top = waveFrame.P.t + 4 + "px";
  cross.style.height = waveFrame.h + "px";
}
$("wave").onpointermove = waveHover;
$("wave").onpointerleave = hideWaveHover;
$("wave").onclick = (e) => {
  if (!waveFrame) return;
  const [x, y] = position(e, $("wave"));
  if (waveFrame.inside(x, y)) {
    pinnedZ = waveFrame.xvalue(x);
    drawWave();
    saveSceneLink();
  }
};
function updateProbe() {
  if (pinnedZ === null) {
    $("probe-input").value = "";
    $("probe-value").textContent = "—";
    return;
  }
  const pt = effectiveFrame().points[selected],
    v = pt ? WaveMath.at(pt.profile, pinnedZ, preferredProfileMode) : null;
  $("probe-input").value = fmt(pinnedZ);
  $("probe-value").textContent = v
    ? pt.profile.variable + " = " + C4.complexText(v)
    : "Pole / no finite value";
}
$("wave").onkeydown = (e) => {
  if (e.key === "Enter") {
    $("probe-input").focus();
    e.preventDefault();
  } else if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
    e.preventDefault();
    pinnedZ =
      (pinnedZ ?? 0) +
      ((e.key === "ArrowLeft" ? -1 : 1) * (lastProfile?.width || 1)) / 100;
    drawWave();
    saveSceneLink();
  } else if (e.key === "Escape") {
    $("probe-clear").click();
  }
};
function download(name, text) {
  const url = URL.createObjectURL(
      new Blob([text], { type: "text/csv;charset=utf-8" }),
    ),
    link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
$("download-wave").onclick = () => {
  if (!lastProfile) return;
  download(
    model().slug + "-wave.csv",
    "z,real,imaginary\n" +
      lastProfile.x
        .map((x, i) =>
          [x, lastProfile.re[i] ?? "", lastProfile.im[i] ?? ""].join(","),
        )
        .join("\n") +
      "\n",
  );
};
$("download-sweep").onclick = () => {
  if (model().kind === "c4") {
    download("c4-" + C4Panel.state().sweep + ".csv", C4.csv(model()));
    return;
  }
  const link = document.createElement("a");
  link.href = "data/" + model().slug + ".csv";
  link.download = model().slug + ".csv";
  link.click();
};
$("cases").innerHTML = CASES.map(
  (c, i) =>
    `<button class="case-card" data-case="${c.id}" aria-pressed="false"><span class="case-number">0${i + 1}</span><span><span class="case-title">${c.name}</span><span class="case-caption">${c.caption}</span></span><span class="case-arrow">↗</span></button>`,
).join("");
$("cases")
  .querySelectorAll("button")
  .forEach((b) => (b.onclick = () => selectCase(b.dataset.case)));
C4Panel.register(ATLAS);
ATLAS.filter((m) => m.kind !== "c4").forEach((m) => {
  const option = document.createElement("option");
  option.value = ATLAS.indexOf(m);
  option.textContent = `(${m.n},${m.p}) ${m.title}`;
  family.appendChild(option);
});
family.onchange = () => {
  setMode("discover", true);
  choose(+family.value, { keepMode: true });
};
$("discover-tab").onclick = () => setMode("discover");
$("lab-tab").onclick = () => setMode("lab");
document.querySelectorAll("[role=tab]").forEach(
  (tab) =>
    (tab.onkeydown = (e) => {
      if (e.key === "ArrowLeft" || e.key === "ArrowRight") {
        e.preventDefault();
        const name = tab.id === "discover-tab" ? "lab" : "discover";
        setMode(name);
        $(name + "-tab").focus();
      }
    }),
);
$("lab-family").onchange = () => {
  labMemory = $("lab-family").value;
  choose(
    ATLAS.findIndex((m) => m.slug === labMemory),
    { keepMode: true },
  );
};
slider.oninput = () => {
  exitJourney();
  $("parameter-status").textContent = "";
  update();
};
$("parameter-input").onchange = () => {
  const value = C4.parseNumber($("parameter-input").value);
  if (value === null) {
    $("parameter-status").textContent = "Enter a finite number or fraction.";
    return;
  }
  if (model().kind === "c4") {
    const p = new URLSearchParams();
    C4Panel.toParams(p);
    p.set("parameter", String(value));
    C4Panel.restore(p);
    C4Panel.rebuild();
  } else {
    const out = value < model().interval[0] || value > model().interval[1];
    toParameter(value);
    if (out || Math.abs(effectiveFrame().t - value) > 1e-9)
      $("parameter-status").textContent =
        "Nearest sampled parameter: " + fmt(effectiveFrame().t);
  }
};
$("profileMode").onchange = () => {
  preferredProfileMode = $("profileMode").value;
  drawWave();
  saveSceneLink();
};
for (const id of ["window-mode", "pulse-reference"])
  $(id).onchange = () => {
    drawWave();
    saveSceneLink();
  };
$("window-span").onchange = () => {
  const v = C4.parseNumber($("window-span").value);
  if (!(v > 0)) {
    $("window-span").value = 60;
  }
  drawWave();
  saveSceneLink();
};
$("probe-input").onchange = () => {
  const x = C4.parseNumber($("probe-input").value);
  if (x === null) {
    $("probe-value").textContent = "Enter a finite coordinate.";
    return;
  }
  pinnedZ = x;
  drawWave();
  saveSceneLink();
};
$("probe-clear").onclick = () => {
  pinnedZ = null;
  drawWave();
  saveSceneLink();
};
$("zoom-in").onclick = () => zoom(0.75);
$("zoom-out").onclick = () => zoom(1 / 0.75);
$("zoom-fit").onclick = fitLattice;
$("journey-position").max = JOURNEY.records.length - 1;
$("journey-position").oninput = () =>
  selectJourney(+$("journey-position").value);
$("journey-play").onclick = playJourney;
$("reset-scene").onclick = () => {
  stopJourney();
  latticeView = null;
  pinnedZ = null;
  preferredProfileMode = "regular";
  $("window-mode").value = journeyIndex === null ? "auto" : "fixed";
  if (journeyIndex !== null) $("window-span").value = JOURNEY.width;
  $("pulse-reference").checked = journeyIndex !== null;
  draw();
  drawWave();
  saveSceneLink();
};
$("home-link").onclick = (e) => {
  e.preventDefault();
  selectCase("kawahara");
};
$("share-scene").onclick = async () => {
  saveSceneLink(true);
  try {
    await navigator.clipboard.writeText(location.href);
    $("share-scene").textContent = "Link copied ✓";
    $("share-status").textContent = "Scene link copied.";
    setTimeout(() => ($("share-scene").textContent = "Copy link ↗"), 1800);
  } catch (error) {
    let input = $("share-fallback");
    if (!input) {
      input = document.createElement("input");
      input.id = "share-fallback";
      input.readOnly = true;
      input.setAttribute("aria-label", "Scene link to copy");
      input.style.cssText =
        "width:100%;margin:10px 0;padding:8px;font-size:11px";
      $("explorer").prepend(input);
    }
    input.value = location.href;
    input.focus();
    input.select();
    $("share-status").textContent = "Select and copy the scene link.";
  }
};
let resizePending = false;
new ResizeObserver(() => {
  if (resizePending) return;
  resizePending = true;
  requestAnimationFrame(() => {
    resizePending = false;
    draw();
    drawWave();
  });
}).observe(canvas);
new ResizeObserver(() => {
  requestAnimationFrame(drawWave);
}).observe($("wave"));
window.addEventListener("hashchange", restoreSceneLink);
document.addEventListener("visibilitychange", () => {
  if (document.hidden) stopJourney();
});
let table =
  "<tr><th>n \\ p</th>" +
  [2, 3, 4, 5, 6].map((p) => "<th>" + p + "</th>").join("") +
  "</tr>";
for (let n = 2; n <= 7; n++) {
  table += "<tr><th>" + n + "</th>";
  for (let p = 2; p <= 6; p++) {
    const m = ATLAS.find((m) => m.n === n && m.p === p && m.primary);
    table += m
      ? `<td class="yes"><button data-family="${ATLAS.indexOf(m)}" aria-label="Explore n ${n}, p ${p}, pole order ${p / (n - 1)}">${p / (n - 1)}</button></td>`
      : "<td>—</td>";
  }
  table += "</tr>";
}
$("coverage").innerHTML = table;
$("coverage")
  .querySelectorAll("button")
  .forEach(
    (b) =>
      (b.onclick = () => {
        setMode("discover", true);
        choose(+b.dataset.family, { keepMode: true });
        $("explorer").scrollIntoView({ behavior: "smooth" });
      }),
  );
if (!restoreSceneLink()) selectCase("kawahara");
