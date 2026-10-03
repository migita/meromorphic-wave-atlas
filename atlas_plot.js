/* Canvas plots with shared scientific exponents and screen-space picking.
 * The selected frame is always an actual solved frame, never an interpolated
 * lattice falsely presented as a solution of an equation.
 */
const AtlasPlot = (() => {
  const palette = [
    [56, 133, 126],
    [87, 139, 172],
    [155, 120, 163],
    [201, 146, 100],
  ];
  const ink = "#253c3e",
    muted = "#77837d",
    grid = "#e8ebe2";
  const fmt = (x) =>
    x == null
      ? "—"
      : Number(x)
          .toPrecision(6)
          .replace(/\.?0+(e|$)/, "$1");
  const superText = (x) =>
    String(x).replace(/[0-9-]/g, (v) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[+v] || "⁻");
  function color(t, m) {
    const f =
        Math.max(
          0,
          Math.min(
            0.999999,
            (t - m.interval[0]) / (m.interval[1] - m.interval[0]),
          ),
        ) * 3,
      i = Math.floor(f);
    return `rgb(${palette[i].map((v, j) => Math.round(v * (1 - (f - i)) + palette[i + 1][j] * (f - i))).join(",")})`;
  }
  function ticks(a, b, n = 4) {
    const raw = (b - a) / n,
      p = 10 ** Math.floor(Math.log10(raw)),
      f = raw / p,
      step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p,
      out = [];
    for (let v = Math.ceil(a / step) * step; v <= b + step * 1e-7; v += step) {
      out.push(Math.abs(v) < step * 1e-7 ? 0 : v);
      if (out.length > 20) break;
    }
    return out;
  }
  function scale(a, b) {
    const v = Math.max(Math.abs(a), Math.abs(b));
    if (!v) return 0;
    const e = Math.floor(Math.log10(v));
    return e <= -3 || e >= 4 ? e : 0;
  }
  function frame(canvas, limits, xlabel, ylabel) {
    const rect = canvas.getBoundingClientRect(),
      W = rect.width,
      H = rect.height,
      ratio = devicePixelRatio || 1;
    canvas.width = Math.max(1, Math.round(W * ratio));
    canvas.height = Math.max(1, Math.round(H * ratio));
    const c = canvas.getContext("2d");
    c.setTransform(ratio, 0, 0, ratio, 0, 0);
    c.clearRect(0, 0, W, H);
    const P = { l: 52, r: 14, t: 26, b: 48 },
      w = Math.max(1, W - P.l - P.r),
      h = Math.max(1, H - P.t - P.b),
      [[x0, x1], [y0, y1]] = limits;
    const px = (x) => P.l + ((x - x0) / (x1 - x0)) * w,
      py = (y) => P.t + ((y1 - y) / (y1 - y0)) * h;
    const tx = ticks(x0, x1, W < 340 ? 3 : 4),
      ty = ticks(y0, y1, 4),
      ex = scale(x0, x1),
      ey = scale(y0, y1);
    c.strokeStyle = grid;
    c.lineWidth = 0.7;
    c.beginPath();
    for (const x of tx) {
      c.moveTo(px(x), P.t);
      c.lineTo(px(x), P.t + h);
    }
    for (const y of ty) {
      c.moveTo(P.l, py(y));
      c.lineTo(P.l + w, py(y));
    }
    c.stroke();
    return {
      c,
      W,
      H,
      P,
      w,
      h,
      limits,
      px,
      py,
      tx,
      ty,
      ex,
      ey,
      xlabel,
      ylabel,
      inside: (x, y) => x >= P.l && x <= P.l + w && y >= P.t && y <= P.t + h,
      xvalue: (x) => x0 + ((x - P.l) / w) * (x1 - x0),
      yvalue: (y) => y1 - ((y - P.t) / h) * (y1 - y0),
    };
  }
  function axes(f) {
    const { c, W, H, P, w, h, px, py, tx, ty, ex, ey, xlabel, ylabel } = f;
    c.strokeStyle = "#cbd3c6";
    c.lineWidth = 1;
    c.beginPath();
    c.moveTo(P.l, P.t);
    c.lineTo(P.l, P.t + h);
    c.lineTo(P.l + w, P.t + h);
    c.stroke();
    c.font = "10px system-ui";
    c.fillStyle = muted;
    c.textAlign = "center";
    for (const v of tx) c.fillText(fmt(v / 10 ** ex), px(v), P.t + h + 18);
    c.textAlign = "right";
    for (const v of ty) c.fillText(fmt(v / 10 ** ey), P.l - 8, py(v) + 3);
    c.font = "12px Georgia";
    c.fillStyle = ink;
    c.textAlign = "center";
    c.fillText(xlabel, P.l + w / 2, H - 5);
    c.save();
    c.translate(13, P.t + h / 2);
    c.rotate(-Math.PI / 2);
    c.fillText(ylabel, 0, 0);
    c.restore();
    c.font = "9px system-ui";
    c.fillStyle = muted;
    if (ex) {
      c.textAlign = "right";
      c.fillText("×10" + superText(ex), P.l + w, H - 5);
    }
    if (ey) {
      c.textAlign = "left";
      c.fillText("×10" + superText(ey), P.l, P.t - 10);
    }
  }
  function clip(f) {
    const { c, P, w, h } = f;
    c.save();
    c.beginPath();
    c.rect(P.l, P.t, w, h);
    c.clip();
  }
  function closestFrame(m, t) {
    let lo = 0,
      hi = m.frames.length - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (m.frames[mid].t < t) lo = mid + 1;
      else hi = mid;
    }
    return lo > 0 &&
      Math.abs(m.frames[lo - 1].t - t) < Math.abs(m.frames[lo].t - t)
      ? lo - 1
      : lo;
  }
  function lattice(canvas, m, points, selected, limits) {
    const f = frame(canvas, limits, "g₂", "g₃"),
      { c, P, w, h, px, py } = f;
    clip(f);
    const [[xmin, xmax], [ymin, ymax]] = limits;
    if (xmax > 0) {
      const xs = Array.from(
        { length: 201 },
        (_, i) => Math.max(0, xmin) + ((xmax - Math.max(0, xmin)) * i) / 200,
      );
      c.fillStyle = "#eaf0e580";
      c.beginPath();
      xs.forEach((x, i) =>
        i
          ? c.lineTo(px(x), py(Math.sqrt(x ** 3 / 27)))
          : c.moveTo(px(x), py(Math.sqrt(x ** 3 / 27))),
      );
      xs.slice()
        .reverse()
        .forEach((x) => c.lineTo(px(x), py(-Math.sqrt(x ** 3 / 27))));
      c.closePath();
      c.fill();
      for (const sign of [-1, 1]) {
        c.strokeStyle = "#a9b4a9";
        c.lineWidth = 1;
        c.setLineDash([4, 3]);
        c.beginPath();
        xs.forEach((x, i) =>
          i
            ? c.lineTo(px(x), py(sign * Math.sqrt(x ** 3 / 27)))
            : c.moveTo(px(x), py(sign * Math.sqrt(x ** 3 / 27))),
        );
        c.stroke();
        c.setLineDash([]);
      }
    }
    const hit = [];
    m.tracks.forEach((tr, k) => {
      let previous = null;
      for (const p of tr.points) {
        if (!p) {
          previous = null;
          continue;
        }
        const fi = closestFrame(m, p[2]),
          pt =
            m.frames[fi].points.find((v) => v.branch === (tr.branch ?? k)) ||
            m.frames[fi].points.find(
              (v) =>
                Math.abs(v.g2 - p[0]) <
                  1e-11 * Math.max(Math.abs(p[0]), 1e-30) &&
                Math.abs(v.g3 - p[1]) < 1e-11 * Math.max(Math.abs(p[1]), 1e-30),
            );
        if (!pt) {
          previous = null;
          continue;
        }
        const cur = {
          x: px(p[0]),
          y: py(p[1]),
          frameIndex: fi,
          branch: pt.branch,
          point: pt,
          t: p[2],
        };
        const broken =
          previous &&
          (tr.breaks.some((v) => previous.t < v && cur.t > v) ||
            Math.max(
              Math.abs(cur.x - previous.x) / w,
              Math.abs(cur.y - previous.y) / h,
            ) > 0.35);
        if (previous && !broken) {
          c.strokeStyle = color((previous.t + cur.t) / 2, m);
          c.lineWidth = 2.4;
          c.lineCap = "round";
          c.beginPath();
          c.moveTo(previous.x, previous.y);
          c.lineTo(cur.x, cur.y);
          c.stroke();
          hit.push({ a: previous, b: cur });
        }
        if (f.inside(cur.x, cur.y)) hit.push({ a: cur, b: cur });
        previous = cur;
      }
    });
    points.forEach((p, i) => {
      const x = px(p.g2),
        y = py(p.g3);
      if (i === selected) {
        c.fillStyle = "#267d7818";
        c.beginPath();
        c.arc(x, y, 11, 0, 2 * Math.PI);
        c.fill();
      }
      c.fillStyle = i === selected ? ink : "#c27d54";
      c.strokeStyle = "#fffefa";
      c.lineWidth = 2;
      c.beginPath();
      c.arc(x, y, i === selected ? 5 : 4, 0, Math.PI * 2);
      c.fill();
      c.stroke();
    });
    c.restore();
    axes(f);
    f.hit = hit;
    const off = points.filter((p) => !f.inside(px(p.g2), py(p.g3))).length;
    if (off) {
      c.font = "9px system-ui";
      c.fillStyle = muted;
      c.textAlign = "right";
      c.fillText(
        `${off} selected ${off === 1 ? "point" : "points"} outside view`,
        P.l + w,
        P.t + 12,
      );
    }
    return f;
  }
  function pick(f, x, y, max = 18) {
    if (!f || !f.inside(x, y)) return null;
    let best = null,
      dist = max;
    for (const seg of f.hit) {
      const { a, b } = seg,
        dx = b.x - a.x,
        dy = b.y - a.y,
        d2 = dx * dx + dy * dy,
        t = d2
          ? Math.max(0, Math.min(1, ((x - a.x) * dx + (y - a.y) * dy) / d2))
          : 0,
        d = Math.hypot(x - a.x - t * dx, y - a.y - t * dy);
      if (d < dist) {
        best = t < 0.5 ? a : b;
        dist = d;
      }
    }
    return best;
  }
  function wave(canvas, p, spec, reference, pinned) {
    const refs = reference ? p.x.map(reference) : [],
      low = Math.min(p.range[0], ...refs),
      high = Math.max(p.range[1], ...refs);
    const margin = reference ? (high - low) * 0.04 : 0,
      f = frame(
        canvas,
        [
          [-p.width / 2, p.width / 2],
          [low - margin, high + margin],
        ],
        "Travelling coordinate z",
        spec.variable + "(z)",
      ),
      { c, P, w, h, px, py } = f;
    clip(f);
    if (p.mean !== null) {
      c.strokeStyle = "#afbbac";
      c.lineWidth = 1;
      c.setLineDash([3, 4]);
      c.beginPath();
      c.moveTo(P.l, py(p.mean));
      c.lineTo(P.l + w, py(p.mean));
      c.stroke();
      c.setLineDash([]);
    }
    const lines = [
      [p.re, "#267d78", []],
      ...(p.anyImag ? [[p.im, "#917ab3", [5, 3]]] : []),
      ...(reference ? [[refs, "#c27d54", [5, 4]]] : []),
    ];
    for (const [values, col, dash] of lines) {
      c.strokeStyle = col;
      c.lineWidth = col === "#267d78" ? 2.5 : 1.8;
      c.lineJoin = "round";
      c.setLineDash(dash);
      c.beginPath();
      let previous = null;
      for (let i = 0; i < values.length; i++) {
        const v = values[i];
        if (v === null || !Number.isFinite(v)) {
          previous = null;
          continue;
        }
        if (previous === null || Math.abs(v - previous) > 4 * (high - low))
          c.moveTo(px(p.x[i]), py(v));
        else c.lineTo(px(p.x[i]), py(v));
        previous = v;
      }
      c.stroke();
      c.setLineDash([]);
    }
    c.strokeStyle = "#bc8f77";
    c.lineWidth = 1;
    c.setLineDash([3, 4]);
    for (const x of p.poles) {
      c.beginPath();
      c.moveTo(px(x), P.t);
      c.lineTo(px(x), P.t + h);
      c.stroke();
    }
    c.setLineDash([]);
    if (pinned !== null) {
      const val = WaveMath.at(spec, pinned, p.mode);
      if (val) {
        c.fillStyle = ink;
        c.strokeStyle = "#fffefa";
        c.lineWidth = 2;
        c.beginPath();
        c.arc(px(pinned), py(val[0]), 4.5, 0, 2 * Math.PI);
        c.fill();
        c.stroke();
      }
    }
    c.restore();
    axes(f);
    return f;
  }
  return { color, fmt, ticks, lattice, wave, pick, closestFrame };
})();
