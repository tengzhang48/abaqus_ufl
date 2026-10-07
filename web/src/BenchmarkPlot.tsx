import { useEffect, useMemo, useRef, useState } from "react";
import type { BenchmarkPlot, CurvePlot, MeshPlot } from "./livebench-report";

const colors = ["#157d79", "#b5783f", "#a4b5bf", "#8c729e"];
const format = (n: number) =>
  n === 0
    ? "0"
    : Math.abs(n) < 0.001 || Math.abs(n) >= 10000
      ? n.toExponential(2)
      : Number(n.toPrecision(4)).toString();

function niceTicks(min: number, max: number, target: number): number[] {
  if (min === max) return [min];
  const base = 10 ** Math.floor(Math.log10((max - min) / target));
  let best: number[] = [];
  for (const mag of [base / 10, base, base * 10]) {
    for (const multiplier of [1, 2, 2.5, 5]) {
      const step = multiplier * mag;
      const ticks: number[] = [];
      for (
        let value = Math.ceil(min / step - 1e-9) * step;
        value <= max + step * 1e-9;
        value += step
      ) {
        ticks.push(Number(value.toPrecision(12)));
      }
      if (
        ticks.length >= 2 &&
        (best.length < 2 ||
          Math.abs(ticks.length - target) < Math.abs(best.length - target))
      ) {
        best = ticks;
      }
    }
  }
  return best;
}

// Ticks in data units; on a log axis, min and max are log10 bounds.
function axisTicks(min: number, max: number, log: boolean, target: number) {
  if (!log) return niceTicks(min, max, target);
  const decades = Array.from(
    { length: Math.max(0, Math.floor(max) - Math.ceil(min) + 1) },
    (_, i) => 10 ** (Math.ceil(min) + i),
  );
  return decades.length >= 2 ? decades : niceTicks(10 ** min, 10 ** max, target);
}

const formatTick = (value: number) =>
  value !== 0 && (Math.abs(value) < 0.001 || Math.abs(value) >= 10000)
    ? value.toExponential()
    : String(value);

// Open on the largest recorded field, so load-unload histories do not
// start on their returned, nearly uniform final frame.
function peakFrame(plot: MeshPlot) {
  let best = 0,
    peak = -Infinity;
  plot.frames?.forEach((frame, index) => {
    const size = Math.max(...frame.values.map(Math.abs));
    if (size >= peak) {
      peak = size;
      best = index;
    }
  });
  return best;
}

function useWidth() {
  const ref = useRef<HTMLDivElement>(null);
  const [width, setWidth] = useState(720);
  useEffect(() => {
    if (!ref.current) return;
    const observer = new ResizeObserver(([entry]) =>
      setWidth(Math.max(240, entry.contentRect.width)),
    );
    observer.observe(ref.current);
    return () => observer.disconnect();
  }, []);
  return { ref, width };
}

function download(name: string, body: string, type: string) {
  const url = URL.createObjectURL(new Blob([body], { type }));
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function Curve({ plot, caseId }: { plot: CurvePlot; caseId: string }) {
  const { ref, width } = useWidth();
  const [point, setPoint] = useState(Math.floor(plot.x.length / 2));
  const [difference, setDifference] = useState(false);
  const [hidden, setHidden] = useState<number[]>([]);
  const computed = plot.series.find((s) => s.role === "computed")!;
  const reference = plot.series.find((s) => s.role === "reference");
  const errors = reference
    ? computed.values.map((v, i) => Math.abs(v - reference.values[i]))
    : null;
  const series: CurvePlot["series"] =
    difference && errors
      ? [{ label: "Absolute discrepancy", role: "computed", values: errors }]
      : plot.series.filter((_, i) => !hidden.includes(i));
  const xLog = plot.x_scale === "log";
  const yLog = !difference && plot.y_scale === "log";
  const tx = (x: number) => (xLog ? Math.log10(x) : x);
  const ty = (y: number) => (yLog ? Math.log10(y) : y);
  const height = width < 500 ? 300 : 370;
  const left = width < 500 ? 58 : 76,
    right = width - 20,
    top = 26,
    bottom = height - 58;
  const xs = plot.x.map(tx);
  const ys = series.flatMap((s) => s.values.map(ty));
  const xmin = Math.min(...xs),
    xmax = Math.max(...xs);
  let ymin = Math.min(...ys),
    ymax = Math.max(...ys);
  const nonnegative = ymin >= 0;
  if (ymin === ymax) {
    ymin -= Math.max(Math.abs(ymin) * 0.05, difference ? 1e-16 : 0.1);
    ymax += Math.max(Math.abs(ymax) * 0.05, difference ? 1e-16 : 0.1);
  }
  const pad = (ymax - ymin) * 0.08;
  ymin = nonnegative && !yLog ? Math.max(0, ymin - pad) : ymin - pad;
  ymax += pad;
  const px = (x: number) =>
    left + ((tx(x) - xmin) / (xmax - xmin)) * (right - left);
  const py = (y: number) =>
    bottom - ((ty(y) - ymin) / (ymax - ymin)) * (bottom - top);
  const tickCount = width < 500 ? 3 : 5;
  const xTicks = axisTicks(xmin, xmax, xLog, tickCount);
  const yTicks = axisTicks(ymin, ymax, yLog, tickCount);
  const color = (index: number) =>
    difference
      ? colors[0]
      : colors[plot.series.indexOf(series[index]) % colors.length];

  function csv() {
    const quote = (s: string) => `"${s.replaceAll('"', '""')}"`;
    const header = [
      plot.x_label,
      ...plot.series.map((s) => `${s.label}: ${plot.y_label}`),
      ...(errors ? [`Absolute discrepancy: ${plot.y_label}`] : []),
    ]
      .map(quote)
      .join(",");
    const rows = plot.x.map((x, i) =>
      [
        x,
        ...plot.series.map((s) => s.values[i]),
        ...(errors ? [errors[i]] : []),
      ].join(","),
    );
    download(
      `${caseId}-${plot.id}.csv`,
      [header, ...rows].join("\n") + "\n",
      "text/csv;charset=utf-8",
    );
  }

  return (
    <div className="result-plot">
      <div className="plot-toolbar">
        <div className="plot-modes">
          {errors && (
            <>
              <button
                type="button"
                aria-pressed={!difference}
                onClick={() => setDifference(false)}
              >
                Response
              </button>
              <button
                type="button"
                aria-pressed={difference}
                onClick={() => setDifference(true)}
              >
                Difference
              </button>
            </>
          )}
        </div>
        <div className="plot-downloads">
          <button type="button" onClick={csv}>
            CSV ↓
          </button>
          {(difference ? plot.difference_figure : plot.figure) && (
            <a
              href={`${import.meta.env.BASE_URL}livebench/${difference ? plot.difference_figure : plot.figure}`}
              download
            >
              Figure SVG ↓
            </a>
          )}
        </div>
      </div>
      <div className="plot-canvas" ref={ref}>
        <svg
          xmlns="http://www.w3.org/2000/svg"
          viewBox={`0 0 ${width} ${height}`}
          role="img"
          aria-label={
            difference ? `${plot.title}: absolute discrepancy` : plot.title
          }
          style={{ fontFamily: "system-ui, sans-serif", background: "white" }}
          onPointerMove={(event) => {
            const box = event.currentTarget.getBoundingClientRect();
            const x = ((event.clientX - box.left) * width) / box.width;
            let nearest = 0;
            plot.x.forEach((v, i) => {
              if (Math.abs(px(v) - x) < Math.abs(px(plot.x[nearest]) - x))
                nearest = i;
            });
            setPoint(nearest);
          }}
        >
          <title>{plot.title}</title>
          <desc>{plot.description}</desc>
          {yTicks.map((value) => {
            const y = py(value);
            return (
              <g key={`y${value}`}>
                <line x1={left} x2={right} y1={y} y2={y} stroke="#e5eaeb" />
                <text
                  x={left - 10}
                  y={y + 4}
                  textAnchor="end"
                  fontSize="11"
                  fill="#5d6976"
                >
                  {formatTick(value)}
                </text>
              </g>
            );
          })}
          {xTicks.map((value) => {
            const x = px(value);
            return (
              <g key={`x${value}`}>
                <line x1={x} x2={x} y1={top} y2={bottom} stroke="#f0f2f2" />
                <text
                  x={x}
                  y={bottom + 22}
                  textAnchor="middle"
                  fontSize="11"
                  fill="#5d6976"
                >
                  {formatTick(value)}
                </text>
              </g>
            );
          })}
          <line x1={left} x2={right} y1={bottom} y2={bottom} stroke="#93a3ae" />
          <text
            x={(left + right) / 2}
            y={height - 8}
            textAnchor="middle"
            fontSize={width < 500 ? 10 : 12}
            fill="#344455"
          >
            {plot.x_label}
            {xLog ? " · log scale" : ""}
          </text>
          <text
            transform={`translate(14 ${(top + bottom) / 2}) rotate(-90)`}
            textAnchor="middle"
            fontSize={width < 500 ? 10 : 12}
            fill="#344455"
          >
            {difference ? "Absolute discrepancy" : plot.y_label}
            {yLog ? " · log scale" : ""}
          </text>
          {[...series].reverse().map((s) => {
            const index = series.indexOf(s);
            return (
              <g key={s.label}>
                <polyline
                  points={plot.x
                    .map((x, i) => `${px(x)},${py(s.values[i])}`)
                    .join(" ")}
                  fill="none"
                  stroke={color(index)}
                  strokeWidth={s.role === "computed" ? 2 : 2.5}
                  strokeDasharray={
                    s.role === "reference"
                      ? "7 5"
                      : s.role === "guide"
                        ? "3 4"
                        : undefined
                  }
                />
                {s.role === "computed" &&
                  plot.x.map(
                    (x, i) =>
                      i % Math.max(1, Math.floor(plot.x.length / 22)) === 0 && (
                        <circle
                          key={i}
                          cx={px(x)}
                          cy={py(s.values[i])}
                          r="3"
                          fill="white"
                          stroke={color(index)}
                          strokeWidth="1.6"
                        />
                      ),
                  )}
              </g>
            );
          })}
          <line
            x1={px(plot.x[point])}
            x2={px(plot.x[point])}
            y1={top}
            y2={bottom}
            stroke="#aebac2"
            strokeDasharray="3 4"
          />
          {series.map((s, i) => (
            <circle
              key={s.label}
              cx={px(plot.x[point])}
              cy={py(s.values[point])}
              r="4"
              fill={color(i)}
              stroke="white"
              strokeWidth="1.5"
            />
          ))}
        </svg>
      </div>
      <div className="plot-legend">
        {plot.series.map((s, i) => (
          <button
            type="button"
            key={s.label}
            aria-pressed={!hidden.includes(i)}
            disabled={
              difference ||
              (!hidden.includes(i) && hidden.length === plot.series.length - 1)
            }
            onClick={() =>
              setHidden((list) =>
                list.includes(i) ? list.filter((v) => v !== i) : [...list, i],
              )
            }
          >
            <span style={{ background: colors[i % colors.length] }} />
            {s.label}
          </button>
        ))}
      </div>
      <div className="plot-inspector">
        <label htmlFor={`point-${plot.id}`}>
          Inspect sample {point + 1} / {plot.x.length}
        </label>
        <input
          id={`point-${plot.id}`}
          type="range"
          min="0"
          max={plot.x.length - 1}
          value={point}
          onChange={(event) => setPoint(Number(event.target.value))}
        />
        <dl>
          <div>
            <dt>{plot.x_label}</dt>
            <dd>{format(plot.x[point])}</dd>
          </div>
          <div>
            <dt>Generated Fortran</dt>
            <dd>{format(computed.values[point])}</dd>
          </div>
          {reference && (
            <>
              <div>
                <dt>{reference.label}</dt>
                <dd>{format(reference.values[point])}</dd>
              </div>
              <div>
                <dt>Absolute difference</dt>
                <dd>{format(errors![point])}</dd>
              </div>
            </>
          )}
        </dl>
      </div>
      {plot.comparison && (
        <p className="plot-agreement">
          <strong>
            Max absolute difference{" "}
            {plot.comparison.max_abs_error.toExponential(2)}
          </strong>
          <span>
            Every sample checked with atol{" "}
            {plot.comparison.atol.toExponential(0)}, rtol{" "}
            {plot.comparison.rtol.toExponential(0)}.
          </span>
        </p>
      )}
    </div>
  );
}

function Mesh({ plot, caseId }: { plot: MeshPlot; caseId: string }) {
  const { ref, width } = useWidth();
  const [node, setNode] = useState(Math.floor(plot.nodes.length / 2));
  const [deformed, setDeformed] = useState(Boolean(plot.displacements));
  const [lines, setLines] = useState(true);
  const [boundaries, setBoundaries] = useState(true);
  const [inspecting, setInspecting] = useState(false);
  const [frameIndex, setFrameIndex] = useState(() => peakFrame(plot));
  const [playing, setPlaying] = useState(false);
  const frame = plot.frames?.[frameIndex];
  const values = frame?.values ?? plot.values;
  const displacements = frame?.displacements ?? plot.displacements;
  const deformationScale = plot.deformation_scale ?? 10;
  const frameLabel = caseId === "ogden_bvp" ? "Load parameter s" : caseId === "plasticity_bvp" ? "Cycle pseudo-time" : "Time";
  useEffect(() => {
    if (!playing || !plot.frames) return;
    const timer = window.setInterval(() => {
      setFrameIndex((index) => (index + 1) % plot.frames!.length);
    }, 350);
    return () => window.clearInterval(timer);
  }, [playing, plot.frames]);
  const shown = plot.nodes.map((p, i) =>
    p.map((v, j) => v + (deformed ? deformationScale * displacements![i][j] : 0)),
  );
  // Fixed geometry and color scales across the recorded history keep changes
  // in deformation and temperature visible during playback.
  const [xmin, xmax, ymin, ymax] = useMemo(() => {
    const box = [Infinity, -Infinity, Infinity, -Infinity];
    const history = deformed
      ? (plot.frames?.map((f) => f.displacements) ?? [plot.displacements])
      : [undefined];
    for (const displacement of history) for (const [i, p] of plot.nodes.entries()) {
      const x = p[0] + (displacement ? deformationScale * displacement[i][0] : 0);
      const y = p[1] + (displacement ? deformationScale * displacement[i][1] : 0);
      box[0] = Math.min(box[0], x); box[1] = Math.max(box[1], x);
      box[2] = Math.min(box[2], y); box[3] = Math.max(box[3], y);
    }
    return box;
  }, [deformed, plot]);
  const height = Math.min(width < 500 ? 300 : 370,
    Math.max(width < 500 ? 150 : 220, 90 + (width - 60) * (ymax - ymin) / (xmax - xmin)));
  const scale = Math.min(
    (width - 60) / (xmax - xmin),
    (height - 90) / (ymax - ymin),
  );
  const ox = (width - scale * (xmax - xmin)) / 2,
    oy = (height - scale * (ymax - ymin)) / 2 - 14;
  const position = (p: number[]) => [
    ox + scale * (p[0] - xmin),
    oy + scale * (ymax - p[1]),
  ];
  const [low, high] = useMemo(() => {
    let minimum = Infinity, maximum = -Infinity;
    for (const v of plot.frames?.map((f) => f.values) ?? [plot.values]) for (const value of v) {
      minimum = Math.min(minimum, value); maximum = Math.max(maximum, value);
    }
    return [minimum, maximum];
  }, [plot]);
  const palette = [
    [68, 1, 84],
    [59, 82, 139],
    [33, 145, 140],
    [94, 201, 98],
    [253, 231, 37],
  ];
  const color = (v: number) => {
    const t = 4 * (high === low ? 0.5 : (v - low) / (high - low));
    const i = Math.min(3, Math.max(0, Math.floor(t))),
      f = t - i;
    return `rgb(${palette[i].map((c, j) => Math.round(c + f * (palette[i + 1][j] - c))).join(",")})`;
  };
  const boundaryColor = (index: number) => ["#915ca5", "#c96c37", "#3674a5", "#458473"][index % 4];
  function csv(history: boolean) {
    const selectedFrames = history ? plot.frames! : [{ time: frame?.time ?? 0, values, displacements }];
    const hasTime = Boolean(plot.frames);
    const quote = (text: string) => `"${text.replaceAll('"', '""')}"`;
    const header = [...(hasTime ? ["time"] : []), "node", "X", "Y", plot.value_label ?? "value",
                    ...(displacements ? ["u1", "u2"] : [])].map(quote).join(",");
    const rows = selectedFrames.flatMap((f) => plot.nodes.map((p, i) =>
      [...(hasTime ? [f.time] : []), i, ...p, f.values[i], ...(f.displacements?.[i] ?? [])].join(",")));
    download(`${caseId}-${plot.id}${history ? "-history" : hasTime ? `-frame-${frameIndex}` : ""}.csv`,
             [header, ...rows].join("\n") + "\n", "text/csv;charset=utf-8");
  }
  return (
    <div className="result-plot">
      <div className="plot-toolbar">
        <div className="plot-modes">
          <button
            type="button"
            aria-pressed={lines}
            onClick={() => setLines((v) => !v)}
          >
            Mesh lines
          </button>
          {plot.displacements && (
            <button
              type="button"
              aria-pressed={deformed}
              onClick={() => setDeformed((v) => !v)}
            >
              Deformation ×{deformationScale}
            </button>
          )}
          {plot.boundaries?.length ? (
            <button type="button" aria-pressed={boundaries} onClick={() => setBoundaries((v) => !v)}>
              Boundary conditions
            </button>
          ) : null}
        </div>
        <div className="plot-downloads">
          <button type="button" onClick={() => csv(false)}>{plot.frames ? "Frame CSV ↓" : "CSV ↓"}</button>
          {plot.frames && <button type="button" onClick={() => csv(true)}>History CSV ↓</button>}
          {plot.figure && (
            <a
              href={`${import.meta.env.BASE_URL}livebench/${plot.figure}`}
              download
            >
              {plot.frames ? "Final frame SVG ↓" : "Figure SVG ↓"}
            </a>
          )}
        </div>
      </div>
      {plot.frames && (
        <div className="field-timeline">
          <button type="button" onClick={() => {
            if (!playing && frameIndex === plot.frames!.length - 1) setFrameIndex(0);
            setPlaying((v) => !v);
          }}>{playing ? "Pause" : "Play"}</button>
          <label htmlFor={`time-${plot.id}`}>{frameLabel} {format(frame!.time)}</label>
          <input id={`time-${plot.id}`} type="range" min="0" max={plot.frames.length - 1}
            value={frameIndex} onChange={(event) => { setPlaying(false); setFrameIndex(Number(event.target.value)); }}
            aria-valuetext={`${frameLabel} ${format(frame!.time)}, step ${frameIndex} of ${plot.frames.length - 1}`} />
          <span>{frameIndex} / {plot.frames.length - 1}</span>
        </div>
      )}
      <div className="plot-canvas" ref={ref}>
        <svg
          viewBox={`0 0 ${width} ${height}`}
          role="img"
          aria-label={plot.title + (frame ? `, ${frameLabel.toLowerCase()} ${format(frame.time)}` : "")}
        >
          <title>{plot.title}</title>
          <desc>{plot.description}</desc>
          {plot.elements.map((element, index) => {
            const nodes = element.map((i) => shown[i]);
            // Bilinear interpolation of the actual nodal values, displayed on
            // a small subdivision of each Q1 cell; the solve is not rerun here.
            const cells = [];
            const subdivisions = plot.elements.length > 200 ? 2 : plot.elements.length > 100 ? 3 : 6;
            function interpolate(x: number, y: number) {
              const weights = [
                (1 - x) * (1 - y),
                x * (1 - y),
                x * y,
                (1 - x) * y,
              ];
              return {
                p: [0, 1].map((d) =>
                  weights.reduce((s, w, i) => s + w * nodes[i][d], 0),
                ),
                value: weights.reduce(
                  (s, w, i) => s + w * values[element[i]],
                  0,
                ),
              };
            }
            for (let y = 0; y < subdivisions; y++)
              for (let x = 0; x < subdivisions; x++) {
                const corners = [
                  [x, y],
                  [x + 1, y],
                  [x + 1, y + 1],
                  [x, y + 1],
                ].map(([a, b]) =>
                  position(
                    interpolate(a / subdivisions, b / subdivisions).p,
                  ).join(","),
                );
                const fill = color(
                  interpolate(
                    (x + 0.5) / subdivisions,
                    (y + 0.5) / subdivisions,
                  ).value,
                );
                cells.push(
                  <polygon
                    key={`${x}-${y}`}
                    points={corners.join(" ")}
                    fill={fill}
                    stroke={fill}
                    strokeWidth="0.4"
                  />,
                );
              }
            return (
              <g key={index}>
                {cells}
                {lines && (
                  <polygon
                    points={nodes.map((p) => position(p).join(",")).join(" ")}
                    fill="none"
                    stroke="#ffffffaa"
                    strokeWidth="1"
                  />
                )}
              </g>
            );
          })}
          {boundaries && plot.boundaries?.map((boundary, index) => (
            <polyline key={`${boundary.field}-${index}`} points={boundary.nodes.map((i) => position(shown[i]).join(",")).join(" ")}
              fill="none" stroke={boundaryColor(index)} strokeWidth="3" strokeDasharray={boundary.kind === "natural" ? "5 3" : undefined}>
              <title>{boundary.label}</title>
            </polyline>
          ))}
          {inspecting && (
            <path
              d={`M ${position(shown[node])[0] - 3} ${position(shown[node])[1]} h 6 M ${position(shown[node])[0]} ${position(shown[node])[1] - 3} v 6`}
              stroke="#17283b"
              strokeWidth="1.2"
              fill="none"
              aria-label={`Inspected node ${node + 1}`}
            />
          )}
          <text
            x={width / 2}
            y={height - (deformed && width < 500 ? 28 : 12)}
            textAnchor="middle"
            fontSize="12"
            fill="#5d6976"
          >
            <tspan>{plot.nodes.length} nodes · {plot.elements.length} Quad4 elements</tspan>
            {deformed && (
              <tspan x={width < 500 ? width / 2 : undefined} dy={width < 500 ? 16 : 0}>
                {`${width < 500 ? "" : " · "}deformation ×${deformationScale}`}
              </tspan>
            )}
          </text>
        </svg>
      </div>
      {boundaries && plot.boundaries?.length ? (
        <div className="boundary-legend">{plot.boundaries.map((boundary, index) => (
          <span key={`${boundary.field}-${index}`}><i style={{background: boundaryColor(index)}} />{boundary.label}</span>
        ))}</div>
      ) : null}
      <div className="field-scale">
        <span>{format(low)}</span>
        <i />
        <span>{format(high)}</span>
      </div>
      {plot.value_label && <p className="field-value-label">{plot.value_label} · scale fixed across all time steps</p>}
      <details className="plot-inspector" onToggle={(event) => setInspecting(event.currentTarget.open)}>
        <summary>Inspect nodal values</summary>
        <label htmlFor={`node-${plot.id}`}>
          Inspect node {node + 1} / {plot.nodes.length}
        </label>
        <input
          id={`node-${plot.id}`}
          type="range"
          min="0"
          max={plot.nodes.length - 1}
          value={node}
          onChange={(event) => setNode(Number(event.target.value))}
        />
        <dl>
          <div>
            <dt>Reference X</dt>
            <dd>{format(plot.nodes[node][0])}</dd>
          </div>
          <div>
            <dt>Reference Y</dt>
            <dd>{format(plot.nodes[node][1])}</dd>
          </div>
          <div>
            <dt>{plot.value_label ?? "Nodal value"}</dt>
            <dd>{format(values[node])}</dd>
          </div>
          {displacements && <>
            <div><dt>u₁</dt><dd>{format(displacements[node][0])}</dd></div>
            <div><dt>u₂</dt><dd>{format(displacements[node][1])}</dd></div>
          </>}
        </dl>
      </details>
    </div>
  );
}

export default function Plot({
  plot,
  caseId,
}: {
  plot: BenchmarkPlot;
  caseId: string;
}) {
  return plot.kind === "curve" ? (
    <Curve plot={plot} caseId={caseId} />
  ) : (
    <Mesh plot={plot} caseId={caseId} />
  );
}
