import { useEffect, useRef, useState } from "react";
import type { BenchmarkPlot, CurvePlot, MeshPlot } from "./livebench-report";

const colors = ["#157d79", "#b5783f", "#a4b5bf", "#8c729e"];
const format = (n: number) =>
  n === 0
    ? "0"
    : Math.abs(n) < 0.001 || Math.abs(n) >= 10000
      ? n.toExponential(2)
      : Number(n.toPrecision(4)).toString();

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
  const ticks = Array.from(
    { length: tickCount },
    (_, i) => i / (tickCount - 1),
  );
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
          {ticks.map((t) => {
            const y = bottom - t * (bottom - top),
              value = ymin + t * (ymax - ymin);
            return (
              <g key={`y${t}`}>
                <line x1={left} x2={right} y1={y} y2={y} stroke="#e5eaeb" />
                <text
                  x={left - 10}
                  y={y + 4}
                  textAnchor="end"
                  fontSize="11"
                  fill="#5d6976"
                >
                  {format(yLog ? 10 ** value : value)}
                </text>
              </g>
            );
          })}
          {ticks.map((t) => {
            const x = left + t * (right - left),
              value = xmin + t * (xmax - xmin);
            return (
              <g key={`x${t}`}>
                <line x1={x} x2={x} y1={top} y2={bottom} stroke="#f0f2f2" />
                <text
                  x={x}
                  y={bottom + 22}
                  textAnchor="middle"
                  fontSize="11"
                  fill="#5d6976"
                >
                  {format(xLog ? 10 ** value : value)}
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
  const height = width < 500 ? 300 : 370;
  const shown = plot.nodes.map((p, i) =>
    p.map((v, j) => v + (deformed ? 10 * plot.displacements![i][j] : 0)),
  );
  const xmin = Math.min(...shown.map((p) => p[0])),
    xmax = Math.max(...shown.map((p) => p[0]));
  const ymin = Math.min(...shown.map((p) => p[1])),
    ymax = Math.max(...shown.map((p) => p[1]));
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
  const low = Math.min(...plot.values),
    high = Math.max(...plot.values);
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
              Deformation ×10
            </button>
          )}
        </div>
        <div className="plot-downloads">
          <button
            type="button"
            onClick={() =>
              download(
                `${caseId}-${plot.id}.csv`,
                `node,X,Y,value${plot.displacements ? ",u1,u2" : ""}\n` +
                  plot.nodes
                    .map((p, i) =>
                      [
                        i,
                        ...p,
                        plot.values[i],
                        ...(plot.displacements?.[i] ?? []),
                      ].join(","),
                    )
                    .join("\n") +
                  "\n",
                "text/csv;charset=utf-8",
              )
            }
          >
            CSV ↓
          </button>
          {plot.figure && (
            <a
              href={`${import.meta.env.BASE_URL}livebench/${plot.figure}`}
              download
            >
              Figure SVG ↓
            </a>
          )}
        </div>
      </div>
      <div className="plot-canvas" ref={ref}>
        <svg
          viewBox={`0 0 ${width} ${height}`}
          role="img"
          aria-label={plot.title}
        >
          <title>{plot.title}</title>
          <desc>{plot.description}</desc>
          {plot.elements.map((element, index) => {
            const nodes = element.map((i) => shown[i]);
            // Bilinear interpolation of the actual nodal values, displayed on
            // a small subdivision of each Q1 cell; the solve is not rerun here.
            const cells = [];
            const subdivisions = 6;
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
                  (s, w, i) => s + w * plot.values[element[i]],
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
          {shown.map((p, i) => (
            <circle
              key={i}
              cx={position(p)[0]}
              cy={position(p)[1]}
              r={i === node ? 5 : 2.5}
              fill={i === node ? "white" : color(plot.values[i])}
              stroke={i === node ? "#17283b" : "#ffffff88"}
              onPointerEnter={() => setNode(i)}
              onClick={() => setNode(i)}
            />
          ))}
          <text
            x={width / 2}
            y={height - 12}
            textAnchor="middle"
            fontSize="12"
            fill="#5d6976"
          >
            {plot.nodes.length} nodes · {plot.elements.length} Quad4 elements
            {deformed ? " · deformation ×10" : ""}
          </text>
        </svg>
      </div>
      <div className="field-scale">
        <span>{format(low)}</span>
        <i />
        <span>{format(high)}</span>
      </div>
      <div className="plot-inspector">
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
            <dt>Nodal value</dt>
            <dd>{format(plot.values[node])}</dd>
          </div>
        </dl>
      </div>
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
