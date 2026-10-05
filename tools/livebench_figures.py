"""Standalone scientific figures from the recorded f2py data, using Matplotlib."""

import textwrap

import numpy as np


def write_figures(plots, output, case_id, case_title, source):
    if not plots:
        return
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.tri import Triangulation

    directory = output / "figures" / case_id
    directory.mkdir(parents=True, exist_ok=True)
    colors = ["#157d79", "#b5783f", "#a4b5bf", "#8c729e"]
    revision = (source["commit"] or "unversioned")[:8] + (" + local changes" if source["dirty"] else "")

    with plt.rc_context({"svg.fonttype": "none", "svg.hashsalt": "abaqus_ufl_livebench", "font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.labelcolor": "#344455", "text.color": "#344455"}):
        for plot in plots:
            modes = [False, True] if plot.get("comparison") else [False]
            for difference in modes:
                fig, ax = plt.subplots(figsize=(8, 5.4))
                try:
                    fig.subplots_adjust(left=0.13, right=0.98, top=0.83, bottom=0.31)
                    fig.suptitle(case_title, x=0.13, y=0.96, ha="left", fontsize=10, color="#657482")
                    ax.set_title(plot["title"] + (" · discrepancy" if difference else ""), loc="left", fontsize=12, pad=12)
                    if plot["kind"] == "curve":
                        x = plot["x"]
                        if difference:
                            computed = next(s for s in plot["series"] if s["role"] == "computed")
                            reference = next(s for s in plot["series"] if s["role"] == "reference")
                            values = abs(np.asarray(computed["values"]) - reference["values"])
                            ax.plot(x, values, color=colors[0], marker="o", markersize=3, linewidth=1.3)
                            ax.set_ylabel("Absolute discrepancy")
                            ax.set_ylim(bottom=0)
                        else:
                            for i in reversed(range(len(plot["series"]))):
                                series = plot["series"][i]
                                ax.plot(x, series["values"], label=series["label"], color=colors[i % len(colors)],
                                        linestyle={"computed": "-", "reference": "--", "guide": ":"}[series["role"]],
                                        marker="o" if series["role"] == "computed" else None,
                                        markersize=3, markevery=max(1, len(x) // 24), linewidth=1.5)
                            handles, labels = ax.get_legend_handles_labels()
                            ax.legend(handles[::-1], labels[::-1], loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, frameon=False, fontsize=8)
                            ax.set_ylabel(plot["y_label"])
                            ax.set_yscale(plot["y_scale"])
                        ax.set_xlabel(plot["x_label"])
                        ax.set_xscale(plot["x_scale"])
                        ax.grid(alpha=0.2)
                    elif plot["kind"] == "mesh":
                        nodes = np.asarray(plot["nodes"], dtype=float)
                        if "displacements" in plot:
                            nodes = nodes + 10 * np.asarray(plot["displacements"])
                        triangles = [triangle for a, b, c, d in plot["elements"] for triangle in ((a, b, c), (a, c, d))]
                        mesh = Triangulation(nodes[:, 0], nodes[:, 1], triangles)
                        field = ax.tripcolor(mesh, plot["values"], shading="gouraud", cmap="viridis")
                        for element in plot["elements"]:
                            loop = nodes[element + [element[0]]]
                            ax.plot(loop[:, 0], loop[:, 1], color="white", linewidth=0.4, alpha=0.7)
                        fig.colorbar(field, ax=ax, label="Nodal value", fraction=0.045, pad=0.04)
                        ax.set_aspect("equal")
                        ax.set_xlabel("Displayed x" if "displacements" in plot else "Reference X")
                        ax.set_ylabel("Displayed y" if "displacements" in plot else "Reference Y")
                    else:
                        raise ValueError("unknown plot kind")
                    description = plot["description"].replace("⁻¹ᐟ²", "^(-1/2)")
                    fig.text(0.13, 0.04, textwrap.fill(description, width=115), fontsize=7, linespacing=1.3)
                    fig.text(0.13, 0.012, "abaqus_ufl · recorded f2py results · source " + revision, fontsize=6.5, color="#657482")
                    suffix = "-difference" if difference else ""
                    name = plot["id"] + suffix + ".svg"
                    path = directory / name
                    fig.savefig(path, format="svg", metadata={"Title": case_title + ": " + plot["title"],
                                                              "Description": plot["description"], "Creator": "abaqus_ufl livebench", "Date": None})
                    plot["difference_figure" if difference else "figure"] = path.relative_to(output).as_posix()
                finally:
                    plt.close(fig)
