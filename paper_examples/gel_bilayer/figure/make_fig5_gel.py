"""Build Figure 5: mixed-order plane-strain gel-bilayer swelling.

Panels (a) and (b) summarize the generated mixed ``u-p-mu`` Quad8 field
supports and the gel/rubber bilayer in this package's accepted deck,
``../abaqus/SwellInducedBending_upmu_quad8.inp``.
The public package does not retain a full-run solver log or machine-readable
comparison record; see the package README for the evidence boundary.

Panels (c)--(f) reproduce the uploaded Abaqus/CAE exports at 0, 30 min, 1 h,
and 6 h.  The colored contours are maximum in-plane principal logarithmic
strain (``LE``) requested from the PHYSICAL CPE8H rubber set; the gel region
appears gray through its negligible-stiffness CPE8R visualization overlay
(E = 1e-20, sharing the UEL nodal motion).  For each nonzero frame, the supplied Abaqus
color ramp and exact endpoint range are retained, while the original 13 long
tick labels are replaced by three rounded values.  The limits differ by frame
and must not be read as a common contour scale.  No solver field or contour is
reconstructed by this script.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from PIL import Image, ImageChops


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
ABAQUS = FIGURES / "abaqus"
OUTPUT = FIGURES / "fig5_gel.pdf"

SNAPSHOTS = (
    ("c", r"Initial, $t=0$", ABAQUS / "gel_bilayer_00min.png", None, None),
    (
        "d",
        r"$t=30$ min",
        ABAQUS / "gel_bilayer_30min.png",
        ABAQUS / "gel_bilayer_30min_legend.png",
        (4.460e-3, 3.407e-2),
    ),
    (
        "e",
        r"$t=1$ h",
        ABAQUS / "gel_bilayer_1h.png",
        ABAQUS / "gel_bilayer_1h_legend.png",
        (8.313e-3, 5.761e-2),
    ),
    (
        "f",
        r"$t=6$ h",
        ABAQUS / "gel_bilayer_6h.png",
        ABAQUS / "gel_bilayer_6h_legend.png",
        (6.999e-3, 9.678e-2),
    ),
)

INK = "#303030"
MUTED = "#666666"
BLUE = "#2f6f9f"
BLUE_LIGHT = "#eaf2f8"
GREEN = "#4d7f4f"
GREEN_LIGHT = "#e8f2e7"
GOLD = "#b7791f"
GOLD_LIGHT = "#fff0c2"
RUBBER = "#d7d7d7"


def prepare_axes(ax: plt.Axes) -> None:
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")


def panel_header(ax: plt.Axes, label: str, title: str) -> None:
    prepare_axes(ax)
    ax.text(
        0.01,
        0.985,
        f"({label}) {title}",
        ha="left",
        va="top",
        fontsize=7.8,
        color=INK,
        fontweight="bold",
    )


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = INK,
    linewidth: float = 0.8,
    mutation: float = 6.0,
) -> None:
    ax.add_patch(FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=mutation,
        linewidth=linewidth,
        color=color,
        shrinkA=0.0,
        shrinkB=0.0,
        clip_on=False,
    ))


def draw_all_node_marker(ax: plt.Axes, points: list[tuple[float, float]]) -> None:
    """Draw the common quadratic displacement/chemical-potential support."""
    ax.scatter(
        [point[0] for point in points],
        [point[1] for point in points],
        s=28,
        marker="o",
        facecolor=BLUE,
        edgecolor="white",
        linewidth=0.45,
        zorder=6,
        clip_on=False,
    )


def draw_pressure_ring(ax: plt.Axes, points: list[tuple[float, float]]) -> None:
    """Draw the bilinear pressure support at the four corner nodes."""
    ax.scatter(
        [point[0] for point in points],
        [point[1] for point in points],
        s=94,
        marker="o",
        facecolor=GOLD_LIGHT,
        edgecolor=GOLD,
        linewidth=1.0,
        zorder=5,
        clip_on=False,
    )


def mixed_order_panel(ax: plt.Axes) -> None:
    prepare_axes(ax)
    ax.text(0.01, 0.985, "(a)", ha="left", va="top", fontsize=9.0,
            color=INK, fontweight="bold")

    x0, y0, width, height = 0.24, 0.475, 0.52, 0.40
    corners = [
        (x0, y0),
        (x0 + width, y0),
        (x0 + width, y0 + height),
        (x0, y0 + height),
    ]
    midsides = [
        (x0 + 0.5 * width, y0),
        (x0 + width, y0 + 0.5 * height),
        (x0 + 0.5 * width, y0 + height),
        (x0, y0 + 0.5 * height),
    ]
    ax.add_patch(Rectangle(
        (x0, y0),
        width,
        height,
        facecolor=BLUE_LIGHT,
        edgecolor=INK,
        linewidth=0.8,
    ))
    draw_pressure_ring(ax, corners)
    draw_all_node_marker(ax, corners + midsides)

    ax.scatter(
        [0.075],
        [0.360],
        s=28,
        facecolor=BLUE,
        edgecolor="white",
        linewidth=0.45,
        zorder=6,
    )
    ax.text(
        0.125,
        0.360,
        r"$u_x,u_y,\mu$: all eight nodes",
        ha="left",
        va="center",
        fontsize=6.4,
        color=BLUE,
    )
    ax.scatter(
        [0.075],
        [0.270],
        s=94,
        facecolor=GOLD_LIGHT,
        edgecolor=GOLD,
        linewidth=1.0,
        zorder=5,
    )
    ax.text(
        0.125,
        0.270,
        r"$p$: four corner nodes",
        ha="left",
        va="center",
        fontsize=6.4,
        color=GOLD,
    )

    ax.text(
        0.055,
        0.175,
        r"corner DOFs: $(1,2,11,12)$",
        ha="left",
        va="center",
        fontsize=6.4,
        color="#7a5210",
    )
    ax.text(
        0.055,
        0.085,
        r"midside DOFs: $(1,2,12)$",
        ha="left",
        va="center",
        fontsize=6.4,
        color=BLUE,
    )


def bilayer_panel(ax: plt.Axes) -> None:
    prepare_axes(ax)
    ax.text(0.01, 0.985, "(b)", ha="left", va="top", fontsize=9.0,
            color=INK, fontweight="bold")

    x0, y0, width, layer_height = 0.14, 0.40, 0.72, 0.165
    ax.add_patch(Rectangle(
        (x0, y0 + layer_height),
        width,
        layer_height,
        facecolor=GREEN_LIGHT,
        edgecolor=INK,
        linewidth=0.8,
    ))
    ax.add_patch(Rectangle(
        (x0, y0),
        width,
        layer_height,
        facecolor=RUBBER,
        edgecolor=INK,
        linewidth=0.8,
    ))
    # Emphasize the conforming gel/rubber interface.
    ax.plot(
        [x0, x0 + width],
        [y0 + layer_height, y0 + layer_height],
        color=INK,
        linewidth=1.4,
        zorder=4,
    )
    ax.text(
        x0 + 0.5 * width,
        y0 + 1.5 * layer_height,
        r"gel: mixed $u$--$p$--$\mu$ UEL",
        ha="center",
        va="center",
        fontsize=6.4,
        color=GREEN,
    )
    ax.text(
        x0 + 0.5 * width,
        y0 + 0.5 * layer_height,
        "rubber: CPE8H",
        ha="center",
        va="center",
        fontsize=6.4,
        color=INK,
    )
    # Layer thicknesses (equal layers), rotated right of the slider guide.
    ax.text(
        x0 + width + 0.078,
        y0 + layer_height,
        r"$2\times2.5$ mm",
        ha="center",
        va="center",
        fontsize=5.9,
        color=MUTED,
        rotation=90,
    )

    # Chemical-potential Dirichlet condition on the exposed gel top.
    ax.plot(
        [x0, x0 + width],
        [y0 + 2.0 * layer_height, y0 + 2.0 * layer_height],
        color=GREEN,
        linewidth=1.8,
    )
    for fraction in (0.14, 0.38, 0.62, 0.86):
        x = x0 + fraction * width
        arrow(
            ax,
            (x, y0 + 2.0 * layer_height + 0.095),
            (x, y0 + 2.0 * layer_height + 0.016),
            color=GREEN,
            mutation=5.2,
        )
    ax.text(
        x0 + 0.5 * width,
        y0 + 2.0 * layer_height + 0.108,
        r"prescribed $\bar{\mu}(t)$ on gel top",
        ha="center",
        va="bottom",
        fontsize=6.4,
        color=GREEN,
    )

    # Left symmetry edge: ux=0 (hatched wall); the lower-left corner
    # additionally pins uy=0.
    ax.plot(
        [x0, x0],
        [y0 - 0.015, y0 + 2.0 * layer_height + 0.015],
        color=BLUE,
        linewidth=1.6,
    )
    for y in (0.05, 0.14, 0.23, 0.31):
        ax.plot(
            [x0 - 0.002, x0 - 0.038],
            [y0 + y, y0 + y - 0.038],
            color=BLUE,
            linewidth=0.6,
        )
    ax.text(
        x0 - 0.052,
        y0 + layer_height,
        r"$u_x=0$",
        ha="center",
        va="center",
        fontsize=6.2,
        color=BLUE,
        rotation=90,
    )
    ax.add_patch(plt.Polygon(
        [(x0, y0), (x0 - 0.020, y0 - 0.042), (x0 + 0.020, y0 - 0.042)],
        closed=True,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=0.8,
        zorder=5,
    ))
    ax.text(
        x0 + 0.030,
        y0 - 0.052,
        r"corner $u_y=0$",
        ha="left",
        va="top",
        fontsize=6.0,
        color=BLUE,
    )

    # Right edge: straight-slider MPC (kept straight while free to
    # translate/rotate); hard contact activates after curling.
    xr = x0 + width
    ax.plot(
        [xr + 0.022, xr + 0.022],
        [y0 - 0.035, y0 + 2.0 * layer_height + 0.035],
        color=MUTED,
        linewidth=0.8,
        linestyle=(0, (4.0, 2.4)),
    )
    ax.add_patch(FancyArrowPatch(
        (xr + 0.022, y0 + 0.055),
        (xr + 0.022, y0 + 2.0 * layer_height - 0.055),
        arrowstyle="<|-|>",
        mutation_scale=5.0,
        linewidth=0.6,
        color=MUTED,
    ))

    # Half-length dimension.
    ax.plot([x0, x0 + width], [y0 - 0.115, y0 - 0.115],
            color=MUTED, linewidth=0.55)
    ax.plot([x0, x0], [y0 - 0.135, y0 - 0.095], color=MUTED, linewidth=0.55)
    ax.plot([xr, xr], [y0 - 0.135, y0 - 0.095], color=MUTED, linewidth=0.55)
    ax.text(
        x0 + 0.62 * width,
        y0 - 0.130,
        r"$L_{\rm half}=50$ mm (not to scale)",
        ha="center",
        va="top",
        fontsize=6.0,
        color=MUTED,
    )



def crop_white(path: Path, margin: int = 14) -> Image.Image:
    """Crop white viewport margins without altering image pixels."""
    image = Image.open(path).convert("RGB")
    background = Image.new("RGB", image.size, "white")
    difference = ImageChops.difference(image, background).convert("L")
    bbox = difference.point(lambda value: 255 if value > 8 else 0).getbbox()
    if bbox is None:
        raise ValueError(f"Snapshot is blank: {path}")
    left, upper, right, lower = bbox
    return image.crop((
        max(0, left - margin),
        max(0, upper - margin),
        min(image.width, right + margin),
        min(image.height, lower + margin),
    ))


def crop_legend_ramp(path: Path, margin: int = 3) -> Image.Image:
    """Retain the supplied Abaqus color ramp while removing its long labels."""
    image = Image.open(path).convert("RGBA")
    # All supplied legend exports place the ramp at the left.  Restricting the
    # crop before finding the alpha bounds removes the original 13 tick strings
    # without sampling or reconstructing any contour color.
    ramp_region = image.crop((0, 0, round(0.19 * image.width), image.height))
    bbox = ramp_region.getchannel("A").getbbox()
    if bbox is None:
        raise ValueError(f"Legend color ramp is blank: {path}")
    left, upper, right, lower = bbox
    cropped = ramp_region.crop((
        max(0, left - margin),
        max(0, upper - margin),
        min(ramp_region.width, right + margin),
        min(ramp_region.height, lower + margin),
    ))
    background = Image.new("RGBA", cropped.size, "white")
    background.alpha_composite(cropped)
    return background.convert("RGB")


def format_legend_value(value: float) -> str:
    """Format three significant digits without Abaqus' leading plus sign."""
    if value < 0.01:
        return f"{value:.4f}"
    return f"{value:.3f}"


def snapshot_panel(
    ax: plt.Axes,
    label: str,
    title: str,
    path: Path,
    legend_path: Path | None,
    legend_range: tuple[float, float] | None,
) -> None:
    panel_header(ax, label, title)

    if legend_path is None:
        result = ax.inset_axes([0.015, 0.030, 0.970, 0.900])
    else:
        result = ax.inset_axes([0.005, 0.030, 0.790, 0.900])
    result.imshow(crop_white(path))
    result.axis("off")

    if legend_path is not None:
        if legend_range is None:
            raise ValueError(f"Missing legend range for {path}")
        legend = ax.inset_axes([0.810, 0.220, 0.040, 0.560])
        legend.imshow(crop_legend_ramp(legend_path), aspect="auto")
        legend.axis("off")
        low, high = legend_range
        values = (high, 0.5 * (low + high), low)
        for y, value in zip((0.780, 0.500, 0.220), values):
            ax.text(
                0.862,
                y,
                format_legend_value(value),
                ha="left",
                va="center",
                fontsize=5.75,
                color=INK,
            )


def required_sources() -> tuple[Path, ...]:
    sources: list[Path] = []
    for _, _, snapshot, legend, _ in SNAPSHOTS:
        sources.append(snapshot)
        if legend is not None:
            sources.append(legend)
    return tuple(sources)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    missing = [path for path in required_sources() if not path.is_file()]
    if missing:
        missing_text = ", ".join(str(path.relative_to(ROOT)) for path in missing)
        raise FileNotFoundError(f"Missing Figure 5 source export(s): {missing_text}")

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "serif",
        "font.size": 8,
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.7,
    })

    # Use the manuscript's natural text-block width so LaTeX does not shrink
    # the already compact schematic and legend typography.
    fig = plt.figure(figsize=(6.50, 4.90))
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=(0.95, 1.30, 1.30),
        height_ratios=(1.0, 1.0),
        left=0.010,
        right=0.995,
        top=0.990,
        bottom=0.020,
        wspace=0.060,
        hspace=0.085,
    )
    mixed_order_panel(fig.add_subplot(grid[0, 0]))
    bilayer_panel(fig.add_subplot(grid[1, 0]))
    for spec, snapshot in zip(
        (grid[0, 1], grid[0, 2], grid[1, 1], grid[1, 2]),
        SNAPSHOTS,
    ):
        snapshot_panel(fig.add_subplot(spec), *snapshot)

    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
