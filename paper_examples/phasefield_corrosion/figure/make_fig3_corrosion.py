"""Create Figure 3 (phase-field stress corrosion) for the EML manuscript.

Panel (a) is a conventional vector schematic of the generated-UEL comparison
deck (SCC_abaqus_ufl_Tdummy_diag.inp): the 0.30 x 0.15 mm semi-elliptical-pit
domain, +/-0.20 micrometre horizontal loading, and the Step-4 scalar boundary
ramp.  The bottom-corner u_y pins removing the vertical rigid-body mode are an
addition of that adapted deck; the Cui-supplied SCC deck does not include
them.  Everything else matches the released benchmark setup.

Panel (b) is a code-to-code comparison, not digitized journal-article data.
The black line is the 747-frame extraction from the original Cui-supplied UEL;
the blue markers are all 58 raw output frames from the generated fully
block-diagonal ``abaqus_ufl`` UEL.  Only the displayed error statistics use
linear interpolation of the generated history to the original frame times.
This deck matches Section 2.2/Figure 5 of the supplied code documentation.
The 2021 article's semi-elliptical-pit curves use different loads and are not
silently substituted.

Panel (c) embeds the one supplied Abaqus/CAE rendering of the generated-UEL
nodal field NT11 = phi.  The supplied color ramp is retained with three short
labels.  PNG metadata identify Abaqus/CAE 6.24 but not the ODB, job, frame, or
display scale, so those details remain an archival provenance item rather than
an inferred claim in the figure.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle
from PIL import Image, ImageChops


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
ABAQUS_DIR = FIG_DIR / "abaqus"
CUI = ROOT.parent / "examples" / "phasefield_corrosion_cui"
ORIGINAL_CURVE_CSV = (
    CUI / "UELCorrosion" / "extracted_results_original"
    / "SCC_original_pit_depth_curve.csv"
)
GENERATED_CURVE_CSV = (
    CUI / "abaqus_test_from_cui"
    / "extracted_results_diag_fixed_order_20260611"
    / "SCC_abaqus_ufl_Tdummy_diag_pit_depth_curve.csv"
)
ALIGNED_CURVE_CSV = (
    CUI / "UELCorrosion" / "extracted_results_original"
    / "cui_original_vs_abaqus_ufl_pit_depth.csv"
)
FIELD_IMAGE = ABAQUS_DIR / "abaqus_ufl_cui_corrison_nt11.png"
FIELD_LEGEND = ABAQUS_DIR / "cui_corrison_legend.png"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def load_curve():
    original_rows = read_rows(ORIGINAL_CURVE_CSV)
    generated_rows = read_rows(GENERATED_CURVE_CSV)
    aligned_rows = read_rows(ALIGNED_CURVE_CSV)
    original_t = np.array([float(row["time_s"]) for row in original_rows])
    original = np.array([float(row["pit_depth_um"]) for row in original_rows])
    generated_t = np.array([float(row["time_s"]) for row in generated_rows])
    generated = np.array([float(row["pit_depth_um"]) for row in generated_rows])
    diff = np.array([
        float(row["abaqus_ufl_minus_original_um"])
        for row in aligned_rows
    ])
    return original_t, original, generated_t, generated, diff


def setup_panel(ax):
    """Draw the geometry and boundary schedule in the released SCC deck."""
    ink = "#303030"
    muted = "#666666"
    blue = "#2f6f9f"
    red = "#a33a35"
    metal = "#e9edf1"

    ax.set_xlim(-0.075, 0.375)
    ax.set_ylim(-0.052, 0.318)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.text(-0.02, 1.045, "(a)", transform=ax.transAxes, ha="left",
            va="top", fontsize=9.0, fontweight="bold", clip_on=False)

    ax.add_patch(Rectangle(
        (0.0, 0.0),
        0.30,
        0.15,
        facecolor=metal,
        edgecolor=ink,
        linewidth=0.9,
    ))

    # Initial pit: 0.02 mm wide and 0.02 mm deep in the released
    # documentation schematic.  The white polygon hides the metal fill.
    angle = np.linspace(np.pi, 0.0, 81)
    pit_x = 0.15 + 0.01 * np.cos(angle)
    pit_y = 0.15 - 0.02 * np.sin(angle)
    pit = np.column_stack([
        np.concatenate(([0.14], pit_x, [0.16])),
        np.concatenate(([0.17], pit_y, [0.17])),
    ])
    ax.add_patch(Polygon(
        pit,
        closed=True,
        facecolor="white",
        edgecolor="none",
        zorder=3,
    ))
    ax.plot(pit_x, pit_y, color=blue, linewidth=1.1, zorder=4)
    ax.plot([0.0, 0.14], [0.15, 0.15], color=ink, linewidth=0.9, zorder=4)
    ax.plot([0.16, 0.30], [0.15, 0.15], color=ink, linewidth=0.9, zorder=4)

    ax.text(
        0.15,
        0.062,
        r"Quad8 UEL: $(u_x,u_y,\phi,c)$",
        ha="center",
        va="center",
        fontsize=7.2,
        color=ink,
    )

    # Horizontal remote displacement applied in Step 3; short labels sit
    # above and below each arrow so nothing collides with the detail view.
    for start, end in (
        ((0.0, 0.075), (-0.055, 0.075)),
        ((0.30, 0.075), (0.355, 0.075)),
    ):
        ax.add_patch(FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=1.1,
            color=red,
            clip_on=False,
        ))
    ax.text(-0.033, 0.086, r"$u_x$", ha="center", va="bottom",
            fontsize=6.6, color=red)
    ax.text(-0.033, 0.063, r"$-0.20\,\mu$m", ha="center", va="top",
            fontsize=6.3, color=red)
    ax.text(0.333, 0.086, r"$u_x$", ha="center", va="bottom",
            fontsize=6.6, color=red)
    ax.text(0.333, 0.063, r"$+0.20\,\mu$m", ha="center", va="top",
            fontsize=6.3, color=red)

    # Bottom-corner pins remove the vertical rigid-body mode.  They are a
    # BC of the adapted generated-UEL deck only; the Cui-supplied deck has
    # no DOF-2 pins (disclosed in the manuscript caption).
    for x_pin in (0.0, 0.30):
        ax.add_patch(Polygon(
            [(x_pin, 0.0), (x_pin - 0.011, -0.020), (x_pin + 0.011, -0.020)],
            closed=True,
            facecolor="white",
            edgecolor=ink,
            linewidth=0.8,
            zorder=5,
        ))
    ax.text(
        0.15,
        -0.014,
        r"corner $u_y=0$",
        ha="center",
        va="top",
        fontsize=6.4,
        color=muted,
    )

    # Overall dimensions along the bottom.
    ax.add_patch(FancyArrowPatch(
        (0.0, -0.036), (0.30, -0.036),
        arrowstyle="<|-|>", mutation_scale=6, linewidth=0.6, color=muted,
    ))
    ax.text(
        0.15,
        -0.041,
        r"$0.30\times0.15$ mm",
        ha="center",
        va="top",
        fontsize=6.4,
        color=muted,
    )

    # --- Magnified pit detail --------------------------------------------
    # A circular marker on the true-scale pit connects to an enlarged view;
    # the scalar-surface ramp annotation attaches to the enlarged view.
    ax.add_patch(Circle(
        (0.15, 0.147), 0.030,
        facecolor="none", edgecolor=muted, linewidth=0.7, zorder=6,
    ))

    detail = ax.inset_axes([0.020, 0.585, 0.400, 0.395])
    detail.set_xlim(0.108, 0.192)
    detail.set_ylim(0.114, 0.170)
    detail.set_aspect("equal", adjustable="box")
    detail.set_xticks([])
    detail.set_yticks([])
    for side in detail.spines.values():
        side.set_color(muted)
        side.set_linewidth(0.7)

    detail.add_patch(Rectangle(
        (0.108, 0.114), 0.084, 0.036,
        facecolor=metal, edgecolor="none",
    ))
    detail.add_patch(Polygon(
        pit, closed=True, facecolor="white", edgecolor="none", zorder=3,
    ))
    detail.plot(pit_x, pit_y, color=blue, linewidth=1.6, zorder=4)
    detail.plot([0.108, 0.14], [0.15, 0.15], color=ink, linewidth=1.0, zorder=4)
    detail.plot([0.16, 0.192], [0.15, 0.15], color=ink, linewidth=1.0, zorder=4)

    # Pit dimensions inside the detail: 20 um wide, 20 um deep.
    detail.add_patch(FancyArrowPatch(
        (0.14, 0.124), (0.16, 0.124),
        arrowstyle="<|-|>", mutation_scale=4.5, linewidth=0.55, color=muted,
    ))
    detail.text(0.15, 0.1215, r"$20\,\mu$m", ha="center", va="top",
                fontsize=6.0, color=muted)
    detail.add_patch(FancyArrowPatch(
        (0.170, 0.150), (0.170, 0.130),
        arrowstyle="<|-|>", mutation_scale=4.5, linewidth=0.55, color=muted,
    ))
    detail.text(0.1735, 0.140, r"$20\,\mu$m", ha="left", va="center",
                fontsize=6.0, color=muted)

    # Dashed leader lines from the true-scale marker circle to the two
    # bottom corners of the detail view (classic magnifier callout).
    ax.plot([0.124, -0.056], [0.158, 0.163], color=muted, linewidth=0.6,
            linestyle=(0, (2.5, 1.8)), zorder=6)
    ax.plot([0.152, 0.112], [0.177, 0.163], color=muted, linewidth=0.6,
            linestyle=(0, (2.5, 1.8)), zorder=6)

    # Scalar boundary ramp on the pit surface, tied to the detail view.
    ax.text(
        0.150,
        0.300,
        "pit surface:\n"
        r"$\phi,c:\ 1\!\rightarrow\!0$ over $0.005$ s;"
        "\nhold",
        ha="left",
        va="top",
        fontsize=6.5,
        color=blue,
        linespacing=1.30,
    )


def crop_white(path: Path, margin: int = 12) -> Image.Image:
    """Crop a supplied Abaqus viewport without changing its pixels."""
    image = Image.open(path).convert("RGB")
    background = Image.new("RGB", image.size, "white")
    difference = ImageChops.difference(image, background).convert("L")
    bbox = difference.point(lambda value: 255 if value > 8 else 0).getbbox()
    if bbox is None:
        raise ValueError(f"Abaqus snapshot is blank: {path}")
    left, upper, right, lower = bbox
    return image.crop((
        max(0, left - margin),
        max(0, upper - margin),
        min(image.width, right + margin),
        min(image.height, lower + margin),
    ))


def crop_color_ramp(path: Path, margin: int = 3) -> Image.Image:
    """Crop the colored ramp from the supplied transparent legend export."""
    rgba = Image.open(path).convert("RGBA")
    background = Image.new("RGBA", rgba.size, "white")
    background.alpha_composite(rgba)
    rgb = background.convert("RGB")
    left_region = rgb.crop((0, 0, round(0.30 * rgb.width), rgb.height))
    pixels = np.asarray(left_region)
    colored = np.ptp(pixels.astype(np.int16), axis=2) > 35
    rows, cols = np.nonzero(colored)
    if not rows.size:
        raise ValueError(f"Abaqus legend has no color ramp: {path}")
    box = (
        max(0, int(cols.min()) - margin),
        max(0, int(rows.min()) - margin),
        min(left_region.width, int(cols.max()) + margin + 1),
        min(left_region.height, int(rows.max()) + margin + 1),
    )
    return left_region.crop(box)


def field_panel(result_ax: plt.Axes) -> None:
    """Place the supplied Abaqus field image and its supplied color ramp."""
    result_ax.set_xlim(0.0, 1.0)
    result_ax.set_ylim(0.0, 1.0)
    result_ax.axis("off")
    result_ax.text(
        0.075, 0.97, "(c)", ha="left", va="top",
        fontsize=9.0, fontweight="bold",
    )

    image_ax = result_ax.inset_axes([0.13, 0.00, 0.64, 0.97])
    image_ax.imshow(crop_white(FIELD_IMAGE))
    image_ax.axis("off")

    legend_ax = result_ax.inset_axes([0.795, 0.02, 0.10, 0.92])
    legend_ax.set_xlim(0.0, 1.0)
    legend_ax.set_ylim(0.0, 1.0)
    legend_ax.axis("off")
    legend_ax.add_patch(Rectangle(
        (0.0, 0.0),
        1.0,
        1.0,
        facecolor="white",
        edgecolor="none",
        zorder=-1,
    ))
    legend_ax.text(
        0.48,
        0.96,
        r"NT11 $=\phi$",
        ha="center",
        va="top",
        fontsize=7.0,
    )
    ramp = legend_ax.inset_axes([0.10, 0.16, 0.26, 0.68])
    ramp.imshow(crop_color_ramp(FIELD_LEGEND))
    ramp.axis("off")
    for y, value in ((0.84, "1.0"), (0.50, "0.5"), (0.16, "0.0")):
        legend_ax.text(
            0.46,
            y,
            value,
            ha="left",
            va="center",
            fontsize=6.8,
            color="#303030",
        )


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    missing = [path for path in (FIELD_IMAGE, FIELD_LEGEND) if not path.is_file()]
    if missing:
        names = ", ".join(str(path.relative_to(ROOT)) for path in missing)
        raise FileNotFoundError(f"Missing supplied Figure 3 Abaqus image(s): {names}")

    original_t, original, generated_t, generated, diff = load_curve()
    if len(original) != 747 or len(generated) != 58 or len(diff) != 747:
        raise ValueError(
            "Figure 3 provenance changed: expected 747 original frames, "
            "58 raw generated frames, and 747 aligned comparison rows"
        )
    if not (
        np.isclose(original[-1], 46.64082080125809)
        and np.isclose(generated[-1], 46.640821562223685)
    ):
        raise ValueError("Figure 3 endpoint values no longer match the audit")

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "serif",
        "font.size": 8.5,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
    })

    fig = plt.figure(figsize=(6.50, 4.35), constrained_layout=True)
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.10, 1.18),
        height_ratios=(1.00, 1.04),
    )

    setup_panel(fig.add_subplot(grid[0, 0]))

    ax = fig.add_subplot(grid[0, 1])
    ax.plot(
        original_t,
        original,
        color="#333333",
        linewidth=1.35,
        label="reference UEL",
        zorder=1,
    )
    ax.plot(
        generated_t,
        generated,
        linestyle="none",
        marker="o",
        markersize=3.0,
        markerfacecolor="none",
        markeredgewidth=0.85,
        markeredgecolor="#2f6f9f",
        label="generated UEL",
        zorder=2,
    )
    ax.set_xlabel("total analysis time (s)", fontsize=8.3)
    ax.set_ylabel(r"SCC-region length ($\mu$m)", fontsize=8.3)
    ax.text(-0.145, 1.045, "(b)", transform=ax.transAxes, ha="left",
            va="top", fontsize=9.0, fontweight="bold", clip_on=False)
    ax.set_xlim(0.0, 910.0)
    ax.set_ylim(0.0, 50.0)
    ax.set_xticks([0.0, 300.0, 600.0, 900.0])
    ax.legend(frameon=False, fontsize=7.0, loc="upper left")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=7.7, width=0.7)

    field_panel(fig.add_subplot(grid[1, :]))

    out = FIG_DIR / "fig3_corrosion.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
