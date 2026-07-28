"""Assemble the mixed-Tet4 design, benchmark setup, and Abaqus snapshots.

Panels (a) and (b) are vector schematics of the generated equal-order Tet4
element and the n=16 quarter-block compression problem, respectively.  Panels
(c) and (d) are direct Abaqus/CAE renders of the completed n=16 ODB.  No
simulated field data are synthesized or interpolated in this script.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Polygon
import numpy as np
from PIL import Image, ImageChops


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
ABAQUS_DIR = FIG_DIR / "abaqus"

U_IMAGE = ABAQUS_DIR / "Tet4_u_mag.png"
THETAT_IMAGE = ABAQUS_DIR / "Tet4_NT11.png"

INK = "#333333"
MUTED = "#666666"
BLUE = "#2f6f9f"
BLUE_LIGHT = "#eaf2f8"
GOLD = "#b7791f"
GOLD_LIGHT = "#fff4d6"
RED = "#a33a35"
RED_LIGHT = "#f6dedb"
NEUTRAL_LIGHT = "#f2f2f2"
def crop_white_image(image: Image.Image, margin: int = 18) -> Image.Image:
    """Crop white viewport margins without changing supplied image pixels."""
    image = image.convert("RGB")
    background = Image.new("RGB", image.size, "white")
    difference = ImageChops.difference(image, background).convert("L")
    bbox = difference.point(lambda value: 255 if value > 8 else 0).getbbox()
    if bbox is None:
        raise RuntimeError(f"snapshot is blank: {path}")
    left, upper, right, lower = bbox
    box = (
        max(0, left - margin), max(0, upper - margin),
        min(image.width, right + margin), min(image.height, lower + margin),
    )
    return image.crop(box)


def crop_model(path: Path, margin: int = 18) -> Image.Image:
    """Separate the largest colored model region from the embedded legend."""
    image = Image.open(path).convert("RGB")
    pixels = np.asarray(image)
    colored = np.ptp(pixels.astype(np.int16), axis=2) > 35
    active = np.count_nonzero(colored, axis=0) > max(12, image.height // 150)

    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, value in enumerate(active):
        if value and start is None:
            start = index
        elif not value and start is not None:
            runs.append((start, index))
            start = None
    if start is not None:
        runs.append((start, image.width))
    if not runs:
        raise RuntimeError(f"no colored Abaqus model found: {path}")

    left, right = max(
        runs,
        key=lambda run: int(np.count_nonzero(colored[:, run[0]:run[1]])),
    )
    region = image.crop((
        max(0, left - margin),
        0,
        min(image.width, right + margin),
        image.height,
    ))
    return crop_white_image(region, margin=margin)


def crop_color_ramp(path: Path, margin: int = 3) -> Image.Image:
    """Crop the supplied Abaqus ramp while discarding its long tick strings."""
    image = Image.open(path).convert("RGB")
    left_region = image.crop((0, 0, round(0.18 * image.width), image.height))
    pixels = np.asarray(left_region)
    colored = np.ptp(pixels.astype(np.int16), axis=2) > 35
    rows, cols = np.nonzero(colored)
    if not rows.size:
        raise RuntimeError(f"no colored Abaqus legend found: {path}")
    return left_region.crop((
        max(0, int(cols.min()) - margin),
        max(0, int(rows.min()) - margin),
        min(left_region.width, int(cols.max()) + margin + 1),
        min(left_region.height, int(rows.max()) + margin + 1),
    ))


def element_panel(ax: plt.Axes) -> None:
    # Hand-tuned 2D embedding that reads as a solid tetrahedron: base nodes
    # 1-2 in front, node 3 behind (its base edges dashed), apex 4 on top.
    verts = {
        1: np.array([0.13, 0.285]),
        2: np.array([0.84, 0.240]),
        3: np.array([0.65, 0.465]),
        4: np.array([0.42, 0.800]),
    }

    for face, fill in (((1, 2, 4), BLUE_LIGHT), ((2, 3, 4), "#d9e6ef")):
        ax.add_patch(Polygon(
            [verts[i] for i in face], closed=True, facecolor=fill,
            edgecolor="none", zorder=2,
        ))
    for a, b in ((1, 2), (2, 3), (1, 4), (2, 4), (3, 4)):
        ax.plot(*zip(verts[a], verts[b]), color=INK, linewidth=0.9, zorder=3)
    ax.plot(*zip(verts[1], verts[3]), color=MUTED, linewidth=0.8,
            linestyle=(0, (3.2, 2.2)), zorder=4)

    label_offsets = {
        1: np.array([-0.052, -0.020]),
        2: np.array([0.050, -0.020]),
        3: np.array([0.055, 0.018]),
        4: np.array([-0.055, 0.018]),
    }
    for index, xy in verts.items():
        ax.add_patch(Circle(
            xy, 0.028, facecolor=GOLD_LIGHT, edgecolor=GOLD,
            linewidth=1.2, zorder=7,
        ))
        ax.add_patch(Circle(
            xy, 0.013, facecolor=BLUE, edgecolor="white",
            linewidth=0.5, zorder=8,
        ))
        ax.text(
            *(xy + label_offsets[index]), str(index), ha="center",
            va="center", fontsize=7, color=INK, zorder=9,
        )

    # One representative displacement-vector glyph; every node carries it.
    anchor = verts[4]
    ax.add_patch(FancyArrowPatch(
        anchor + np.array([0.024, 0.016]),
        anchor + np.array([0.150, 0.062]),
        arrowstyle="-|>", mutation_scale=7, linewidth=1.0,
        color=BLUE, zorder=9,
    ))
    ax.text(
        anchor[0] + 0.165, anchor[1] + 0.068, r"$\mathbf{u}_a$",
        fontsize=7.2, color=BLUE, ha="left", va="center",
    )

    ax.text(0.02, 0.98, "(a)", transform=ax.transAxes,
            ha="left", va="top", fontsize=9.0, fontweight="bold")

    # Legend and formulation summary: three aligned lines.
    ax.scatter([0.065], [0.105], s=25, color=BLUE, edgecolor="white",
               linewidth=0.5, zorder=5)
    ax.text(0.105, 0.105, r"nodal $\mathbf{u}_a\in\mathbb{R}^3$",
            fontsize=6.7, ha="left", va="center", color=INK)
    ax.scatter([0.525], [0.105], s=46, facecolor=GOLD_LIGHT,
               edgecolor=GOLD, linewidth=1.0, zorder=4)
    ax.text(0.565, 0.105, r"nodal $\widetilde{\vartheta}_a=\vartheta_a-1$",
            fontsize=6.7, ha="left", va="center", color=INK)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")


def bvp_panel(ax: plt.Axes) -> None:
    """Draw the exact quarter-block geometry, load patch, and supports."""
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")
    ax.text(
        0.02,
        0.98,
        "(b)",
        ha="left",
        va="top",
        fontsize=9.0,
        color=INK,
        fontweight="bold",
    )

    # Isometric projection viewed from the (-x,-y) octant so BOTH symmetry
    # planes (x=0 and y=0) are the visible side faces of the quarter block.
    def iso(x: float, y: float, z: float) -> np.ndarray:
        return np.array([
            0.50 + 0.262 * (y - x),
            0.120 + 0.330 * z + 0.138 * (x + y),
        ])

    c000, c100, c010, c110 = iso(0, 0, 0), iso(1, 0, 0), iso(0, 1, 0), iso(1, 1, 0)
    c001, c101, c011, c111 = iso(0, 0, 1), iso(1, 0, 1), iso(0, 1, 1), iso(1, 1, 1)

    # Visible faces: y=0 symmetry (gold, left), x=0 symmetry (blue, right),
    # loaded top surface (white).
    ax.add_patch(Polygon(
        [c000, c100, c101, c001], closed=True,
        facecolor=GOLD_LIGHT, edgecolor=GOLD, linewidth=0.9, zorder=2,
    ))
    ax.add_patch(Polygon(
        [c000, c010, c011, c001], closed=True,
        facecolor=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.9, zorder=2,
    ))
    ax.add_patch(Polygon(
        [c001, c101, c111, c011], closed=True,
        facecolor="white", edgecolor=INK, linewidth=0.9, zorder=3,
    ))

    # Hidden edges of the far corner, dashed.
    for a, b in ((c100, c110), (c010, c110), (c110, c111)):
        ax.plot(*zip(a, b), color=MUTED, linewidth=0.6,
                linestyle=(0, (3.0, 2.2)), zorder=1)

    # Top-face partition at x=0.5 and y=0.5; the 0.5 x 0.5 follower patch
    # sits adjacent to both symmetry planes (nearest the viewer).
    ax.plot(*zip(iso(0.5, 0, 1), iso(0.5, 1, 1)), color=MUTED,
            linewidth=0.55, zorder=4)
    ax.plot(*zip(iso(0, 0.5, 1), iso(1, 0.5, 1)), color=MUTED,
            linewidth=0.55, zorder=4)
    ax.add_patch(Polygon(
        [iso(0, 0, 1), iso(0.5, 0, 1), iso(0.5, 0.5, 1), iso(0, 0.5, 1)],
        closed=True, facecolor=RED_LIGHT, edgecolor=RED,
        linewidth=1.0, zorder=5,
    ))

    # Follower-pressure arrows onto the patch.
    for xq, yq in ((0.14, 0.14), (0.38, 0.12), (0.12, 0.38)):
        tip = iso(xq, yq, 1)
        ax.add_patch(FancyArrowPatch(
            tip + np.array([0.0, 0.145]), tip + np.array([0.0, 0.020]),
            arrowstyle="-|>", mutation_scale=6.0, linewidth=0.95,
            color=RED, zorder=8,
        ))

    # Symmetry-face labels, rotated to lie along each face.
    ax.text(
        *(0.25 * (c000 + c100 + c101 + c001) + np.array([0.0, -0.01])),
        r"$y{=}0$: $u_2{=}0$",
        ha="center", va="center", fontsize=6.3, color="#7a5210",
        rotation=-27, zorder=6,
    )
    ax.text(
        *(0.25 * (c000 + c010 + c011 + c001) + np.array([0.0, -0.01])),
        r"$x{=}0$: $u_1{=}0$",
        ha="center", va="center", fontsize=6.3, color=BLUE,
        rotation=27, zorder=6,
    )

    # Base support: z=0 with u3=0, marked by hatch ticks under the two
    # visible bottom edges.
    for edge_a, edge_b in ((c000, c100), (c000, c010)):
        for frac in (0.22, 0.46, 0.70, 0.94):
            base = edge_a + frac * (edge_b - edge_a)
            ax.plot(
                [base[0], base[0] - 0.016],
                [base[1] - 0.004, base[1] - 0.030],
                color=MUTED, linewidth=0.6, zorder=2,
            )
    ax.text(
        0.98, 0.150, r"$z{=}0$: $u_3{=}0$",
        ha="right", va="center", fontsize=6.2, color=MUTED,
    )

    # In-plane restraint of the loaded top face.
    ax.text(
        0.98, 0.720, r"$z{=}1$: $u_1{=}u_2{=}0$",
        ha="right", va="center", fontsize=6.2, color=INK,
    )
    ax.plot(
        [0.795, c011[0] + 0.008], [0.712, c011[1] + 0.028],
        color=MUTED, linewidth=0.5, zorder=4,
    )

    ax.text(
        0.50,
        0.010,
        r"$1\times1\times1$ mm quarter block",
        ha="center",
        va="bottom",
        fontsize=6.4,
        color=MUTED,
    )


def result_panel(
    ax: plt.Axes,
    image: Image.Image,
    ramp: Image.Image,
    title: str,
    tick_labels: tuple[str, str, str],
    colorbar_label: str,
) -> None:
    """Place one supplied Abaqus model and its shortened supplied legend."""
    ax.axis("off")
    ax.text(0.02, 0.99, title, transform=ax.transAxes, ha="left",
            va="top", fontsize=9.0, fontweight="bold")

    model_ax = ax.inset_axes([0.005, 0.14, 0.990, 0.83])
    model_ax.imshow(image)
    model_ax.axis("off")

    ramp_ax = ax.inset_axes([0.12, 0.095, 0.76, 0.035])
    ramp_ax.imshow(ramp.rotate(-90, expand=True), aspect="auto")
    ramp_ax.axis("off")
    for x, label in zip((0.12, 0.50, 0.88), tick_labels):
        ax.text(
            x,
            0.073,
            label,
            ha="center",
            va="top",
            fontsize=6.5,
            color=INK,
        )
    ax.text(
        0.50,
        0.005,
        colorbar_label,
        ha="center",
        va="bottom",
        fontsize=7.0,
        color=INK,
    )


def main() -> None:
    missing = [path for path in (U_IMAGE, THETAT_IMAGE) if not path.exists()]
    if missing:
        names = ", ".join(str(path) for path in missing)
        raise FileNotFoundError(f"render Abaqus snapshots first; missing: {names}")

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    u_image = crop_model(U_IMAGE)
    u_ramp = crop_color_ramp(U_IMAGE)
    thetat_image = crop_model(THETAT_IMAGE)
    thetat_ramp = crop_color_ramp(THETAT_IMAGE)

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "serif",
        "font.size": 8,
        "axes.linewidth": 0.7,
        "mathtext.fontset": "dejavuserif",
    })
    # Match the natural PDF width to the manuscript text block so labels are
    # not needlessly down-scaled when LaTeX places the figure.
    fig = plt.figure(figsize=(6.50, 3.70))
    grid = fig.add_gridspec(
        1, 3, width_ratios=(1.02, 1.0, 1.0),
        left=0.015, right=0.995, top=0.985, bottom=0.135, wspace=0.025,
    )

    design = grid[0].subgridspec(
        2,
        1,
        height_ratios=(0.92, 1.08),
        hspace=0.035,
    )
    element_panel(fig.add_subplot(design[0]))
    bvp_panel(fig.add_subplot(design[1]))
    result_panel(
        fig.add_subplot(grid[1]), u_image, u_ramp,
        "(c)",
        ("0", "0.348", "0.696"),
        r"$\|\mathbf{u}\|$ (mm)",
    )
    result_panel(
        fig.add_subplot(grid[2]), thetat_image, thetat_ramp,
        "(d)",
        ("$-7.5$", "$-3.25$", "$+1.0$"),
        r"$\widetilde{\vartheta}=\vartheta-1$ ($\times 10^{-4}$)",
    )

    pdf = FIG_DIR / "fig4_utheta_block.pdf"
    fig.savefig(pdf, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"Wrote {pdf}")


if __name__ == "__main__":
    main()
