"""Build Figure 6 from the six supplied Abaqus/CAE pasta snapshots.

The exports show ``UVARM1``, which the Hex8 local-pressure implementation
stores as polymer volume fraction ``phi``.  This script does not reconstruct a
solver field.  It crops the same viewport from every supplied image, retains
each native Abaqus color ramp, and replaces the embedded 13-value legend with
three short values.  The contour limits differ by frame and therefore must not
be interpreted as a common scale.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
ABAQUS = FIGURES / "abaqus"
OUTPUT = FIGURES / "fig6_pasta.pdf"

# The supplied 4000 x 2895 exports use an identical Abaqus viewport.  A common
# crop removes the embedded legend and triad while preserving the relative
# scale and position of the deforming model across all six frames.
VIEWPORT_BOX = (985, 805, 3650, 2395)

# Exact colored portion of the embedded Abaqus 13-band legend.  Rotating it
# clockwise places the low value at left and high value at right.
RAMP_BOX = (164, 418, 293, 1440)

SNAPSHOTS = (
    ("a", 0, ABAQUS / "pasta_t0.png", (0.9900, 0.9900, 0.9900)),
    ("b", 45, ABAQUS / "pasta_t45.png", (0.5679, 0.7743, 0.9807)),
    ("c", 90, ABAQUS / "pasta_t90.png", (0.2729, 0.5931, 0.9132)),
    ("d", 150, ABAQUS / "pasta_t150.png", (0.2039, 0.5123, 0.8207)),
    ("e", 200, ABAQUS / "pasta_t200.png", (0.1898, 0.4725, 0.7552)),
    ("f", 360, ABAQUS / "pasta_t360.png", (0.1764, 0.4026, 0.6288)),
)

INK = "#303030"
MUTED = "#666666"


def crop_model(path: Path) -> Image.Image:
    """Return the unchanged model pixels within the common viewport."""
    cropped = Image.open(path).convert("RGB").crop(VIEWPORT_BOX)
    # The right edge of the original legend box enters the common viewport at
    # its upper-left corner.  Remove only that empty-viewport annotation; this
    # rectangle is disjoint from the model in every supplied frame.
    cropped.paste("white", (0, 0, 50, 730))
    return cropped


def crop_ramp(path: Path) -> Image.Image:
    """Return the supplied Abaqus ramp without its long numeric labels."""
    return (
        Image.open(path)
        .convert("RGB")
        .crop(RAMP_BOX)
        .rotate(-90, expand=True)
    )


def snapshot_panel(
    ax: plt.Axes,
    label: str,
    time: int,
    path: Path,
    values: tuple[float, float, float],
) -> None:
    """Place one native viewport and its compact panel-specific legend."""
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis("off")
    ax.text(
        0.01,
        0.985,
        rf"({label})",
        ha="left",
        va="top",
        fontsize=8.0,
        color=INK,
        fontweight="bold",
    )

    result = ax.inset_axes([0.005, 0.185, 0.990, 0.755])
    result.imshow(crop_model(path), interpolation="none")
    result.axis("off")

    low, middle, high = values
    if low == high:
        ax.text(
            0.5,
            0.075,
            rf"$\phi={low:.3f}$",
            ha="center",
            va="center",
            fontsize=6.5,
            color=INK,
        )
        return

    legend = ax.inset_axes([0.205, 0.095, 0.590, 0.040])
    legend.imshow(crop_ramp(path), aspect="auto", interpolation="none")
    legend.axis("off")
    for x, value, alignment in (
        (0.205, low, "left"),
        (0.500, middle, "center"),
        (0.795, high, "right"),
    ):
        ax.text(
            x,
            0.060,
            f"{value:.3f}",
            ha=alignment,
            va="center",
            fontsize=5.9,
            color=INK,
        )


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    missing = [path for _, _, path, _ in SNAPSHOTS if not path.is_file()]
    if missing:
        missing_text = ", ".join(
            str(path.relative_to(ROOT)) for path in missing
        )
        raise FileNotFoundError(
            f"Missing Figure 6 source export(s): {missing_text}"
        )

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "serif",
        "font.size": 8,
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.7,
    })
    # The formulation and field identity (Hex8, UVARM1 = phi, panel-specific
    # ranges) are stated in the manuscript caption, so the figure itself
    # carries no duplicate title.
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(6.50, 4.25),
        gridspec_kw={
            "left": 0.008,
            "right": 0.995,
            "top": 0.985,
            "bottom": 0.050,
            "wspace": 0.020,
            "hspace": 0.080,
        },
    )

    for ax, snapshot in zip(axes.flat, SNAPSHOTS):
        snapshot_panel(ax, *snapshot)

    fig.savefig(OUTPUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
