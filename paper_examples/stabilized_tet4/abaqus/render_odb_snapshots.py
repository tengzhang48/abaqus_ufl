"""Render publication snapshots from the completed Figure 14 Abaqus ODB.

Run with Abaqus/CAE rather than standard Python. On installations that do not
forward arguments through the Abaqus wrapper, use the environment variables::

    export SCOVAZZI_ODB=/scratch/user/job/model.odb
    export SCOVAZZI_SNAPSHOT_DIR=/path/to/figures/source
    abaqus cae noGUI=render_odb_snapshots.py

The script isolates the C3D4T visualization overlay, plots the final deformed
configuration at true scale, and writes undecorated high-resolution PNGs. The
manuscript assembly script adds consistent titles and color bars.
"""

from __future__ import print_function

import os
import sys

from abaqus import session
from abaqusConstants import (
    ALL, COLOR, CONTOURS_ON_DEF, FILLED, INVARIANT, NODAL, OFF, PNG,
    SOLID, UNIFORM,
)
from caeModules import dgo
import visualization


U_MIN = 0.0
U_MAX = 0.8
THETAT_MIN = -7.5e-4
THETAT_MAX = 1.0e-4
IMAGE_SIZE = (1800, 1500)


def command_arguments():
    if "--" in sys.argv:
        args = sys.argv[sys.argv.index("--") + 1:]
    else:
        args = sys.argv[1:]
    if len(args) != 2:
        odb_path = os.environ.get("SCOVAZZI_ODB")
        output_directory = os.environ.get("SCOVAZZI_SNAPSHOT_DIR")
        if odb_path and output_directory:
            return os.path.abspath(odb_path), os.path.abspath(output_directory)
        raise RuntimeError(
            "expected ODB_PATH OUTPUT_DIRECTORY; sys.argv={!r}".format(
                sys.argv))
    return os.path.abspath(args[0]), os.path.abspath(args[1])


def configure_viewport(odb):
    name = "Figure 14 publication render"
    if name in session.viewports:
        del session.viewports[name]
    viewport = session.Viewport(name=name, origin=(0, 0), width=180, height=145)
    viewport.setValues(displayedObject=odb)
    viewport.odbDisplay.setFrame(step=0, frame=-1)

    leaf = dgo.LeafFromElementSets(
        elementSets=("PART-1-1.VISUALIZATION",))
    viewport.odbDisplay.displayGroup.replace(leaf=leaf)

    viewport.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF,))
    viewport.odbDisplay.commonOptions.setValues(
        renderStyle=FILLED,
        visibleEdges=ALL,
        deformationScaling=UNIFORM,
        uniformScaleFactor=1.0,
    )

    # Look toward the loaded symmetry corner x=y=0, matching the orientation
    # of the unstructured-mesh panel in the source paper.
    viewport.view.setValues(
        cameraPosition=(-3.0, -3.0, 2.35),
        cameraTarget=(0.5, 0.5, 0.35),
        cameraUpVector=(0.0, 0.0, 1.0),
    )
    viewport.view.fitView()
    viewport.viewportAnnotationOptions.setValues(
        triad=OFF,
        legend=OFF,
        title=OFF,
        state=OFF,
        annotations=OFF,
        compass=OFF,
    )
    return viewport


def configure_spectrum():
    name = "EML viridis"
    if name not in session.spectrums:
        session.Spectrum(
            name=name,
            colors=(
                "#440154", "#482878", "#3E4989", "#31688E",
                "#26828E", "#1F9E89", "#35B779", "#6DCD59",
                "#B4DE2C", "#FDE725",
            ),
        )
    return name


def render(viewport, output_base, variable, lower, upper, spectrum,
           invariant=None):
    if invariant is None:
        viewport.odbDisplay.setPrimaryVariable(
            variableLabel=variable, outputPosition=NODAL)
    else:
        viewport.odbDisplay.setPrimaryVariable(
            variableLabel=variable,
            outputPosition=NODAL,
            refinement=(INVARIANT, invariant),
        )
    viewport.odbDisplay.contourOptions.setValues(
        spectrum=spectrum,
        intervalType=UNIFORM,
        numIntervals=12,
        minAutoCompute=OFF,
        maxAutoCompute=OFF,
        minValue=lower,
        maxValue=upper,
        outsideLimitsAboveColor="#FDE725",
        outsideLimitsBelowColor="#440154",
    )
    session.printToFile(
        fileName=output_base,
        format=PNG,
        canvasObjects=(viewport,),
    )


def main():
    odb_path, output_directory = command_arguments()
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    session.graphicsOptions.setValues(
        backgroundStyle=SOLID, backgroundColor="#FFFFFF")
    session.printOptions.setValues(
        rendition=COLOR, vpDecorations=OFF, vpBackground=OFF)
    session.pngOptions.setValues(imageSize=IMAGE_SIZE)

    odb = session.openOdb(name=odb_path, readOnly=True)
    try:
        viewport = configure_viewport(odb)
        spectrum = configure_spectrum()
        render(
            viewport,
            os.path.join(output_directory, "scovazzi_fig14_u_magnitude"),
            "U", U_MIN, U_MAX, spectrum, invariant="Magnitude",
        )
        render(
            viewport,
            os.path.join(output_directory, "scovazzi_fig14_thetat"),
            "NT11", THETAT_MIN, THETAT_MAX, spectrum,
        )
    finally:
        odb.close()

    print("Wrote Abaqus snapshots to {}".format(output_directory))


if __name__ == "__main__":
    main()
