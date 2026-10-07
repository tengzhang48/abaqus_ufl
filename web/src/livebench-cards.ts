export const modelCards = [
  {
    id: "heated_plate_bvp",
    title: "Heat entering a plate",
    kind: "Heat transfer",
    summary: "Follow heat conservation into a transient temperature field.",
    plot: "temperature",
    alt: "Computed temperature field in a square plate heated along its top edge",
    caption: "Temperature · final recorded step",
  },
  {
    id: "thermal_bending_bvp",
    title: "A strip that bends as it heats",
    kind: "Thermomechanics",
    summary: "Connect heat diffusion and thermal stress to a bending strip.",
    plot: "deformation",
    alt: "Computed deformation of a thermally loaded strip with its left end clamped",
    caption: "Displacement magnitude · deformation ×10",
  },
  {
    id: "ogden_bvp",
    title: "Rubber under finite shear",
    kind: "Ogden elasticity",
    summary:
      "Turn a strain-energy function into a finite-strain mesh simulation.",
    plot: "deformation",
    alt: "Computed finite-strain deformation of an Ogden block under shear",
    caption: "Displacement magnitude · deformation ×1",
  },
  {
    id: "plasticity_bvp",
    title: "Yield, harden, then unload",
    kind: "J2 plasticity",
    summary:
      "See how a stress update and material history govern a loading cycle.",
    plot: "plastic-strain",
    alt: "Nodal projection of computed Gauss-point equivalent plastic strain in a sheared block",
    caption: "Equivalent plastic strain · nodal projection",
  },
] as const;
