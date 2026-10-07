import heatSource from "../../examples/scalar_diffusion_uel/build.py?raw";
import plateSource from "../../tools/check_heated_plate.py?raw";
import bendingSource from "../../tools/check_thermal_bending.py?raw";
import compiledSource from "../../abaqus_ufl/fe/compiled_element.py?raw";
import fortranSource from "../../examples/scalar_diffusion_uel/scalar_diffusion_uel.for?raw";

export const repository = "https://github.com/tengzhang48/abaqus_ufl";
export const steps = [
  {
    id: "equations",
    title: "Equations",
    detail: "Physics, fields & boundaries",
  },
  { id: "weak-form", title: "Weak form", detail: "Test functions & residuals" },
  { id: "python", title: "Python", detail: "Declare the model" },
  { id: "fortran", title: "Generated Fortran", detail: "Build the Abaqus UEL" },
  {
    id: "simulation",
    title: "Simulation",
    detail: "Results & independent checks",
  },
] as const;
export type StepId = (typeof steps)[number]["id"];

export type Equation = {
  title: string;
  expression: string;
  explanation: string;
};
export type Walkthrough = {
  title: string;
  summary: string;
  kind: string;
  domain: string;
  fields: string;
  boundaries: string[];
  properties: string;
  equations: Equation[];
  weakForms: Equation[];
  testFunctions: string;
  scope: string;
  reference: string;
  sourcePath: string;
};

const heatEquation: Equation = {
  title: "Conservation of heat",
  expression: "ρcₚ ∂T/∂t + Div_X q = 0,     q = −k Grad_X T",
  explanation:
    "Temperature changes because heat flows down its gradient. The public model uses reference coordinates X and has no volumetric heat source.",
};
const heatResidual: Equation = {
  title: "Discrete thermal residual",
  expression:
    "Rₜ(θ) = ∫Ω₀ [θ ρcₚ(Tⁿ − Tⁿ⁻¹)/Δt + k Grad_X θ · Grad_X Tⁿ] dΩ₀ = 0",
  explanation:
    "Multiply the heat balance by θ and integrate the flux divergence by parts over the reference domain Ω₀. Insert q = −k Grad_X T and use backward Euler for storage. The boundary term vanishes where θ = 0 on prescribed-temperature edges, or where q · N = 0 on insulated edges.",
};

export const walkthroughs: Record<string, Walkthrough> = {
  heated_plate_bvp: {
    title: "A plate heated from one edge",
    summary:
      "Start with heat conservation. Build its weak form and Python declaration, generate the user element, then follow heat into the plate.",
    kind: "Heat transfer · Quad4 UEL",
    domain: "Unit square, 0 ≤ X ≤ 1 and 0 ≤ Y ≤ 1. Consistent model units.",
    fields:
      "T(X,t): temperature rise. Displacement is fixed throughout the mesh for this heat-transfer example.",
    boundaries: [
      "Initially T = 0.",
      "Top: T = sin(πX)(1 − exp(−2t)).",
      "Bottom, left and right: T = 0.",
    ],
    properties:
      "k = 0.5 and ρcₚ = 1. The reused coupled declaration has G = 1, K = 10 and α = 0.001; its displacement DOFs are prescribed here.",
    equations: [heatEquation],
    weakForms: [heatResidual],
    testFunctions:
      "Choose θ = 0 on all four prescribed-temperature edges. The interior test functions are free. The fixed displacement field contributes no free mechanical equations in this example.",
    scope:
      "A transient heat-transfer boundary-value problem on a 24 × 24 mesh, with 529 interior temperatures solved over 20 increments of Δt = 0.05.",
    reference:
      "Independent Fourier series check temperature and stored heat. Mesh refinement checks second-order spatial accuracy; step refinement checks first-order time accuracy. Prescribed values, residuals, supplied heat and thermal dissipation are checked at every increment.",
    sourcePath: "tools/check_heated_plate.py",
  },
  thermal_bending_bvp: {
    title: "A strip that bends as it heats",
    summary:
      "A temperature difference produces thermal stress. Follow the coupled equations into a generated user element, then see the clamped strip bend.",
    kind: "Thermomechanics · Quad4 UEL",
    domain:
      "Plane-strain strip, 0 ≤ X ≤ 8 and 0 ≤ Y ≤ 1. Consistent model units.",
    fields:
      "u(X,t): displacement. T(X,t): temperature rise. F = I + Grad_X u and J = det F.",
    boundaries: [
      "Initially u = 0 and T = 0.",
      "Left end: u₁ = u₂ = 0. Remaining mechanical edges: zero traction.",
      "Top/bottom: T = ±(1 − exp(−4t)). Left/right: q · N = 0 (insulated).",
    ],
    properties:
      "G = 1, K = 2, α = 0.001, k = 0.5 and ρcₚ = 1. K is the logarithmic volumetric coefficient; α controls thermal pressure in the stated law.",
    equations: [
      {
        title: "Mechanical equilibrium",
        expression: "Div_X P = 0",
        explanation:
          "Mechanical equilibrium is quasi-static: each temperature increment has a converged displacement field. Div_X acts in the reference coordinates.",
      },
      {
        title: "Stress law",
        expression: "P = G(F − F⁻ᵀ) + K ln(J) F⁻ᵀ − KαT F⁻ᵀ",
        explanation:
          "The first two terms supply finite-strain elasticity. The last term couples temperature to the first Piola stress.",
      },
      heatEquation,
    ],
    weakForms: [
      {
        title: "Mechanical residual",
        expression: "Rᵤ(v) = ∫Ω₀ Grad_X v : P(F,T) dΩ₀ = 0",
        explanation:
          "Multiply equilibrium by v and integrate by parts over the reference domain Ω₀. The traction term vanishes on the free edges; v vanishes at the clamp.",
      },
      heatResidual,
    ],
    testFunctions:
      "Choose v = 0 at the clamped left end and θ = 0 on the top and bottom baths. The insulated end faces contribute no thermal-boundary term. Solve both residuals for every free displacement and temperature DOF.",
    scope:
      "A 48 × 8 mesh solves 343 thermal and 864 displacement DOFs over 20 increments. Heat drives finite-strain deformation through one-way coupling; this material has no mechanical-heating feedback or deformation-dependent transport.",
    reference:
      "An independent Fourier series checks temperature. A slender-beam approximation uses the independent thermal moment to check tip deflection (3% final tolerance; it omits clamp end effects and finite strain). Refinement, thermal dissipation, supplied heat and clamp force/moment balance provide separate checks.",
    sourcePath: "tools/check_thermal_bending.py",
  },
};

export const supportCases: Record<
  string,
  { title: string; summary: string; reference: string }
> = {
  serial_fe: {
    title: "Heat diffusion on a mesh",
    summary:
      "A one-step sinusoidal diffusion check and affine mechanical patch.",
    reference:
      "A closed-form Q1 consistent-mass eigenmode checks the discrete solve. Refinement approaches the continuum backward-Euler solution at second order.",
  },
  neo_hookean_umat: {
    title: "Neo-Hookean elasticity",
    summary:
      "Stretch and shear a material point to exercise elastic stress and its tangent.",
    reference:
      "Hand-derived Cauchy stress formulas and an independent finite-difference check of the compiled material tangent.",
  },
  ogden_umat: {
    title: "Ogden elasticity",
    summary:
      "An isochoric stretch sweep exercises the spectral material response.",
    reference:
      "Principal-stretch powers supply an independent stress reference. Separate checks cover rotated states and repeated-spectrum tangents.",
  },
  small_strain_j2_umat: {
    title: "Yield and plastic hardening",
    summary:
      "A material point follows shear loading through elastic yield into plastic flow.",
    reference:
      "The exact proportional-loading J2 consistency solution with linear isotropic hardening checks every increment and returned state.",
  },
  small_strain_viscoelastic_umat: {
    title: "Viscoelastic relaxation",
    summary:
      "Apply a step in shear, hold strain, and follow stress relaxation at a material point.",
    reference:
      "The exact backward-Euler recurrence checks the compiled history. A continuous curve shows the time-step approximation.",
  },
  scalar_diffusion_uel: {
    title: "Thermal expansion and storage",
    summary:
      "A constrained Quad4 checks thermal forces and uniform heat storage at element level.",
    reference:
      "Closed-form shape-gradient integrals check thermal loading; conservation gives the exact signed storage rate.",
  },
  thermo_mechanics_quad8: {
    title: "Heat flux under deformation",
    summary:
      "A stretched mixed-order Quad8 checks the pull-back of a linear thermal gradient.",
    reference:
      "The C⁻¹ pull-back gives a thermal-load ratio of 1/λ², checked against the compiled element.",
  },
};

function section(source: string, start: string, end: string) {
  const first = source.indexOf(start),
    last = source.indexOf(end, first + start.length);
  if (first < 0 || last < 0)
    throw new Error("A walkthrough source excerpt is missing.");
  return source.slice(first, last).trim();
}

// Display excerpts from the shipped source, rather than maintaining example code twice.
export const sourceExcerpts = {
  material: section(
    heatSource,
    "class HeatDiffusionMaterial",
    "class HeatDiffusionProblem",
  ),
  weakForm: section(
    heatSource,
    "class HeatDiffusionProblem",
    "def verification_state",
  ),
  generation: section(
    compiledSource,
    "    source = build /",
    "    result = subprocess.run",
  ),
  interface: fortranSource.split("\n").slice(0, 22).join("\n"),
  boundaries: {
    heated_plate_bvp: section(
      plateSource,
      "def conditions",
      "def assert_plate_accuracy",
    ),
    thermal_bending_bvp: section(
      bendingSource,
      "def conditions",
      "def tip_values",
    ),
  } as Record<string, string>,
};
