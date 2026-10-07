import heatSource from "../../examples/scalar_diffusion_uel/build.py?raw";
import plateSource from "../../tools/check_heated_plate.py?raw";
import bendingSource from "../../tools/check_thermal_bending.py?raw";
import fortranSource from "../../examples/scalar_diffusion_uel/scalar_diffusion_uel.for?raw";
import ogdenMaterialSource from "../../examples/ogden_umat/build.py?raw";
import ogdenSource from "../../tools/check_ogden_bvp.py?raw";
import ogdenFortran from "../../examples/ogden_umat/ogden_plane_strain_uel.for?raw";
import plasticityMaterialSource from "../../examples/small_strain_j2_umat/build.py?raw";
import plasticityHostSource from "../../tools/plasticity_bvp_common.py?raw";
import plasticityFortran from "../../examples/small_strain_j2_umat/small_strain_j2.for?raw";

export const repository = "https://github.com/tengzhang48/abaqus_ufl";
export const steps = [
  {
    id: "equations",
    title: "Equations",
    detail: "Physics, fields & boundaries",
  },
  { id: "weak-form", title: "Weak form", detail: "Test functions & residuals" },
  { id: "python", title: "Python", detail: "Declare the model" },
  {
    id: "fortran",
    title: "Generated Fortran",
    detail: "Build the Abaqus subroutine",
  },
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
  constitutive?: "ogden" | "plasticity";
};

const heatEquation: Equation = {
  title: "Conservation of heat",
  expression: "ρcₚ ∂T/∂t + Div q = 0,     q = −k Grad T",
  explanation:
    "Temperature changes because heat flows down its gradient. Div and Grad act on the reference coordinates X; there is no volumetric heat source.",
};
const heatResidual: Equation = {
  title: "Discrete thermal residual",
  expression:
    "Rₜ(θ) = ∫Ω₀ [θ ρcₚ(Tⁿ − Tⁿ⁻¹)/Δt + k Grad θ · Grad Tⁿ] dΩ₀ = 0",
  explanation:
    "Multiply the heat balance by θ and integrate the flux divergence by parts over the reference domain Ω₀. Insert q = −k Grad T and use backward Euler for storage. The boundary term vanishes where θ = 0 on prescribed-temperature edges, or where q · N = 0 on insulated edges.",
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
      "u(X,t): displacement. T(X,t): temperature rise. F = I + Grad u and J = det F.",
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
        expression: "Div P = 0",
        explanation:
          "Mechanical equilibrium is quasi-static: each temperature increment has a converged displacement field. Div and Grad act on the reference coordinates X.",
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
        expression: "Rᵤ(v) = ∫Ω₀ Grad v : P(F,T) dΩ₀ = 0",
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
  ogden_bvp: {
    title: "Rubber under finite shear",
    summary:
      "Start with Ogden strain energy. Derive stress and the mechanical weak form, generate a displacement-only user element, then shear a block with free sides.",
    kind: "Ogden elasticity · finite-strain Quad4 UEL",
    constitutive: "ogden",
    domain:
      "Unit-square reference domain Ω₀, plane strain and unit out-of-plane thickness. Consistent model units.",
    fields:
      "u(X,s): displacement at quasi-static load parameter s. F = I + Grad u, F₃₃ = 1, and J = det F; Grad acts on the reference coordinates X. The full three-dimensional stretch spectrum enters the material law.",
    boundaries: [
      "Initially u = 0.",
      "Bottom: u₁ = u₂ = 0. Top: u₁ = 0.6s and u₂ = 0, for 0 ≤ s ≤ 1.",
      "Left and right: P N = 0. No body force.",
    ],
    properties:
      "μ = 1, α = 3.5 and K = 10. This simulation explicitly overrides the shipped Ogden material's default K = 100; it demonstrates moderate compressibility with a displacement-only Q1 element.",
    equations: [
      {
        title: "Strain energy",
        expression:
          "W(F) = (2μ/α²)[Σᵢ λ̄ᵢᵅ − 3] + (K/2)[ln J]²,     λ̄ᵢ = J⁻¹ᐟ³λᵢ",
        explanation:
          "The principal stretches λᵢ come from F. Isochoric stretch powers describe distortional elasticity; the logarithmic term penalizes volume change. All three principal stretches contribute, including the plane-strain direction.",
      },
      {
        title: "Stress and equilibrium",
        expression: "P = ∂W/∂F,     Div P = 0",
        explanation:
          "Differentiate energy to obtain first Piola stress. Each load increment satisfies equilibrium in the reference domain. This model is stateless: unloading would retrace its elastic response.",
      },
    ],
    weakForms: [
      {
        title: "Mechanical residual",
        expression: "Rᵤ(v) = ∫Ω₀ Grad v : P(I + Grad u) dΩ₀ = 0",
        explanation:
          "Multiply reference equilibrium by v and integrate by parts. The boundary term vanishes because the sides are traction-free and v = 0 at both prescribed grips.",
      },
    ],
    testFunctions:
      "Choose v = 0 on the top and bottom grips. Side and interior displacement test functions remain free. No temperature field or material history is solved.",
    scope:
      "A non-affine, plane-strain finite-shear mesh solve. The shipped Ogden constitutive law is reused in a displacement-only verification UEL; the material example also ships a UMAT. The load parameter is not physical time. These are compiled local verification runs.",
    reference:
      "Independent principal-stretch stress and energy, integrated with separately written quadrature, check equilibrium, reactions and virtual-grip work. A closed-form shear material check and an α = 2 eigenvalue-free limit exercise the constitutive response. Positive Jacobians, force/current-coordinate moment balance and mesh refinement supply separate gates.",
    sourcePath: "tools/check_ogden_bvp.py",
  },
  plasticity_bvp: {
    title: "Yield, harden, then unload",
    summary:
      "Follow an incremental J2 stress update into a mesh loading cycle. Track elastic loading, plastic flow, unloading, reverse yielding, and accepted material history.",
    kind: "J2 plasticity · compiled UMAT in a Quad4 host",
    constitutive: "plasticity",
    domain:
      "Unit-square small-strain domain, plane strain and unit out-of-plane thickness. Consistent model units.",
    fields:
      "u(X,t): displacement. ε = sym Grad u with ε₃₃ = 0. Every Gauss point retains all six stress components and equivalent plastic strain eₚ in STATEV(1).",
    boundaries: [
      "Initially u = 0, σ = 0 and eₚ = 0.",
      "Bottom: u₁ = u₂ = 0. Top: u₁ = 0.02t for 0 ≤ t ≤ 1, then 0.02(2 − t) for 1 < t ≤ 2; u₂ = 0. t is load-cycle pseudo-time.",
      "Left and right: σ N = 0. No body force.",
    ],
    properties:
      "G = 10, λ = 20, σᵧ = 0.1 and H = 5. The host uses the UMAT's tension-positive Abaqus stresses and engineering shear γ₁₂ = 2ε₁₂.",
    equations: [
      {
        title: "Small-strain equilibrium",
        expression: "Div σ = 0,     ε = ½(Grad u + Grad uᵀ)",
        explanation:
          "The case-local host imposes plane strain and evaluates the full three-dimensional material update. Stress σ₃₃ remains part of the deviatoric response.",
      },
      {
        title: "Elastic predictor and yield",
        expression:
          "σᵗʳ = σⁿ + 2G Δε + λ tr(Δε) I,     qᵗʳ = √(3/2 sᵗʳ:sᵗʳ),     fᵗʳ = qᵗʳ − (σᵧ + H eₚⁿ)",
        explanation:
          "The strain increment is Δε = ε(u) − εⁿ, using the accepted old strain. Trial stress uses the accepted old stress. The deviator sᵗʳ removes its hydrostatic part. A nonpositive yield function gives an elastic step.",
      },
      {
        title: "Plastic radial return",
        expression:
          "Δγ = max(0, fᵗʳ)/(3G + H),     n = 3sᵗʳ/(2qᵗʳ),     σⁿ⁺¹ = σᵗʳ − 2G Δγ n,     eₚⁿ⁺¹ = eₚⁿ + Δγ",
        explanation:
          "Evaluate n only on the plastic branch, where qᵗʳ > 0. Plastic flow returns the trial stress to the expanded yield surface. Equivalent plastic strain can grow during loading and reverse yielding. Unconverged trials never become accepted history.",
      },
    ],
    weakForms: [
      {
        title: "Incremental mechanical residual",
        expression: "Rᵤ(v) = ∫Ω ε(v) : σⁿ⁺¹(u; σⁿ, εⁿ, eₚⁿ) dΩ = 0",
        explanation:
          "Integrate equilibrium against admissible test displacement v. Free sides supply no traction term. The constitutive update is evaluated from fixed, accepted Gauss-point history throughout the global Newton solve.",
      },
      {
        title: "Consistent element tangent",
        expression: "Kₑ = Σᵍ wᵍ det Jᵍ Bᵍᵀ Dᵍ Bᵍ",
        explanation:
          "Sum over the four Gauss points g with weights wᵍ. Dᵍ is the algorithmic tangent DDSDDE returned by the UMAT, and B maps nodal displacements to the active engineering-strain components. The UMAT supplies its full 6 × 6 DDSDDE; the plane-strain host selects indices xx, yy and xy without a plane-stress condensation.",
      },
    ],
    testFunctions:
      "Choose v = 0 at both grips. Solve the free side and interior displacement DOFs. State variables remain local to integration points; they are not independent nodal unknowns.",
    scope:
      "A case-local small-strain Quad4 verification host calls the shipped generated 3D UMAT at every Gauss point. Accepted stress, strain and STATEV are committed only after equilibrium. Plastic-strain images are a volume-weighted nodal projection of Gauss-point state. This case does not add a generic stateful FE runtime or claim an Abaqus job.",
    reference:
      "An independent radial-return implementation checks stresses, state and global residuals. Closed-form elastic and proportional-shear checks cover conventions. The exact discrete endpoint-work balance checks elastic/hardening energy, physical plastic dissipation and numerical dissipation. Yield consistency, unloading, reverse plasticity, state rollback, force/moment balance and mesh refinement are checked separately.",
    sourcePath: "tools/check_plasticity_bvp.py",
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

// Python excerpts drop their common indent. Fortran keeps its fixed-form
// columns, so continuation markers stay in column 6.
function section(source: string, start: string, end: string, dedent = true) {
  const marker = source.indexOf(start),
    last = source.indexOf(end, marker + start.length);
  if (marker < 0 || last < 0)
    throw new Error("A walkthrough source excerpt is missing.");
  const lines = source
    .slice(source.lastIndexOf("\n", marker) + 1, last)
    .trimEnd()
    .split("\n");
  const indent = dedent
    ? Math.min(
        ...lines
          .filter((line) => line.trim())
          .map((line) => line.length - line.trimStart().length),
      )
    : 0;
  return lines.map((line) => line.slice(indent)).join("\n");
}

const gaussPointBlock = (source: string) =>
  section(
    source,
    "C       Real material evaluation (for RHS)",
    "C       AMATRX assembly",
    false,
  );

// Display excerpts from the shipped source, rather than maintaining example code twice.
export type CodeBlock = { title: string; path: string; code: string };

const uelInterface = (source: string) =>
  section(
    source,
    "      SUBROUTINE UEL(",
    "C     --- Local parameters ---",
    false,
  );

const heatPython: CodeBlock[] = [
  {
    title: "Material law",
    path: "examples/scalar_diffusion_uel/build.py",
    code: section(
      heatSource,
      "class HeatDiffusionMaterial",
      "class HeatDiffusionProblem",
    ),
  },
  {
    title: "Weak-form declaration",
    path: "examples/scalar_diffusion_uel/build.py",
    code: section(heatSource, "class HeatDiffusionProblem", "def verification_state"),
  },
];

const thermalCode = {
  pythonIntro:
    "The material returns stress, heat storage and heat flux. The weak-form declaration names the fields and their interpolation; the generator builds both residuals and their tangents from these methods.",
  python: heatPython,
  generationTitle: "Generate the Abaqus UEL",
  generationSummary:
    "The declaration generates a standard Quad4 user element. Each of its four nodes carries u₁, u₂ and T, giving 12 element DOFs. At every Gauss point the element evaluates the material law, assembles both residuals, and differentiates them by complex step for the tangent blocks. The simulation compiles this shipped file and checks that it matches byte for byte; material values arrive as element properties.",
  fortranPath: "examples/scalar_diffusion_uel/scalar_diffusion_uel.for",
  fortranExcerpt: gaussPointBlock(fortranSource),
  fortranCaption: "Gauss-point residual assembly",
  compilePath: "abaqus_ufl/fe/compiled_element.py",
  interface: uelInterface(fortranSource),
  interfaceTitle: "Inspect the generated UEL interface",
  runtimeNote:
    "The UEL returns RHS = −R and AMATRX = ∂R/∂U. The local driver preserves accepted temperature history throughout each Newton solve. Mesh, loads, and steps belong to the analysis setup.",
};

export const modelCode: Record<
  string,
  {
    pythonIntro: string;
    python: CodeBlock[];
    pythonMore?: CodeBlock;
    pythonNote?: string;
    boundaries: CodeBlock;
    generationTitle: string;
    generationSummary: string;
    fortranPath: string;
    fortranExcerpt: string;
    fortranCaption: string;
    compilePath: string;
    interface: string;
    interfaceTitle: string;
    runtimeNote: string;
  }
> = {
  heated_plate_bvp: {
    ...thermalCode,
    boundaries: {
      title: "Boundary values in the simulation script",
      path: "tools/check_heated_plate.py",
      code: section(plateSource, "def conditions", "def assert_plate_accuracy"),
    },
  },
  thermal_bending_bvp: {
    ...thermalCode,
    pythonNote:
      "This case constructs HeatDiffusionProblem(K=2.0). The material excerpt shows the shared declaration's defaults; the simulation overrides K.",
    boundaries: {
      title: "Boundary values in the simulation script",
      path: "tools/check_thermal_bending.py",
      code: section(bendingSource, "def conditions", "def tip_values"),
    },
  },
  ogden_bvp: {
    pythonIntro:
      "The material returns first Piola stress from the Ogden strain energy, through the principal stretches of C. The weak-form declaration reuses that law for a displacement-only plane-strain element; the generator performs element integration and tangent construction.",
    python: [
      {
        title: "Material law",
        path: "examples/ogden_umat/build.py",
        code: section(
          ogdenMaterialSource,
          "class OgdenOneTerm",
          "class OgdenPlaneStrainProblem",
        ),
      },
      {
        title: "Weak-form declaration",
        path: "examples/ogden_umat/build.py",
        code: section(
          ogdenMaterialSource,
          "class OgdenPlaneStrainProblem",
          "def generate(",
        ),
      },
    ],
    pythonNote:
      "The shipped material default is K = 100. This simulation constructs OgdenPlaneStrainProblem(mu=1, alpha=3.5, K=10).",
    boundaries: {
      title: "Boundary values in the simulation script",
      path: "tools/check_ogden_bvp.py",
      code: section(ogdenSource, "def conditions", "def solve_history"),
    },
    generationTitle: "Generate the Abaqus UEL",
    generationSummary:
      "The plane-strain declaration generates a standard Quad4 UEL with eight element DOFs. Each Gauss point evaluates the spectral Ogden stress and its complex-step tangent, then assembles the residual and stiffness. The simulation compiles this shipped file and checks that it matches byte for byte; μ, α and K arrive as element properties.",
    fortranPath: "examples/ogden_umat/ogden_plane_strain_uel.for",
    fortranExcerpt: gaussPointBlock(ogdenFortran),
    fortranCaption: "Gauss-point residual assembly",
    compilePath: "abaqus_ufl/fe/compiled_element.py",
    interface: uelInterface(ogdenFortran),
    interfaceTitle: "Inspect the generated UEL interface",
    runtimeNote:
      "The UEL returns RHS = −R and AMATRX = ∂R/∂U. Its only nodal field is displacement, with no material state. The same law also ships as a three-dimensional UMAT, examples/ogden_umat/ogden_umat.for.",
  },
  plasticity_bvp: {
    pythonIntro:
      "The material declares an incremental stress update with one state variable. The generator writes it as a UMAT and differentiates the update for DDSDDE. A case-local Quad4 host assembles element equilibrium; it is separate from the generator.",
    python: [
      {
        title: "Material law",
        path: "examples/small_strain_j2_umat/build.py",
        code: section(plasticityMaterialSource, "class SmallStrainJ2", "def generate"),
      },
    ],
    pythonMore: {
      title: "Read the element assembly and accepted-state commit",
      path: "tools/plasticity_bvp_common.py",
      code: section(
        plasticityHostSource,
        "    def evaluate_element",
        "def oracle_assembly",
      ),
    },
    boundaries: {
      title: "Boundary values in the simulation script",
      path: "tools/plasticity_bvp_common.py",
      code: section(plasticityHostSource, "def top_displacement", "def assert_trial_purity"),
    },
    generationTitle: "Generate the UMAT, then assemble mesh equilibrium",
    generationSummary:
      "generate_small_strain_umat writes the stress update below. DDSDDE comes from complex-step differentiation of that same update, one perturbed call per strain component. The check requires byte parity with the shipped source before compiling it. At each of four Gauss points per Quad4, the host calls the full 3D interface (NDI = 3, NSHR = 3, NTENS = 6) while imposing zero out-of-plane strain, and integrates Bᵀσ and BᵀDDSDDE B for its eight nodal DOFs.",
    fortranPath: "examples/small_strain_j2_umat/small_strain_j2.for",
    fortranExcerpt: section(
      plasticityFortran,
      "      cutback_flag = DCMPLX(0.0d0, 0.0d0)",
      "      DO ii = 1, 3",
      false,
    ),
    fortranCaption: "Elastic predictor and radial return",
    compilePath: "tools/plasticity_bvp_common.py",
    interface: plasticityFortran.split("\n").slice(0, 39).join("\n"),
    interfaceTitle: "Inspect the generated 3D UMAT interface",
    runtimeNote:
      "The UMAT supplies tension-positive Abaqus stresses and an engineering-shear DDSDDE. Every trial starts from copies of accepted stress, strain, and STATEV. Only a converged global equilibrium can commit all Gauss-point states. A rejected trial or cutback leaves accepted history untouched. This state ownership belongs to this verification host.",
  },
};
