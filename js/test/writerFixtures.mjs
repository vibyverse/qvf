// Shared by the writer tests: one archive that uses every typed builder.

import { Writer } from "../dist/index.js";

export const SOURCE = { program: "qvf-js-test", version: "1.0", calculation: "h2" };
export const FIXED = { timestamp: "2026-01-02T03:04:05Z", uuid: "0f1e2d3c-4b5a-4968-8776-655443322110" };
export const H2 = [
  { symbol: "H", position: [0, 0, 0] },
  { symbol: "H", position: [0, 0, 0.74] },
];
export const EYE = [
  [1, 0, 0],
  [0, 1, 0],
  [0, 0, 1],
];
// Two normalized s primitives: the overlap diagonal is 1 (spec Appendix A.1).
export const SHELLS = [
  { center: 0, l: 0, exponents: [1.0], coefficients: [1.0] },
  { center: 1, l: 0, exponents: [1.0], coefficients: [1.0] },
];
export const OVERLAP = [
  [1, 0.6],
  [0.6, 1],
];

/**
 * An archive with every builder: the kinds the Python writer's names cover.
 * @param {import("../dist/index.js").WriterOptions} [options]
 */
export async function fullArchive(options = {}) {
  const w = new Writer(SOURCE, { ...FIXED, ...options });
  w.setMetadata({ title: "H2 test archive", license: "CC0-1.0" });
  w.setProvenance({ method: "rhf", basis: "minimal" });
  await w.addStructure(H2, { charge: 0, multiplicity: 1, name: "dihydrogen", label: "H2" });
  await w.addStructure(H2, { sectionId: "structure_p", charge: 0, multiplicity: 1 });
  await w.addStructure(
    [{ symbol: "Na", position: [0, 0, 0] }, { atomic_number: 17, position: [1.4, 1.4, 1.4] }],
    { sectionId: "nacl", pbc: [true, true, true], latticeVectors: [[2.8, 0, 0], [0, 2.8, 0], [0, 0, 2.8]], dimensionality: 3 },
  );
  await w.addTrajectory({ symbols: ["H", "H"], kind: "md", pbc: [false, false, false], energies: [-1.1, -1.2, -1.15] }, [
    H2.map((a) => a.position),
    H2.map((a) => [a.position[0], a.position[1], a.position[2] * 1.01]),
    H2.map((a) => [a.position[0], a.position[1], a.position[2] * 0.99]),
  ], {
    times: new Float64Array([0, 0.5, 1.0]),
    energies: [-1.1, -1.2, -1.15],
    velocities: new Float64Array(18),
    precision: "float32",
  });
  await w.addWavefunctionGto(SHELLS, {
    moCoefficients: [
      [0.55, 0.55],
      [1.2, -1.2],
    ],
    energies: [-0.58, 0.67],
    occupations: [2, 0],
    nElectrons: 2,
    overlap: OVERLAP,
    densityMatrix: new Float64Array([0.6, 0.6, 0.6, 0.6]),
  });
  await w.addWavefunctionGto(SHELLS, {
    sectionId: "wf_u",
    moCoefficientsAlpha: { data: new Float64Array([0.55, 0.55, 1.2, -1.2]), shape: [2, 2] },
    moCoefficientsBeta: [[0.55, 0.55]],
    occupationsAlpha: [1, 0],
    occupationsBeta: [1],
    energiesAlpha: [-0.6, 0.6],
    energiesBeta: [-0.55],
    nElectrons: 2,
    overlap: OVERLAP,
    symmetryLabels: { alpha: ["sg", "su"], beta: ["sg"] },
  });
  const grid = { data: new Float32Array(8).fill(0.25), shape: [2, 2, 2] };
  await w.addDensity(grid, [0, 0, 0], EYE, { sectionId: "rho", structureRef: "structure", units: { origin: "bohr" } });
  await w.addOrbital([[[1, -1], [0.5, 0]], [[0, 0], [0, 0]]], [0, 0, 0], EYE, { component: "homo", dtype: "float64" });
  await w.addVolume("volume.vector", { data: new Float32Array(24), shape: [2, 2, 2, 3] }, [0, 0, 0], EYE, { sectionId: "grad_rho" });
  await w.addVolume("volume.complex", { data: new Float32Array(16), shape: [2, 2, 2, 2] }, [0, 0, 0], EYE, { sectionId: "psi_c" });
  await w.addVolume("volume.series", { data: new Float32Array(24), shape: [3, 2, 2, 2] }, [0, 0, 0], EYE, {
    sectionId: "rho_t",
    axis: { kind: "frame", values: [0, 1, 2], trajectory_ref: "trajectory" },
  });
  await w.addVolume("volume.difference", grid, [0, 0, 0], EYE, { sectionId: "drho", operandA: "rho", operandB: "rho" });
  await w.addVolume("volume.transition_density", grid, [0, 0, 0], EYE, { sectionId: "tdens", excitedStatesRef: "excited_states", stateIndex: 1 });
  await w.addExcitedStates(
    [
      {
        index: 1,
        excitation_energy: 12.1,
        multiplicity: 1,
        oscillator_strength: 0.5,
        transition_dipole: [0, 0, 1.1],
        contributions: [{ from_mo: 0, to_mo: 1, weight: 0.98 }],
        density_ref: "tdens",
      },
      { index: 2, excitation_energy: { value: 0.6, units: "hartree" }, multiplicity: 3 },
    ],
    { wavefunctionRef: "wf", meta: { method: "tddft", n_states: 2, tda: true } },
  );
  await w.addSpectrum("spectra.ir", { frequencies: [4400], intensities: [0.0], xUnits: "cm-1", peaks: [{ x: 4400 }] });
  await w.addSpectrum("spectra.uvvis", { payload: { frequencies: [12.1], intensities: [0.5], x_units: "ev", states_ref: "excited_states" } });
  await w.addSpectrum("spectra.nmr", { payload: { isotope: "1H", reference: "TMS", chemical_shifts: [{ atom_index: 0, isotropic_shift_ppm: 0.5 }] } });
  await w.addReactionScheme(
    [
      { id: "s1", reactants: ["structure"], products: [{ structure_ref: "structure_p", coefficient: 1 }], arrow: "equilibrium", conditions: ["300 K"], yield_percent: 50 },
    ],
    { label: "toy scheme" },
  );
  w.addSection({ id: "notes", kind: "x_test.notes", label: "free text" });
  await w.addOpaque("notes", "log", "line one\nline two\n", { mediaType: "text/plain" });
  await w.addJson("notes", "extra", { answer: 42, values: new Float64Array([1, 2]) });
  await w.addBinary("notes", "ids", [1n, -2n, 2n ** 62n], { dtype: "int64" });
  w.annotateSection("rho", { group: "densities" });
  return w;
}
