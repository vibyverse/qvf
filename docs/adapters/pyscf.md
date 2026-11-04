# PySCF adapter

`qvf.adapters.pyscf` writes molecular PySCF SCF/DFT and closed-shell
excited-state results, and periodic Γ-point SCF/DFT results, as QVF
archives (WP-62). It is the second independent producer of
`wavefunction.gto` beside vibe-qc, and it uses only existing QVF v1 kinds and
members: no format change. The adapter is unreleased and is not included in
toolkit v0.6.0. From a checkout containing it, install the optional dependency
and validator support:

```sh
pip install -e 'python/[pyscf,validate]'  # PySCF 2.14 or later, and NumPy
```

`import qvf.adapters.pyscf` works without PySCF; calling an export function
then raises an `ImportError` that names the extra. Pure convention helpers
such as `real_solid_harmonics` do not need the chemistry engine.

For a runnable first calculation and installed-package checks, follow the
[PySCF-to-QVF tutorial](../tutorials/pyscf-to-qvf.md). The adapter is included
in the prepared Python package candidate `0.7.0rc1`; that candidate is separate
from the existing toolkit v0.6.0 tag.

## SCF example

```python
from pyscf import dft, gto, scf
from qvf import Reader, validate_qvf
from qvf.adapters.pyscf import write_scf

mol = gto.M(atom="O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692",
            basis="def2-svp")
mf = dft.RKS(mol, xc="pbe0").run()
write_scf(mf, "water.qvf")
assert validate_qvf("water.qvf")["findings"] == []

oh = gto.M(atom="O 0 0 0; H 0 0 0.97", basis="def2-svp", spin=1, cart=True)
write_scf(scf.UHF(oh).run(), "oh-uhf-cartesian.qvf")

with Reader("water.qvf") as reader:
    wf = reader.get("wf")
    C = wf.coefficients().numpy()              # [n_mo, n_ao], QVF AO order
    S = reader.tensor("wf", "overlap").numpy()
    print(abs(C @ S @ C.T).max())              # orthonormal orbitals
```

To add the sections to a writer you already have, use
`add_scf(writer, mf)`, or the single steps `add_structure`, `add_wavefunction`,
`add_populations`, `add_job_spec` and `set_scf_provenance`; the
[Python API reference](../api/python.md) lists their options.

## What is written

| PySCF | QVF |
|---|---|
| `mf.mol` geometry, `charge`, `spin` | `structure` (positions in Å from the bohr coordinates of the calculation, `charge`, `multiplicity = spin + 1`) |
| `RHF`/`RKS`, `ROHF`/`ROKS`, `UHF`/`UKS` orbitals | `wavefunction.gto`, `spin` `restricted`, `restricted_open` or `unrestricted` |
| AO basis (`mol.cart` false or true) | `basis` payload, `pure` true or false, on `N_i`-normalized primitives |
| `mo_coeff`, `mo_energy`, `mo_occ` | coefficients and `mo_metadata` (energies in hartree, occupations; point-group labels with `symmetry=True`) |
| `mol.intor('int1e_ovlp')` | `overlap`, in the QVF AO basis |
| `mf.make_rdm1()` | `density_matrix`, or `density_matrix_alpha` and `density_matrix_beta` |
| ECP (`mol._ecpbas`) | `basis.ecp` entries with `n_core_electrons`, the ECP name and the scalar `projectors` ([below](#effective-core-potentials)). `ecp_projectors=False` writes the name and core count only |
| Mulliken and Löwdin analysis | `atom_properties`: `mulliken_charge`, `loewdin_charge`; open shells add `spin_population` (Mulliken α − β) |
| method, `xc`, basis name | `job.spec` (`functional` is `null` for Hartree–Fock), `provenance` |
| `e_tot`, `energy_nuc()` | `provenance.scf_energy`, root `energies` (`total`, `nuclear_repulsion`, `electronic`) |

Source program and version are PySCF's. A non-converged SCF is written with
`scf_converged: false` and `run_status: "failed"`.

The charges are computed from the archive's own density and overlap. For a
spherical basis the QVF AOs are PySCF's AOs up to order, so the values are
PySCF's own. For a Cartesian basis the Löwdin partition uses the AO
normalization of Appendix A.1 (axial components of unit norm).

## Excited states, spectra and NTOs

`write_td(td, path)` writes a completed molecular response calculation beside
its SCF reference. `add_td(writer, td)` adds the same sections to an open
writer. Both accept closed-shell RHF/RKS references and use existing QVF
`excited_states`, `spectra.uvvis` and `wavefunction.gto` kinds; there is no
new schema or archive version.

```python
from pyscf import dft, gto
from qvf import Reader
from qvf.adapters.pyscf import write_td

mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="6-31g", verbose=0)
mf = dft.RKS(mol, xc="pbe0").run()
tda = mf.TDA().run(nstates=2)
ids = write_td(tda, "h2-tda.qvf")

with Reader("h2-tda.qvf") as reader:
    assert reader.validate()["ok"]
    payload = reader.json(ids["excited_states"], "states")
    for state in payload["states"]:
        print(state["index"], state["excitation_energy"], state["oscillator_strength"])
    hole_meta = reader.json(ids["nto_holes"][0], "mo_metadata")
    assert hole_meta["occupation_semantics"] == "transition_weight"

full = mf.TDDFT().run(nstates=2)
full_ids = write_td(full, "h2-tddft.qvf")
assert len(full_ids["nto_holes"]) == len(full.e)
with Reader("h2-tddft.qvf") as reader:
    assert reader.validate()["ok"]
    weights = reader.json(full_ids["nto_holes"][0], "mo_metadata")["occupations"]
    print(sum(weights))  # full-response weights are not rescaled to sum to one
```

| Source result | State table and UV-visible sticks | Paired NTOs |
|---|---|---|
| Singlet RHF/RKS TDA | Yes; RHF TDA is labelled `cis`, RKS TDA `tddft` with `tda: true` | Written by default |
| Singlet TDHF/TDDFT, including pure-functional Casida response | Yes; `tda: false` | Written by default, retaining X and Y |
| Triplet response from a closed-shell RHF/RKS reference | Yes; multiplicity 3, with PySCF's zero electric transition strengths | Omitted |
| Unrestricted or restricted-open reference, spin-flip, complex or periodic response | Refused | Refused |
| Solvent response or direct-RPA variants | Refused | Refused |

The state table preserves the source root order and numbers it from 1.
Excitation energies are converted from hartree to eV with QVF's conversion
constant; state totals are the SCF total plus each vertical excitation.
Transition dipoles are in atomic units and oscillator strengths use the
length gauge. Every occupied-to-virtual contribution is included, with
normalized `X²−Y²` weight; negative terms are retained. A frozen-orbital
mask maps these contributions back to rows of the **full** reference
wavefunction, including when frozen virtual orbitals make that map
non-contiguous. Spectrum peaks refer to these same state indices.

Root `provenance`, `energies` and the optional `job.spec` describe the **SCF
reference**. The response method and TDA flag belong to the excited-state
metadata. The reference job is not a complete request to repeat the response
calculation; `job_spec=False` omits it. Export does not run SCF/response
iterations or alter the molecule, orbitals, energies or X/Y arrays. It
requires a converged reference and converged roots, real finite data and
PySCF's original normalization `Σ(X²−Y²) = 1/2`. Invalid source results
are rejected before `write_td` opens its target.

For singlet TDA, the adapter factorizes the normalized transition matrix
`X / ||X||` by SVD. Squared singular values are the pair weights; the hole
and particle sections have corresponding rows, identical descending
weights, `state_ref` / `state_index` and zero orbital energies. The AO map
and overlap are the reference's. A pair reconstructs the source transition
in the occupied/virtual space, including its relative signs. See
R. L. Martin, *J. Chem. Phys.* **118**, 4775–4777 (2003),
[doi:10.1063/1.1558471](https://doi.org/10.1063/1.1558471), eqs. 4–10,
and spec § 5.3.

For **singlet full response**, the matrix is instead
`T = sqrt(2) * (X + Y)`, in PySCF's restricted convention. Both amplitudes
contribute, and the squared singular-value weights retain their original
sum `2 * sum((X + Y)**2)`. That sum may be below or above one: do not
rescale it or interpret the weights as electron occupations or normalized
percentages. This uses QVF's existing allowance for non-unit full-response
weights (spec § 5.3), without a new field or validator rule.

This is the spin-adapted **spatial charge-transition** convention: it
represents the real, spin-independent Hermitian transition properties.
For a real symmetric one-electron operator `O`, the pair contraction is
`sqrt(2) * sum_k(sqrt(weight_k) * <hole_k|O|particle_k>)`. In particular,
it reconstructs PySCF's length-gauge transition dipole `2 * sum_ia((X+Y)_ia
* O_ia)`. It does not retain the separate X/Y blocks needed for
anti-Hermitian operators, such as the velocity-gauge momentum operator,
or define an excited-state difference density. The SVD of X+Y is described
by J. M. Herbert, *Density-functional theory for electronic excited states*,
in *Theoretical and Computational Photochemistry* (2023), eq. 3.51,
[doi:10.1016/B978-0-323-91738-4.00005-1](https://doi.org/10.1016/B978-0-323-91738-4.00005-1).
The factor `sqrt(2)` converts PySCF's half-normalized restricted amplitudes
to the adapter's spatial convention. The accepted normalized TDA path is
unchanged. NTO signs within each pair retain the source transition's phase.

The adapter does **not** call PySCF 2.14's `get_nto`: that helper normalizes
the stored X array in place and ignores Y in full TDHF/TDDFT. `ntos=False`
omits NTOs; `ntos=True` requires a supported singlet result and otherwise
raises `NotImplementedError`. Triplet NTOs remain deferred: the charge
transition cancels between spins, and a spin-transition representation
must be mapped explicitly. Unrestricted and restricted-open response
also remain unsupported.

The returned IDs include the `add_scf` entries plus `excited_states`,
`spectrum`, `nto_holes` and `nto_particles`; the last two are ordered lists
and are empty when omitted. `states_id`, `spectrum_id` and `nto_prefix`
customize response IDs, while `add_scf` options customize reference IDs and
optional sections. All references use those selected IDs.

These APIs are unreleased. The output uses existing v0.6.0 schema fields,
but consumer display support remains specific to each application; see the
[family workflow](../family_workflow.md). Frozen consumer pins stay v0.6.0.

## Periodic Γ-point calculations

`write_scf` and `add_scf` also accept a PySCF periodic mean-field object of a
`pbc.gto.Cell` at the Γ point: `pbc.scf.RHF`, `pbc.scf.UHF`, `pbc.dft.RKS`
and `pbc.dft.UKS`, with real orbitals. The archive uses the existing
`structure` (spec § 5.1) and `wavefunction.gto` (spec § 5.3) contracts, so
there is no format change.

```python
from pyscf.pbc import gto, scf
from qvf import validate_qvf
from qvf.adapters.pyscf import write_scf

cell = gto.M(atom="C 0 0 0; C 0.8917 0.8917 0.8917",
             a="0 1.7834 1.7834; 1.7834 0 1.7834; 1.7834 1.7834 0",
             basis="gth-szv", pseudo="gth-pade", mesh=[15, 15, 15])
mf = scf.RHF(cell).run()
write_scf(mf, "diamond-gamma.qvf")
assert validate_qvf("diamond-gamma.qvf")["findings"] == []
```

| PySCF | QVF |
|---|---|
| `cell.lattice_vectors()` (bohr, rows) | `structure.lattice_vectors`, rows in Å through `Lattice` and the format's bohr factor |
| `cell.dimension` | `pbc`: the first `dimension` axes are periodic, as PySCF treats them |
| atom positions | Cartesian positions in Å, as for a molecule, not wrapped into the cell |
| `mf.kpt` (zero) | `mo_metadata.k_point` `[0, 0, 0]` |
| orbitals, AO map | as for a molecule (Appendix A), one AO set per cell |
| `cell.intor('int1e_ovlp')` | `overlap`: the AO overlap of spec § 5.3, without images |
| `mf.make_rdm1()` | `density_matrix`, or the α and β densities |
| GTH pseudopotential | `basis.ecp` entries with `n_core_electrons` and the potential name as `label`, also for a zero core (H) |
| `mf.energy_nuc()` | `energies.components.nuclear_repulsion`: the Ewald energy of the point nuclei per cell |
| request | `job.spec` with `job_type: "periodic"`, `kpoints: [1, 1, 1]` and `options` `mesh`, `integrals` (`fftdf`, `gdf`, …), `pseudo` and `exxdiv` |

**The overlap and the Γ metric.** Spec § 5.3 defines `overlap` as
`⟨χ_μ|χ_ν⟩` of the contracted AOs, and validators check its diagonal against
Appendix A.1 (`E_NORMALIZATION`). The adapter therefore writes the AO overlap
of one cell without its images. The Γ-point orbitals are Bloch sums at k = 0
(spec § 5.3.1 with N = 1), orthonormal in the lattice-summed overlap
`S_Γ[μ, ν] = Σ_T ⟨χ_μ(r)|χ_ν(r − T)⟩`, which is PySCF's
`cell.pbc_intor('int1e_ovlp')`. Its diagonal exceeds one wherever an AO
overlaps its own images, 2.37 for the diffuse carbon s function of GTH-SZV
in diamond, so writing it as `overlap` would fail the Appendix A.1 check.
A consumer that needs `Cᵀ S_Γ C = 1`, the electron count `tr(P S_Γ)` or a
population builds `S_Γ` from the basis and the lattice. With the stored
`overlap`, `tr(P S)` is not the electron count of a cell. The Mulliken and
Löwdin charges in `atom_properties` are computed with `S_Γ`, so they equal
PySCF's periodic Mulliken analysis.

**Pseudopotentials.** A GTH pseudopotential has a local part with polynomial
prefactors and separable nonlocal projectors. It is not the semi-local form
of the spec's `projectors`, so the adapter records the core count and the
name only. Such an entry cannot be evaluated: rebuild the potential by name
in a program that has it. `n_electrons` is the valence count, the sum of the
atomic numbers minus the core counts minus the charge, and the populations use
the effective nuclear charges `Z − n_core_electrons`.

**Energies.** `energies.total` and `provenance.scf_energy` are `mf.e_tot`,
including the exchange-divergence (Madelung) correction that `exxdiv`
selects for HF and hybrid functionals. `nuclear_repulsion` is
PySCF's Ewald energy per cell of the point nuclei with their effective
charges, in a neutralizing background for a charged cell. `electronic` is
the rest.

**Refused.** k-point sampling objects (`KRHF`, `KUHF`, `KRKS` and the like:
their orbitals belong in `wavefunction.bloch`, a later tranche), a nonzero
`kpt`, complex orbitals, periodic ROHF, ROKS and GHF, a cell of dimension
zero, periodic response and grids of a cell (`add_density_grid`,
`add_orbital_grid`, `molecular_grid`), which would need lattice-summed AOs.
Space-group symmetry labels are not written.

## Effective core potentials

Every ECP atom gets a `basis.ecp` entry (spec
[§ 5.3](../../spec/qvf-format-spec.md#effective-core-potentials)) with its
core count, the ECP name as `label`, and by default the radial functions of
the potential as `projectors`, so a consumer can evaluate the potential from
the archive. `write_scf(mf, path, ecp_projectors=False)` (also accepted by
`add_scf`, `add_wavefunction`, `add_td` and `write_td`) writes the label and
core count only, as earlier versions of the adapter did.

The projectors are read from `mol._ecpbas` and `mol._env`, the input of the
`ECPscalar` integrals that built the core Hamiltonian
(`gto.mole.make_ecp_env`), not from the library file. They therefore hold
exactly the potential of the calculation, also for an ECP given inline. The
PySCF conventions map onto the spec without conversion:

| PySCF | QVF | Source |
|---|---|---|
| `ANG_OF = -1` | the local part, `local: true`, `l` one above the highest semi-local `l` of the centre (0 without channels) | `make_ecp_env`, spec § 5.3 |
| `RADI_POWER = n` (the `rorder` index of `mol._ecp`) | `powers` entry `n` of `c r^(n−2) exp(−α r²)` | NWChem ECP format (NWChem documentation, "Effective Core Potentials"), kept by `gto.basis.parse_nwchem.parse_ecp` |
| `ANG_OF = l ≥ 0`: `V = U_L + Σ_l \|l⟩ U_l ⟨l\|` | the semi-local channel, which already is `U_l − U_L` | Gaussian/NWChem/Molpro tabulation (the def2 tables contain `−U_L` in every channel) |
| exponents, coefficients | bohr⁻², hartree · bohr^(2−n) | atomic units of the basis-set libraries |
| `SO_TYPE_OF = 1` | not written, and the entry carries `x_pyscf_spin_orbit_omitted: true` | `ECPscalar` contains no spin-orbit term |

A potential without a local part (Stuttgart RSC for Ag, for example, whose
`U_L` is empty) is written without a `local: true` function, which the spec
reads as `U_L = 0`. The terms of one channel appear in PySCF's integral
order: by radial power, then by decreasing exponent.

Spin-orbit ECP terms (CRENBL, for example) enter only the spin-orbit
integrals of a two-component or GHF calculation, which the adapter refuses.
The scalar mean field the adapter writes never used them, so the written
projectors are still the complete potential of that calculation. The
`x_pyscf_spin_orbit_omitted` marker records that the named ECP has more.

## Densities and orbitals on a grid

`add_density_grid` and `add_orbital_grid` evaluate the SCF result on a
regular grid and add it as volume sections (spec § 5.2). Add the SCF
sections first, because the volumes refer to them:

```python
import pyscf
from pyscf import gto, scf
from qvf import Writer, validate_qvf
from qvf.adapters.pyscf import add_density_grid, add_orbital_grid, add_scf

mol = gto.M(atom="O 0 0 0; H 0 0 0.97", basis="def2-svp", spin=1)
mf = scf.UHF(mol).run()
with Writer("oh-grid.qvf", program="PySCF", version=pyscf.__version__) as writer:
    add_scf(writer, mf)
    add_density_grid(writer, mf)                     # density and spin density
    add_orbital_grid(writer, mf, ["homo", "lumo"])   # alpha set
    add_orbital_grid(writer, mf, ["homo"], spin="beta")
assert validate_qvf("oh-grid.qvf")["findings"] == []
```

| Function | QVF kind | Holds | Unit |
|---|---|---|---|
| `add_density_grid` | `volume.density` | ρ(r) of `mf.make_rdm1()`, both spins | e/bohr³ |
| `add_density_grid` (open shell) | `volume.spin` | ρ_α − ρ_β | e/bohr³ |
| `add_orbital_grid` | `volume.orbital`, one per orbital | ψ_i(r), PySCF's sign | bohr^(−3/2) |

`volume.spin` is written by default for UHF, UKS, ROHF and ROKS references
and never for RHF or RKS (`spin=True` then raises, `spin=False` omits it).
`add_orbital_grid` takes 0-based orbital rows, the rows of the coefficient
member of `wavefunction.gto`, or the names `"homo"` and `"lumo"`. An
unrestricted reference selects its set with `spin="alpha"` or `"beta"`.
Section ids are `mo_<row>`, or `mo_<spin>_<row>` for an unrestricted set,
and the labels name the row, the spin and HOMO or LUMO.

**Grid.** By default `molecular_grid(mol)` builds an axis-aligned box that
reaches `DEFAULT_GRID_MARGIN_BOHR = 4.0` bohr beyond the outermost nuclei on
every side, sampled with `DEFAULT_GRID_SPACING_BOHR = 0.2` bohr steps. The
box edge is rounded up to whole steps and centred on the nuclei. Pass
`spacing=` and `margin=` to change it, or `shape=[nx, ny, nz]` to sample the
box exactly as PySCF's `tools.cubegen.Cube(mol, nx, ny, nz, margin=...)`
does: origin at the lower corner and `shape[a] − 1` steps across each edge.
Any other grid, including non-orthogonal voxel vectors, is given as
`grid={"origin": ..., "voxel_vectors": ..., "shape": ...}` in bohr, the
fields of the grid payload. Point `(i, j, k)` is
`origin + i·v_0 + j·v_1 + k·v_2`, stored in row-major `[i][j][k]` order. The
grids are finite boxes, so `periodic` and `includes_endpoint` are not written.

**Evaluation.** Values come from `pyscf.dft.numint.eval_ao` and `eval_rho`
(orbitals as `ao @ C`) in blocks of `chunk` points. The default block holds
2²³ AO values (64 MiB), so memory does not grow with the grid. All orbitals
of one call share one AO evaluation. Data are stored as float32, the format
default, or as float64 with `dtype="float64"`.

**Links.** Every volume carries `structure_ref` at section and grid level.
The wavefunction id (`wavefunction_ref`, default `"wf"`, which must already
be in the writer) is written as the section field `wavefunction_ref`
(§ 5.2, *Volumes of a wavefunction*), and each `volume.orbital` also names
its coefficient row as `mo_index` and, for an unrestricted reference, its
channel as `spin` (`"alpha"` or `"beta"`). A viewer reads the orbital's
energy and occupation from that row of the wavefunction's `mo_metadata`.
The same id is also written as the section `group` (§ 4.1), so viewers that
group by it keep working. Pass `wavefunction_ref=None` to write neither.
The section fields are part of a proposed additive format change (see
`GOVERNANCE.md`, 2026-11-04), and toolkits that predate it ignore them.

**Integration accuracy.** The grid sum `Σ f(r) ΔV` is a check, not an
accurate quadrature, for an all-electron density. A uniform grid aliases the
tight core primitives (α up to 2266 bohr⁻² for O in def2-SVP), and the error
changes sign with the grid's offset from the nuclei. Measured with a
4 bohr margin:

| System | Spacing / bohr | ∫ρ − N / e | ∫(ρ_α − ρ_β) − 2S | ∫ψ²_HOMO − 1 |
|---|---|---|---|---|
| H₂O RHF/def2-SVP | 0.10 | −0.015 | | −1.0e-4 |
| H₂O RHF/def2-SVP | 0.15 | −0.015 | | −7.5e-5 |
| H₂O RHF/def2-SVP | 0.20 | −0.061 | | −1.3e-4 |
| H₂O RHF/def2-SVP | 0.30 | −0.71 | | +1.5e-3 |
| OH UHF/def2-SVP | 0.15 | +0.012 | −1.4e-5 | −5.9e-5 |
| OH UHF/def2-SVP | 0.20 | −0.22 | −3.9e-4 | +1.6e-6 |

The density error is about twice that of the O 1s orbital. Valence
orbitals and the spin density are smooth and converge to about 1e-4 at
0.2 bohr. Diffuse virtual orbitals are limited by the box instead: the
def2-SVP LUMO of water loses 2 % outside a 4 bohr margin and 12 % outside
3 bohr. Without a tight core the sum converges in the spacing and only the
box limits it: H₂ RHF and H₃ UHF in def2-SVP integrate to within 5e-5 e
with 0.25 bohr steps and a 6 bohr margin, against 2e-3 to 9e-3 e with a
4 bohr margin. Use 0.1 to 0.15 bohr for quantitative integration of a
density with heavy-atom cores, and a margin of 6 bohr or more for diffuse
orbitals.

## How the AO conventions map

These are the three places where a PySCF archive could silently disagree with
another producer (spec [Appendix A](../../spec/qvf-format-spec.md#appendix-a-ao-ordering-primitive-normalization-and-component-naming)).

**Primitive normalization (A.1).** libcint stores each contraction
coefficient multiplied by the radial norm `gto_norm(l, α)` and the contraction
norm, and gives its Cartesian s and p functions the common factors `1/(2√π)`
and `√(3/(4π))`. The adapter takes the stored coefficient
(`mol.bas_ctr_coeff(ib) * gto.gto_norm(l, α)`), applies that factor, divides
by the Appendix A.1 `N_i`, and scales the contraction to a unit axial
self-overlap. Because `gto_norm(l, α) = N_i · √(4π/(2l+1))` for every
exponent, the result is PySCF's own contraction on `N_i`-normalized
primitives. A PySCF Cartesian AO with `l ≥ 2` has the norm `√(4π/(2l+1))`
(its `d_xx` has self-overlap 4π/5); in the archive it is divided by that, so
`⟨xx|xx⟩ = 1` and `⟨xy|xy⟩ = 1/3`, exactly as Appendix A.1 requires.

**In-shell order and phases (A.2).** PySCF's spherical `p` shell is
`x, y, z` and is permuted to `m = −1, 0, +1` (`y, z, x`). For `l ≥ 2` the
spec's harmonics are the conventional real solid harmonics: `cos(mφ)` for
`m > 0`, `sin(|m|φ)` for `m < 0`, no Condon–Shortley phase, each normalized
like the axial `x^l` (Helgaker, Jørgensen and Olsen, eqs. 6.4.47–6.4.50;
the convention of libint). The adapter does not assume that PySCF uses
them: it solves for the map from PySCF's `gto.cart2sph(l)` table onto
`real_solid_harmonics(l)`, refuses a basis whose map is not exact, and
writes the coefficients through that map. With PySCF 2.14 the table is
`√((2l+1)/(4π))` times the same harmonics in the same `m = −l … +l` order,
so the map is the identity for every `l ≥ 2`. Cartesian shells use libint's
lexicographic order in both programs.

**General contractions.** A PySCF shell with several contracted functions
(cc-pVDZ's s shells, for example) becomes one QVF shell per function, in
PySCF's AO order; primitives with a zero coefficient in that function are
left out.

All matrices are transformed with one block-diagonal AO map `T`
(`qvf.adapters.pyscf.ao_map(mol)`): `C_qvf = T C`, `P_qvf = T P Tᵀ` and
`S_qvf = T⁻ᵀ S T⁻¹`. The orbitals are therefore PySCF's orbitals, sign
included, and `C S Cᵀ = 1` holds in the archive.

## Checks

`python/tests/test_adapter_pyscf.py` (skipped without PySCF) writes water
RHF and RKS (PBE0, B3LYP) in def2-SVP and cc-pVDZ, spherical and Cartesian,
and the OH radical in UHF, UKS and ROHF. Every archive validates with no
finding, `E_NORMALIZATION` included; `C S Cᵀ = 1`; the populations add up to
the electron count; the stored densities equal PySCF's. An AO evaluator
written from Appendix A alone, with SciPy's spherical harmonics, reproduces
every PySCF orbital on random points, which also catches a sign or order
error inside a spherical shell that the overlap diagonal cannot see.
Candidate CI installs PySCF 2.14.0 on GitHub's
Python 3.14 job so these tests run; the Python 3.9 job keeps PySCF optional.

`scripts/pyscf_cross_producer.py REFERENCE.qvf` rebuilds the molecule of any
`wavefunction.gto` archive in PySCF, writes it with the adapter and compares
the SCF energy and the AO density matrix in the QVF AO convention. Against
vibe-qc archives (libint) of water RHF/def2-SVP, formaldehyde
RHF/def2-TZVP (with f shells) and the OH radical UHF/6-31G\*, the energies
agree to 3e-12 hartree and the densities to `max |ΔP|` ≈ 1e-8, the SCF
convergence level. The OH comparison starts PySCF from the reference density,
because the UHF π hole has two solutions of the same energy.

`python/tests/test_adapter_pyscf_ecp.py` checks the ECP projectors against
PySCF in three independent ways. The archive's radial functions equal those
evaluated from `gto.basis.load_ecp`, channel by channel, for def2
(Stuttgart-Köln ECP28MWB on I), LANL2DZ (with `r⁻²` and `r⁻¹` terms),
Stuttgart large-core I, Stuttgart RSC Ag and CRENBL (with spin-orbit terms).
Adaptive radial quadrature of the archive's potential reproduces the
one-centre `ECPscalar` integrals of s, p, d and f primitives to 1e-10
relative, which fixes the `r^(n−2)` power and shows that a channel acts in
addition to the local part. Finally a molecule rebuilt in PySCF from the
archive's `basis` payload alone (shells and ECP) reproduces the
`ECPscalar` matrix to 5e-14 and the RHF energy to 2e-13 hartree for HI in
def2-SVP, I₂ in LANL2DZ, AgH in def2-SVP with the Stuttgart RSC ECP and HI
in CRENBL. Every archive validates with no finding. The response sections'
NTO bases carry the same entries, and `ecp_projectors=False` keeps the
label-only entries.

`python/tests/test_adapter_pyscf_td.py` checks response conversion separately:
RHF/RKS TDA and full response, the Casida class, triplet state tables,
spherical/Cartesian bases, rectangular NTO spaces and frozen MO indices.
It reads every archive back, checks validation, energies, transition
properties and references, reconstructs the TDA or full-response transition
and length-gauge moments from the stored NTOs, and verifies that source
arrays are unchanged. Explicit-Y controls catch X-only, X−Y and
unit-normalized mappings; cases retain weights both above and below a unit
sum. Error cases include
unconverged roots and X amplitudes changed by PySCF's NTO helper. These are
conversion checks, not a new scientific accuracy claim for either engine.

`python/tests/test_adapter_pyscf_grid.py` checks the grid output for water
RHF (spherical) and RKS (Cartesian) and the OH radical in UHF and ROHF on
0.15 bohr grids. Every archive validates with no finding. The density, the
spin density and the HOMO, LUMO, lowest and β HOMO orbitals integrate within
the bounds above. On sampled points, including the densest ones near the
nuclei, the stored values equal ρ, ρ_α − ρ_β and ψ evaluated from the
archive's own `basis`, density matrices and coefficients with the Appendix A
evaluator of `test_adapter_pyscf.py`, to float32 precision. With `shape` and
`margin` the points equal `tools.cubegen.Cube.get_coords()` to 1e-12 bohr,
the values equal `cubegen.density` and `cubegen.orbital`, and the written
cube file, imported through `qvf.convert.import_cube`, has the archive's
grid. A skewed explicit grid and a small `chunk` give the same values.

`python/tests/test_adapter_pyscf_pbc.py` checks the Γ-point path on diamond
(two-atom primitive cell, RHF and RKS/PBE) and a hydrogen atom in a 4 Å cubic
box (UHF and UKS), all with GTH-SZV, GTH-Pade and a 15³ mesh, in about
30 s. Every archive validates with no finding. The lattice vectors and the
fractional coordinates round-trip to 1e-12 bohr. `S_Γ` is computed twice:
by quadrature over the cell of the lattice-summed AOs from the Appendix A
evaluator, and by PySCF from a cell rebuilt from the archive alone (the
structure, the basis shells, the pseudopotential labels and the mesh). The two
agree to 2e-11 for diamond and 1e-7 for the hydrogen box, and with either the
stored orbitals satisfy `C S_Γ Cᵀ = 1` and `tr(P S_Γ) = N`. The stored
`overlap` equals the AO overlap without images, with unit diagonal. An Ewald
sum written in the test reproduces `nuclear_repulsion` to 1e-13 hartree,
and `total = electronic + nuclear_repulsion`. The rebuilt cell reproduces
`e_tot` from the stored density to 1e-8 hartree and by a fresh SCF to 1e-7
hartree (observed: 2e-15). Charges equal PySCF's Mulliken analysis with
`S_Γ`.

## Not yet covered

Later WP-62 tranches: unrestricted/open-shell response, triplet NTOs and
k-point sampled periodic calculations (`wavefunction.bloch`). GTH projectors
(named only, see above), spin-orbit ECP terms (omitted, see above), solvent
models (refused), GHF and relativistic two-component wavefunctions (refused)
and ghost atoms (refused) are not mapped either.
