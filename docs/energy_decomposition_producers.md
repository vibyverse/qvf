# Energy decompositions from common programs

This page is informative. It shows how the output of programs that compute
energy decompositions maps onto the `energy_decomposition` kind (spec § 5.16).
Each mapping was checked by converting a real output and validating the
archive, including every declared sum rule. The program versions, systems and
methods are listed with each table; numbers are rounded for display.

The rule of thumb for every scheme is the same:

- Write the terms as printed, with the program's own name in `label`.
- Nest sub-terms under the term they are part of (`parent`).
- Declare every additivity the program guarantees as a sum rule, at the
  precision it prints.
- Use an `x_` name for a term that has no general meaning outside its scheme.

## AMS ADF: EDA and ETS-NOCV

ADF 2026.1, PBE-D3(BJ)/TZ2P all-electron; water dimer, water trimer (three
fragments), and ethane from two methyl radicals with one spin flipped.

| ADF output | Term | Notes |
|---|---|---|
| Total Bonding Energy | `interaction` | |
| Total Pauli Repulsion | `pauli` | parts: `kinetic` (ΔT⁰), `coulomb` (ΔV^Pauli Coulomb), `x_lda_xc`, `x_gga_exchange`, `x_gga_correlation` |
| Electrostatic Interaction | `electrostatic` | |
| Total Steric Interaction | `frozen` | sum rule: `pauli` + `electrostatic` |
| Total Orbital Interactions | `orbital` | parts per irreducible representation (`x_irrep_A`, …); a second sum rule over the alternative split `kinetic`, `coulomb`, `exchange_correlation` (with `order: "alternative"`) |
| Residu | `residual` | needed for the total sum rule |
| Dispersion Energy | `dispersion` | |
| ETS-NOCV eigenvalue pairs, ΔE_k (kcal/mol) | `orbital_pairs`, `nocv`, `per_spin` | unrestricted fragments: `alpha` and `beta` pairs |
| `fragoccupations` with β > α | `spin_flip: true`, `multiplicity: 2` | |

The hartree column has 15 decimals, so the default tolerance (10⁻⁶ hartree)
applies. Water dimer: ΔE_int −4.59, Pauli 4.74, electrostatic −6.02, orbital
−2.93, dispersion −0.38 kcal/mol; 7 NOCV pairs.

## AMS BAND: periodic EDA and pEDA-NOCV

BAND 2026.1, PBE/TZP; H₂O on h-BN, CO on MgO(001), H₂O on H-terminated Si(001)
at Γ and with KSpace Quality Normal.

| BAND output | Term | Notes |
|---|---|---|
| E_int | `interaction` | |
| E_Pauli, E_elstat, E_orb | `pauli`, `electrostatic`, `orbital` | |
| E^0 = T^0 + Elst^0 + XC^0 | `frozen` with parts `kinetic`, `electrostatic`, `exchange_correlation` | |
| T_orb + Elst_orb + XC_orb | parts of `orbital` | |
| PEDANOCV table (kJ/mol) | `orbital_pairs`, `nocv`, `per_spin` | T and V per pair as `x_` keys; BAND computes NOCVs at Γ only (`orbital_pairs_gamma_only`) |

BAND prints the energies to 5 decimals in hartree. The sum rules therefore carry
`tolerance: {value: 3e-5, units: hartree}`.

## ORCA: EDA-NOCV

ORCA 6.1.1, PBE-D3BJ/def2-TZVP; water dimer and ethane (spin flip with `FRAG1_SF`).

| ORCA output | Term | Notes |
|---|---|---|
| Bond Energy | `interaction` | |
| Electrostatic Energy | `electrostatic` | |
| Pauli Energy + Delta E^0(XC) | `pauli` | ORCA prints the Pauli term without the exchange–correlation change of the frozen state; § 5.16's `pauli` is the sum, written with the two printed values as parts (`x_pauli_without_xc`, `exchange_correlation`) |
| Orbital Energy | `orbital` | |
| Delta Dispersion | `dispersion` | |
| NOCV analysis (DE_k, kcal/mol, 2 decimals) | `orbital_pairs`, `nocv` | DT_k and DV_k as `x_` keys; open shell: `alpha` and `beta` tables |

In the open-shell spin-flip case, ORCA 6.1.1's printed terms miss its own
Bond Energy by 5.9·10⁻⁵ hartree (0.04 kcal/mol). The sum rule then needs a
tolerance, and the payload says why in `notes`. This is the case the
declared sum rules exist for.

## ORCA: local energy decomposition (DLPNO-CCSD(T))

ORCA 6.1.1, DLPNO-CCSD(T)/def2-TZVP, TightPNO; water dimer and cyclic water trimer.

| ORCA output | Term | Notes |
|---|---|---|
| Sum of INTER-fragment total energies | `interaction` | without monomer calculations this excludes the electronic preparation of the fragments |
| Interaction of Fragments i and j: Total | `interaction`, `fragments: [i, j]` | the pairs add up to the interaction |
| Interfragment reference / correlation | `reference`, `correlation` (pair, parent = pair total) | |
| Electrostatics (REF.), Exchange (REF.) | `electrostatic`, `exchange` (parent = reference) | |
| Inter strong pairs / triples / weak pairs | `x_strong_pairs`, `x_triples`, `x_weak_pairs` (parent = correlation) | |
| Dispersion (strong pairs), Charge Transfer i to j | `dispersion`, `charge_transfer` with `directed: true` (parent = strong pairs) | not additive to the strong pairs, so in no sum rule |
| Intra fragment REF. energy | `x_intra_reference`, `fragments: [i]` | absolute energies, not interaction terms |

## Psi4: SAPT0 and F-SAPT

Psi4 1.9.1, jun-cc-pVDZ; water dimer. One archive holds two decompositions.

| Psi4 output | Term | Notes |
|---|---|---|
| Total SAPT0 | `interaction` | |
| Electrostatics, Exchange, Induction, Dispersion | `electrostatic`, `exchange`, `induction`, `dispersion` | |
| Elst10,r; Exch10 | `x_elst10r`, `x_exch10` (with `order`) | |
| Exch10(S^2) | `x_exch10_s2` under `exchange`, in no sum rule | an alternative estimate, not a part |
| Ind20,r; Exch-Ind20,r; delta HF,r (2) | `x_ind20r`, `exchange_induction`, `delta_hf` | sum rule to `induction` |
| Disp20; Exch-Disp20, each with (SS) and (OS) | `x_disp20`, `exchange_dispersion`; the same-spin and opposite-spin parts nested below them | |

F-SAPT (`fisapt0`, links by charge) partitions each monomer into atoms and
"link" fragments, the electrons of the bonds that the partition cuts. Atoms and
links are sub-fragments (`parent`). A link is written as the atom pair of its
bond with an `x_` note, because a fragment is a set of atoms. Every non-zero
pair of the full analysis is an `interaction` with `fragments: [a, b]`, and its
six columns (Elst, Exch, IndAB, IndBA, Disp, EDisp) are parts of it. The
values are printed in kcal/mol to three decimals, so the pair sum rules carry a
tolerance of 0.005 kcal/mol, and the sum over pairs one of 0.02 kcal/mol.

## CP2K: ALMO-EDA

CP2K 2023.1, PBE-D3(BJ)/DZVP-MOLOPT-GTH, `ALMO_SCF` with `ANALYSIS` and XALMO
delocalization; water dimer and water trimer in a 15 Å box. The fragments are
the molecules: set `&TOPOLOGY CONN_FILE_FORMAT GENERATE`, because the default
makes every atom a molecule.

| CP2K output | Term | Notes |
|---|---|---|
| total − Σ isolated molecules | `interaction` | a difference of printed energies |
| block-diagonal ALMOs − Σ isolated | `x_frozen_polarization` | CP2K does not split frozen and polarization in this analysis |
| ENERGY LOWERING | `charge_transfer` | |
| `ct_energy_terms` file, i ≠ j | `charge_transfer`, `directed: true`, `fragments: [i, j]` | the transferred charge from `charge_terms` as an `x_` key; the diagonal (intra-fragment) entries as `x_` terms, so the parts sum to the total |

## vibe-qc: periodic EDA, LMO-EDA, Kitaura–Morokuma

vibe-qc writes `energy_decomposition` directly. Its periodic EDA gives NOCV
eigenvalues and energies per k-point (`kpoint_eigenvalues`,
`kpoint_energies`) and the NOCV deformation densities as `volume.difference`
sections. Its molecular LMO-EDA and Kitaura–Morokuma analyses of one set of
calculations are two `energy_decomposition` sections of one archive.
