"""PySCF adapter: molecular SCF and DFT results as QVF archives (WP-62).

Writes a converged (or settled) PySCF mean-field object through the QVF
writer API: the molecule as a ``structure`` section, the orbitals as a
``wavefunction.gto`` section with its ``overlap`` and density members,
Mulliken and Löwdin populations as an ``atom_properties`` section, the
request as a ``job.spec`` section, and the SCF energy in the root
``provenance`` and ``energies`` blocks. Supported are restricted,
unrestricted and restricted open-shell Hartree–Fock and Kohn–Sham objects of
a molecular ``gto.Mole`` (``RHF``, ``UHF``, ``ROHF``, ``RKS``, ``UKS``,
``ROKS``, including their density-fitted and second-order variants), in a
spherical or a Cartesian basis.

``add_td`` and ``write_td`` add molecular closed-shell TDA/TDHF/TDDFT
states and UV-visible sticks beside that SCF reference. Singlet TDA and
full response can also carry paired natural transition orbitals. Full
response retains both X and Y; triplet NTOs are not mapped (see ``add_td``).

PySCF is imported lazily: ``import qvf.adapters.pyscf`` works without it, and
calling an export function then raises ``ImportError``. From a checkout
containing this unreleased adapter, install the optional dependencies with
``python -m pip install 'python/[pyscf,validate]'``. The format is unchanged:
everything is written with existing QVF v1 kinds and members.

The AO conventions (spec Appendix A) are mapped as follows; the
[PySCF adapter guide](../adapters/pyscf.md) has the derivation.

* **Primitive normalization (A.1).** libcint stores each contraction
  coefficient multiplied by the radial norm ``gto_norm(l, α)`` and by the
  contraction norm. The adapter takes those raw coefficients
  (``mol.bas_ctr_coeff(ib) * gto.gto_norm(l, α)``), applies libcint's common
  s and p factor (``1/(2√π)`` and ``√(3/(4π))``, 1 for ``l ≥ 2``), divides by
  the Appendix A.1 ``N_i`` and scales the contraction so that its axial
  self-overlap is 1. The written ``overlap`` therefore has a diagonal of 1
  for every spherical AO and every axial Cartesian AO, and 1/3 for ``d_xy``.
* **In-shell order (A.2).** PySCF's spherical ``p`` shell is ``x, y, z``; it
  is permuted to ``m = −1, 0, +1`` (``y, z, x``). For ``l ≥ 2`` the order
  ``m = −l … +l`` and the real-solid-harmonic phases are not assumed: the
  map is solved from PySCF's ``gto.cart2sph`` table against the spec's real
  solid harmonics (``real_solid_harmonics``) and checked to be exact.
  With PySCF 2.14 it is the identity. Cartesian shells share libint's
  lexicographic order and need no permutation.
* **Effective core potentials** are written as ``basis.ecp`` entries with
  their scalar ``projectors`` (spec § 5.3), read from ``mol._ecpbas``. PySCF's
  radial power and difference channels are the spec's; spin-orbit terms are
  omitted (``_ecp_projectors`` has the mapping).
* **General contractions** (one PySCF shell with several contracted
  functions, as in cc-pVDZ) become one QVF shell per contracted function, in
  PySCF's AO order. Primitives whose coefficient is exactly zero in that
  function are left out.

Coefficients, densities and the overlap are transformed with the same AO map
(``AOMap``), so the orbitals are exactly PySCF's orbitals and
``Cᵀ S C = 1`` holds in the archive.
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Any

from .._legacy_writer import BOHR_PER_ANGSTROM, primitive_norm

__all__ = [
    "AOMap", "ao_map", "add_job_spec", "add_populations", "add_scf", "add_td",
    "add_structure", "add_wavefunction", "basis_name", "cartesian_powers",
    "method_name", "real_solid_harmonics", "set_scf_provenance", "write_scf", "write_td",
]

_INSTALL_HINT = ("qvf.adapters.pyscf needs the optional qvf-writer[pyscf] dependencies; "
                 "from a toolkit checkout run python -m pip install 'python/[pyscf,validate]'")


def _pyscf():
    """Import PySCF (and NumPy) on first use, with an actionable error."""
    try:
        import pyscf
        from pyscf import dft, gto, scf  # noqa: F401  (submodules used below)
    except ImportError as exc:  # pragma: no cover - exercised without PySCF
        raise ImportError(_INSTALL_HINT) from exc
    return pyscf


def _np():
    try:
        import numpy
    except ImportError as exc:  # pragma: no cover - PySCF depends on NumPy
        raise ImportError(_INSTALL_HINT) from exc
    return numpy


# -- conventions --------------------------------------------------------------

def cartesian_powers(l: int) -> list[tuple[int, int, int]]:
    """The ``(l_x, l_y, l_z)`` of a Cartesian shell in spec Appendix A.2 order.

    Decreasing ``l_x`` first, then decreasing ``l_y`` (libint's order, which
    is also libcint's): ``xx, xy, xz, yy, yz, zz`` for a d shell.
    """
    return [(lx, ly, l - lx - ly) for lx in range(l, -1, -1)
            for ly in range(l - lx, -1, -1)]


def real_solid_harmonics(l: int) -> list[list[float]]:
    """The real solid harmonics of a QVF spherical shell on its Cartesians.

    Returns a ``(l+1)(l+2)/2`` by ``2l+1`` table: row ``k`` is the Cartesian
    component ``cartesian_powers(l)[k]``, column ``j`` the harmonic
    ``m = j − l`` (spec Appendix A.2 order). The harmonics are the
    conventional real solid harmonics ``S_lm`` (T. Helgaker, P. Jørgensen,
    J. Olsen, *Molecular Electronic-Structure Theory*, Wiley 2000,
    eqs. 6.4.47–6.4.50; H. B. Schlegel, M. J. Frisch, *Int. J. Quantum Chem.*
    **54**, 83 (1995)): ``cos(mφ)`` for ``m > 0``, ``sin(|m|φ)`` for
    ``m < 0``, without the Condon–Shortley phase, so ``S_1,−1 = y``,
    ``S_10 = z``, ``S_11 = x`` and ``S_2,−2 = √3 xy``. Each has the norm of
    the axial component ``x^l``, so on the ``N_i``-normalized primitives of
    Appendix A.1 every spherical AO has the self-overlap ``S_axial``.
    """
    if not isinstance(l, int) or l < 0:
        raise ValueError("l must be an integer >= 0")
    index = {p: k for k, p in enumerate(cartesian_powers(l))}
    table = [[0.0] * (2 * l + 1) for _ in index]
    for m in range(-l, l + 1):
        am = abs(m)
        wm = 0 if m >= 0 else 1  # 2 v_m
        norm = math.sqrt(2 * math.factorial(l + am) * math.factorial(l - am)
                         / (2 if m == 0 else 1)) / (2 ** am * math.factorial(l))
        for t in range((l - am) // 2 + 1):
            for u in range(t + 1):
                for w in range(wm, am + 1, 2):  # w = 2 v
                    c = (Fraction(-1) ** (t + (w - wm) // 2) * Fraction(1, 4 ** t)
                         * math.comb(l, t) * math.comb(l - t, am + t)
                         * math.comb(t, u) * math.comb(am, w))
                    powers = (2 * t + am - 2 * u - w, 2 * u + w, l - 2 * t - am)
                    table[index[powers]][m + l] += norm * float(c)
    return table


def _axial_self_overlap(l, exponents, coefficients):
    """``S_axial`` of spec Appendix A.1 for one contraction (ratio form)."""
    total = 0.0
    for a, ca in zip(exponents, coefficients):
        for b, cb in zip(exponents, coefficients):
            ratio = min(a, b) / max(a, b)
            total += ca * cb * (2.0 * math.sqrt(ratio) / (1.0 + ratio)) ** (l + 1.5)
    return total


def _libcint_sp_factor(l):
    """libcint's common factor of Cartesian s and p functions (1 for l >= 2)."""
    return math.sqrt((2 * l + 1) / (4.0 * math.pi)) if l < 2 else 1.0


def _snap(block):
    """Round a block that is a signed permutation to within 1e-9 to exactly one."""
    np = _np()
    rounded = np.rint(block)
    if (np.abs(block - rounded).max() <= 1e-9
            and set(np.unique(rounded).tolist()) <= {-1.0, 0.0, 1.0}
            and (np.count_nonzero(rounded, axis=0) == 1).all()
            and (np.count_nonzero(rounded, axis=1) == 1).all()):
        return rounded
    return block


# -- the AO map ---------------------------------------------------------------

class AOMap:
    """The map from a PySCF AO basis to the QVF ``basis`` payload it writes.

    ``shells`` is the QVF shell list (spec § 5.3), ``pure`` the payload's
    angular type and ``n_ao`` its AO count. The map is block diagonal: a
    PySCF AO ``ν`` is ``Σ_μ T[μ, ν] χ_μ`` in the QVF AOs ``χ_μ``, one block
    per contracted shell. Its methods transform PySCF matrices into the QVF
    AO basis: MO coefficients ``C_qvf = T C`` (returned with one row per MO,
    as the ``mo_coefficients`` member stores them), densities
    ``P_qvf = T P Tᵀ`` and the overlap ``S_qvf = T⁻ᵀ S T⁻¹``.
    """

    def __init__(self, shells, blocks, n_ao, pure):
        self.shells = shells
        self.pure = pure
        self.n_ao = n_ao
        self._blocks = blocks  # (start, T, T^-1)

    @property
    def matrix(self):
        """The dense ``[n_ao, n_ao]`` AO map ``T``."""
        np = _np()
        out = np.zeros((self.n_ao, self.n_ao))
        for start, block, _ in self._blocks:
            n = block.shape[0]
            out[start:start + n, start:start + n] = block
        return out

    def _left(self, matrix, inverse_transpose=False):
        np = _np()
        matrix = np.asarray(matrix, dtype=float)
        if matrix.shape[0] != self.n_ao:
            raise ValueError(f"expected {self.n_ao} AO rows, got {matrix.shape[0]}")
        out = np.empty_like(matrix)
        for start, block, inverse in self._blocks:
            n = block.shape[0]
            op = inverse.T if inverse_transpose else block
            out[start:start + n] = op @ matrix[start:start + n]
        return out

    def coefficients(self, mo_coeff):
        """PySCF ``mo_coeff`` ``[n_ao, n_mo]`` as QVF rows ``[n_mo, n_ao]``."""
        return _np().ascontiguousarray(self._left(mo_coeff).T)

    def density(self, dm):
        """A PySCF AO density matrix in the QVF AO basis, ``T P Tᵀ``."""
        return _np().ascontiguousarray(self._left(self._left(dm).T).T)

    def overlap(self, s):
        """A PySCF AO overlap matrix in the QVF AO basis, ``T⁻ᵀ S T⁻¹``."""
        once = self._left(s, inverse_transpose=True)
        return _np().ascontiguousarray(self._left(once.T, inverse_transpose=True).T)

    def ao_atoms(self):
        """The atom index of every QVF AO, in AO order."""
        out = []
        for shell in self.shells:
            l = shell["l"]
            out.extend([shell["center"]] * ((2 * l + 1) if shell["pure"]
                                            else (l + 1) * (l + 2) // 2))
        return out


def ao_map(mol) -> AOMap:
    """Build the ``AOMap`` of a PySCF ``gto.Mole`` (spec Appendix A).

    Raises ``ValueError`` for a basis this tranche cannot express, and
    ``RuntimeError`` if PySCF's spherical transform is not a linear map onto
    the spec's real solid harmonics (an unknown convention).
    """
    _pyscf()
    from pyscf import gto
    np = _np()
    _check_molecule(mol)
    pure = not mol.cart
    ao_loc = mol.ao_loc_nr()
    harmonics = {}
    shells, blocks = [], []
    for ib in range(mol.nbas):
        l = int(mol.bas_angular(ib))
        if mol.bas_kappa(ib) != 0:
            raise ValueError("spinor (kappa != 0) shells are not supported")
        if l > 7:
            raise ValueError(f"QVF shells have l <= 7, PySCF shell {ib} has l = {l}")
        center = int(mol.bas_atom(ib))
        exps = np.asarray(mol.bas_exp(ib), dtype=float)
        raw = np.asarray(mol.bas_ctr_coeff(ib), dtype=float) * gto.gto_norm(l, exps)[:, None]
        n_cart = (l + 1) * (l + 2) // 2
        n_comp = 2 * l + 1 if pure else n_cart
        nctr = raw.shape[1]
        if ao_loc[ib + 1] - ao_loc[ib] != nctr * n_comp:
            raise RuntimeError(f"unexpected AO count of PySCF shell {ib}")
        for ic in range(nctr):
            d = [_libcint_sp_factor(l) * raw[i, ic] / primitive_norm(float(exps[i]), l)
                 for i in range(len(exps))]
            keep = [i for i in range(len(exps)) if d[i] != 0.0]
            if not keep:
                raise ValueError(f"PySCF shell {ib} has an all-zero contraction")
            e_keep = [float(exps[i]) for i in keep]
            d_keep = [d[i] for i in keep]
            sigma = math.sqrt(_axial_self_overlap(l, e_keep, d_keep))
            if not sigma > 0.0 or not math.isfinite(sigma):
                raise ValueError(f"PySCF shell {ib} has a contraction without a norm")
            shells.append({"center": center, "l": l, "exponents": e_keep,
                           "coefficients": [c / sigma for c in d_keep], "pure": pure})
            if pure:
                if l not in harmonics:
                    h = np.asarray(real_solid_harmonics(l))
                    c2s = np.asarray(gto.cart2sph(l, normalized="sp"), dtype=float)
                    # libcint: cart2sph acts on its Cartesians, which carry the s/p factor.
                    m, *_ = np.linalg.lstsq(h, c2s, rcond=None)
                    if np.abs(h @ m - c2s).max() > 1e-10 * max(1.0, np.abs(c2s).max()):
                        raise RuntimeError(
                            f"PySCF's l = {l} spherical functions are not the spec's "
                            f"real solid harmonics; refusing to guess the convention")
                    harmonics[l] = m
                block = _snap(sigma * harmonics[l])
            else:
                block = _snap(sigma * np.eye(n_cart))
            blocks.append((int(ao_loc[ib]) + ic * n_comp, block, np.linalg.inv(block)))
    n_ao = int(ao_loc[-1])
    return AOMap(shells, blocks, n_ao, pure)


# -- names --------------------------------------------------------------------

def _is_ks(mf):
    from pyscf import dft
    return isinstance(mf, dft.rks.KohnShamDFT)


def _spin_kind(mf):
    """'restricted', 'restricted_open' or 'unrestricted' for a supported SCF."""
    from pyscf import scf
    if isinstance(mf, scf.uhf.UHF):
        return "unrestricted"
    if isinstance(mf, scf.rohf.ROHF):
        return "restricted_open"
    if isinstance(mf, scf.hf.RHF):
        return "restricted"
    raise NotImplementedError(
        f"{type(mf).__name__}: the PySCF adapter writes RHF, UHF, ROHF, RKS, UKS "
        f"and ROKS results; GHF, Dirac and periodic objects are later tranches")


def method_name(mf) -> str:
    """The lowercase method name of a supported SCF object, e.g. ``"rks"``."""
    _pyscf()
    prefix = {"restricted": "r", "unrestricted": "u", "restricted_open": "ro"}[_spin_kind(mf)]
    return prefix + ("ks" if _is_ks(mf) else "hf")


def _functional(mf):
    xc = getattr(mf, "xc", None)
    if not isinstance(xc, str) or not xc.strip():
        raise ValueError("a Kohn-Sham object needs its xc functional name")
    return xc.strip().lower()


def basis_name(mol) -> str | None:
    """The basis-set name of a molecule, or ``None`` when it has none.

    A name given as one string, or the same string for every element, is
    returned as written. A per-element mixture or an explicit basis has no
    single name; the archive's ``basis`` payload is then the only record.
    """
    basis = mol.basis
    if isinstance(basis, str):
        return basis.strip() or None
    if isinstance(basis, dict) and basis:
        values = {v.strip() if isinstance(v, str) else None for v in basis.values()}
        if len(values) == 1:
            (value,) = values
            return value or None
    return None


def _check_molecule(mol):
    if hasattr(mol, "lattice_vectors"):
        raise NotImplementedError("periodic cells are a later WP-62 tranche")
    for i in range(mol.natm):
        if mol.atom_charge(i) == 0:
            raise ValueError(f"atom {i} ({mol.atom_symbol(i)}) is a ghost or dummy "
                             f"atom; the adapter writes real atoms only")


def _check_scf(mf):
    _pyscf()
    if getattr(mf, "mo_coeff", None) is None or getattr(mf, "mo_occ", None) is None:
        raise ValueError("the SCF object has no orbitals; run it first")
    if getattr(mf, "with_solvent", None) is not None:
        raise NotImplementedError("solvent models are not mapped by this tranche")
    _spin_kind(mf)
    _check_molecule(mf.mol)


# -- prepared data shared by the sections --------------------------------------

class _Prepared:
    def __init__(self, mf):
        _check_scf(mf)
        np = _np()
        self.mf = mf
        self.mol = mol = mf.mol
        self.spin = _spin_kind(mf)
        self.map = ao_map(mol)
        self.overlap = self.map.overlap(mol.intor_symmetric("int1e_ovlp"))
        dm = np.asarray(mf.make_rdm1(), dtype=float)
        if dm.ndim == 3:  # UHF and ROHF return the alpha and beta densities
            self.dm_alpha = self.map.density(dm[0])
            self.dm_beta = self.map.density(dm[1])
            self.dm_total = self.dm_alpha + self.dm_beta
        else:
            self.dm_alpha = self.dm_beta = None
            self.dm_total = self.map.density(dm)


def _symmetry_labels(mol, mo_coeff):
    """Irreducible-representation labels of one set of orbitals, or None."""
    if not getattr(mol, "symmetry", False):
        return None
    from pyscf import symm
    np = _np()
    n_mo = np.shape(mo_coeff)[-1]
    try:
        orbsym = getattr(mo_coeff, "orbsym", None)
        if orbsym is None or np.shape(orbsym) != (n_mo,):
            # An unrestricted set: label each spin channel by projection.
            orbsym = symm.label_orb_symm(mol, mol.irrep_id, mol.symm_orb,
                                         np.asarray(mo_coeff), check=False)
        return [str(symm.irrep_id2name(mol.groupname, int(ir))) for ir in orbsym]
    except Exception:  # symmetry-broken orbitals carry no clean labels
        return None


# -- sections -------------------------------------------------------------------

def add_structure(writer, mol, *, section_id: str = "structure",
                  label: str | None = None) -> str:
    """Add the molecule as a ``structure`` section and return its id.

    Positions are those of the calculation (``mol.atom_coords()``, bohr,
    including any reorientation PySCF applied for point-group symmetry),
    written in Å with the format's conversion factor (spec § 4.3), so a
    consumer recovers the bohr positions the integrals used. The section
    carries ``charge`` and ``multiplicity`` (``mol.spin + 1``).
    """
    _pyscf()
    from pyscf.data import elements
    _check_molecule(mol)
    coords = mol.atom_coords(unit="Bohr")
    atoms = []
    for i in range(mol.natm):
        symbol = mol.atom_pure_symbol(i)
        atoms.append({"symbol": symbol, "atomic_number": int(elements.charge(symbol)),
                      "position": [float(x) / BOHR_PER_ANGSTROM for x in coords[i]]})
    return writer.add_structure(atoms, charge=int(mol.charge),
                                multiplicity=int(mol.spin) + 1,
                                section_id=section_id, label=label)


def _ecp_projectors(mol, atom):
    """The scalar ECP radial functions of one atom, as spec § 5.3 projectors.

    Read from the integral input ``mol._ecpbas`` / ``mol._env`` that PySCF's
    ``ECPscalar`` integrals use (``gto.mole.make_ecp_env``), so the archive
    holds exactly the potential of the calculation. The conventions map
    without conversion:

    * ``ANG_OF = -1`` is the local part. It is labelled one above the
      highest semi-local ``l`` of the centre (0 when there is none).
    * ``RADI_POWER`` is the NWChem power ``n`` of ``r^(n-2) exp(-α r²)``:
      PySCF's NWChem parser (``gto.basis.parse_nwchem.parse_ecp``) keeps the
      tabulated ``n`` (NWChem documentation, "Effective Core Potentials":
      ``U_l(r) = Σ d r^(n-2) exp(-ζ r²)``), and the ``ECPscalar`` integrals
      use ``r^(n-2)`` (checked by radial quadrature in
      ``test_adapter_pyscf_ecp.py``). It is the spec's ``n``.
    * A semi-local channel multiplies the projector onto ``l`` in addition
      to the local part (``V = U_L + Σ_l |l⟩ U_l ⟨l|``), so it already is the
      difference ``U_l − U_L`` the spec stores. The def2 tables, for one,
      carry ``−U_L`` inside every channel.
    * Exponents and coefficients are atomic units (bohr⁻², hartree ·
      bohr^(2−n)), as in the basis-set libraries.

    Spin-orbit rows (``SO_TYPE_OF = 1``) are skipped: the scalar ``ECPscalar``
    Hamiltonian of a non-relativistic mean field does not contain them.
    Returns ``(projectors, has_spin_orbit)``.
    """
    from pyscf.gto import mole
    np = _np()
    channels: dict[int, list[tuple[float, float, int]]] = {}
    spin_orbit = False
    for row in np.asarray(mol._ecpbas):
        if int(row[mole.ATOM_OF]) != atom:
            continue
        if int(row[mole.SO_TYPE_OF]) != 0:
            spin_orbit = True
            continue
        n_prim, power = int(row[mole.NPRIM_OF]), int(row[mole.RADI_POWER])
        exps = mol._env[row[mole.PTR_EXP]:row[mole.PTR_EXP] + n_prim]
        coefs = mol._env[row[mole.PTR_COEFF]:row[mole.PTR_COEFF] + n_prim]
        channels.setdefault(int(row[mole.ANG_OF]), []).extend(
            (float(e), float(c), power) for e, c in zip(exps, coefs))
    semi = sorted(l for l in channels if l >= 0)
    projectors = []
    for l in ([-1] if -1 in channels else []) + semi:
        terms = channels[l]
        p: dict[str, Any] = {"l": (semi[-1] + 1 if semi else 0) if l < 0 else l,
                             "exponents": [t[0] for t in terms],
                             "coefficients": [t[1] for t in terms],
                             "powers": [t[2] for t in terms]}
        if l < 0:
            p["local"] = True
        projectors.append(p)
    return projectors, spin_orbit


def _ecp_entries(mol, projectors: bool = True):
    """``basis.ecp`` entries, with the radial functions when ``projectors``."""
    entries = []
    for i in range(mol.natm):
        core = int(mol.atom_nelec_core(i))
        radial, spin_orbit = _ecp_projectors(mol, i) if projectors else ([], False)
        if core <= 0 and not radial:
            continue
        entry: dict[str, Any] = {"center": i, "n_core_electrons": core}
        ecp = mol.ecp
        name = ecp if isinstance(ecp, str) else (
            ecp.get(mol.atom_symbol(i), ecp.get(mol.atom_pure_symbol(i)))
            if isinstance(ecp, dict) else None)
        if isinstance(name, str) and name.strip():
            entry["label"] = name.strip()
        if radial:
            entry["projectors"] = radial
        if spin_orbit:
            entry["x_pyscf_spin_orbit_omitted"] = True
        entries.append(entry)
    return entries or None


def _add_wavefunction(writer, p, *, structure_ref, section_id, density, label,
                      ecp_projectors=True):
    np = _np()
    mf = p.mf
    kwargs: dict[str, Any] = {}
    if p.spin == "unrestricted":
        ca, cb = mf.mo_coeff
        kwargs.update(
            mo_coefficients_alpha=p.map.coefficients(ca),
            mo_coefficients_beta=p.map.coefficients(cb),
            energies_alpha=np.asarray(mf.mo_energy[0], dtype=float).tolist(),
            energies_beta=np.asarray(mf.mo_energy[1], dtype=float).tolist(),
            occupations_alpha=np.asarray(mf.mo_occ[0], dtype=float).tolist(),
            occupations_beta=np.asarray(mf.mo_occ[1], dtype=float).tolist())
        labels = {spin: _symmetry_labels(p.mol, c) for spin, c in (("alpha", ca), ("beta", cb))}
        labels = {k: v for k, v in labels.items() if v is not None} or None
        if density:
            kwargs.update(density_matrix_alpha=p.dm_alpha, density_matrix_beta=p.dm_beta)
    else:
        kwargs.update(
            mo_coefficients=p.map.coefficients(mf.mo_coeff),
            energies=np.asarray(mf.mo_energy, dtype=float).tolist(),
            occupations=np.asarray(mf.mo_occ, dtype=float).tolist())
        labels = _symmetry_labels(p.mol, mf.mo_coeff)
        if density:
            kwargs["density_matrix"] = p.dm_total
    sid = writer.add_wavefunction_gto(
        p.map.shells, spin=p.spin, pure=p.map.pure, structure_ref=structure_ref,
        overlap=p.overlap, n_electrons=int(p.mol.nelectron),
        ecp=_ecp_entries(p.mol, ecp_projectors),
        symmetry_labels=labels, section_id=section_id, label=label, **kwargs)
    writer.annotate_section(sid, structure_ref=structure_ref)
    return sid


def add_wavefunction(writer, mf, *, structure_ref: str = "structure",
                     section_id: str = "wf", density: bool = True,
                     label: str | None = None, ecp_projectors: bool = True) -> str:
    """Add the orbitals of ``mf`` as a ``wavefunction.gto`` section.

    Writes the ``basis`` payload (spherical or Cartesian as ``mol.cart``
    says), ``mo_metadata`` (spin shape, orbital energies in hartree,
    occupations, point-group labels when the molecule has symmetry), the
    coefficients in the QVF AO basis, the AO ``overlap`` of PySCF's integral
    code in that basis, and with ``density`` the one-particle density
    matrices (``density_matrix``, or ``density_matrix_alpha`` and
    ``density_matrix_beta``). ``n_electrons`` and a ``basis.ecp`` entry per
    ECP atom (core count and ECP name) are written too. With
    ``ecp_projectors`` (the default) each entry also carries the scalar
    radial functions of the potential (``_ecp_projectors``), and
    ``ecp_projectors=False`` keeps only the label and core count. Spin-orbit
    ECP terms are never written. An entry whose PySCF ECP has them is marked
    ``x_pyscf_spin_orbit_omitted: true``. The structure ``structure_ref`` names must already be
    in the writer. Returns the section id.
    """
    return _add_wavefunction(writer, _Prepared(mf), structure_ref=structure_ref,
                             section_id=section_id, density=density, label=label,
                             ecp_projectors=ecp_projectors)


def _populations(p):
    np = _np()
    s = p.overlap
    eigenvalues, vectors = np.linalg.eigh(s)
    if eigenvalues.min() <= 0.0:
        raise ValueError("the AO overlap matrix is not positive definite")
    s_half = (vectors * np.sqrt(eigenvalues)) @ vectors.T
    atoms = np.asarray(p.map.ao_atoms())
    natm = p.mol.natm
    z_eff = np.asarray(p.mol.atom_charges(), dtype=float)

    def per_atom(values):
        return np.bincount(atoms, weights=values, minlength=natm)

    def mulliken(dm):
        return per_atom(np.einsum("ij,ji->i", dm, s))

    def loewdin(dm):
        return per_atom(np.einsum("ij,jk,ki->i", s_half, dm, s_half))

    out = {"mulliken_charge": z_eff - mulliken(p.dm_total),
           "loewdin_charge": z_eff - loewdin(p.dm_total)}
    if p.dm_alpha is not None:
        out["spin_population"] = mulliken(p.dm_alpha) - mulliken(p.dm_beta)
    return out


def _add_populations(writer, p, *, structure_ref, section_id, label):
    sid = writer.add_atom_properties(section_id=section_id, label=label, **_populations(p))
    writer.annotate_section(sid, structure_ref=structure_ref)
    return sid


def add_populations(writer, mf, *, structure_ref: str = "structure",
                    section_id: str = "atom_properties",
                    label: str | None = None) -> str:
    """Add Mulliken and Löwdin charges as an ``atom_properties`` section.

    ``mulliken_charge`` and ``loewdin_charge`` are net atomic charges, the
    effective nuclear charge (``Z`` minus any ECP core) minus the atom's
    population, computed from the archive's own density and overlap in the
    QVF AO basis. For a spherical basis that basis is PySCF's up to signs and
    order, so the values are PySCF's; for a Cartesian basis the Löwdin
    partition uses the Appendix A.1 AO normalization. An open-shell
    (``UHF``/``ROHF``-type) result adds ``spin_population``, the Mulliken α
    minus β population. Returns the section id.
    """
    return _add_populations(writer, _Prepared(mf), structure_ref=structure_ref,
                            section_id=section_id, label=label)


def add_job_spec(writer, mf, *, structure_ref: str = "structure",
                 section_id: str = "job_spec", label: str | None = None) -> str:
    """Add the request that produced ``mf`` as a ``job.spec`` section.

    ``job_type`` is ``"molecular"``; ``method`` is ``method_name``;
    ``basis`` is ``basis_name`` (omitted when the basis has no single
    name); ``functional`` is the lowercase ``xc`` string of a Kohn–Sham
    object and JSON ``null`` for Hartree–Fock, which takes none (spec § 5.9);
    ``charge`` and ``multiplicity`` are the molecule's; ``tasks`` is
    ``["single_point"]``; ``scf`` carries ``conv_tol``, ``max_cycle`` and the
    initial guess name. A ``disp`` setting becomes ``dispersion`` and an
    X2C Hamiltonian ``relativistic: "x2c"``. Returns the section id.
    """
    _check_scf(mf)
    mol = mf.mol
    scf_block: dict[str, Any] = {"convergence": float(mf.conv_tol),
                                 "max_iterations": int(mf.max_cycle)}
    guess = getattr(mf, "init_guess", None)
    if isinstance(guess, str) and guess.strip():
        scf_block["guess"] = guess.strip().lower()
    kwargs: dict[str, Any] = {
        "job_type": "molecular", "method": method_name(mf),
        "functional": _functional(mf) if _is_ks(mf) else None,
        "charge": int(mol.charge), "multiplicity": int(mol.spin) + 1,
        "tasks": ["single_point"], "structure_ref": structure_ref, "scf": scf_block,
        "section_id": section_id, "label": label}
    name = basis_name(mol)
    if name is not None:
        kwargs["basis"] = name
    disp = getattr(mf, "disp", None)
    if isinstance(disp, str) and disp.strip():
        kwargs["dispersion"] = disp.strip().lower()
    if getattr(mf, "with_x2c", None) is not None:
        kwargs["relativistic"] = "x2c"
    return writer.add_job_spec(**kwargs)


def _software(mf):
    pyscf = _pyscf()
    software = [{"name": "PySCF", "version": str(pyscf.__version__)}]
    if _is_ks(mf):
        library = getattr(getattr(mf, "_numint", None), "libxc", None)
        name = getattr(library, "__name__", "")
        entry = {"name": "XCFun" if name.endswith("xcfun") else "Libxc"}
        version = getattr(library, "__version__", None)
        if isinstance(version, str) and version:
            entry["version"] = version
        software.append(entry)
    return software


def set_scf_provenance(writer, mf) -> None:
    """Set the root ``provenance`` and ``energies`` blocks for ``mf``.

    ``provenance`` gets the method, functional (Kohn–Sham), basis name,
    charge, multiplicity, electron count, convergence flag and iteration
    count, ``scf_energy`` (``mf.e_tot``, hartree), ``run_status``
    (``"converged"``, or ``"failed"`` for an unconverged SCF) and
    ``software``. ``energies`` restates ``e_tot`` as the total of that
    method, with the nuclear repulsion and the electronic energy (and a
    dispersion correction, when PySCF reports one) as components.
    """
    _check_scf(mf)
    mol = mf.mol
    method = method_name(mf)
    e_tot = float(mf.e_tot)
    fields: dict[str, Any] = {
        "method": method, "charge": int(mol.charge),
        "multiplicity": int(mol.spin) + 1, "n_electrons": int(mol.nelectron),
        "scf_converged": bool(mf.converged), "scf_energy_eh": e_tot,
        "run_status": "converged" if mf.converged else "failed",
        "software": _software(mf)}
    if _is_ks(mf):
        fields["functional"] = _functional(mf)
    name = basis_name(mol)
    if name is not None:
        fields["basis"] = name
    cycles = getattr(mf, "cycles", None)
    if isinstance(cycles, int) and not isinstance(cycles, bool) and cycles >= 0:
        fields["n_scf_iterations"] = cycles
    writer.set_provenance(**fields)
    e_nuc = float(mf.energy_nuc())
    components = {"nuclear_repulsion": e_nuc}
    if not getattr(mf, "disp", None):
        components["electronic"] = e_tot - e_nuc
    else:
        # e_tot includes the dispersion correction; without PySCF's share
        # the electronic energy is unknown and is left out.
        e_disp = (getattr(mf, "scf_summary", None) or {}).get("dispersion")
        if e_disp is not None:
            components["dispersion"] = float(e_disp)
            components["electronic"] = e_tot - e_nuc - float(e_disp)
    writer.set_energies(e_tot, method=method, components=components)


def add_scf(writer, mf, *, structure_id: str = "structure", wavefunction_id: str = "wf",
            properties_id: str = "atom_properties", job_id: str = "job_spec",
            populations: bool = True, job_spec: bool = True, density: bool = True,
            provenance: bool = True, ecp_projectors: bool = True) -> dict[str, str]:
    """Add everything the adapter writes for ``mf`` to an open QVF writer.

    ``writer`` is a ``qvf.Writer`` (or the compatibility ``QvfWriter``).
    Adds the structure, the wavefunction, the populations (``populations``),
    the job specification (``job_spec``) and the provenance and energies
    blocks (``provenance``). ``ecp_projectors`` is that of
    ``add_wavefunction``. Returns the ids of the sections added, keyed
    ``structure``, ``wavefunction``, ``atom_properties`` and ``job_spec``.
    """
    p = _Prepared(mf)
    ids = {"structure": add_structure(writer, p.mol, section_id=structure_id)}
    ids["wavefunction"] = _add_wavefunction(
        writer, p, structure_ref=ids["structure"], section_id=wavefunction_id,
        density=density, label=None, ecp_projectors=ecp_projectors)
    if populations:
        ids["atom_properties"] = _add_populations(
            writer, p, structure_ref=ids["structure"], section_id=properties_id, label=None)
    if job_spec:
        ids["job_spec"] = add_job_spec(writer, mf, structure_ref=ids["structure"],
                                       section_id=job_id)
    if provenance:
        set_scf_provenance(writer, mf)
    return ids


def write_scf(mf, target, *, calculation: str | None = None, **options) -> dict[str, str]:
    """Write a complete archive for ``mf`` to ``target`` (a path or binary sink).

    The archive's ``source`` is PySCF and its version; ``calculation``
    defaults to ``method/functional/basis``. ``options`` are those of
    ``add_scf``. Returns the section ids that ``add_scf`` returns.

    ```python
    from pyscf import gto, dft
    from qvf.adapters.pyscf import write_scf

    mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.469; H 0 -0.757 -0.469",
                basis="def2-svp")
    mf = dft.RKS(mol, xc="pbe0").run()
    write_scf(mf, "water.qvf")
    ```
    """
    pyscf = _pyscf()
    from ..writer import Writer
    _check_scf(mf)
    if calculation is None:
        parts = [method_name(mf), _functional(mf) if _is_ks(mf) else None,
                 basis_name(mf.mol)]
        calculation = "/".join(x for x in parts if x)
    with Writer(target, program="PySCF", version=str(pyscf.__version__),
                calculation=calculation) as writer:
        ids = add_scf(writer, mf, **options)
    return ids


# -- molecular closed-shell excited states ------------------------------------

def _td_real(value, name, shape=None):
    np = _np()
    array = np.asarray(value)
    if np.iscomplexobj(array):
        raise NotImplementedError(f"complex {name} are not supported by the TD adapter")
    if array.dtype.kind not in "fiu" or not np.isfinite(array).all():
        raise ValueError(f"{name} must contain finite real numbers")
    if shape is not None and array.shape != shape:
        raise ValueError(f"{name} has shape {array.shape}, expected {shape}")
    return np.asarray(array, dtype=float)


class _PreparedTD:
    """Check source results before writing or opening a destination file."""

    def __init__(self, td, ntos):
        _pyscf()
        from pyscf.tdscf import rhf, rks
        from ..units import HARTREE_IN_EV
        np = _np()
        if not isinstance(td, (rhf.TDA, rhf.TDHF)) or isinstance(td, (rks.dRPA, rks.dTDA)):
            raise NotImplementedError("the TD adapter supports molecular RHF/RKS TDA and TDHF/TDDFT")
        if getattr(td, "with_solvent", None) is not None:
            raise NotImplementedError("solvent response is not mapped by the TD adapter")
        mf = td._scf
        _check_scf(mf)
        if _spin_kind(mf) != "restricted" or mf.mol.spin != 0:
            raise NotImplementedError("the TD adapter needs a closed-shell RHF/RKS reference")
        if not mf.converged:
            raise ValueError("the SCF reference must be converged")
        if td.mol is not mf.mol:
            raise ValueError("TD and SCF results must refer to the same molecule")
        occ = _td_real(mf.mo_occ, "MO occupations")
        if occ.ndim != 1 or not np.isin(occ, [0, 2]).all():
            raise NotImplementedError("fractional or open-shell occupations are not supported")
        coeff = _td_real(mf.mo_coeff, "MO coefficients", (mf.mol.nao_nr(), len(occ)))
        _td_real(mf.mo_energy, "MO energies", occ.shape)
        reference_energy = float(_td_real(mf.e_tot, "SCF energy", ()))
        mask = np.asarray(td.get_frozen_mask())
        if mask.dtype.kind != "b" or mask.shape != occ.shape:
            raise ValueError("the frozen-orbital mask must match the reference MO rows")
        occupied = np.flatnonzero(mask & (occ == 2))
        virtual = np.flatnonzero(mask & (occ == 0))
        if not len(occupied) or not len(virtual):
            raise ValueError("the TD result needs active occupied and virtual orbitals")
        energies = _td_real(td.e, "excitation energies")
        if energies.ndim != 1 or not len(energies) or (energies <= 0).any():
            raise ValueError("the TD result needs a non-empty array of positive excitation energies")
        n = len(energies)
        if not isinstance(td.xy, (list, tuple)) or len(td.xy) != n:
            raise ValueError("X/Y amplitudes must have one pair per excited state")
        converged = np.asarray(td.converged)
        if converged.dtype.kind != "b" or converged.shape != (n,) or not converged.all():
            raise ValueError("every exported TD root must be converged")
        if not isinstance(td.singlet, (bool, np.bool_)):
            raise ValueError("the restricted TD result must identify singlet or triplet states")
        self.singlet = bool(td.singlet)
        # CasidaTDDFT inherits both TDHF and TDA; it is a full-response method.
        self.tda = isinstance(td, rhf.TDA) and not isinstance(td, rhf.TDHF)
        eligible_ntos = self.singlet
        if ntos is not None and not isinstance(ntos, bool):
            raise TypeError("ntos must be True, False or None")
        if ntos is True and not eligible_ntos:
            raise NotImplementedError("paired NTO export currently supports singlet RHF/RKS response only")
        self.with_ntos = eligible_ntos if ntos is None else ntos
        self.method = "tddft" if _is_ks(mf) else ("cis" if self.tda else "tdhf")
        self.mf = mf
        self.states = []
        self.ntos = []
        shape = (len(occupied), len(virtual))
        for k, (energy, xy) in enumerate(zip(energies, td.xy), 1):
            if not isinstance(xy, (tuple, list)) or len(xy) != 2:
                raise ValueError(f"state {k}: expected an X/Y pair")
            x = _td_real(xy[0], f"state {k} X amplitudes", shape)
            y = _td_real(xy[1], f"state {k} Y amplitudes")
            if y.shape == () and y == 0:
                y = np.zeros(shape)
            elif y.shape != shape:
                raise ValueError(f"state {k}: Y amplitudes must match X or be scalar zero")
            if self.tda and np.any(y != 0):
                raise ValueError(f"state {k}: TDA must have zero Y amplitudes")
            weights = x * x - y * y
            norm = float(weights.sum())
            if not math.isclose(norm, 0.5, rel_tol=1e-6, abs_tol=1e-8):
                raise ValueError(f"state {k}: X/Y amplitudes must retain PySCF normalization "
                                 "sum(X^2-Y^2)=1/2; do not renormalize them before export")
            state = {
                "index": k,
                "excitation_energy": {"value": float(energy * HARTREE_IN_EV), "units": "ev"},
                "total_energy": {"value": reference_energy + float(energy), "units": "hartree"},
                "multiplicity": 1 if self.singlet else 3,
                "contributions": [
                    {"from_mo": int(i), "to_mo": int(a), "weight": float(weights[u, v] / norm)}
                    for u, i in enumerate(occupied) for v, a in enumerate(virtual)],
            }
            self.states.append(state)
            if self.with_ntos:
                # Martin, JCP 118, 4775 (2003), eqs. 4-10; Herbert (2023),
                # doi:10.1016/B978-0-323-91738-4.00005-1, eq. 3.51.
                # PySCF uses half-normalized restricted X/Y. For full
                # response, retain the spin-adapted charge transition
                # sqrt(2)(X+Y), including its non-unit Frobenius norm.
                # Keep the accepted TDA convention; neither path mutates X/Y.
                transition = (x / math.sqrt(norm) if self.tda
                              else math.sqrt(2) * (x + y))
                u, singular, vh = np.linalg.svd(transition, full_matrices=False)
                self.ntos.append((coeff[:, occupied] @ u, coeff[:, virtual] @ vh.T,
                                  (singular * singular).tolist()))
        dipoles = _td_real(td.transition_dipole(), "transition dipoles", (n, 3))
        strengths = _td_real(td.oscillator_strength(gauge="length"), "oscillator strengths", (n,))
        if (strengths < 0).any():
            raise ValueError("oscillator strengths must be nonnegative")
        for state, dipole, strength in zip(self.states, dipoles, strengths):
            state["transition_dipole"] = dipole.tolist()
            state["oscillator_strength"] = float(strength)
        self.reference = _Prepared(mf)


def _add_td(writer, p, *, states_id, spectrum_id, nto_prefix, scf_options):
    ids = add_scf(writer, p.mf, **scf_options)
    ref = p.reference
    ecp = _ecp_entries(ref.mol, scf_options.get("ecp_projectors", True))
    holes, particles = [], []
    for state, (hole, particle, weights) in zip(p.states, p.ntos):
        for role, coefficients, output in (("hole", hole, holes),
                                            ("particle", particle, particles)):
            sid = writer.add_wavefunction_gto(
                ref.map.shells, mo_coefficients=ref.map.coefficients(coefficients),
                occupations=weights, energies=[0.0] * len(weights), spin="restricted",
                pure=ref.map.pure, structure_ref=ids["structure"],
                orbital_kind="nto", state_ref=states_id, state_index=state["index"],
                nto_role=role, overlap=ref.overlap, ecp=ecp,
                section_id=f"{nto_prefix}_{state['index']}_{role}")
            writer.annotate_section(sid, structure_ref=ids["structure"])
            output.append(sid)
        state["nto_ref"] = holes[-1]
    ids["excited_states"] = writer.add_excited_states(
        p.states, wavefunction_ref=ids["wavefunction"], section_id=states_id,
        meta={"method": p.method, "reference_state": "SCF reference", "n_states": len(p.states),
              "tda": p.tda, "spin_flip": False})
    writer.annotate_section(states_id, structure_ref=ids["structure"])
    ids["spectrum"] = writer.add_spectrum("spectra.uvvis", section_id=spectrum_id, payload={
        "frequencies": [s["excitation_energy"]["value"] for s in p.states],
        "intensities": [s["oscillator_strength"] for s in p.states],
        "x_units": "ev", "y_units": "dimensionless", "states_ref": states_id,
        "peaks": [{"x": s["excitation_energy"]["value"], "y": s["oscillator_strength"],
                   "assignment": {"state_index": s["index"]}} for s in p.states]})
    writer.annotate_section(spectrum_id, structure_ref=ids["structure"])
    ids.update(nto_holes=holes, nto_particles=particles)
    return ids


def add_td(writer, td, *, states_id: str = "excited_states", spectrum_id: str = "uvvis",
           nto_prefix: str = "nto", ntos: bool | None = None, **scf_options) -> dict[str, Any]:
    """Add a closed-shell molecular PySCF TDA/TDHF/TDDFT result.

    Adds the SCF reference through ``add_scf`` (``scf_options`` are its
    arguments), an ``excited_states`` table, and length-gauge UV-visible
    sticks with state assignments. All roots must be converged and retain
    PySCF's real X/Y normalization. Frozen orbital indices are mapped back
    to rows of the full reference wavefunction. Contributions are the
    normalized X²−Y² weights, including negative terms when present.

    ``ntos=None`` writes paired NTOs for singlet TDA and full response;
    ``False`` omits them, and ``True`` requires a supported NTO result.
    Triplet NTOs, unrestricted/open-shell references, complex amplitudes,
    solvent and direct-RPA variants are not mapped. Source arrays are not
    changed and no SCF or response iterations are run. NTO energies are
    zero and their occupations are squared singular-value transition weights.

    TDA uses the normalized X matrix, with weights summing to one. Full
    response uses ``sqrt(2) * (X + Y)`` in PySCF's restricted convention,
    with no rescaling of the resulting weights. Their sum may differ from
    one. These are spatial charge-transition NTOs for real spin-independent
    Hermitian operators; their orbital-pair contraction times ``sqrt(2)``
    recovers the source length-gauge transition dipole. They do not encode
    the separate X/Y blocks needed for anti-Hermitian operators.

    Root provenance/energies and the optional ``job.spec`` describe the
    **SCF reference**, not the excitation calculation. The excitation
    method and TDA flag are in the states' metadata; state totals are the
    reference total plus the vertical excitation. The job is therefore
    not a complete request for repeating TD response. Pass ``job_spec=False``
    to omit the reference request.

    Returns the ``add_scf`` IDs plus ``excited_states``, ``spectrum`` and
    ordered lists ``nto_holes`` / ``nto_particles`` (empty when omitted).
    """
    return _add_td(writer, _PreparedTD(td, ntos), states_id=states_id,
                   spectrum_id=spectrum_id, nto_prefix=nto_prefix, scf_options=scf_options)


def write_td(td, target, *, calculation: str | None = None,
             states_id: str = "excited_states", spectrum_id: str = "uvvis",
             nto_prefix: str = "nto", ntos: bool | None = None,
             **scf_options) -> dict[str, Any]:
    """Write a complete archive for a PySCF TD result, with ``add_td`` options.

    Source-result checks finish before the target is opened. The source
    program/version is PySCF's; the default calculation label includes the
    excitation method, reference method, functional and basis. The target
    is a path or binary sink, as for ``write_scf``.
    """
    p = _PreparedTD(td, ntos)
    from ..writer import Writer
    if calculation is None:
        method = p.method + ("-tda" if p.tda and p.method != "cis" else "")
        parts = [method, method_name(p.mf), _functional(p.mf) if _is_ks(p.mf) else None,
                 basis_name(p.mf.mol)]
        calculation = "/".join(x for x in parts if x)
    with Writer(target, program="PySCF", version=str(_pyscf().__version__),
                calculation=calculation) as writer:
        return _add_td(writer, p, states_id=states_id, spectrum_id=spectrum_id,
                       nto_prefix=nto_prefix, scf_options=scf_options)
