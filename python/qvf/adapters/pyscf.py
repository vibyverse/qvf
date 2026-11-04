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

Periodic Γ-point (k = 0) calculations of a ``pbc.gto.Cell`` are written the
same way (``pbc.scf.RHF``, ``UHF`` and ``pbc.dft.RKS``, ``UKS``, real
orbitals): the cell as a periodic ``structure`` with ``pbc`` and
``lattice_vectors``, the orbitals as a ``wavefunction.gto`` with
``k_point`` ``[0, 0, 0]``. Its ``overlap`` member is the AO overlap of the
spec, ``⟨χ_μ|χ_ν⟩`` of the cell's AOs without images. The Γ-point orbitals
are orthonormal in the lattice-summed overlap ``Σ_T ⟨χ_μ|χ_ν(r − T)⟩``,
which a consumer rebuilds from the basis and the lattice. The populations
use it. GTH pseudopotentials are recorded by name and core count only.
k-point sampling (``KRHF`` and the like) and a nonzero ``kpt`` are refused.

``add_td`` and ``write_td`` add molecular closed-shell TDA/TDHF/TDDFT
states and UV-visible sticks beside that SCF reference. Singlet TDA and
full response can also carry paired natural transition orbitals. Full
response retains both X and Y; triplet NTOs are not mapped (see ``add_td``).

``add_density_grid`` and ``add_orbital_grid`` evaluate the SCF density, the
spin density of an open-shell reference and selected orbitals on a regular
grid (``volume.density``, ``volume.spin``, ``volume.orbital``);
``molecular_grid`` builds the default box around the nuclei.

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
    "AOMap", "DEFAULT_GRID_MARGIN_BOHR", "DEFAULT_GRID_SPACING_BOHR", "ao_map",
    "add_density_grid", "add_job_spec", "add_orbital_grid", "add_populations", "add_scf",
    "add_td", "add_structure", "add_wavefunction", "basis_name", "cartesian_powers",
    "method_name", "molecular_grid", "real_solid_harmonics", "set_scf_provenance",
    "write_scf", "write_td",
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


def _is_cell(mol):
    """Whether ``mol`` is a periodic ``pbc.gto.Cell``."""
    return hasattr(mol, "lattice_vectors")


def _periodic_spin_kind(mf):
    """'restricted' or 'unrestricted' for a supported Γ-point SCF of a cell."""
    from pyscf.pbc.scf import ghf, hf, khf, rohf, uhf
    if isinstance(mf, khf.KSCF):
        raise NotImplementedError(
            f"{type(mf).__name__}: k-point sampling is not mapped. The adapter writes "
            f"Γ-point pbc.scf.RHF/UHF and pbc.dft.RKS/UKS objects (k-resolved "
            f"orbitals belong in wavefunction.bloch, a later tranche)")
    if isinstance(mf, uhf.UHF):
        return "unrestricted"
    if isinstance(mf, hf.RHF) and not isinstance(mf, (rohf.ROHF, ghf.GHF)):
        return "restricted"
    raise NotImplementedError(
        f"{type(mf).__name__}: for periodic cells the PySCF adapter writes Γ-point "
        f"RHF, UHF, RKS and UKS results only")


def _spin_kind(mf):
    """'restricted', 'restricted_open' or 'unrestricted' for a supported SCF."""
    from pyscf import scf
    if _is_cell(getattr(mf, "mol", None)):
        return _periodic_spin_kind(mf)
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
    if _is_cell(mol) and int(mol.dimension) not in (1, 2, 3):
        raise NotImplementedError(f"a cell of dimension {mol.dimension} is not mapped, "
                                  f"use a molecular gto.Mole")
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
    if _is_cell(mf.mol):
        _check_gamma(mf)


def _check_gamma(mf):
    """Refuse a periodic SCF that is not a real Γ-point result."""
    np = _np()
    kpt = np.asarray(getattr(mf, "kpt", np.zeros(3)), dtype=float).reshape(-1)
    if kpt.shape != (3,) or np.abs(kpt).max() > 1e-12:
        raise NotImplementedError(
            f"kpt {kpt.tolist()} is not the Γ point. Orbitals at other k-points "
            f"belong in wavefunction.bloch, a later tranche")
    for name in ("mo_coeff", "mo_energy", "mo_occ"):
        if np.iscomplexobj(np.asarray(getattr(mf, name))):
            raise NotImplementedError(f"complex {name}: the adapter writes real "
                                      f"Γ-point orbitals only")


# -- prepared data shared by the sections --------------------------------------

class _Prepared:
    def __init__(self, mf):
        _check_scf(mf)
        np = _np()
        self.mf = mf
        self.mol = mol = mf.mol
        self.spin = _spin_kind(mf)
        self.map = ao_map(mol)
        # For a cell this is the AO overlap of spec § 5.3, without images.
        # The Γ-point orbitals are orthonormal in the lattice sum, ``metric``.
        self.overlap = self.map.overlap(mol.intor_symmetric("int1e_ovlp"))
        self.periodic = _is_cell(mol)
        self.metric = (self.map.overlap(mol.pbc_intor("int1e_ovlp", hermi=1))
                       if self.periodic else self.overlap)
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

    For a ``pbc.gto.Cell`` the section is periodic: ``lattice_vectors`` are
    the rows of ``cell.lattice_vectors()`` (bohr) in Å through ``Lattice``,
    and ``pbc`` marks the first ``cell.dimension`` axes periodic, as PySCF
    treats them. The positions are those of the cell, not wrapped.
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
    periodic: dict[str, Any] = {}
    if _is_cell(mol):
        from ..model import Lattice
        rows_bohr = _np().asarray(mol.lattice_vectors(), dtype=float)
        lattice = Lattice.from_columns_bohr(rows_bohr.T.tolist())
        periodic = {"pbc": [i < int(mol.dimension) for i in range(3)],
                    "lattice_vectors": lattice.rows_angstrom()}
    return writer.add_structure(atoms, charge=int(mol.charge),
                                multiplicity=int(mol.spin) + 1,
                                section_id=section_id, label=label, **periodic)


def _pseudo_name(mol, i):
    """The name of atom ``i``'s GTH pseudopotential in a cell, else None."""
    if not _is_cell(mol) or mol.atom_symbol(i) not in (getattr(mol, "_pseudo", None) or {}):
        return None
    pseudo = mol.pseudo
    name = pseudo if isinstance(pseudo, str) else (
        pseudo.get(mol.atom_symbol(i), pseudo.get(mol.atom_pure_symbol(i)))
        if isinstance(pseudo, dict) else None)
    return name.strip() if isinstance(name, str) and name.strip() else "gth"


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
        gth = _pseudo_name(mol, i)
        if gth is not None:
            # A GTH potential is not the semi-local form of spec § 5.3: name
            # and core count only, also for a zero core (hydrogen).
            entries.append({"center": i, "n_core_electrons": core, "label": gth})
            continue
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
    if p.periodic:
        kwargs["k_point"] = [0.0, 0.0, 0.0]
        labels = None  # space-group labels are not point-group irreps
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

    A Γ-point cell result adds ``k_point`` ``[0, 0, 0]``. Its ``overlap`` is
    the AO overlap without images. The orbitals are orthonormal in the
    lattice-summed overlap instead (see the module notes).
    """
    return _add_wavefunction(writer, _Prepared(mf), structure_ref=structure_ref,
                             section_id=section_id, density=density, label=label,
                             ecp_projectors=ecp_projectors)


def _populations(p):
    np = _np()
    s = p.metric
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
    X2C Hamiltonian ``relativistic: "x2c"``. A Γ-point cell result is a
    ``"periodic"`` job with ``kpoints`` ``[1, 1, 1]`` and ``options`` naming
    the FFT ``mesh``, the integral scheme (``fftdf``, ``gdf``, …), the
    pseudopotential and the exact-exchange treatment ``exxdiv`` (``null``
    for none).
    Returns the section id.
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
    if _is_cell(mol):
        kwargs.update(job_type="periodic", kpoints=[1, 1, 1], options=_periodic_options(mf))
    return writer.add_job_spec(**kwargs)


def _periodic_options(mf):
    """PySCF settings of a Γ-point job that change its energy."""
    cell = mf.mol
    options: dict[str, Any] = {"mesh": [int(n) for n in cell.mesh],
                               "integrals": type(mf.with_df).__name__.lower()}
    pseudo = cell.pseudo
    if isinstance(pseudo, str) and pseudo.strip():
        options["pseudo"] = pseudo.strip()
    elif isinstance(pseudo, dict) and pseudo:
        options["pseudo"] = {str(k): str(v) for k, v in pseudo.items()}
    exxdiv = getattr(mf, "exxdiv", None)  # exact-exchange divergence treatment
    options["exxdiv"] = exxdiv if isinstance(exxdiv, str) else None
    return options


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
    dispersion correction, when PySCF reports one) as components. For a
    cell, ``nuclear_repulsion`` is PySCF's Ewald energy of the point nuclei
    (effective charges ``Z`` minus the pseudopotential core) per cell, and
    ``total`` includes any exchange-divergence (Madelung) correction.
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
        if _is_cell(getattr(td, "mol", None)) or _is_cell(getattr(getattr(td, "_scf", None), "mol", None)):
            raise NotImplementedError("periodic response is not mapped by the TD adapter")
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


# -- densities and orbitals on a regular grid ----------------------------------

#: Default distance, in bohr, between the outermost nuclei and the grid box faces.
DEFAULT_GRID_MARGIN_BOHR = 4.0
#: Default point spacing, in bohr, of ``molecular_grid``.
DEFAULT_GRID_SPACING_BOHR = 0.2
# AO values held per evaluation block: 2**23 doubles, 64 MiB.
_AO_BLOCK_VALUES = 1 << 23


def _refuse_periodic_grid(mol):
    if _is_cell(mol):
        raise NotImplementedError("grids of a periodic cell are not mapped: they need "
                                  "lattice-summed AOs on a periodic grid")


def molecular_grid(mol, *, spacing: float | None = DEFAULT_GRID_SPACING_BOHR,
                   shape=None, margin: float = DEFAULT_GRID_MARGIN_BOHR) -> dict[str, Any]:
    """A finite, axis-aligned grid around the nuclei of ``mol`` (spec § 5.2).

    The box spans the nuclear coordinates (bohr) widened by ``margin`` bohr on
    every side. With ``shape`` (``[nx, ny, nz]``, each at least 2) the box
    is sampled as PySCF's ``tools.cubegen.Cube`` does with the same ``margin``:
    the origin is the lower box corner and the step along axis ``a`` is the
    box edge over ``shape[a] − 1``, so both ends lie on the box faces. Without
    ``shape``, ``spacing`` bohr is used on every axis and the box edge is
    rounded up to a whole number of steps, centred on the nuclei. Returns the
    ``origin``, ``voxel_vectors`` and ``shape`` of a grid payload, in bohr.
    """
    _pyscf()
    np = _np()
    _check_molecule(mol)
    _refuse_periodic_grid(mol)
    margin = float(margin)
    if not math.isfinite(margin) or margin < 0.0:
        raise ValueError("margin must be a finite number of bohr >= 0")
    coords = np.asarray(mol.atom_coords(unit="Bohr"), dtype=float)
    lo, hi = coords.min(axis=0) - margin, coords.max(axis=0) + margin
    box = hi - lo
    if shape is not None:
        shape = [int(n) for n in shape]
        if len(shape) != 3 or min(shape) < 2:
            raise ValueError("shape must be three point counts >= 2")
        if (box <= 0.0).any():
            raise ValueError("the box is flat along an axis; give margin > 0")
        steps = [box[a] / (shape[a] - 1) for a in range(3)]
        origin = lo
    else:
        h = float(spacing) if spacing is not None else math.nan
        if not math.isfinite(h) or h <= 0.0:
            raise ValueError("spacing must be a positive number of bohr")
        counts = np.maximum(np.ceil(box / h - 1e-9), 1).astype(int)
        shape = [int(n) + 1 for n in counts]
        steps = [h] * 3
        origin = (lo + hi) / 2 - counts * h / 2
    return {"origin": [float(x) for x in origin],
            "voxel_vectors": [[float(steps[a]) if b == a else 0.0 for b in range(3)]
                              for a in range(3)],
            "shape": shape}


def _grid_spec(mol, grid, spacing, shape, margin):
    """Normalize a grid argument to (origin, voxel vectors, shape) arrays."""
    np = _np()
    if grid is None:
        grid = molecular_grid(mol, spacing=spacing, shape=shape, margin=margin)
    elif shape is not None:
        raise ValueError("pass either grid or shape, not both")
    origin = np.asarray(grid["origin"], dtype=float)
    vectors = np.asarray(grid["voxel_vectors"], dtype=float)
    dims = tuple(int(n) for n in grid["shape"])
    if origin.shape != (3,) or vectors.shape != (3, 3) or len(dims) != 3 or min(dims) < 1:
        raise ValueError("grid needs origin [3], voxel_vectors [3][3] and shape of three counts >= 1")
    if not (np.isfinite(origin).all() and np.isfinite(vectors).all()):
        raise ValueError("grid origin and voxel_vectors must be finite")
    if abs(np.linalg.det(vectors)) <= 1e-12 * max(1.0, np.abs(vectors).max()) ** 3:
        raise ValueError("grid voxel_vectors are linearly dependent")
    return origin, vectors, dims


def _on_grid(mol, origin, vectors, dims, fields, n_fields, chunk):
    """Evaluate ``fields(ao) -> [n_fields, n_points]`` block by block.

    The points are those of spec § 5.2, ``origin + i v0 + j v1 + k v2``, in
    the row-major ``[i][j][k]`` order of the ``data`` member. Only one block
    of AO values is held at a time.
    """
    from pyscf.dft import numint
    np = _np()
    n_points = math.prod(dims)
    if chunk is None:
        chunk = max(256, _AO_BLOCK_VALUES // max(1, mol.nao_nr()))
    chunk = int(chunk)
    if chunk < 1:
        raise ValueError("chunk must be a positive number of points")
    out = np.empty((n_fields, n_points))
    for start in range(0, n_points, chunk):
        stop = min(n_points, start + chunk)
        index = np.stack(np.unravel_index(np.arange(start, stop), dims), axis=1)
        points = origin + index @ vectors
        out[:, start:stop] = fields(numint.eval_ao(mol, points, deriv=0))
    return out.reshape((n_fields,) + dims)


def _volume_links(writer, structure_ref, wavefunction_ref):
    sections = {s["id"]: s["kind"] for s in writer._sections}
    if sections.get(structure_ref) != "structure":
        raise ValueError(f"structure_ref {structure_ref!r} names no structure section "
                         f"added so far; add the SCF result first (add_scf)")
    if wavefunction_ref is not None and sections.get(wavefunction_ref) != "wavefunction.gto":
        raise ValueError(f"wavefunction_ref {wavefunction_ref!r} names no "
                         f"wavefunction.gto section added so far")


def _add_grid_volume(writer, kind, data, origin, vectors, *, section_id, label,
                     structure_ref, wavefunction_ref, dtype, mo_index=None, spin=None):
    if wavefunction_ref is None:
        mo_index = spin = None
    sid = writer.add_volume(kind, data, origin.tolist(), vectors.tolist(), label=label,
                            dtype=dtype, section_id=section_id, structure_ref=structure_ref,
                            wavefunction_ref=wavefunction_ref, mo_index=mo_index, spin=spin)
    writer.annotate_section(sid, structure_ref=structure_ref, group=wavefunction_ref)
    return sid


def add_density_grid(writer, mf, *, grid=None, spacing: float | None = DEFAULT_GRID_SPACING_BOHR,
                     shape=None, margin: float = DEFAULT_GRID_MARGIN_BOHR,
                     spin: bool | None = None, structure_ref: str = "structure",
                     wavefunction_ref: str | None = "wf", density_id: str = "density",
                     spin_id: str = "spin_density", dtype: str = "float32",
                     chunk: int | None = None) -> dict[str, str]:
    """Add the SCF electron density, and the spin density, on a regular grid.

    ``volume.density`` holds ρ(r) = Σ_μν P_μν χ_μ(r) χ_ν(r) of
    ``mf.make_rdm1()`` (both spins) in e/bohr³; ``volume.spin`` holds
    ρ_α − ρ_β in the same unit. ``spin=None`` writes the spin density for an
    open-shell reference (UHF, UKS, ROHF, ROKS) and not for a closed-shell
    RHF/RKS; ``True`` requires an open-shell reference and ``False`` omits it.

    The grid is ``grid`` (a mapping with the ``origin``, ``voxel_vectors``
    and ``shape`` of a spec § 5.2 grid payload, in bohr; any non-singular
    voxel vectors) or else ``molecular_grid(mf.mol, spacing=, shape=,
    margin=)``: by default 0.2 bohr steps and a 4 bohr margin around the
    nuclei. Values come from ``pyscf.dft.numint.eval_ao``/``eval_rho`` in
    blocks of ``chunk`` points (by default a 64 MiB AO block), stored as
    ``dtype`` (``"float32"``, the format's default, or ``"float64"``).

    The sections name ``structure_ref`` at section and grid level.
    ``wavefunction_ref`` (the ``wavefunction.gto`` id, which must exist
    already, or ``None``) is written as the section field
    ``wavefunction_ref`` (spec § 5.2, Volumes of a wavefunction) and, for
    grouping in viewers, as the section ``group`` (spec § 4.1). Returns the
    ids, keyed ``density`` and, when written, ``spin``.
    """
    _check_scf(mf)
    _refuse_periodic_grid(mf.mol)
    np = _np()
    from pyscf.dft import numint
    mol = mf.mol
    open_shell = _spin_kind(mf) != "restricted"
    if spin is None:
        spin = open_shell
    elif spin and not open_shell:
        raise ValueError("a closed-shell RHF/RKS reference has no spin density")
    _volume_links(writer, structure_ref, wavefunction_ref)
    origin, vectors, dims = _grid_spec(mol, grid, spacing, shape, margin)
    dm = np.asarray(mf.make_rdm1(), dtype=float)
    if dm.ndim == 3:
        matrices = [dm[0] + dm[1]] + ([dm[0] - dm[1]] if spin else [])
    else:
        matrices = [dm]

    def fields(ao):
        return [numint.eval_rho(mol, ao, p, xctype="LDA", hermi=1) for p in matrices]

    values = _on_grid(mol, origin, vectors, dims, fields, len(matrices), chunk)
    ids = {"density": _add_grid_volume(
        writer, "volume.density", values[0], origin, vectors, section_id=density_id,
        label="SCF electron density", structure_ref=structure_ref,
        wavefunction_ref=wavefunction_ref, dtype=dtype)}
    if spin:
        ids["spin"] = _add_grid_volume(
            writer, "volume.spin", values[1], origin, vectors, section_id=spin_id,
            label="SCF spin density (alpha - beta)", structure_ref=structure_ref,
            wavefunction_ref=wavefunction_ref, dtype=dtype)
    return ids


def _resolve_mo(index, occupations, energies):
    """A 0-based MO row from an integer or ``"homo"``/``"lumo"``."""
    np = _np()
    n = len(occupations)
    if isinstance(index, str):
        name = index.strip().lower()
        occupied = np.flatnonzero(np.asarray(occupations) > 0)
        virtual = np.flatnonzero(np.asarray(occupations) == 0)
        if name == "homo" and len(occupied):
            return int(max(occupied, key=lambda i: (energies[i], i)))
        if name == "lumo" and len(virtual):
            return int(min(virtual, key=lambda i: (energies[i], i)))
        raise ValueError(f"no {index!r} orbital in this orbital set")
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)):
        raise TypeError("MO indices are integers or 'homo'/'lumo'")
    if not 0 <= int(index) < n:
        raise ValueError(f"MO index {index} is outside 0..{n - 1}")
    return int(index)


def add_orbital_grid(writer, mf, mo_indices, *, spin: str = "alpha", grid=None,
                     spacing: float | None = DEFAULT_GRID_SPACING_BOHR, shape=None,
                     margin: float = DEFAULT_GRID_MARGIN_BOHR,
                     structure_ref: str = "structure", wavefunction_ref: str | None = "wf",
                     section_prefix: str = "mo", dtype: str = "float32",
                     chunk: int | None = None) -> list[str]:
    """Add selected molecular orbitals as ``volume.orbital`` sections.

    ``mo_indices`` lists 0-based orbital rows, the rows of the coefficient
    member of the ``wavefunction.gto`` section (spec § 5.3), or the names
    ``"homo"`` and ``"lumo"`` (the highest occupied and lowest unoccupied
    orbital by energy). For an unrestricted reference ``spin`` picks the
    ``"alpha"`` or ``"beta"`` set; restricted and restricted open-shell
    orbitals have one set and accept only ``"alpha"``. Each orbital is the
    real amplitude ψ_i(r) = Σ_μ C_μi χ_μ(r) in bohr^(−3/2) with PySCF's
    sign, which is also the archive's (``C_qvf = T C``).

    The grid, ``chunk``, ``dtype`` and links are those of
    ``add_density_grid``; every orbital in one call shares the grid and one
    AO evaluation. With a ``wavefunction_ref`` each section also names its
    row as ``mo_index`` and, for an unrestricted set, its ``spin`` (spec
    § 5.2), so a viewer finds the orbital's energy and occupation. Section ids are ``{section_prefix}_{row}``, with the spin
    inserted for an unrestricted set (``mo_beta_4``), and the label names
    the row, the spin and HOMO/LUMO. Returns the ids in ``mo_indices`` order.
    """
    _check_scf(mf)
    _refuse_periodic_grid(mf.mol)
    np = _np()
    mol = mf.mol
    if spin not in ("alpha", "beta"):
        raise ValueError("spin must be 'alpha' or 'beta'")
    unrestricted = _spin_kind(mf) == "unrestricted"
    if not unrestricted and spin != "alpha":
        raise ValueError("restricted orbitals have one set; use spin='alpha'")
    channel = ("alpha", "beta").index(spin)
    coeff = np.asarray(mf.mo_coeff[channel] if unrestricted else mf.mo_coeff, dtype=float)
    occ = np.asarray(mf.mo_occ[channel] if unrestricted else mf.mo_occ, dtype=float)
    energies = np.asarray(mf.mo_energy[channel] if unrestricted else mf.mo_energy, dtype=float)
    if isinstance(mo_indices, (int, str, np.integer)):
        mo_indices = [mo_indices]
    rows = [_resolve_mo(i, occ, energies) for i in mo_indices]
    if not rows:
        raise ValueError("no MO indices given")
    if len(set(rows)) != len(rows):
        raise ValueError("MO indices name the same orbital twice")
    named = {}
    for name in ("homo", "lumo"):
        try:
            named[_resolve_mo(name, occ, energies)] = name.upper()
        except ValueError:
            pass
    _volume_links(writer, structure_ref, wavefunction_ref)
    origin, vectors, dims = _grid_spec(mol, grid, spacing, shape, margin)
    selected = coeff[:, rows]
    values = _on_grid(mol, origin, vectors, dims, lambda ao: (ao @ selected).T, len(rows), chunk)
    ids = []
    for slab, row in zip(values, rows):
        sid = f"{section_prefix}_{spin}_{row}" if unrestricted else f"{section_prefix}_{row}"
        label = f"MO {row}" + (f" {spin}" if unrestricted else "")
        if row in named:
            label += f" ({named[row]})"
        ids.append(_add_grid_volume(
            writer, "volume.orbital", slab, origin, vectors, section_id=sid, label=label,
            structure_ref=structure_ref, wavefunction_ref=wavefunction_ref, dtype=dtype,
            mo_index=row, spin=spin if unrestricted else None))
    return ids
