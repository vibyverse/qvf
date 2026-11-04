"""qvf_writer — a standalone reference writer for the QVF format.

QVF (``.qvf``, *Quantum Visualization Format*) is a ZIP container with a
mandatory ``manifest.json`` plus typed JSON / binary members. This module is a
dependency-light reference implementation of the QVF **producer** contract: it
depends only on the Python standard library and NumPy, so any code can vendor it
without pulling in vibe-qc.

It is the executable companion to ``spec/qvf-format-spec.md``. Every method below
maps onto a canonical section kind from § 5 of that spec; the class enforces the
byte-level rules (little-endian binary, ``sha256`` over stored bytes,
``itemsize·∏(shape)`` sizing, Appendix-A.1 GTO coefficient normalization).

Typical use::

    from qvf_writer import QvfWriter
    w = QvfWriter(program="my-code", version="1.0", calculation="h2o/rhf")
    w.add_structure(atoms=[{"symbol": "O", "position": [0, 0, 0.117],
                            "atomic_number": 8}, ...])
    w.add_spectrum("spectra.ir", frequencies=[1595.0], intensities=[67.0])
    w.write("h2o.qvf")

License: Apache-2.0.
"""

from __future__ import annotations

import hashlib
import json
import math
import numbers
import re
import sys
import uuid as _uuid
import zipfile
from datetime import datetime as _datetime, timezone as _timezone
from typing import Any, Iterable, Mapping, Sequence

class _OptionalNumPy:
    """Load numeric convenience support only when it is used."""
    def __getattr__(self, name):
        import importlib
        try:
            module = importlib.import_module("numpy")
        except ImportError as exc:
            raise ImportError("This numeric convenience needs NumPy; install qvf-writer[numpy]") from exc
        globals()["np"] = module
        return getattr(module, name)


np = _OptionalNumPy()

__all__ = ["QvfWriter", "Lattice", "BOHR_PER_ANGSTROM", "primitive_norm",
           "QVF_VERSION", "TOOLKIT_VERSION"]

QVF_VERSION = 1
TOOLKIT_VERSION = "0.7.0rc1"

# Bohr per Ångström: the reciprocal of CODATA 2018's bohr radius,
# 1 bohr = 0.529177210903 Å (1 / 0.529177210903 = 1.88972612462577...), rounded
# to ten decimals. The QVF spec (§ 4.3, "Lattice conventions") uses this value.
BOHR_PER_ANGSTROM = 1.8897261246

_SCHEMA_URI = "https://vibe-qc.org/spec/qvf/1/manifest.schema.json"
# Byte order of every binary member, declared at the manifest root (spec § 2.2):
# "little" is the only v1 value.
_BYTE_ORDER = "little"

# NumPy dtype names permitted for binary members by the schema (NumpyDtype enum).
_ALLOWED_DTYPES = frozenset(
    {
        "int8", "int16", "int32", "int64",
        "uint8", "uint16", "uint32", "uint64",
        "float32", "float64",
    }
)

_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
_VENDOR_RE = re.compile(r"^x_[A-Za-z0-9_]+(\.[A-Za-z0-9_.-]+)?$")
# viewer_defaults keys that are not section ids (spec § 3.1, § 3.4).
_RESERVED_IDS = frozenset({"auto_open", "bookmarks", "camera", "crossfade"})
_DRIVE_RE = re.compile(r"^[A-Za-z]:")

# Root `metadata` block (spec § 3.5, $defs/Metadata): its keys, the forms of
# its time stamps and uuid, and the keys of an author object.
_METADATA_STRINGS = ("title", "description", "license", "doi")
_METADATA_LISTS = ("keywords", "derived_from", "tags")
_METADATA_KEYS = frozenset(_METADATA_STRINGS + _METADATA_LISTS + (
    "created_utc", "modified_utc", "uuid", "authors"))
_AUTHOR_KEYS = frozenset({"name", "orcid", "affiliation", "email"})
_UTC_RE = re.compile(r"^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])"
                     r"T([01][0-9]|2[0-3]):[0-5][0-9]:([0-5][0-9]|60)"
                     r"(\.[0-9]+)?Z$")
_UUID_RE = re.compile(r"^[0-9a-f]{8}-[0-9a-f]{4}-[1-8][0-9a-f]{3}-[89ab][0-9a-f]{3}"
                      r"-[0-9a-f]{12}$")
# Generic section references (spec § 4.1): field -> the kind it must name.
_SECTION_REF_KINDS = {"structure_ref": "structure", "run_ref": "run.record"}

# Excited states (spec § 5.12, $defs/ExcitedState and
# $defs/ExcitedStateContribution): the fields of a state and of a
# contribution, and the kinds a state's density_ref may name.
_STATE_FIELDS = frozenset({
    "index", "excitation_energy", "total_energy", "multiplicity", "symmetry",
    "oscillator_strength", "rotatory_strength_cgs", "transition_dipole",
    "transition_magnetic_dipole", "contributions", "nto_ref", "density_ref",
})
_CONTRIBUTION_FIELDS = frozenset({"from_mo", "to_mo", "spin", "weight"})
_STATE_DENSITY_KINDS = ("volume.transition_density", "volume.difference")
# Typed optional fields of the root `provenance` block that the writer checks
# (spec § 3.2, the schema's root `provenance` property). The block is open, so
# every other key passes through unchecked.
_PROVENANCE_INT_MINIMUM = {"n_threads": 1, "n_mpi_ranks": 1, "n_gpus": 0}
_PROVENANCE_UTC_KEYS = ("started_utc", "finished_utc")
_PROVENANCE_TYPED = frozenset(
    ("software", "git_commit") + _PROVENANCE_UTC_KEYS
    + tuple(_PROVENANCE_INT_MINIMUM))
_GIT_COMMIT_RE = re.compile(r"^[0-9a-fA-F]{7,64}$")
# Keys of one `provenance.software` item ($defs/ProvenanceSoftware).
_SOFTWARE_KEYS = frozenset({"name", "version", "url"})

# The geometry optimiser's coordinate systems (spec § 5.9, $defs/JobSpecGeometry).
_GEOMETRY_COORDINATES = ("cartesian", "internal", "redundant")

# Per-atom biomolecular fields a structure payload may carry (spec § 5.1,
# $defs/StructureAtom). Passed through verbatim so a PDB-derived structure
# keeps the identity a cartoon renderer needs; insertion_code, alt_loc,
# occupancy and hetero complete the PDB/mmCIF atom-site identity.
_ATOM_BIOMOLECULE_FIELDS = (
    "atom_name", "residue_name", "residue_seq", "chain_id", "b_factor",
    "insertion_code", "alt_loc", "occupancy", "hetero",
)

# The atom_properties member roles, in schema order ($defs/SectionAtomProperties,
# spec § 5.7). The schema, the spec, this writer, the C++ writer and the C ABI
# share this one list. The scalar roles are float64 [n_atoms]; the last two are
# the vector roles, float64 [n_atoms, 3] ($defs/AtomVectorBinary).
_ATOM_PROPERTY_ROLES = (
    "mulliken_charge", "loewdin_charge", "hirshfeld_charge", "iao_charge",
    "spin_population",
    "cm5_charge", "chelpg_charge", "resp_charge", "esp_charge", "nbo_charge",
    "mbis_charge", "bader_charge", "ddec_charge",
    "valence", "bonded_valence", "fukui_plus", "fukui_minus", "fukui_zero",
    "dual_descriptor", "atomic_energy_eh", "atomic_volume_bohr3",
    "nmr_shielding_iso_ppm",
    "magnetic_moment", "force_eh_bohr",
)
_ATOM_VECTOR_ROLES = ("magnetic_moment", "force_eh_bohr")

# A custom atom_properties role is ``custom.<name>`` (schema patternProperties).
_CUSTOM_NAME_RE = re.compile(r"^[A-Za-z0-9_]+$")
_CUSTOM_META_KEYS = ("units", "description", "kind")

# The named components of a root `energies` block ($defs/EnergyComponents,
# spec § 3.6). A component key outside this list must begin with `x_`.
_ENERGY_COMPONENTS = (
    "nuclear_repulsion", "electronic", "one_electron", "two_electron",
    "exchange_correlation", "exchange", "correlation", "dispersion",
    "solvation", "zpve", "thermal", "embedding",
)

# The property kinds (spec § 5.13): the keys each payload and block allows
# beyond keys that begin with `x_`, and the enumerations of the schema
# ($defs/ElectricPropertiesPayload and its siblings).
_ELECTRIC_KEYS = frozenset({"origin", "dipole", "quadrupole", "octupole",
                            "polarizability", "hyperpolarizability", "units"})
_ELECTRIC_BLOCKS = ("dipole", "quadrupole", "octupole", "polarizability",
                    "hyperpolarizability")
_MULTIPOLE_CONVENTIONS = ("buckingham", "cartesian")
_ELASTIC_MODULI = ("bulk_modulus_gpa", "shear_modulus_gpa", "youngs_modulus_gpa",
                   "poisson_ratio")
_ELASTIC_AVERAGES = ("voigt", "reuss", "hill")
# Largest disagreement the validators accept between two statements of one
# property: matrix symmetry (relative to the largest entry), static against
# electronic + ionic.
_PROPERTY_TOLERANCE = 1e-6

_ELEMENT_SYMBOLS = (
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)


def primitive_norm(alpha: float, l: int) -> float:
    """Norm of the *axial* Cartesian Gaussian primitive of a shell.

    ``N = (2α/π)^{3/4} · (4α)^{l/2} / √((2l−1)!!)`` where ``l`` is the total
    angular momentum. Engines that store contraction coefficients already
    multiplied by this factor (libint / libcint / PySCF, and GBW-style storage)
    must divide by ``primitive_norm`` before writing them into ``basis.json``.
    See spec Appendix A.

    Note the factor depends only on ``l``, not on the individual
    ``(l_x, l_y, l_z)``: it is the norm of the axial component ``(l, 0, 0)``
    and is applied uniformly across the shell. So for ``pure: false`` shells
    the primitives are *not* individually unit-normalized -- ``⟨xy|xy⟩ = 1/3``
    for a d shell, against ``⟨xx|xx⟩ = 1``. That is the convention the format
    defines; consumers must not add a per-component correction on top of it.
    See spec Appendix A.1.
    """
    radial = (2.0 * alpha / math.pi) ** 0.75
    angular = (4.0 * alpha) ** (l / 2.0)
    df = 1.0
    for k in range(1, 2 * l, 2):  # (2l-1)!!
        df *= k
    return radial * angular / math.sqrt(df)


def _basis_shells(shells: Sequence[Mapping[str, Any]], pure: bool,
                  coeffs_are_libint_normalized: bool
                  ) -> tuple[list[dict[str, Any]], int]:
    """The ``shells`` of a ``basis`` payload (spec § 5.3, Appendix A) and its
    AO count.

    Each input shell is ``{"center", "l", "exponents", "coefficients",
    "pure"?}``; a shell without ``pure`` takes ``pure``. With
    ``coeffs_are_libint_normalized`` every contraction coefficient is divided
    by the primitive norm ``N_i`` (Appendix A.1), so that it applies to
    ``N_i``-normalized primitives. Shared by ``add_wavefunction_gto`` and
    ``add_wavefunction_bloch``.
    """
    out_shells = []
    for sh in shells:
        l = int(sh["l"])
        exps = [float(x) for x in sh["exponents"]]
        coeffs = [float(c) for c in sh["coefficients"]]
        if coeffs_are_libint_normalized:
            coeffs = [c / primitive_norm(a, l) for a, c in zip(exps, coeffs)]
        out_shells.append({
            "center": int(sh["center"]), "l": l,
            "exponents": exps, "coefficients": coeffs,
            "pure": bool(sh.get("pure", pure)),
        })
    n_ao = 0
    for sh in out_shells:
        l = sh["l"]
        n_ao += (2 * l + 1) if sh["pure"] else ((l + 1) * (l + 2) // 2)
    return out_shells, n_ao


def _symbol(z: int) -> str:
    return _ELEMENT_SYMBOLS[z] if 0 <= z < len(_ELEMENT_SYMBOLS) else "X"


def _sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(obj: Any) -> bytes:
    """Serialize ``obj`` as canonical UTF-8 JSON bytes for a QVF member."""
    return json.dumps(
        obj, ensure_ascii=False, indent=2, allow_nan=False
    ).encode("utf-8")


def _le_bytes(data: Any, dtype: str) -> tuple[bytes, str, list[int]]:
    """Return (little-endian raw bytes, dtype-name, shape) for a binary member.

    The stored bytes are C-contiguous and little-endian regardless of host
    byte order, satisfying the spec § 2.2 contract.
    """
    if dtype not in _ALLOWED_DTYPES:
        raise ValueError(f"dtype {dtype!r} is not a permitted QVF dtype")
    le = np.dtype(dtype).newbyteorder("<")
    arr = np.ascontiguousarray(np.asarray(data).astype(le, copy=False))
    return arr.tobytes(), dtype, list(arr.shape)


def _real_array(data: Any, what: str) -> np.ndarray:
    """``data`` as an array of real numbers; anything else is a ValueError."""
    arr = np.asarray(data)
    if not np.issubdtype(arr.dtype, np.number) or np.iscomplexobj(arr):
        raise ValueError(f"{what} must hold real numbers, got dtype {arr.dtype}")
    return arr


def _energy(value: Any, units: str, what: str) -> dict[str, Any]:
    """A ``{value, units}`` energy: a number is wrapped with ``units``, a
    mapping with ``value`` and ``units`` is copied (spec § 3.6)."""
    if isinstance(value, Mapping):
        if set(value) != {"value", "units"}:
            raise ValueError(f"{what} must be a number or a {{value, units}} "
                             f"object, got keys {sorted(value)}")
        number, unit = value["value"], value["units"]
    else:
        number, unit = value, units
    # Python scalars must not load the optional NumPy convenience dependency.
    from numbers import Real
    if isinstance(number, bool) or not isinstance(number, Real):
        raise ValueError(f"{what} must be a number, got {number!r}")
    if not math.isfinite(float(number)):
        raise ValueError(f"{what} must be finite, got {number!r}")
    if not isinstance(unit, str) or not unit:
        raise ValueError(f"{what}: units must be a non-empty string")
    return {"value": float(number), "units": unit}


def _energy_components(components: Mapping[str, Any], units: str,
                       what: str) -> dict[str, Any]:
    """An ``energies`` components object (``$defs/EnergyComponents``)."""
    out: dict[str, Any] = {}
    for name, value in components.items():
        if name not in _ENERGY_COMPONENTS and not str(name).startswith("x_"):
            raise ValueError(
                f"{what}: unknown energy component {name!r}; use one of "
                f"{', '.join(_ENERGY_COMPONENTS)} or an x_ key (spec § 3.6)")
        out[str(name)] = _energy(value, units, f"{what}.{name}")
    return out


def _atom_property_array(role: str, data: Any, *, vector: bool) -> np.ndarray:
    """One ``atom_properties`` role as a checked real array: ``[n_atoms]`` for
    a scalar role, ``[n_atoms, 3]`` for a vector role (spec § 5.7)."""
    arr = _real_array(data, f"atom_properties {role}")
    want = "[n_atoms, 3]" if vector else "[n_atoms]"
    if arr.ndim != (2 if vector else 1) or (vector and arr.shape[1] != 3):
        raise ValueError(f"atom_properties {role} must be {want}, "
                         f"got shape {list(arr.shape)}")
    if arr.shape[0] == 0:
        raise ValueError(f"atom_properties {role} has no values")
    return arr


def _custom_role_name(name: Any) -> str:
    """The role ``custom.<name>`` for a ``custom=`` / ``meta=`` key; a leading
    ``custom.`` is accepted. The name must match ``[A-Za-z0-9_]+``."""
    if not isinstance(name, str):
        raise ValueError(f"atom_properties custom name must be a string, "
                         f"got {name!r}")
    bare = name[len("custom."):] if name.startswith("custom.") else name
    if not _CUSTOM_NAME_RE.match(bare):
        raise ValueError(f"atom_properties custom name {bare!r} must match "
                         f"[A-Za-z0-9_]+ (role custom.<name>, spec § 5.7)")
    return f"custom.{bare}"


def _custom_atom_roles(custom: Mapping[str, Any] | None,
                       meta: Mapping[str, Mapping[str, Any]] | None
                       ) -> dict[str, tuple[Any, Mapping[str, Any]]]:
    """Resolve ``custom=`` and ``meta=`` into ``{role: (array, entry)}``.

    Each custom role gets its entry from the ``(array, entry)`` tuple or, for a
    bare array, from ``meta``; one source only. ``meta`` may not describe a
    role that ``custom`` does not carry.
    """
    explicit: dict[str, Mapping[str, Any]] = {}
    for name, entry in (meta or {}).items():
        role = _custom_role_name(name)
        if role in explicit:
            raise ValueError(f"atom_properties meta names {role!r} twice")
        explicit[role] = entry
    out: dict[str, tuple[Any, Mapping[str, Any]]] = {}
    for name, value in (custom or {}).items():
        role = _custom_role_name(name)
        if role in out:
            raise ValueError(f"atom_properties custom names {role!r} twice")
        if isinstance(value, tuple):
            if len(value) != 2:
                raise ValueError(f"atom_properties {role} must be an "
                                 f"(array, meta entry) pair")
            if role in explicit:
                raise ValueError(f"atom_properties {role} has a meta entry "
                                 f"both in custom and in meta")
            out[role] = (value[0], value[1])
        elif role in explicit:
            out[role] = (value, explicit[role])
        else:
            raise ValueError(f"atom_properties {role} needs a meta entry "
                             f"(units, description): pass (array, entry) or "
                             f"meta={{...}}")
    unused = sorted(set(explicit) - set(out))
    if unused:
        raise ValueError(f"atom_properties meta describes {unused}, which "
                         f"custom does not carry")
    return out


def _custom_meta_entry(role: str, entry: Any, vector: bool) -> dict[str, Any]:
    """The ``meta`` entry of a custom role (``$defs/AtomPropertyMetaEntry``):
    ``units`` and ``description`` non-empty strings, ``kind`` derived from the
    array shape and, when given, in agreement with it; ``x_*`` keys pass."""
    if not isinstance(entry, Mapping):
        raise ValueError(f"atom_properties {role} meta entry must be a mapping, "
                         f"got {type(entry).__name__}")
    unknown = [k for k in entry
               if k not in _CUSTOM_META_KEYS
               and not (isinstance(k, str) and k.startswith("x_"))]
    if unknown:
        raise ValueError(f"atom_properties {role} meta entry has unknown "
                         f"key(s) {sorted(map(str, unknown))}; allowed: "
                         f"{list(_CUSTOM_META_KEYS)} and x_* keys")
    for key in ("units", "description"):
        value = entry.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError(f"atom_properties {role} meta entry needs a "
                             f"non-empty string {key!r}")
    kind = "vector" if vector else "scalar"
    if entry.get("kind", kind) != kind:
        raise ValueError(f"atom_properties {role} meta kind "
                         f"{entry['kind']!r} disagrees with the array, which "
                         f"is a {kind}")
    out = {"units": entry["units"], "description": entry["description"],
           "kind": kind}
    out.update({k: v for k, v in entry.items() if k not in out})
    return out


def _phonon_eigenvectors(eigenvectors: Any, frequencies: Any) -> np.ndarray:
    """Check ``phonon_bands`` eigenvectors and return them as float64.

    Real input must be ``[n_q, n_modes, n_atoms, 3]``; complex input is either
    a complex array of that shape, which is split into ``[..., 2]`` = ``[re,
    im]``, or already the real ``[n_q, n_modes, n_atoms, 3, 2]`` form (spec
    § 5.6). The leading ``[n_q, n_modes]`` must match ``frequencies``.
    """
    arr = np.asarray(eigenvectors)
    if np.iscomplexobj(arr):
        arr = np.stack([arr.real, arr.imag], axis=-1)
    elif not np.issubdtype(arr.dtype, np.number):
        raise ValueError("phonon_bands eigenvectors must hold numbers, "
                         f"got dtype {arr.dtype}")
    shape = arr.shape
    if not ((arr.ndim == 4 and shape[3] == 3)
            or (arr.ndim == 5 and shape[3:] == (3, 2))):
        raise ValueError(
            "phonon_bands eigenvectors must be [n_q, n_modes, n_atoms, 3] "
            "(real) or [n_q, n_modes, n_atoms, 3, 2] (complex, last axis "
            f"[re, im]); got shape {list(shape)}")
    freq_shape = tuple(np.shape(frequencies))
    if len(freq_shape) != 2 or shape[:2] != freq_shape:
        raise ValueError(
            f"phonon_bands eigenvectors [n_q, n_modes] = {list(shape[:2])} "
            f"must match frequencies shape {list(freq_shape)}")
    return arr.astype(np.float64)


# --- wavefunction.gto (spec § 5.3) ---------------------------------------------

_WF_SPINS = frozenset({"restricted", "unrestricted", "restricted_open",
                       "generalized"})
_ECP_KEYS = frozenset({"center", "n_core_electrons", "label", "projectors"})
_ECP_PROJECTOR_KEYS = frozenset({"l", "exponents", "coefficients", "powers",
                                 "local"})


def _is_int(x: Any) -> bool:
    return isinstance(x, (int, np.integer)) and not isinstance(x, bool)


def _wf_coefficients(data: Any, what: str, complex_coefficients: bool,
                     ao_extent: int) -> np.ndarray:
    """Check one MO coefficient matrix and return it as float64 (spec § 5.3).

    Real coefficients are ``[n_mo, ao_extent]``. With ``complex_coefficients``
    a complex ``[n_mo, ao_extent]`` array is split into ``[..., 2]`` = ``[re,
    im]``, and a real ``[n_mo, ao_extent, 2]`` array is taken as that form.
    """
    arr = np.asarray(data)
    if np.iscomplexobj(arr):
        if not complex_coefficients:
            raise ValueError(f"{what} is complex; pass complex_coefficients=True")
        if arr.ndim != 2:
            raise ValueError(f"{what} must be [n_mo, n_ao], got shape "
                             f"{list(arr.shape)}")
        arr = np.stack([arr.real, arr.imag], axis=-1)
    elif not np.issubdtype(arr.dtype, np.number):
        raise ValueError(f"{what} must hold numbers, got dtype {arr.dtype}")
    want = 3 if complex_coefficients else 2
    if arr.ndim != want or (complex_coefficients and arr.shape[2] != 2):
        form = ("[n_mo, n_ao, 2] (or a complex [n_mo, n_ao])"
                if complex_coefficients else "[n_mo, n_ao]")
        raise ValueError(f"{what} must be {form}, got shape {list(arr.shape)}")
    if arr.shape[1] != ao_extent:
        raise ValueError(f"{what} has {arr.shape[1]} AO columns but the basis "
                         f"gives {ao_extent}")
    return arr.astype(np.float64)


def _wf_ao_matrix(data: Any, what: str, n_ao: int) -> np.ndarray:
    """An ``overlap`` or density matrix: real float64 ``[n_ao, n_ao]``."""
    arr = _real_array(data, what)
    if arr.shape != (n_ao, n_ao):
        raise ValueError(f"{what} must be [n_ao, n_ao] = [{n_ao}, {n_ao}], got "
                         f"shape {list(arr.shape)}")
    return arr.astype(np.float64)


def _wf_numbers(values: Any, what: str, *, positive: bool = False) -> list[float]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Iterable):
        raise ValueError(f"{what} must be a list of numbers")
    out = []
    for x in values:
        if isinstance(x, bool) or not isinstance(x, (int, float, np.number)):
            raise ValueError(f"{what} must be a list of numbers, got {x!r}")
        if not math.isfinite(float(x)) or (positive and not float(x) > 0):
            raise ValueError(f"{what} must hold finite numbers"
                             + (" above 0" if positive else "") + f", got {x!r}")
        out.append(float(x))
    if not out:
        raise ValueError(f"{what} must not be empty")
    return out


def _wf_ecp(ecp: Any) -> list[dict[str, Any]]:
    """Check the ``ecp`` entries of a ``basis`` payload (spec § 5.3) and return
    them in their on-disk form: at most one per atom; parallel projector
    arrays; at most one local part per centre, with the centre's highest
    ``l``."""
    if isinstance(ecp, Mapping) or not isinstance(ecp, Sequence):
        raise ValueError("ecp must be a list of {center, n_core_electrons, "
                         "label?, projectors?} entries")
    out, seen = [], set()
    for k, entry in enumerate(ecp):
        what = f"ecp[{k}]"
        if not isinstance(entry, Mapping):
            raise ValueError(f"{what} must be a mapping")
        unknown = {key for key in entry
                   if key not in _ECP_KEYS and not str(key).startswith("x_")}
        if unknown:
            raise ValueError(f"{what}: unknown keys {sorted(unknown)}")
        center, n_core = entry.get("center"), entry.get("n_core_electrons")
        if not _is_int(center) or center < 0:
            raise ValueError(f"{what}: center must be an integer >= 0")
        if not _is_int(n_core) or n_core < 0:
            raise ValueError(f"{what}: n_core_electrons must be an integer >= 0")
        if int(center) in seen:
            raise ValueError(f"{what}: atom {int(center)} already has an ecp "
                             f"entry; there is at most one per atom")
        seen.add(int(center))
        item: dict[str, Any] = {"center": int(center),
                                "n_core_electrons": int(n_core)}
        if entry.get("label") is not None:
            if not isinstance(entry["label"], str):
                raise ValueError(f"{what}: label must be a string")
            item["label"] = entry["label"]
        if entry.get("projectors") is not None:
            projectors, local_l, semi_l = [], [], []
            for j, proj in enumerate(entry["projectors"]):
                pwhat = f"{what}.projectors[{j}]"
                if not isinstance(proj, Mapping):
                    raise ValueError(f"{pwhat} must be a mapping")
                unknown = {key for key in proj if key not in _ECP_PROJECTOR_KEYS
                           and not str(key).startswith("x_")}
                if unknown:
                    raise ValueError(f"{pwhat}: unknown keys {sorted(unknown)}")
                l = proj.get("l")
                if not _is_int(l) or l < 0:
                    raise ValueError(f"{pwhat}: l must be an integer >= 0")
                exps = _wf_numbers(proj.get("exponents"), f"{pwhat}.exponents",
                                   positive=True)
                coefs = _wf_numbers(proj.get("coefficients"),
                                    f"{pwhat}.coefficients")
                powers = proj.get("powers")
                if (isinstance(powers, (str, bytes))
                        or not isinstance(powers, Iterable)
                        or not all(_is_int(n) for n in powers)):
                    raise ValueError(f"{pwhat}.powers must be a list of integers")
                powers = [int(n) for n in powers]
                if not len(exps) == len(coefs) == len(powers):
                    raise ValueError(f"{pwhat}: exponents, coefficients and "
                                     f"powers must have the same length")
                local = proj.get("local", False)
                if not isinstance(local, bool):
                    raise ValueError(f"{pwhat}: local must be a boolean")
                (local_l if local else semi_l).append(int(l))
                p: dict[str, Any] = {"l": int(l), "exponents": exps,
                                     "coefficients": coefs, "powers": powers}
                if local:
                    p["local"] = True
                p.update({key: value for key, value in proj.items()
                          if str(key).startswith("x_")})
                projectors.append(p)
            if len(local_l) > 1:
                raise ValueError(f"{what}: more than one local projector")
            if local_l and semi_l and local_l[0] <= max(semi_l):
                raise ValueError(f"{what}: the local projector's l "
                                 f"({local_l[0]}) must exceed the l of every "
                                 f"semi-local channel")
            item["projectors"] = projectors
        item.update({key: value for key, value in entry.items()
                     if str(key).startswith("x_")})
        out.append(item)
    return out


def _wf_labels(labels: Any, what: str, n: int) -> list[str]:
    if (isinstance(labels, (str, bytes)) or not isinstance(labels, Sequence)
            or not all(isinstance(s, str) for s in labels)):
        raise ValueError(f"{what} must be a list of strings")
    if len(labels) != n:
        raise ValueError(f"{what} has {len(labels)} entries but there are {n} "
                         f"MOs; there is one label per MO")
    return list(labels)


def _finite_number(value: Any, what: str) -> float:
    """``value`` as a finite float; a bool or anything else is a ValueError."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(
            value, (int, float, np.integer, np.floating)):
        raise ValueError(f"{what} must be a number, got {value!r}")
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{what} must be finite, got {value!r}")
    return out


def _real_matrix(value: Any, rows: int, cols: int, what: str) -> list[list[float]]:
    """``value`` as a ``rows`` x ``cols`` nested list of finite floats."""
    arr = _real_array(value, what)
    if arr.shape != (rows, cols):
        raise ValueError(f"{what} must be {rows} x {cols}, got shape "
                         f"{list(arr.shape)}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{what} must be finite")
    return arr.astype(np.float64).tolist()


def _real_vector(value: Any, n: int, what: str) -> list[float]:
    """``value`` as a list of ``n`` finite floats."""
    arr = _real_array(value, what)
    if arr.shape != (n,):
        raise ValueError(f"{what} must have {n} numbers, got shape "
                         f"{list(arr.shape)}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{what} must be finite")
    return arr.astype(np.float64).tolist()


def _real_tensor3(value: Any, what: str) -> list[float]:
    """A rank-3 tensor of shape ``[3, 3, 3]``, given as that array or as the 27
    numbers in row-major order, as the flat list the payload stores."""
    arr = _real_array(value, what)
    if arr.shape not in ((3, 3, 3), (27,)):
        raise ValueError(f"{what} must be [3, 3, 3] or 27 numbers in "
                         f"row-major order, got shape {list(arr.shape)}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{what} must be finite")
    return arr.astype(np.float64).reshape(27).tolist()


def _require_symmetric(matrix: list[list[float]], what: str) -> None:
    """Raise unless ``matrix`` is symmetric to within 1e-6 of its largest
    absolute entry, the rule the validators apply (spec § 5.13)."""
    scale = max(abs(x) for row in matrix for x in row)
    n = len(matrix)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(matrix[i][j] - matrix[j][i]) > _PROPERTY_TOLERANCE * scale:
                raise ValueError(
                    f"{what} must be symmetric: entry ({i}, {j}) is "
                    f"{matrix[i][j]!r} but entry ({j}, {i}) is "
                    f"{matrix[j][i]!r} (spec § 5.13)")


def _check_keys(obj: Mapping[str, Any], allowed: Iterable[str], what: str) -> None:
    """Raise unless ``obj`` is a mapping whose keys are in ``allowed`` or begin
    with ``x_`` (a strict payload object, spec § 5)."""
    if not isinstance(obj, Mapping):
        raise ValueError(f"{what} must be a mapping, got {type(obj).__name__}")
    allowed = set(allowed)
    unknown = sorted(str(k) for k in obj
                     if k not in allowed and not str(k).startswith("x_"))
    if unknown:
        raise ValueError(f"{what}: unknown keys {unknown}; allowed are "
                         f"{sorted(allowed)} and keys beginning x_")


def _choice(value: Any, options: Sequence[str], what: str) -> str:
    if not isinstance(value, str) or value not in options:
        raise ValueError(f"{what} must be one of {', '.join(options)}, "
                         f"got {value!r}")
    return value


class Lattice:
    """A periodic cell that converts between QVF's two lattice conventions.

    QVF v1 carries a lattice in two members that differ in layout and unit
    (spec § 4.3, "Lattice conventions"):

    * ``structure.lattice_vectors``: three **row** vectors in **Ångström**.
    * ``reaction.path.lattice``: a 3x3 matrix whose **columns** are a, b, c, in
      **bohr** (or one such matrix per frame).

    A ``Lattice`` stores the rows in Ångström internally and converts at the
    edges, so the conversion is written once. Build it from either convention
    and read it back in either::

        cell = Lattice.from_rows_angstrom([[2.46, 0, 0], [1.23, 2.1304, 0],
                                           [0, 0, 15.875]])
        cell.columns_bohr()        # the reaction.path.lattice matrix
        w.add_structure(atoms, pbc=[True, True, False], lattice_vectors=cell)
        w.add_reaction_path(meta, coords, lattice=cell)

    The columns of ``columns_bohr()`` are the rows of ``rows_angstrom()`` times
    ``BOHR_PER_ANGSTROM``. Instances are immutable by convention (the stored
    array is read-only). Only the shape (3x3) and finiteness are validated; a lattice with a zero-length or non-physical row
    (an aperiodic axis, spec § 5.1) is accepted.
    """

    __slots__ = ("_rows",)

    def __init__(self, rows_angstrom: Any) -> None:
        """Use :meth:`from_rows_angstrom` or :meth:`from_columns_bohr`.

        ``rows_angstrom`` is a 3x3 array-like of row vectors in Ångström.
        """
        self._rows = self._checked(rows_angstrom, "rows_angstrom")

    @staticmethod
    def _checked(matrix: Any, what: str) -> np.ndarray:
        try:
            arr = np.array(matrix, dtype=np.float64)  # copy
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Lattice {what} must be a 3x3 array of numbers: "
                             f"{exc}") from exc
        if arr.shape != (3, 3):
            raise ValueError(f"Lattice {what} must have shape (3, 3); "
                             f"got {tuple(arr.shape)}")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"Lattice {what} must be finite (no NaN or inf)")
        arr.setflags(write=False)
        return arr

    @classmethod
    def from_rows_angstrom(cls, rows: Any) -> "Lattice":
        """From ``structure.lattice_vectors``: row ``i`` is axis ``i``, in Å."""
        return cls(rows)

    @classmethod
    def from_columns_bohr(cls, matrix: Any) -> "Lattice":
        """From one ``reaction.path.lattice`` matrix: columns a, b, c, in bohr."""
        arr = cls._checked(matrix, "columns_bohr")
        return cls(arr.T / BOHR_PER_ANGSTROM)

    def rows_angstrom(self) -> list[list[float]]:
        """The ``structure.lattice_vectors`` form: 3x3 nested list, rows, Å."""
        return self._rows.tolist()

    def columns_bohr(self) -> list[list[float]]:
        """The ``reaction.path.lattice`` form: 3x3 nested list, columns, bohr."""
        return (self._rows.T * BOHR_PER_ANGSTROM).tolist()

    def volume_angstrom3(self) -> float:
        """Cell volume in Å^3: ``|a . (b x c)|`` (zero for a degenerate cell)."""
        a, b, c = self._rows
        return float(abs(np.dot(a, np.cross(b, c))))

    def volume_bohr3(self) -> float:
        """Cell volume in bohr^3."""
        return self.volume_angstrom3() * BOHR_PER_ANGSTROM ** 3

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Lattice):
            return NotImplemented
        return bool(np.array_equal(self._rows, other._rows))

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"Lattice.from_rows_angstrom({self.rows_angstrom()!r})"


def _lattice_as_rows(lattice_vectors: Any) -> list[list[float]]:
    """``add_structure`` ``lattice_vectors`` as a nested list of floats (Å)."""
    if isinstance(lattice_vectors, Lattice):
        return lattice_vectors.rows_angstrom()
    return [[float(x) for x in row] for row in lattice_vectors]


def _reaction_path_lattice(lattice: Any) -> np.ndarray:
    """``add_reaction_path`` ``lattice`` as float64 columns-in-bohr data.

    A ``Lattice`` becomes ``[3, 3]``, a list/tuple of ``Lattice`` becomes
    ``[n_frames, 3, 3]``; anything else is taken to be columns in bohr already.
    """
    if isinstance(lattice, Lattice):
        return np.asarray(lattice.columns_bohr(), dtype=np.float64)
    if isinstance(lattice, (list, tuple)) and any(
            isinstance(x, Lattice) for x in lattice):
        if not all(isinstance(x, Lattice) for x in lattice):
            raise TypeError(
                "reaction.path lattice: a per-frame sequence must contain only "
                "Lattice instances, not a mix of Lattice and arrays")
        return np.asarray([x.columns_bohr() for x in lattice], dtype=np.float64)
    return np.asarray(lattice, dtype=np.float64)


def _trajectory_lattice(lattice: Any) -> np.ndarray:
    """``add_trajectory`` ``lattice`` as float64 rows-in-Ångström data, the
    structure convention (spec § 4.3, "Lattice conventions").

    A ``Lattice`` becomes ``[3, 3]``, a list/tuple of ``Lattice`` becomes
    ``[n_frames, 3, 3]``; anything else is taken to be rows in Å already.
    """
    if isinstance(lattice, Lattice):
        return np.asarray(lattice.rows_angstrom(), dtype=np.float64)
    if isinstance(lattice, (list, tuple)) and any(
            isinstance(x, Lattice) for x in lattice):
        if not all(isinstance(x, Lattice) for x in lattice):
            raise TypeError(
                "trajectory lattice: a per-frame sequence must contain only "
                "Lattice instances, not a mix of Lattice and arrays")
        return np.asarray([x.rows_angstrom() for x in lattice], dtype=np.float64)
    return _real_array(lattice, "trajectory lattice").astype(np.float64)


def _check_lattice_frames(lat: np.ndarray, n_frames: int, kind: str) -> None:
    """A trajectory-family ``lattice`` is ``[3, 3]`` or ``[n_frames, 3, 3]``."""
    if lat.shape != (3, 3) and not (lat.ndim == 3 and lat.shape[1:] == (3, 3)):
        raise ValueError(f"{kind} lattice must be [3, 3] or [n_frames, 3, 3]; "
                         f"got shape {tuple(lat.shape)}")
    if lat.ndim == 3 and lat.shape[0] != n_frames:
        raise ValueError(f"{kind} per-frame lattice has {lat.shape[0]} frames "
                         f"but coords has {n_frames}")


# The optional per-frame members of the trajectory family (spec § 5.5), in the
# order the writer emits them: True for [n_frames, n_atoms, 3] (velocities in
# Å/fs, the energy gradient in hartree/bohr), False for [n_frames] (times in fs,
# energies in hartree, cell volumes in Å³, temperatures in kelvin).
_TRAJECTORY_FRAME_ROLES = {"velocities": True, "gradient": True, "times": False,
                           "energies": False, "cell_volumes": False,
                           "temperatures": False}
# Largest difference, in hartree, between a binary `energies` member and the
# metadata `energies` array of the same section (spec § 5.5).
_TRAJECTORY_ENERGY_TOLERANCE_EH = 1e-10


def _trajectory_frame_arrays(kind: str, coords: np.ndarray,
                             arrays: Mapping[str, Any],
                             meta: Mapping[str, Any]) -> list[tuple[str, np.ndarray]]:
    """The given per-frame members of a trajectory or reaction path as checked
    float64 arrays, ``(role, array)`` in emission order (spec § 5.5).

    ``velocities`` and ``gradient`` must be ``[n_frames, n_atoms, 3]`` and the
    others ``[n_frames]``, with the extents of ``coords``; a binary
    ``energies`` must agree with ``meta["energies"]`` within 1e-10 hartree
    when both are given. Anything else raises ``ValueError``.
    """
    n_frames, n_atoms = coords.shape[0], coords.shape[1]
    out: list[tuple[str, np.ndarray]] = []
    for role, vector in _TRAJECTORY_FRAME_ROLES.items():
        data = arrays.get(role)
        if data is None:
            continue
        arr = _real_array(data, f"{kind} {role}").astype(np.float64)
        want = (n_frames, n_atoms, 3) if vector else (n_frames,)
        if arr.shape != want:
            form = "[n_frames, n_atoms, 3]" if vector else "[n_frames]"
            raise ValueError(f"{kind} {role} must be {form} = {list(want)} to "
                             f"match coords, got shape {list(arr.shape)}")
        out.append((role, arr))
    binary = dict(out).get("energies")
    listed = meta.get("energies")
    if binary is not None and listed is not None:
        try:
            ref = np.asarray(listed, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{kind} metadata energies must be numbers: "
                             f"{exc}") from exc
        if ref.shape != binary.shape or not np.all(
                np.abs(ref - binary) <= _TRAJECTORY_ENERGY_TOLERANCE_EH):
            raise ValueError(f"{kind} energies and metadata energies must "
                             f"agree within 1e-10 hartree (spec § 5.5)")
    return out


def _with_units(payload: dict[str, Any], units: Mapping[str, str] | None,
                what: str, defname: str) -> dict[str, Any]:
    """``payload`` with its ``units`` declaration checked (spec § 4.3,
    "Declared units"), after adding ``units`` as that declaration.

    Every value must be a unit name of the vocabulary, and a key whose unit
    the schema fixes for ``defname`` must name that unit; anything else
    raises ``ValueError``. Declaring never converts a value.
    """
    if units is not None:
        if "units" in payload:
            raise ValueError(f"{what}: units given both as an argument and in "
                             f"the payload")
        payload["units"] = units
    if "units" in payload:
        if not isinstance(payload["units"], Mapping):
            raise ValueError(f"{what}: units must be a mapping of field or "
                             f"member names to unit names")
        payload["units"] = _checked_units(payload["units"], what, defname)
    return payload


def _checked_units(units: Mapping[str, Any], what: str,
                   defname: str | None = None) -> dict[str, str]:
    # Deferred import: qvf_units imports this module for BOHR_PER_ANGSTROM.
    from .units import _check_declared
    return _check_declared(units, what, defname)


def _checked_unit_name(name: Any, what: str) -> str:
    from .units import is_unit
    if not is_unit(name):
        raise ValueError(f"{what}: {name!r} is not a QVF unit name "
                         f"(spec § 4.3)")
    return name


def _positive_int(value: Any, what: str, minimum: int = 1) -> int:
    """``value`` as an int of at least ``minimum``; a bool, a float or a
    smaller number raises ``ValueError``."""
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{what} must be an integer, got {value!r}")
    if int(value) < minimum:
        raise ValueError(f"{what} must be >= {minimum}, got {value!r}")
    return int(value)


def _real_number(value: Any, what: str, minimum: float | None = None) -> float:
    """``value`` as a finite float (at least ``minimum`` when given)."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"{what} must be a number, got {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{what} must be finite, got {value!r}")
    if minimum is not None and number < minimum:
        raise ValueError(f"{what} must be >= {minimum}, got {value!r}")
    return number


def _vec3(value: Any, what: str) -> list[float]:
    """Three finite numbers as a list of floats."""
    if isinstance(value, (str, bytes)) or not hasattr(value, "__len__") \
            or len(value) != 3:
        raise ValueError(f"{what} must be three numbers, got {value!r}")
    return [_real_number(x, what) for x in value]


def _unit_energy(value: Any, units: str, what: str) -> dict[str, Any]:
    """A ``{value, units}`` energy whose unit is a unit name (spec § 4.3); a
    plain number is written in ``units``."""
    energy = _energy(value, units, what)
    _checked_unit_name(energy["units"], f"{what} units")
    return energy


def _excited_state(state: Any, k: int, rows: dict[str, int] | None,
                   check_ref: Any) -> dict[str, Any]:
    """One state of an ``excited_states`` payload, checked and normalized
    (spec § 5.12). ``rows`` maps a coefficient channel to its row count when
    the reference wavefunction is known; ``check_ref(field, ref, kinds)``
    raises unless ``ref`` names an added section of one of ``kinds``."""
    what = f"excited_states: states[{k}]"
    if not isinstance(state, Mapping):
        raise ValueError(f"{what} must be a mapping")
    unknown = sorted(key for key in state if key not in _STATE_FIELDS
                     and not str(key).startswith("x_"))
    if unknown:
        raise ValueError(f"{what}: unknown field(s) {unknown}; a state carries "
                         f"{', '.join(sorted(_STATE_FIELDS))} or x_ keys "
                         f"(spec § 5.12)")
    if "index" not in state or "excitation_energy" not in state:
        raise ValueError(f"{what} needs index and excitation_energy "
                         f"(spec § 5.12)")
    out: dict[str, Any] = {}
    for key, value in state.items():
        if key == "index":
            out[key] = _positive_int(value, f"{what}.index")
        elif key == "excitation_energy":
            out[key] = _unit_energy(value, "ev", f"{what}.excitation_energy")
        elif key == "total_energy":
            out[key] = _unit_energy(value, "hartree", f"{what}.total_energy")
        elif key == "multiplicity":
            out[key] = _positive_int(value, f"{what}.multiplicity")
        elif key == "symmetry":
            if not isinstance(value, str):
                raise ValueError(f"{what}.symmetry must be a string")
            out[key] = value
        elif key == "oscillator_strength":
            out[key] = _real_number(value, f"{what}.oscillator_strength", 0.0)
        elif key == "rotatory_strength_cgs":
            out[key] = _real_number(value, f"{what}.rotatory_strength_cgs")
        elif key in ("transition_dipole", "transition_magnetic_dipole"):
            out[key] = _vec3(value, f"{what}.{key}")
        elif key == "contributions":
            if isinstance(value, (str, bytes, Mapping)) or not isinstance(
                    value, Iterable):
                raise ValueError(f"{what}.contributions must be a sequence")
            out[key] = [_contribution(c, f"{what}.contributions[{j}]", rows)
                        for j, c in enumerate(value)]
        elif key == "nto_ref":
            check_ref(f"{what}.nto_ref", value, ("wavefunction.gto",))
            out[key] = value
        elif key == "density_ref":
            check_ref(f"{what}.density_ref", value, _STATE_DENSITY_KINDS)
            out[key] = value
        else:
            out[key] = value                                  # an x_ key
    return out


def _contribution(contrib: Any, what: str,
                  rows: dict[str, int] | None) -> dict[str, Any]:
    """One orbital excitation of a state (``$defs/ExcitedStateContribution``);
    with the reference wavefunction's ``rows``, its orbitals must exist."""
    if not isinstance(contrib, Mapping):
        raise ValueError(f"{what} must be a mapping")
    unknown = sorted(key for key in contrib if key not in _CONTRIBUTION_FIELDS
                     and not str(key).startswith("x_"))
    if unknown or not {"from_mo", "to_mo", "weight"} <= set(contrib):
        raise ValueError(f"{what} carries from_mo, to_mo and weight, optionally "
                         f"spin and x_ keys (spec § 5.12); got {sorted(contrib)}")
    out = dict(contrib)
    spin = contrib.get("spin")
    if "spin" in contrib and spin not in ("alpha", "beta"):
        raise ValueError(f"{what}.spin must be 'alpha' or 'beta', got {spin!r}")
    out["weight"] = _real_number(contrib["weight"], f"{what}.weight")
    if rows is not None:
        if "" in rows:
            n = rows[""]
        elif spin in ("alpha", "beta"):
            n = rows.get(spin)
        else:
            known = [rows[c] for c in ("alpha", "beta") if c in rows]
            n = min(known) if known else None
    else:
        n = None
    for key in ("from_mo", "to_mo"):
        out[key] = _positive_int(contrib[key], f"{what}.{key}", 0)
        if n is not None and out[key] >= n:
            raise ValueError(f"{what}.{key} {out[key]} is not one of the {n} "
                             f"orbital rows of the reference wavefunction "
                             f"(spec § 5.12)")
    return out

def _utc_timestamp(value: Any, key: str, *, block: str = "metadata",
                   rule: str = "§ 3.5") -> str:
    """A time stamp of ``block`` (spec ``rule``: ``metadata``, § 3.5, or
    ``provenance``, § 3.2) from a string or an aware datetime."""
    if isinstance(value, _datetime):
        if value.tzinfo is None or value.utcoffset() is None:
            raise ValueError(f"{block} {key}: a naive datetime has no time "
                             f"zone; pass an aware one or a string ending 'Z'")
        value = value.astimezone(_timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    if not isinstance(value, str) or not _UTC_RE.fullmatch(value):
        raise ValueError(f"{block} {key} must be an ISO 8601 UTC time stamp "
                         f"'YYYY-MM-DDThh:mm:ssZ' (spec {rule}), got {value!r}")
    return value


def _checked_text(value: Any, what: str) -> str:
    """``value`` as a non-empty string, else ``ValueError`` naming ``what``."""
    if not isinstance(value, str) or not value:
        raise ValueError(f"{what} must be a non-empty string, got {value!r}")
    return value


def _checked_int(value: Any, what: str, minimum: int) -> int:
    """``value`` as an integer of at least ``minimum`` (a bool is not one)."""
    if isinstance(value, bool) or not isinstance(value, numbers.Integral):
        raise ValueError(f"{what} must be an integer >= {minimum}, got {value!r}")
    if value < minimum:
        raise ValueError(f"{what} must be an integer >= {minimum}, got {value!r}")
    return int(value)


def _checked_number(value: Any, what: str, *, positive: bool = False) -> float:
    """``value`` as a finite number, greater than zero when ``positive``."""
    if (isinstance(value, bool) or not isinstance(value, numbers.Real)
            or not math.isfinite(value)):
        raise ValueError(f"{what} must be a finite number, got {value!r}")
    if positive and not value > 0:
        raise ValueError(f"{what} must be a number > 0, got {value!r}")
    return float(value)


def _software_entry(entry: Any, what: str) -> dict[str, Any]:
    """One item of ``provenance.software`` (spec § 3.2, ``$defs/ProvenanceSoftware``)."""
    if not isinstance(entry, Mapping):
        raise ValueError(f"{what} must be a mapping with a 'name', got {entry!r}")
    extra = {str(k) for k in entry
             if not str(k).startswith("x_")} - _SOFTWARE_KEYS
    if extra:
        raise ValueError(f"{what}: keys {sorted(extra)} are not one of "
                         f"{sorted(_SOFTWARE_KEYS)}; vendor keys must begin "
                         f"with 'x_' (spec § 3.2)")
    _checked_text(entry.get("name"), f"{what} name")
    for key in ("version", "url"):
        if key in entry and not isinstance(entry[key], str):
            raise ValueError(f"{what} {key} must be a string, got {entry[key]!r}")
    return dict(entry)


def _provenance_field(key: str, value: Any) -> Any:
    """Check and normalise one typed field of the root ``provenance`` block
    (spec § 3.2): ``software``, ``git_commit``, ``n_threads``, ``n_mpi_ranks``,
    ``n_gpus``, ``started_utc`` and ``finished_utc``."""
    what = f"provenance {key}"
    if key == "software":
        if isinstance(value, (str, bytes, Mapping)) or not isinstance(
                value, Iterable):
            raise ValueError("provenance software must be a list of mappings "
                             "{name, version?, url?} (spec § 3.2), got "
                             f"{value!r}")
        return [_software_entry(item, f"provenance software[{i}]")
                for i, item in enumerate(value)]
    if key == "git_commit":
        if not isinstance(value, str) or not _GIT_COMMIT_RE.fullmatch(value):
            raise ValueError("provenance git_commit must be a git commit hash "
                             "of 7 to 64 hexadecimal characters (spec § 3.2), "
                             f"got {value!r}")
        return value
    if key in _PROVENANCE_UTC_KEYS:
        return _utc_timestamp(value, key, block="provenance", rule="§ 3.2")
    return _checked_int(value, what, _PROVENANCE_INT_MINIMUM[key])


def _job_spec_coordinates(value: Any, what: str) -> str:
    if value not in _GEOMETRY_COORDINATES:
        raise ValueError(f"{what} must be one of "
                         f"{list(_GEOMETRY_COORDINATES)}, got {value!r}")
    return value


# The typed blocks of a job.spec payload (spec § 5.9, $defs/JobSpecSolvation,
# JobSpecScf, JobSpecGeometry, JobSpecResources): block -> (field -> check,
# required fields). A check is ``(value, what) -> normalised value``.
_JOB_SPEC_BLOCKS: dict[str, tuple[dict[str, Any], tuple[str, ...]]] = {
    "solvation": ({"model": _checked_text, "solvent": _checked_text,
                   "epsilon": lambda v, w: _checked_number(v, w, positive=True)},
                  ("model",)),
    "scf": ({"convergence": lambda v, w: _checked_number(v, w, positive=True),
             "max_iterations": lambda v, w: _checked_int(v, w, 1),
             "guess": _checked_text,
             "damping": _checked_number}, ()),
    "geometry": ({"convergence": lambda v, w: _checked_number(v, w, positive=True),
                  "max_steps": lambda v, w: _checked_int(v, w, 1),
                  "coordinates": _job_spec_coordinates}, ()),
    "resources": ({"n_cores": lambda v, w: _checked_int(v, w, 1),
                   "memory_mb": lambda v, w: _checked_int(v, w, 1),
                   "walltime_s": lambda v, w: _checked_number(v, w, positive=True)},
                  ()),
}


def _job_spec_block(name: str, value: Any) -> dict[str, Any]:
    """One typed block of a ``job.spec`` payload, checked and normalised.

    A key outside the block's fields must begin ``x_``; engine keywords go
    under ``options`` (spec § 5.9).
    """
    if not isinstance(value, Mapping):
        raise ValueError(f"job.spec {name} must be a mapping (spec § 5.9), "
                         f"got {value!r}")
    rules, required = _JOB_SPEC_BLOCKS[name]
    block: dict[str, Any] = {}
    for key, item in value.items():
        if isinstance(key, str) and key.startswith("x_"):
            block[key] = item
        elif key in rules:
            block[key] = rules[key](item, f"job.spec {name}.{key}")
        else:
            raise ValueError(f"job.spec {name}: key {key!r} is not one of "
                             f"{sorted(rules)}; vendor keys must begin with "
                             f"'x_' and engine keywords go under options "
                             f"(spec § 5.9)")
    for key in required:
        if key not in block:
            raise ValueError(f"job.spec {name} requires {key!r} (spec § 5.9)")
    return block


_SMILES_RE = re.compile(r"^\S+$")
_INCHI_RE = re.compile(r"^InChI=\S+$")
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
_MEDIA_TYPE_RE = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]*/[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]*$")
_ARROWS = ("forward", "equilibrium", "retrosynthetic", "resonance")
_PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
_JPEG_SIGNATURE = b"\xff\xd8\xff"


def _structure_identity(**fields: Any) -> dict[str, Any]:
    """The chemical-identity fields of a ``structure`` section (spec § 5.1),
    checked; ``None`` values are left out."""
    out: dict[str, Any] = {}
    for key in ("charge", "multiplicity"):
        value = fields[key]
        if value is None:
            continue
        if not _is_int(value):
            raise ValueError(f"structure {key} must be an integer, got {value!r}")
        if key == "multiplicity" and int(value) < 1:
            raise ValueError(f"structure multiplicity must be >= 1, got {value!r}")
        out[key] = int(value)
    patterns = {"name": None, "smiles": _SMILES_RE, "inchi": _INCHI_RE,
                "inchikey": _INCHIKEY_RE}
    for key, pattern in patterns.items():
        value = fields[key]
        if value is None:
            continue
        if not isinstance(value, str) or not value:
            raise ValueError(f"structure {key} must be a non-empty string, "
                             f"got {value!r}")
        if pattern is not None and not pattern.fullmatch(value):
            raise ValueError(f"structure {key} {value!r} is not of the form "
                             f"spec § 5.1 gives")
        out[key] = value
    return out


def _charge_multiplicity_problem(sum_z: int, charge: Any,
                                 multiplicity: Any) -> str | None:
    """Why ``charge`` and ``multiplicity`` are impossible for atoms whose
    atomic numbers sum to ``sum_z`` (spec § 5.1), or ``None``."""
    if charge is None:
        return None
    electrons = sum_z - charge
    if electrons < 0:
        return (f"charge {charge} leaves {electrons} electrons for atoms with "
                f"atomic numbers summing to {sum_z}")
    if multiplicity is None:
        return None
    unpaired = multiplicity - 1
    if unpaired > electrons or (electrons - unpaired) % 2:
        return (f"multiplicity {multiplicity} is impossible with {electrons} "
                f"electrons (charge {charge}): the unpaired electrons must not "
                f"exceed the electron count and must have its parity")
    return None


def _participants(items: Any, what: str, required: bool) -> list[dict[str, Any]]:
    """A ``reaction.scheme`` role as participant objects (spec § 5.14): each
    item is a structure section id or a ``{structure_ref, coefficient?}``
    mapping (further keys must begin ``x_``)."""
    if items is None:
        if required:
            raise ValueError(f"{what} is required")
        return []
    if isinstance(items, (str, Mapping)) or not isinstance(items, Iterable):
        raise ValueError(f"{what} must be a list of participants")
    out = []
    for k, item in enumerate(items):
        where = f"{what}[{k}]"
        if isinstance(item, str):
            entry: dict[str, Any] = {"structure_ref": item}
        elif isinstance(item, Mapping):
            _check_keys(item, ("structure_ref", "coefficient"), where)
            entry = dict(item)
        else:
            raise ValueError(f"{where} must be a structure section id or a "
                             f"mapping, got {item!r}")
        ref = entry.get("structure_ref")
        if not isinstance(ref, str) or not ref:
            raise ValueError(f"{where}: structure_ref must be a non-empty string")
        if "coefficient" in entry:
            c = entry["coefficient"]
            _checked_number(c, f"{where} coefficient", positive=True)
            entry["coefficient"] = int(c) if _is_int(c) else float(c)
        out.append(entry)
    if required and not out:
        raise ValueError(f"{what} needs at least one participant")
    refs = [e["structure_ref"] for e in out]
    repeated = sorted({r for r in refs if refs.count(r) > 1})
    if repeated:
        raise ValueError(f"{what} names {repeated} more than once; a "
                         f"coefficient says how many (spec § 5.14)")
    return out


def _image_media_type(data: bytes, allowed: Sequence[str], what: str,
                      media_type: str | None) -> str:
    """The media type of preview bytes (spec § 5.15): ``media_type`` when
    given, checked against the bytes, else sniffed from them."""
    if media_type is None:
        if data.startswith(_PNG_SIGNATURE):
            media_type = "image/png"
        elif data.startswith(_JPEG_SIGNATURE):
            media_type = "image/jpeg"
        elif b"<svg" in data[:4096] and "image/svg+xml" in allowed:
            media_type = "image/svg+xml"
        else:
            raise ValueError(f"{what}: not PNG or JPEG bytes; pass "
                             f"{what}_media_type to say what they are")
    if media_type not in allowed:
        raise ValueError(f"{what}: media type {media_type!r} is not one of "
                         f"{list(allowed)} (spec § 5.15)")
    problem = _preview_bytes_problem(data, media_type)
    if problem:
        raise ValueError(f"{what}: {problem}")
    return media_type


def _preview_bytes_problem(data: bytes, media_type: str) -> str | None:
    """Why ``data`` are not what ``media_type`` says (spec § 5.15), or None."""
    if media_type == "image/png" and not data.startswith(_PNG_SIGNATURE):
        return "image/png bytes do not begin with the PNG signature"
    if media_type == "image/jpeg" and not data.startswith(_JPEG_SIGNATURE):
        return "image/jpeg bytes do not begin with FF D8 FF"
    if media_type == "image/svg+xml":
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            return "image/svg+xml bytes do not decode as UTF-8"
        if "<svg" not in text:
            return "image/svg+xml bytes contain no <svg element"
    return None


def _metadata_value(key: str, value: Any) -> Any:
    """Check and normalise one field of the root ``metadata`` block (spec § 3.5)."""
    if key.startswith("x_"):
        return value
    if key not in _METADATA_KEYS:
        raise ValueError(f"metadata key {key!r} is not one of "
                         f"{sorted(_METADATA_KEYS)}; vendor keys must begin "
                         f"with 'x_' (spec § 3.5)")
    if key in ("created_utc", "modified_utc"):
        return _utc_timestamp(value, key)
    if key == "uuid":
        text = str(value).lower() if isinstance(value, (str, _uuid.UUID)) else value
        if not isinstance(text, str) or not _UUID_RE.match(text):
            raise ValueError(f"metadata uuid must be an RFC 4122 UUID in "
                             f"8-4-4-4-12 hex form (spec § 3.5), got {value!r}")
        return text
    if key in _METADATA_STRINGS:
        if not isinstance(value, str):
            raise ValueError(f"metadata {key} must be a string, got {value!r}")
        if not value and key in ("license", "doi"):
            raise ValueError(f"metadata {key} must not be empty")
        return value
    if key in _METADATA_LISTS:
        if isinstance(value, (str, bytes)) or not all(
                isinstance(x, str) for x in value):
            raise ValueError(f"metadata {key} must be a list of strings, "
                             f"got {value!r}")
        return list(value)
    # authors
    if isinstance(value, (str, bytes, Mapping)):
        raise ValueError("metadata authors must be a list of author objects")
    authors = []
    for author in value:
        if not isinstance(author, Mapping):
            raise ValueError(f"metadata author {author!r} must be a mapping")
        extra = {k for k in author if not str(k).startswith("x_")} - _AUTHOR_KEYS
        if extra:
            raise ValueError(f"metadata author keys {sorted(extra)} are not "
                             f"one of {sorted(_AUTHOR_KEYS)} (spec § 3.5)")
        name = author.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError(f"metadata author {author!r} needs a non-empty "
                             f"'name'")
        for k in ("orcid", "affiliation", "email"):
            if k in author and not isinstance(author[k], str):
                raise ValueError(f"metadata author {k} must be a string")
        authors.append(dict(author))
    return authors


class _Unset:
    """A missing keyword, distinct from an explicit JSON null."""


_UNSET = _Unset()



def _check_residue_insertion_codes(section_id: str, residues: Any,
                                   atoms: Sequence[Mapping[str, Any]]) -> None:
    """Raise if a residue and one of its atoms name different insertion codes.

    Residue identity is (chain, seq, insertion_code), absent meaning none
    (spec § 5.1, *Biomolecular atom identity*). Atoms without the per-atom
    field, and out-of-range indices (left to the validator), are not compared.
    """
    if not isinstance(residues, Sequence):
        return
    for r, residue in enumerate(residues):
        if not isinstance(residue, Mapping):
            continue
        code = residue.get("insertion_code")
        for index in residue.get("atom_indices") or ():
            if not isinstance(index, int) or not 0 <= index < len(atoms):
                continue
            atom = atoms[index]
            if "insertion_code" in atom and atom["insertion_code"] != code:
                raise ValueError(
                    f"structure {section_id!r}: residues[{r}] has insertion_code "
                    f"{code!r} but atom {index} has {atom['insertion_code']!r} "
                    f"(spec § 5.1)")

class QvfWriter:
    """Accumulate sections and root metadata, then write a ``.qvf`` archive.

    All ``add_*`` methods return the section ``id`` they created. Payloads are
    accepted as plain Python lists / dicts and NumPy arrays; the writer converts,
    checksums, and lays out the ZIP for you.

    With ``auto_metadata=True`` (the default) the manifest's root ``metadata``
    block gets a ``created_utc`` time stamp and a random (version 4) ``uuid``
    unless :meth:`set_metadata` supplied them (spec § 3.5). Both are made once,
    on the first :meth:`build_manifest` (or ``to_bytes`` / ``write``), and
    reused afterwards. Pass ``auto_metadata=False``, or set the attribute of
    the same name, for a manifest that is the same on every run, e.g. a test
    vector.
    """

    def __init__(self, program: str, version: str, calculation: str = "", *,
                 auto_metadata: bool = True) -> None:
        if not program or not version:
            raise ValueError("source.program and source.version are required")
        self._source = {
            "program": str(program),
            "version": str(version),
            "calculation": str(calculation),
        }
        self._sections: list[dict[str, Any]] = []
        self._files: dict[str, bytes] = {}
        #: Member paths written uncompressed (already compressed images).
        self._stored: set[str] = set()
        self._ids: set[str] = set()
        self._root: dict[str, Any] = {}
        #: Add ``created_utc`` and ``uuid`` to the root ``metadata`` (spec § 3.5).
        self.auto_metadata = bool(auto_metadata)
        self._auto_values: dict[str, str] | None = None

    # -- root metadata ----------------------------------------------------

    def set_provenance(self, *, allow_hostname: bool = False,
                       **fields: Any) -> None:
        """Set the root ``provenance`` block (method, basis, energies, …).

        Repeated calls merge into the existing block. ``scf_energy`` /
        ``fermi_energy`` may be passed as ``scf_energy_eh`` /
        ``fermi_energy_ev`` convenience floats, which are wrapped into the
        ``{value, units}`` shape the schema requires, with the unit names
        ``"hartree"`` and ``"ev"`` (spec § 4.3).

        ``hostname`` is deprecated (spec § 3.2): archives are meant to be
        shared, and machine names are private operational detail. A
        ``hostname`` field passed here is therefore dropped silently (no
        warning) unless the keyword-only ``allow_hostname=True`` is given, in
        which case it is written as passed.

        The typed fields are checked, and a wrong one raises ``ValueError``
        before anything changes: ``software`` (a list of ``{name, version?,
        url?}`` mappings, ``name`` non-empty; further keys must begin
        ``x_``), ``git_commit`` (7 to 64 hexadecimal characters),
        ``n_threads`` and ``n_mpi_ranks`` (integers >= 1), ``n_gpus`` (an
        integer >= 0), and ``started_utc`` and ``finished_utc`` (strings
        ``YYYY-MM-DDThh:mm:ssZ``, or timezone-aware ``datetime`` objects,
        which are converted to UTC). A value of ``None`` removes that key.
        The other fields of the open block are written as passed.
        """
        if not allow_hostname:
            fields.pop("hostname", None)
        typed: dict[str, Any] = {}
        for key in [k for k in fields if k in _PROVENANCE_TYPED]:
            value = fields.pop(key)
            typed[key] = None if value is None else _provenance_field(key, value)
        prov = dict(self._root.get("provenance", {}))
        if "scf_energy_eh" in fields:
            prov["scf_energy"] = {"value": float(fields.pop("scf_energy_eh")),
                                  "units": "hartree"}
        if "fermi_energy_ev" in fields:
            prov["fermi_energy"] = {"value": float(fields.pop("fermi_energy_ev")),
                                    "units": "ev"}
        prov.update(fields)
        for key, value in typed.items():
            if value is None:
                prov.pop(key, None)
            else:
                prov[key] = value
        self._root["provenance"] = prov

    def set_metadata(self, **fields: Any) -> None:
        """Set fields of the root ``metadata`` block (spec § 3.5).

        Keys: ``title``, ``description``, ``license`` (an SPDX identifier is
        recommended), ``doi`` (strings); ``keywords``, ``derived_from``,
        ``tags`` (lists of strings); ``created_utc``, ``modified_utc`` (strings
        ``YYYY-MM-DDThh:mm:ssZ``, or timezone-aware ``datetime`` objects, which
        are converted to UTC); ``uuid`` (an RFC 4122 string or ``uuid.UUID``,
        written in lowercase); ``authors`` (a list of ``{name, orcid?,
        affiliation?, email?}``); and vendor keys beginning ``x_``. Anything
        else, or a value of the wrong form, is a ``ValueError``. ``email`` and
        ``affiliation`` are personal data: pass them only when the person
        concerned wants them in the archive.

        Repeated calls merge into the existing block; a value of ``None``
        removes that key. Fields set here take precedence over the automatic
        ``created_utc`` / ``uuid`` (see the class docstring).
        """
        meta = dict(self._root.get("metadata", {}))
        updates = {}
        for key, value in fields.items():
            if value is None:
                meta.pop(key, None)
            else:
                updates[key] = _metadata_value(key, value)
        meta.update(updates)
        self._root["metadata"] = meta

    def set_thermochemistry(self, **fields: Any) -> None:
        self._root["thermochemistry"] = dict(fields)

    def set_energies(self, total_eh: float | Mapping[str, Any], *,
                     components: Mapping[str, Any] | None = None,
                     by_method: Sequence[Mapping[str, Any]] | None = None,
                     units: str = "hartree", method: str | None = None,
                     units_note: str | None = None) -> None:
        """Set the root ``energies`` block (spec § 3.6, ``$defs/Energies``).

        ``total_eh`` is the energy of the archive's primary result and
        ``method`` its method (e.g. ``"rks"``). ``components`` maps component
        names (``nuclear_repulsion``, ``electronic``, ``one_electron``,
        ``two_electron``, ``exchange_correlation``, ``exchange``,
        ``correlation``, ``dispersion``, ``solvation``, ``zpve``, ``thermal``,
        ``embedding`` or an ``x_`` key) to energies. ``by_method`` lists
        ``{"method", "total", "components"?}`` entries for further methods.

        Every energy may be a number, which is wrapped as ``{"value": x,
        "units": units}``, or an explicit ``{value, units}`` mapping. Nothing
        is converted: ``units`` (default ``"hartree"``, the unit producers
        SHOULD write) only labels the numbers. A later call replaces the
        block. An energy stated for the method ``provenance.method`` names
        must agree with ``provenance.scf_energy`` (the validator checks).
        """
        if not isinstance(units, str) or not units:
            raise ValueError("energies: units must be a non-empty string")
        block: dict[str, Any] = {"total": _energy(total_eh, units, "energies.total")}
        if method is not None:
            if not isinstance(method, str) or not method:
                raise ValueError("energies: method must be a non-empty string")
            block["method"] = method
        if components is not None:
            block["components"] = _energy_components(components, units,
                                                     "energies.components")
        if by_method is not None:
            entries = []
            for k, entry in enumerate(by_method):
                what = f"energies.by_method[{k}]"
                unknown = set(entry) - {"method", "total", "components"}
                if unknown:
                    raise ValueError(f"{what}: unknown keys {sorted(unknown)}")
                name = entry.get("method")
                if not isinstance(name, str) or not name:
                    raise ValueError(f"{what}: method must be a non-empty string")
                if "total" not in entry:
                    raise ValueError(f"{what}: total is required")
                item: dict[str, Any] = {
                    "method": name,
                    "total": _energy(entry["total"], units, f"{what}.total")}
                if entry.get("components") is not None:
                    item["components"] = _energy_components(
                        entry["components"], units, f"{what}.components")
                entries.append(item)
            block["by_method"] = entries
        if units_note is not None:
            block["units_note"] = str(units_note)
        self._root["energies"] = block

    def set_dipole_moment(self, total_debye: float | None = None,
                          vector_debye: Sequence[float] | None = None,
                          origin: str | None = None) -> None:
        d: dict[str, Any] = {}
        if total_debye is not None:
            d["total_debye"] = float(total_debye)
        if vector_debye is not None:
            d["vector_debye"] = [float(x) for x in vector_debye]
        if origin is not None:
            d["origin"] = str(origin)
        self._root["dipole_moment"] = d

    def set_constraints(self, frozen_atoms: Sequence[int] | None = None,
                        distance_constraints: Sequence[Mapping[str, Any]] | None = None
                        ) -> None:
        c: dict[str, Any] = {}
        if frozen_atoms is not None:
            c["frozen_atoms"] = [int(i) for i in frozen_atoms]
        if distance_constraints is not None:
            c["distance_constraints"] = [dict(x) for x in distance_constraints]
        self._root["constraints"] = c

    def set_viewer_defaults(self, defaults: Mapping[str, Any]) -> None:
        self._root["viewer_defaults"] = dict(defaults)

    def set_extensions(self, extensions: Mapping[str, Mapping[str, Any]]) -> None:
        """Declare vendor extension namespaces (spec § 7.4). Each key is the
        namespace with its ``x_`` prefix, as in the section kinds, and each
        value ``{version, schema_uri?, critical?}`` with ``version`` a
        string: ``{"x_orca": {"version": "1.0"}}``."""
        self._root["extensions"] = {k: dict(v) for k, v in extensions.items()}

    # -- internal helpers -------------------------------------------------

    def _check_new_id(self, section_id: str) -> None:
        if not isinstance(section_id, str) or not _ID_RE.fullmatch(section_id):
            raise ValueError(f"section id {section_id!r} must match [A-Za-z0-9_.-]+")
        if section_id in _RESERVED_IDS:
            raise ValueError(f"section id {section_id!r} is reserved "
                             f"(a viewer_defaults key, spec § 3.1, § 3.4)")
        if section_id in self._ids:
            raise ValueError(f"duplicate section id {section_id!r}")

    def _reserve_id(self, section_id: str) -> str:
        self._check_new_id(section_id)
        self._ids.add(section_id)
        return section_id

    def _add_file(self, path: str, data: bytes) -> str:
        """Queue ``data`` under ``path`` and return the path it gets. A path
        that the archive being updated already holds is never rewritten
        (spec § 2.3): the member gets a fresh one, ``<stem>.<n><suffix>``."""
        # Spec § 2.1: relative, '/'-separated, no '..' segment.
        if (".." in path.split("/") or path.startswith("/") or "\\" in path
                or "\0" in path or _DRIVE_RE.match(path)):
            raise ValueError(f"illegal member path {path!r}")
        if path in self._files:
            if self._files[path] is not None:
                raise ValueError(f"duplicate member path {path!r}")
            head, slash, tail = path.rpartition("/")
            stem, dot, suffix = tail.partition(".")
            n = 2
            while True:
                fresh = f"{head}{slash}{stem}.{n}{dot}{suffix}"
                if fresh not in self._files:
                    break
                n += 1
            path = fresh
        self._files[path] = data
        return path

    def _json_member(self, path: str, obj: Any) -> dict[str, Any]:
        data = _json_bytes(obj)
        path = self._add_file(path, data)
        return {"path": path, "format": "json", "sha256": _sha256_hex(data),
                "size": len(data)}

    def _binary_member(self, path: str, data: Any, dtype: str,
                       *, with_shape: bool = True) -> dict[str, Any]:
        raw, dname, shape = _le_bytes(data, dtype)
        path = self._add_file(path, raw)
        member: dict[str, Any] = {"path": path, "format": "binary",
                                  "sha256": _sha256_hex(raw)}
        if with_shape:
            member["dtype"] = dname
            member["shape"] = shape
        member["size"] = len(raw)
        return member

    def _opaque_member(self, path: str, data: bytes,
                       media_type: str | None = None) -> dict[str, Any]:
        """An opaque binary member (no dtype/shape, spec § 4.2) of ``data``,
        with its ``media_type`` when given."""
        if media_type is not None and not (isinstance(media_type, str)
                                           and _MEDIA_TYPE_RE.fullmatch(media_type)):
            raise ValueError(f"member {path!r}: media type {media_type!r} is "
                             f"not type/subtype (spec § 4.2)")
        path = self._add_file(path, data)
        member: dict[str, Any] = {"path": path, "format": "binary",
                                  "sha256": _sha256_hex(data), "size": len(data)}
        if media_type is not None:
            member["media_type"] = media_type
        return member

    def _push(self, section: dict[str, Any], *,
              label: str | None = None) -> str:
        """Append ``section``; a ``label`` (spec § 4.1) is written after its
        ``kind``."""
        if label is not None:
            if not isinstance(label, str):
                raise ValueError(f"section {section['id']!r}: label must be a "
                                 f"string, got {label!r}")
            ordered: dict[str, Any] = {}
            for key, value in section.items():
                if key != "label":
                    ordered[key] = value
                if key == "kind":
                    ordered["label"] = label
            section = ordered
        self._sections.append(section)
        return section["id"]

    # -- section references and grouping ----------------------------------

    def annotate_section(self, section_id: str, *,
                         structure_ref: str | None = None,
                         run_ref: str | None = None,
                         group: str | None = None,
                         label: str | None = None) -> str:
        """Set the generic fields of an already added section (spec § 4.1).

        ``structure_ref`` names the ``structure`` section the section
        describes, ``run_ref`` the ``run.record`` section of the run that
        produced it; both must already have been added, with that kind.
        ``group`` is a free label that groups related sections for display,
        and ``label`` the section's display label. ``None`` leaves a field as
        it is. A section's ``structure_ref`` takes precedence over the one in
        a ``wavefunction.gto`` ``basis`` payload, so keep the two equal.
        Returns ``section_id``; raises ``ValueError`` (and changes nothing) on
        an unknown section, a reference that does not resolve to the right
        kind, or a value that is not a string.
        """
        by_id = {s["id"]: s for s in self._sections}
        section = by_id.get(section_id)
        if section is None:
            raise ValueError(f"annotate_section: no section with id "
                             f"{section_id!r}")
        fields = {"structure_ref": structure_ref, "run_ref": run_ref,
                  "group": group, "label": label}
        updates: dict[str, str] = {}
        for key, value in fields.items():
            if value is None:
                continue
            if not isinstance(value, str):
                raise ValueError(f"annotate_section: {key} must be a string, "
                                 f"got {value!r}")
            target = _SECTION_REF_KINDS.get(key)
            if target is not None:
                ref = by_id.get(value)
                if ref is None or ref["kind"] != target:
                    found = ("no such section" if ref is None
                             else f"a {ref['kind']!r} section")
                    raise ValueError(
                        f"annotate_section: {key} {value!r} must name an "
                        f"existing {target!r} section, found {found} "
                        f"(spec § 4.1)")
            updates[key] = value
        section.update(updates)
        return section_id

    def _check_ref(self, what: str, ref: Any, kind: str) -> None:
        """Raise ``ValueError`` unless ``ref`` names an added ``kind`` section."""
        found = next((s for s in self._sections if s["id"] == ref), None)
        if not isinstance(ref, str) or found is None or found["kind"] != kind:
            got = ("no such section" if found is None
                   else f"a {found['kind']!r} section")
            raise ValueError(f"{what} {ref!r} must name an existing {kind!r} "
                             f"section, found {got}")

    def remove_section(self, section_id: str) -> dict[str, Any]:
        """Remove the section ``section_id`` from the manifest and return it.

        In a new archive its members are dropped. In a writer from
        :meth:`open_for_update` they stay in the archive as superseded
        entries (``W_ORPHAN_MEMBER`` until :func:`repack`, or ``write`` /
        ``to_bytes`` with ``repack=True``), as spec § 2.3 prescribes.

        To **replace** a section, remove it and add the new one under the
        same id; members that would reuse a path the archive already holds
        get a fresh one, because an update never rewrites an entry. Nothing
        else is changed: references to the removed section from other
        sections (``structure_ref``, ``run_ref``, ``job_ref``, a
        ``reaction.scheme`` participant, …) dangle until a section with the
        id is added again, and a validator reports any left dangling.
        Raises ``ValueError`` for an unknown id.
        """
        for index, section in enumerate(self._sections):
            if section.get("id") == section_id:
                break
        else:
            raise ValueError(f"remove_section: no section with id "
                             f"{section_id!r}")
        del self._sections[index]
        self._ids.discard(section_id)
        members = section.get("members")
        for member in (members.values() if isinstance(members, dict) else ()):
            path = member.get("path") if isinstance(member, dict) else None
            if isinstance(path, str) and self._files.get(path) is not None:
                del self._files[path]       # not yet written: drop it
        return section

    # -- structure / bonds -----------------------------------------------

    def add_structure(self, atoms: Sequence[Mapping[str, Any]], *,
                      pbc: Sequence[bool] | None = None,
                      lattice_vectors: Sequence[Sequence[float]] | Lattice | None = None,
                      dimensionality: int | None = None,
                      bonds: Sequence[Mapping[str, Any]] | None = None,
                      biomolecule: Mapping[str, Any] | None = None,
                      charge: int | None = None,
                      multiplicity: int | None = None,
                      name: str | None = None,
                      smiles: str | None = None,
                      inchi: str | None = None,
                      inchikey: str | None = None,
                      label: str | None = None,
                      section_id: str = "structure") -> str:
        """Add a ``structure`` section. ``atoms`` positions are in Ångström.

        Each atom is ``{"symbol", "position": [x, y, z], "atomic_number"}``;
        missing ``symbol``/``atomic_number`` are inferred from each other.
        ``bonds`` (optional) is a list of ``{i, j, order}`` written as the
        section's ``bonds`` member. ``biomolecule`` adds ``chains`` /
        ``residues`` / ``secondary_structure`` / ``b_factors`` peer fields for
        cartoon rendering (spec § 5.1). The per-atom biomolecular fields
        ``atom_name`` / ``residue_name`` / ``residue_seq`` / ``chain_id`` /
        ``b_factor`` / ``insertion_code`` / ``alt_loc`` / ``occupancy`` /
        ``hetero`` are passed through from each atom mapping when present.
        A ``residues`` entry whose ``insertion_code`` disagrees with that of
        one of its atoms raises ``ValueError`` (spec § 5.1, *Biomolecular
        atom identity*).

        ``pbc`` marks which axes are periodic and is the normative carrier of
        that information (spec § 5.1); ``lattice_vectors`` is required whenever
        any axis is periodic, and its row ``i`` pairs with ``pbc[i]``. It is a
        3x3 nested sequence of row vectors in Å (written verbatim) or a
        :class:`Lattice`, which is written as its ``rows_angstrom()``; the
        payload is identical either way (spec § 4.3, "Lattice conventions").
        ``dimensionality`` is derived — pass it only to have it written out, and
        it must equal ``sum(pbc)``. Passing ``dimensionality`` without ``pbc``
        derives the leading axes as the periodic ones, which is this writer's
        convenience, not a rule of the format.

        The chemical identity of the species (spec § 5.1, *Chemical
        identity*) goes on the section object: ``charge`` (an integer),
        ``multiplicity`` (an integer >= 1), ``name``, ``smiles``, ``inchi``
        (beginning ``InChI=``) and ``inchikey``. ``label`` is the section's
        display label, for a species in a scheme its compound number
        (``"1"``, ``"3a"``). A charge and multiplicity impossible for the
        atoms (a negative electron count, or unpaired electrons beyond it or
        of the other parity) raise ``ValueError``.
        """
        self._check_new_id(section_id)
        identity = _structure_identity(charge=charge, multiplicity=multiplicity,
                                       name=name, smiles=smiles, inchi=inchi,
                                       inchikey=inchikey)
        norm_atoms = []
        for a in atoms:
            z = a.get("atomic_number")
            sym = a.get("symbol")
            if z is None and sym is not None:
                z = _ELEMENT_SYMBOLS.index(sym) if sym in _ELEMENT_SYMBOLS else 0
            if sym is None and z is not None:
                sym = _symbol(int(z))
            atom: dict[str, Any] = {
                "symbol": sym,
                "position": [float(x) for x in a["position"]],
                "atomic_number": int(z) if z is not None else 0,
            }
            for extra in _ATOM_BIOMOLECULE_FIELDS:
                if a.get(extra) is not None:
                    atom[extra] = a[extra]
            norm_atoms.append(atom)
        if pbc is not None:
            flags = [bool(x) for x in pbc]
            if len(flags) != 3:
                raise ValueError(f"pbc must have 3 entries, got {len(flags)}")
        elif dimensionality is not None:
            flags = [i < int(dimensionality) for i in range(3)]
        else:
            flags = [False, False, False]

        if dimensionality is not None and int(dimensionality) != sum(flags):
            raise ValueError(
                f"dimensionality {int(dimensionality)} != sum(pbc) {sum(flags)}; "
                f"the QVF invariant is dimensionality == sum(pbc) (spec § 5.1)")

        lattice = (
            _lattice_as_rows(lattice_vectors)
            if lattice_vectors is not None else None
        )
        if any(flags) and lattice is None:
            raise ValueError(
                f"pbc {flags} marks a periodic axis, so lattice_vectors is "
                f"required (spec § 5.1)")

        struct: dict[str, Any] = {
            "atoms": norm_atoms,
            "pbc": flags,
            "lattice_vectors": lattice,
        }
        problem = _charge_multiplicity_problem(
            sum(a["atomic_number"] for a in norm_atoms),
            identity.get("charge"), identity.get("multiplicity"))
        if problem:
            raise ValueError(f"structure {section_id!r}: {problem} (spec § 5.1)")
        if dimensionality is not None:
            struct["dimensionality"] = int(dimensionality)
        members = {"structure": self._json_member(
            f"{section_id}/structure.json", struct)}
        if bonds is not None:
            members["bonds"] = self._json_member(
                f"{section_id}/bonds.json", {"pairs": [dict(b) for b in bonds]})
        section: dict[str, Any] = {
            "id": self._reserve_id(section_id), "kind": "structure",
            "members": members,
        }
        if biomolecule is not None:
            for key in ("chains", "residues", "secondary_structure", "b_factors"):
                if key in biomolecule:
                    section[key] = biomolecule[key]
            _check_residue_insertion_codes(section_id,
                                           biomolecule.get("residues"),
                                           norm_atoms)
        section.update(identity)
        return self._push(section, label=label)

    def add_bonds(self, pairs: Sequence[Mapping[str, Any]], *,
                 structure_ref: str | None = None,
                 label: str | None = None, section_id: str = "bonds") -> str:
        """Add a standalone ``bonds`` connectivity section. Each pair is
        ``{i, j, order?}`` (``$defs/BondsPayload``, spec § 5.7).

        ``structure_ref`` names the ``structure`` section, already added,
        whose atoms the indices address, and is written on the section
        object. Without it the bonds belong to the archive's first structure
        (spec § 5.1), so an archive with several structures should pass it
        (or give each structure its own ``bonds`` with :meth:`add_structure`).
        A ``structure_ref`` that names no added ``structure`` section raises
        ``ValueError`` before anything is added; validators check the indices
        against that structure's atoms.
        """
        self._check_new_id(section_id)
        if structure_ref is not None:
            self._check_ref("bonds structure_ref", structure_ref, "structure")
        sid = self._reserve_id(section_id)
        member = self._json_member(f"{sid}/connectivity.json",
                                   {"pairs": [dict(p) for p in pairs]})
        section: dict[str, Any] = {"id": sid, "kind": "bonds"}
        if structure_ref is not None:
            section["structure_ref"] = structure_ref
        section["members"] = {"bonds": member}
        return self._push(section, label=label)

    def add_bond_orders(self, method: str, pairs: Sequence[Mapping[str, Any]], *,
                       label: str | None = None,
                       section_id: str = "bond_orders") -> str:
        """Add a ``bond_orders`` table (``$defs/BondOrdersPayload``, spec
        § 5.7). ``pairs`` each carry ``i``, ``j``, ``order`` and optionally
        ``distance_ang``, ``symbol_i``, ``symbol_j``."""
        self._check_new_id(section_id)
        payload = {"method": str(method), "pairs": [dict(p) for p in pairs]}
        member = self._json_member(f"{section_id}/orders.json", payload)
        return self._push({"id": self._reserve_id(section_id),
                          "kind": "bond_orders",
                          "members": {"bond_orders": member}}, label=label)

    # -- volumes ----------------------------------------------------------

    # The scalar volume kinds: `data` is [nx, ny, nz] (spec § 5.2).
    _VOLUME_KINDS = frozenset({
        "volume.density", "volume.orbital", "volume.spin", "volume.elf",
        "volume.generic", "volume.potential", "volume.rdg", "volume.difference",
        "volume.transition_density",
        "volume.laplacian",
    })
    # The volume kinds with 4-D `data` (spec § 5.2), each with the extent its
    # trailing axis must have (None: any) and the form of its shape.
    _VOLUME_4D_KINDS = {
        "volume.vector": (3, "[nx, ny, nz, 3]"),
        "volume.complex": (2, "[nx, ny, nz, 2]"),
        "volume.series": (None, "[n, nx, ny, nz]"),
    }
    _SERIES_AXIS_KINDS = ("time", "frame", "state", "custom")

    def add_volume(self, kind: str, data: np.ndarray, origin: Sequence[float],
                  voxel_vectors: Sequence[Sequence[float]], *,
                  label: str | None = None, dtype: str = "float32",
                  component: str | None = None, section_id: str | None = None,
                  operand_a: str | None = None, operand_b: str | None = None,
                  description: str | None = None,
                  units: Mapping[str, str] | None = None,
                  periodic: Sequence[bool] | None = None,
                  includes_endpoint: bool | None = None,
                  structure_ref: str | None = None,
                  axis: Mapping[str, Any] | None = None,
                  excited_states_ref: str | None = None,
                  state_index: int | None = None,
                  wavefunction_ref: str | None = None,
                  mo_index: int | None = None,
                  spin: str | None = None) -> str:
        """Add any ``volume.*`` section; ``origin`` and ``voxel_vectors`` are
        in bohr (spec § 5.2). For ``volume.difference`` you may pass
        ``operand_a``/``operand_b`` section ids.

        ``data`` is a 3-D array ``[nx, ny, nz]`` for the scalar kinds
        (``volume.density``, ``orbital``, ``spin``, ``elf``, ``generic``,
        ``potential``, ``rdg``, ``difference``, ``laplacian``);
        ``[nx, ny, nz, 3]`` for ``volume.vector`` (Cartesian x, y, z on the
        last axis); ``[nx, ny, nz, 2]`` for ``volume.complex`` (``[re, im]``,
        or pass a complex 3-D array and it is split); and ``[n, nx, ny, nz]``
        for ``volume.series``, which also needs ``axis``: a mapping with
        ``kind`` (``"time"``, ``"frame"``, ``"state"`` or ``"custom"``),
        ``values`` (one number per slab) and optionally ``units`` (a unit
        name), ``labels`` (one string per slab) and ``trajectory_ref`` (an
        added ``trajectory`` or ``reaction.path`` section). The grid
        ``shape`` is the three spatial extents of ``data``.

        For ``volume.transition_density`` (the transition density of one
        excitation, e/bohr³) ``excited_states_ref`` names the
        ``excited_states`` section of the state and ``state_index`` its
        ``index`` there, which needs ``excited_states_ref``. The reference is
        not resolved here, because the ``excited_states`` section usually
        names this one through ``density_ref`` and is added after it; the
        validators resolve it (spec § 5.12).

        ``wavefunction_ref`` (``volume.density``, ``volume.spin`` and
        ``volume.orbital``) names the ``wavefunction.gto`` section the field
        was computed from, which must already have been added. For a
        ``volume.orbital``, ``mo_index`` is the 0-based row of the orbital in
        that section's coefficient member and ``spin`` (``"alpha"`` or
        ``"beta"``) the member of an unrestricted wavefunction; ``spin`` is
        required with the alpha and beta coefficients and refused with
        ``mo_coefficients`` (spec § 5.2, Volumes of a wavefunction).

        ``units`` (optional) is written as the grid's ``units`` declaration,
        e.g. ``{"origin": "bohr", "voxel_vectors": "bohr"}``: unit names of
        spec § 4.3, keyed by a grid field or a member role (``data`` for the
        unit of the values). It declares and never converts, so ``origin``
        and ``voxel_vectors`` can only restate ``"bohr"``; an unknown unit
        name raises ``ValueError``.

        ``periodic`` (three booleans), ``includes_endpoint`` and
        ``structure_ref`` (an added ``structure`` section) are written to the
        grid as given (spec § 5.2). Along a periodic axis a the grid must span
        lattice vector a of its structure in bohr: ``voxel_vectors[a]`` times
        ``shape[a] - 1`` when ``includes_endpoint`` is true (the default,
        cube style) or times ``shape[a]`` when it is false (CHGCAR style).
        The validator checks the span; this writer does not."""
        if section_id is not None:
            self._check_new_id(section_id)
        trailing_form = self._VOLUME_4D_KINDS.get(kind)
        if kind not in self._VOLUME_KINDS and trailing_form is None:
            raise ValueError(f"{kind!r} is not a volume kind")
        arr = np.asarray(data)
        if kind == "volume.complex" and np.iscomplexobj(arr) and arr.ndim == 3:
            arr = np.stack([arr.real, arr.imag], axis=-1)
        if trailing_form is None:
            if arr.ndim != 3:
                raise ValueError(f"{kind}: data must be 3-D, got ndim={arr.ndim}")
            spatial = list(arr.shape)
        else:
            trailing, form = trailing_form
            if np.iscomplexobj(arr):
                raise ValueError(f"{kind}: data must be real, got {arr.dtype}")
            if arr.ndim != 4 or (trailing is not None and arr.shape[3] != trailing):
                raise ValueError(f"{kind}: data must be {form}, got shape "
                                 f"{list(arr.shape)}")
            spatial = list(arr.shape[1:] if kind == "volume.series"
                           else arr.shape[:3])
        if (axis is not None) != (kind == "volume.series"):
            raise ValueError("volume.series needs axis=, and only volume.series "
                             f"takes one ({kind})")
        axis_payload = (self._series_axis(axis, arr.shape[0])
                        if axis is not None else None)
        state_fields: dict[str, Any] = {}
        if excited_states_ref is not None or state_index is not None:
            if kind != "volume.transition_density":
                raise ValueError(f"{kind}: excited_states_ref and state_index "
                                 f"apply to volume.transition_density only")
            if excited_states_ref is None:
                raise ValueError("volume.transition_density: state_index needs "
                                 "excited_states_ref (spec § 5.2)")
            if not isinstance(excited_states_ref, str) or not excited_states_ref:
                raise ValueError("volume.transition_density: excited_states_ref "
                                 "must be a section id")
            state_fields["excited_states_ref"] = excited_states_ref
            if state_index is not None:
                state_fields["state_index"] = _positive_int(
                    state_index, "volume.transition_density state_index")
        state_fields.update(self._wavefunction_volume_fields(
            kind, wavefunction_ref, mo_index, spin))
        sid = section_id or (kind.replace(".", "_") + f"_{len(self._sections)}")
        self._check_new_id(sid)
        stem = sid
        grid: dict[str, Any] = {
            "origin": [float(x) for x in origin],
            "voxel_vectors": [[float(x) for x in v] for v in voxel_vectors],
            "shape": spatial,
        }
        if periodic is not None:
            flags = list(periodic)
            if len(flags) != 3 or not all(isinstance(x, (bool, np.bool_))
                                          for x in flags):
                raise ValueError(f"{kind}: periodic must be three booleans, "
                                 f"got {periodic!r}")
            grid["periodic"] = [bool(x) for x in flags]
        if includes_endpoint is not None:
            if not isinstance(includes_endpoint, (bool, np.bool_)):
                raise ValueError(f"{kind}: includes_endpoint must be a boolean")
            grid["includes_endpoint"] = bool(includes_endpoint)
        if structure_ref is not None:
            if not any(s["id"] == structure_ref and s["kind"] == "structure"
                       for s in self._sections):
                raise ValueError(f"{kind}: structure_ref {structure_ref!r} names "
                                 f"no structure section added so far")
            grid["structure_ref"] = structure_ref
        _with_units(grid, units, kind, "VolumeGridPayload")
        members = {
            "grid": self._json_member(f"volumes/{stem}_grid.json", grid),
            "data": self._binary_member(f"volumes/{stem}.dat", arr, dtype),
        }
        if axis_payload is not None:
            members["axis"] = self._json_member(f"volumes/{stem}_axis.json",
                                                axis_payload)
        section: dict[str, Any] = {"id": self._reserve_id(sid), "kind": kind,
                                   "members": members}
        if label is not None:
            section["label"] = label
        if component is not None:
            section["component"] = component
        if kind == "volume.difference" and (operand_a or operand_b):
            if not (operand_a and operand_b):
                raise ValueError("volume.difference needs both operand_a and operand_b")
            section["operand_a"] = operand_a
            section["operand_b"] = operand_b
            if description is not None:
                section["description"] = description
        section.update(state_fields)
        return self._push(section)

    def _wavefunction_volume_fields(self, kind: str, wavefunction_ref: Any,
                                    mo_index: Any, spin: Any) -> dict[str, Any]:
        """The ``wavefunction_ref``, ``mo_index`` and ``spin`` section fields
        of a volume (spec § 5.2, Volumes of a wavefunction), checked against
        the ``wavefunction.gto`` section already added."""
        if wavefunction_ref is None and mo_index is None and spin is None:
            return {}
        if kind not in ("volume.density", "volume.spin", "volume.orbital"):
            raise ValueError(f"{kind}: wavefunction_ref applies to volume.density, "
                             f"volume.spin and volume.orbital only")
        if (mo_index is not None or spin is not None) and kind != "volume.orbital":
            raise ValueError(f"{kind}: mo_index and spin apply to volume.orbital only")
        if wavefunction_ref is None:
            raise ValueError(f"{kind}: mo_index needs wavefunction_ref (spec § 5.2)")
        if spin is not None and mo_index is None:
            raise ValueError(f"{kind}: spin needs mo_index (spec § 5.2)")
        self._check_ref(f"{kind} wavefunction_ref", wavefunction_ref, "wavefunction.gto")
        fields: dict[str, Any] = {"wavefunction_ref": wavefunction_ref}
        if mo_index is None:
            return fields
        if (isinstance(mo_index, (bool, np.bool_))
                or not isinstance(mo_index, (int, np.integer)) or mo_index < 0):
            raise ValueError(f"{kind}: mo_index must be an integer >= 0, got {mo_index!r}")
        members = next(s for s in self._sections
                       if s["id"] == wavefunction_ref)["members"]
        if "mo_coefficients" in members:
            if spin is not None:
                raise ValueError(f"{kind}: wavefunction {wavefunction_ref!r} has one "
                                 f"set of orbitals; spin must be omitted (spec § 5.2)")
            role = "mo_coefficients"
        else:
            if spin not in ("alpha", "beta"):
                raise ValueError(f"{kind}: wavefunction {wavefunction_ref!r} is "
                                 f"unrestricted; spin must be 'alpha' or 'beta' "
                                 f"(spec § 5.2)")
            role = f"mo_coefficients_{spin}"
        shape = (members.get(role) or {}).get("shape")
        if not shape:
            raise ValueError(f"{kind}: wavefunction {wavefunction_ref!r} has no "
                             f"{role} member with a shape")
        rows = shape[0]
        if int(mo_index) >= rows:
            raise ValueError(f"{kind}: mo_index {mo_index} is outside the {rows} "
                             f"rows of {role} of {wavefunction_ref!r}")
        fields["mo_index"] = int(mo_index)
        if spin is not None:
            fields["spin"] = spin
        return fields

    def _series_axis(self, axis: Mapping[str, Any], n: int) -> dict[str, Any]:
        """A ``volume.series`` ``axis`` payload (``$defs/VolumeSeriesAxisPayload``,
        spec § 5.2) checked against the ``n`` slabs of the data."""
        if not isinstance(axis, Mapping):
            raise ValueError("volume.series axis must be a mapping")
        known = {"kind", "values", "units", "labels", "trajectory_ref"}
        unknown = sorted(k for k in axis if k not in known
                         and not str(k).startswith("x_"))
        if unknown:
            raise ValueError(f"volume.series axis: unknown keys {unknown}")
        kind = axis.get("kind")
        if kind not in self._SERIES_AXIS_KINDS:
            raise ValueError(f"volume.series axis kind must be one of "
                             f"{list(self._SERIES_AXIS_KINDS)}, got {kind!r}")
        if "values" not in axis:
            raise ValueError("volume.series axis: values is required")
        values = [float(v) for v in axis["values"]]
        if len(values) != n:
            raise ValueError(f"volume.series axis has {len(values)} values but "
                             f"data has {n} slabs")
        if kind == "frame":
            if not all(v.is_integer() and v >= 0 for v in values):
                raise ValueError("volume.series frame axis values must be "
                                 "0-based frame indices")
            values = [int(v) for v in values]
        out: dict[str, Any] = {"kind": kind, "values": values}
        if axis.get("units") is not None:
            out["units"] = _checked_unit_name(axis["units"],
                                              "volume.series axis units")
        if axis.get("labels") is not None:
            labels = [str(x) for x in axis["labels"]]
            if len(labels) != n:
                raise ValueError(f"volume.series axis has {len(labels)} labels "
                                 f"but data has {n} slabs")
            out["labels"] = labels
        ref = axis.get("trajectory_ref")
        if ref is not None:
            if not any(s["id"] == ref and s["kind"] in ("trajectory", "reaction.path")
                       for s in self._sections):
                raise ValueError(f"volume.series axis trajectory_ref {ref!r} names "
                                 f"no trajectory or reaction.path section added "
                                 f"so far")
            out["trajectory_ref"] = ref
        out.update({k: v for k, v in axis.items() if str(k).startswith("x_")})
        return out

    def add_density(self, data, origin, voxel_vectors, **kw) -> str:
        return self.add_volume("volume.density", data, origin, voxel_vectors, **kw)

    def add_orbital(self, data, origin, voxel_vectors, **kw) -> str:
        return self.add_volume("volume.orbital", data, origin, voxel_vectors, **kw)

    def add_basis_ao(self, ao_metadata: Mapping[str, Any], data: np.ndarray,
                    origin: Sequence[float], voxel_vectors: Sequence[Sequence[float]],
                    *, label: str | None = None, dtype: str = "float32",
                    section_id: str | None = None) -> str:
        """Add a ``basis.ao`` section: one AO evaluated on a grid + quantum
        numbers in the section-level ``ao_metadata`` peer field."""
        if section_id is not None:
            self._check_new_id(section_id)
        arr = np.asarray(data)
        if arr.ndim != 3:
            raise ValueError(f"basis.ao: data must be 3-D, got ndim={arr.ndim}")
        sid = section_id or f"ao_{len(self._sections)}"
        self._check_new_id(sid)
        grid = {
            "origin": [float(x) for x in origin],
            "voxel_vectors": [[float(x) for x in v] for v in voxel_vectors],
            "shape": list(arr.shape),
        }
        members = {
            "grid": self._json_member(f"basis_ao/{sid}_grid.json", grid),
            "data": self._binary_member(f"basis_ao/{sid}.dat", arr, dtype),
        }
        section: dict[str, Any] = {"id": self._reserve_id(sid), "kind": "basis.ao",
                                   "members": members,
                                   "ao_metadata": dict(ao_metadata)}
        if label is not None:
            section["label"] = label
        return self._push(section)

    # -- wavefunction -----------------------------------------------------

    def add_wavefunction_gto(self, shells: Sequence[Mapping[str, Any]], *,
                            mo_coefficients: np.ndarray | None = None,
                            mo_coefficients_alpha: np.ndarray | None = None,
                            mo_coefficients_beta: np.ndarray | None = None,
                            energies: Sequence[float] | None = None,
                            occupations: Sequence[float] | None = None,
                            energies_alpha: Sequence[float] | None = None,
                            occupations_alpha: Sequence[float] | None = None,
                            energies_beta: Sequence[float] | None = None,
                            occupations_beta: Sequence[float] | None = None,
                            structure_ref: str = "structure",
                            orbital_kind: str = "canonical",
                            pure: bool = True,
                            k_point: Sequence[float] | None = None,
                            coeffs_are_libint_normalized: bool = False,
                            spin: str | None = None,
                            ecp: Sequence[Mapping[str, Any]] | None = None,
                            n_electrons: int | None = None,
                            symmetry_labels: Sequence[str] | Mapping[str, Sequence[str]]
                            | None = None,
                            overlap: np.ndarray | None = None,
                            density_matrix: np.ndarray | None = None,
                            density_matrix_alpha: np.ndarray | None = None,
                            density_matrix_beta: np.ndarray | None = None,
                            complex_coefficients: bool = False,
                            label: str | None = None, section_id: str = "wf",
                            state_ref: str | None = None,
                            state_index: int | None = None,
                            nto_role: str | None = None,
                            localization_reference_basis: str | None = None) -> str:
        """Add a ``wavefunction.gto`` section (spec § 5.3 + Appendix A).

        ``shells`` are ``{"center", "l", "exponents": [...],
        "coefficients": [...], "pure"?}``. ``mo_coefficients`` is ``[n_mo, n_ao]``
        (rows are MOs); for unrestricted pass ``*_alpha`` + ``*_beta`` instead.

        If your engine stores contraction coefficients pre-multiplied by the
        primitive norm ``N_i`` (libint / libcint / GBW-style), pass
        ``coeffs_are_libint_normalized=True`` and this writer divides them out so
        they apply to ``N_i``-normalized primitives, as the format requires
        (spec Appendix A.1 -- ``N_i`` is the shell's *axial* norm, one factor
        per shell from the total ``l``, not a per-component unit norm).

        Optional, spec § 5.3:

        * ``spin``: ``"restricted"`` (the default with ``mo_coefficients``),
          ``"restricted_open"`` or ``"generalized"`` for a single set;
          ``"unrestricted"`` (the default with ``*_alpha`` + ``*_beta``). A
          ``generalized`` (two-component) wavefunction needs
          ``complex_coefficients=True`` and ``2 n_ao`` AO columns, the α
          components of the AOs first, then their β components.
        * ``complex_coefficients``: the coefficients are complex. Pass complex
          arrays (split into the ``[..., 2]`` = ``[re, im]`` form) or real
          ``[n_mo, n_ao, 2]`` arrays; ``mo_metadata`` gets ``complex: true``.
        * ``localization_reference_basis``: the non-empty minimal reference
          basis name for localized orbitals (for example ``"minao"``).
        * ``ecp``: effective core potentials, a list of ``{"center",
          "n_core_electrons", "label"?, "projectors"?: [{"l", "exponents",
          "coefficients", "powers", "local"?}]}``, at most one per atom.
        * ``n_electrons``: the number of electrons the wavefunction describes
          (Σ Z − Σ n_core_electrons − charge); the occupations must sum to it.
        * ``symmetry_labels``: one string per MO; for an unrestricted
          wavefunction a mapping ``{"alpha": [...], "beta": [...]}``.
        * ``overlap``: the AO overlap matrix ``[n_ao, n_ao]`` of the AOs as
          written (after the ``N_i`` division), so that a validator can check
          the normalization (spec Appendix A.1); producers SHOULD pass it.
        * ``density_matrix`` (with ``mo_coefficients``) or
          ``density_matrix_alpha`` + ``density_matrix_beta`` (with the α and β
          coefficients): one-particle density matrices ``[n_ao, n_ao]``.

        Everything is checked before anything is written; a problem raises
        ``ValueError``.

        Natural transition orbitals (spec § 5.3) are two sections per state,
        written with ``orbital_kind="nto"``: ``nto_role`` is ``"hole"`` or
        ``"particle"``, ``state_ref`` the id of the ``excited_states`` section
        and ``state_index`` the state's ``index`` there; ``occupations`` are
        the pair weights and ``energies`` may be left out (zeros). The writer
        then also writes ``occupation_semantics: "transition_weight"``.
        Unrestricted NTOs still use two sections: pass alpha/beta coefficients
        and occupations to each. The spin row counts may differ; weights must
        be nonnegative and descending within each spin. Pair consistency is
        checked by the validator after both sections have been written. Both
        spins together sum to one for a normalized state; full
        TDDFT/RPA weights need not. References and missing partners are checked
        by validation after all sections have been written.
        ``state_ref`` is not resolved here, because the ``excited_states``
        section is usually added after its orbitals.
        """
        self._check_new_id(section_id)
        nto_fields: dict[str, Any] = {}
        if localization_reference_basis is not None:
            nto_fields["localization_reference_basis"] = _checked_text(
                localization_reference_basis, "localization_reference_basis")
        if state_ref is not None:
            if not isinstance(state_ref, str) or not state_ref:
                raise ValueError("wavefunction.gto state_ref must be a section id")
            nto_fields["state_ref"] = state_ref
        if state_index is not None:
            if state_ref is None:
                raise ValueError("wavefunction.gto state_index needs state_ref "
                                 "(spec § 5.3)")
            nto_fields["state_index"] = _positive_int(
                state_index, "wavefunction.gto state_index")
        if nto_role is not None:
            if orbital_kind != "nto" or nto_role not in ("hole", "particle"):
                raise ValueError("wavefunction.gto nto_role is 'hole' or "
                                 "'particle', for orbital_kind 'nto' (spec § 5.3)")
            nto_fields["nto_role"] = nto_role
        if orbital_kind == "nto":
            nto_fields["occupation_semantics"] = "transition_weight"
        out_shells, n_ao = _basis_shells(shells, pure,
                                         coeffs_are_libint_normalized)

        unrestricted_given = (mo_coefficients_alpha is not None
                              and mo_coefficients_beta is not None)
        if mo_coefficients is None and not unrestricted_given:
            raise ValueError("provide mo_coefficients or both alpha and beta")
        if spin is None:
            spin = "restricted" if mo_coefficients is not None else "unrestricted"
        if spin not in _WF_SPINS:
            raise ValueError(f"spin {spin!r} is not one of {sorted(_WF_SPINS)}")
        unrestricted = spin == "unrestricted"
        if unrestricted and not unrestricted_given:
            raise ValueError("an unrestricted wavefunction needs "
                             "mo_coefficients_alpha and mo_coefficients_beta")
        if not unrestricted and mo_coefficients is None:
            raise ValueError(f"a {spin} wavefunction needs mo_coefficients")
        if spin == "generalized" and not complex_coefficients:
            raise ValueError("a generalized wavefunction needs "
                             "complex_coefficients=True (spec § 5.3)")
        complex_coefficients = bool(complex_coefficients)
        ao_extent = 2 * n_ao if spin == "generalized" else n_ao

        def numbers(values, n):
            return list(map(float, values)) if values is not None else [0.0] * n

        if not unrestricted:
            C = _wf_coefficients(mo_coefficients, "mo_coefficients",
                                 complex_coefficients, ao_extent)
            n_mo = int(C.shape[0])
            mo_meta: dict[str, Any] = {
                "n_mo": n_mo, "n_ao": n_ao, "spin": spin,
                "orbital_kind": orbital_kind,
                "energies": numbers(energies, n_mo),
                "occupations": numbers(occupations, n_mo),
            }
            if symmetry_labels is not None:
                mo_meta["symmetry_labels"] = _wf_labels(
                    symmetry_labels, "symmetry_labels", n_mo)
            coefficients = {"mo_coefficients": C}
            occupation_lists = [mo_meta["occupations"]]
        else:
            Ca = _wf_coefficients(mo_coefficients_alpha, "mo_coefficients_alpha",
                                  complex_coefficients, ao_extent)
            Cb = _wf_coefficients(mo_coefficients_beta, "mo_coefficients_beta",
                                  complex_coefficients, ao_extent)
            channels = {}
            for name, C, e, o in (("alpha", Ca, energies_alpha, occupations_alpha),
                                  ("beta", Cb, energies_beta, occupations_beta)):
                n = int(C.shape[0])
                channels[name] = {"n_mo": n, "energies": numbers(e, n),
                                  "occupations": numbers(o, n)}
            if symmetry_labels is not None:
                if not isinstance(symmetry_labels, Mapping) or (
                        set(symmetry_labels) - {"alpha", "beta"}):
                    raise ValueError("symmetry_labels of an unrestricted "
                                     "wavefunction must be a mapping "
                                     "{'alpha': [...], 'beta': [...]}")
                for name, labels in symmetry_labels.items():
                    channels[name]["symmetry_labels"] = _wf_labels(
                        labels, f"symmetry_labels[{name!r}]",
                        channels[name]["n_mo"])
            mo_meta = {"n_ao": n_ao, "spin": "unrestricted",
                       "orbital_kind": orbital_kind,
                       "alpha": channels["alpha"], "beta": channels["beta"]}
            coefficients = {"mo_coefficients_alpha": Ca,
                            "mo_coefficients_beta": Cb}
            occupation_lists = [channels["alpha"]["occupations"],
                                channels["beta"]["occupations"]]
        if k_point is not None:
            mo_meta["k_point"] = [float(x) for x in k_point]
        if complex_coefficients:
            mo_meta["complex"] = True
        mo_meta.update(nto_fields)
        if orbital_kind == "nto" and unrestricted:
            if state_ref is None or state_index is None or nto_role is None:
                raise ValueError("unrestricted NTOs require state_ref, state_index and nto_role")
            for ch, C, supplied in (("alpha", Ca, occupations_alpha), ("beta", Cb, occupations_beta)):
                block = mo_meta[ch]
                n = C.shape[0]
                if supplied is None or len(block["occupations"]) != n or len(block["energies"]) != n:
                    raise ValueError(f"NTO {ch} needs one weight and energy per coefficient row")
                weights = block["occupations"]
                if any(not math.isfinite(x) or x < 0 for x in weights) or any(a < b for a, b in zip(weights, weights[1:])):
                    raise ValueError(f"NTO {ch} weights must be finite, nonnegative and non-increasing")

        basis_json: dict[str, Any] = {"structure_ref": structure_ref, "pure": pure,
                                      "n_ao": n_ao, "shells": out_shells}
        if ecp is not None:
            basis_json["ecp"] = _wf_ecp(ecp)
        if n_electrons is not None:
            if not _is_int(n_electrons) or n_electrons < 0:
                raise ValueError("n_electrons must be an integer >= 0")
            total = math.fsum(x for occ in occupation_lists for x in occ)
            if orbital_kind != "nto" and abs(total - n_electrons) > 1e-6:
                raise ValueError(
                    f"n_electrons is {n_electrons} but the occupations sum to "
                    f"{total:.10g}; they must agree (spec § 5.3)")
            basis_json["n_electrons"] = int(n_electrons)

        matrices: dict[str, np.ndarray] = {}
        if overlap is not None:
            matrices["overlap"] = _wf_ao_matrix(overlap, "overlap", n_ao)
        if density_matrix is not None:
            if unrestricted:
                raise ValueError("density_matrix accompanies mo_coefficients; "
                                 "pass density_matrix_alpha and "
                                 "density_matrix_beta for an unrestricted "
                                 "wavefunction")
            matrices["density_matrix"] = _wf_ao_matrix(
                density_matrix, "density_matrix", n_ao)
        if density_matrix_alpha is not None or density_matrix_beta is not None:
            if not unrestricted:
                raise ValueError("density_matrix_alpha and density_matrix_beta "
                                 "accompany the alpha and beta coefficients; "
                                 "pass density_matrix for a single set")
            if density_matrix_alpha is None or density_matrix_beta is None:
                raise ValueError("pass density_matrix_alpha and "
                                 "density_matrix_beta together")
            matrices["density_matrix_alpha"] = _wf_ao_matrix(
                density_matrix_alpha, "density_matrix_alpha", n_ao)
            matrices["density_matrix_beta"] = _wf_ao_matrix(
                density_matrix_beta, "density_matrix_beta", n_ao)

        sid = self._reserve_id(section_id)
        members = {
            "basis": self._json_member(f"{sid}/basis.json", basis_json),
            "mo_metadata": self._json_member(f"{sid}/mo_metadata.json", mo_meta),
        }
        for role, C in coefficients.items():
            members[role] = self._binary_member(f"{sid}/{role}.dat", C, "float64")
        for role, matrix in matrices.items():
            members[role] = self._binary_member(f"{sid}/{role}.dat", matrix,
                                                "float64")
        return self._push({"id": sid, "kind": "wavefunction.gto",
                           "members": members}, label=label)

    def add_wavefunction_bloch(self, shells: Sequence[Mapping[str, Any]],
                               kpoints: Sequence[Sequence[float]],
                               mo_coefficients: np.ndarray,
                               energies: np.ndarray,
                               occupations: np.ndarray, *,
                               weights: Sequence[float] | None = None,
                               spin: str = "restricted",
                               orbital_kind: str = "canonical",
                               structure_ref: str = "structure",
                               coeffs_are_libint_normalized: bool = False,
                               fermi_energy_ev: float | None = None,
                               pure: bool = True,
                               labels: Sequence[str] | None = None,
                               mesh: Sequence[int] | None = None,
                               shift: Sequence[float] | None = None,
                               time_reversal: bool | None = None,
                               n_mo_per_k: Sequence[int] | None = None,
                               label: str | None = None,
                               section_id: str = "wf_bloch") -> str:
        """Add a ``wavefunction.bloch`` section (spec § 5.3.1, provisional).

        The k-resolved Bloch orbitals of a periodic system over the
        atom-centred Gaussian basis of its cell, with the Bloch sum
        ``phi_mu,k(r) = N^(-1/2) sum_T exp(+i k.T) chi_mu(r - T)``.
        ``shells`` and ``coeffs_are_libint_normalized`` are as for
        ``add_wavefunction_gto`` (Appendix A). ``kpoints`` is ``[n_k, 3]`` in
        fractional reciprocal-lattice coordinates. ``mo_coefficients`` is
        complex ``[n_k, n_spin, n_mo, n_ao]`` or real ``[n_k, n_spin, n_mo,
        n_ao, 2]`` with the last axis (real, imaginary); ``energies``
        (hartree) and ``occupations`` (electrons per orbital, not weighted by
        k-point weights) are ``[n_k, n_spin, n_mo]``. ``spin`` is
        ``"restricted"`` (n_spin 1) or ``"unrestricted"`` (n_spin 2).
        ``weights`` (one per k-point, summing to 1 within 1e-6), ``labels``,
        ``mesh``, ``shift`` and ``time_reversal`` fill the ``kpoints``
        payload; ``fermi_energy_ev`` and ``n_mo_per_k`` (the number of
        defined orbitals at each k-point, for a producer whose count varies)
        fill ``mo_metadata``. ``structure_ref`` must name a periodic
        structure. Raises ``ValueError`` before anything is written when the
        shapes disagree.
        """
        self._check_new_id(section_id)
        if spin not in ("restricted", "unrestricted"):
            raise ValueError(f"spin must be 'restricted' or 'unrestricted', "
                             f"got {spin!r}")
        n_spin = 1 if spin == "restricted" else 2
        kpts = _real_array(kpoints, "kpoints").astype(np.float64)
        if kpts.ndim != 2 or kpts.shape[1] != 3 or kpts.shape[0] < 1:
            raise ValueError(f"kpoints must be [n_k, 3] with n_k >= 1, got shape "
                             f"{list(kpts.shape)}")
        n_k = int(kpts.shape[0])
        raw = np.asarray(mo_coefficients)
        if np.iscomplexobj(raw):
            if raw.ndim != 4:
                raise ValueError("complex mo_coefficients must be [n_k, n_spin, "
                                 f"n_mo, n_ao], got shape {list(raw.shape)}")
            coeff = np.stack([raw.real, raw.imag], axis=-1).astype(np.float64)
        else:
            coeff = _real_array(raw, "mo_coefficients").astype(np.float64)
            if coeff.ndim != 5 or coeff.shape[4] != 2:
                raise ValueError("real mo_coefficients must be [n_k, n_spin, "
                                 f"n_mo, n_ao, 2], got shape {list(coeff.shape)}")
        if coeff.shape[0] != n_k or coeff.shape[1] != n_spin:
            raise ValueError(f"mo_coefficients {list(coeff.shape)} must have "
                             f"n_k = {n_k} and n_spin = {n_spin} ({spin})")
        n_mo, n_ao_coeff = int(coeff.shape[2]), int(coeff.shape[3])
        eps = _real_array(energies, "energies").astype(np.float64)
        occ = _real_array(occupations, "occupations").astype(np.float64)
        for name, arr in (("energies", eps), ("occupations", occ)):
            if arr.shape != coeff.shape[:3]:
                raise ValueError(f"{name} {list(arr.shape)} must be [n_k, n_spin, "
                                 f"n_mo] = {list(coeff.shape[:3])}")
        out_shells, n_ao = _basis_shells(shells, pure,
                                         coeffs_are_libint_normalized)
        if n_ao != n_ao_coeff:
            raise ValueError(f"the basis has {n_ao} AOs but mo_coefficients has "
                             f"n_ao = {n_ao_coeff}")
        kp: dict[str, Any] = {"kpoints": kpts.tolist()}
        if weights is not None:
            w = [float(x) for x in weights]
            if len(w) != n_k or any(x < 0 for x in w) or \
                    abs(math.fsum(w) - 1.0) > 1e-6:
                raise ValueError("weights must be n_k non-negative numbers that "
                                 "sum to 1 within 1e-6")
            kp["weights"] = w
        if labels is not None:
            if len(labels) != n_k:
                raise ValueError(f"labels must have {n_k} entries")
            kp["labels"] = [str(x) for x in labels]
        if mesh is not None:
            if len(mesh) != 3 or any(int(n) < 1 for n in mesh):
                raise ValueError("mesh must be three integers >= 1")
            kp["mesh"] = [int(n) for n in mesh]
        if shift is not None:
            if len(shift) != 3:
                raise ValueError("shift must have three entries")
            kp["shift"] = [float(x) for x in shift]
        if time_reversal is not None:
            kp["time_reversal"] = bool(time_reversal)
        meta: dict[str, Any] = {
            "spin": spin, "n_spin": n_spin, "n_mo": n_mo, "n_ao": n_ao,
            "orbital_kind": orbital_kind, "complex": True,
            "k_point_convention": "fractional",
        }
        if n_mo_per_k is not None:
            per_k = [int(n) for n in n_mo_per_k]
            if len(per_k) != n_k or any(n < 0 or n > n_mo for n in per_k):
                raise ValueError(f"n_mo_per_k must have {n_k} entries between 0 "
                                 f"and n_mo = {n_mo}")
            meta["n_mo_per_k"] = per_k
        if fermi_energy_ev is not None:
            meta["fermi_energy"] = {"value": float(fermi_energy_ev), "units": "ev"}
        sid = self._reserve_id(section_id)
        members = {
            "basis": self._json_member(f"{sid}/basis.json", {
                "structure_ref": structure_ref, "pure": pure, "n_ao": n_ao,
                "shells": out_shells}),
            "kpoints": self._json_member(f"{sid}/kpoints.json", kp),
            "mo_metadata": self._json_member(f"{sid}/mo_metadata.json", meta),
            "mo_coefficients": self._binary_member(
                f"{sid}/mo_coefficients.bin", coeff, "float64"),
            "energies": self._binary_member(f"{sid}/energies.bin", eps,
                                            "float64"),
            "occupations": self._binary_member(f"{sid}/occupations.bin", occ,
                                               "float64"),
        }
        return self._push({"id": sid, "kind": "wavefunction.bloch",
                           "members": members}, label=label)

    # -- atom properties --------------------------------------------------

    def add_atom_properties(
            self, *,
            mulliken_charge: Sequence[float] | None = None,
            loewdin_charge: Sequence[float] | None = None,
            hirshfeld_charge: Sequence[float] | None = None,
            iao_charge: Sequence[float] | None = None,
            spin_population: Sequence[float] | None = None,
            cm5_charge: Sequence[float] | None = None,
            chelpg_charge: Sequence[float] | None = None,
            resp_charge: Sequence[float] | None = None,
            esp_charge: Sequence[float] | None = None,
            nbo_charge: Sequence[float] | None = None,
            mbis_charge: Sequence[float] | None = None,
            bader_charge: Sequence[float] | None = None,
            ddec_charge: Sequence[float] | None = None,
            valence: Sequence[float] | None = None,
            bonded_valence: Sequence[float] | None = None,
            fukui_plus: Sequence[float] | None = None,
            fukui_minus: Sequence[float] | None = None,
            fukui_zero: Sequence[float] | None = None,
            dual_descriptor: Sequence[float] | None = None,
            atomic_energy_eh: Sequence[float] | None = None,
            atomic_volume_bohr3: Sequence[float] | None = None,
            nmr_shielding_iso_ppm: Sequence[float] | None = None,
            magnetic_moment: Sequence[Sequence[float]] | None = None,
            force_eh_bohr: Sequence[Sequence[float]] | None = None,
            custom: Mapping[str, Any] | None = None,
            meta: Mapping[str, Mapping[str, Any]] | None = None,
            label: str | None = None, section_id: str = "atom_properties") -> str:
        """Add an ``atom_properties`` section (spec § 5.7).

        Pass one or more of the schema's roles by keyword. The scalar roles
        (``mulliken_charge`` ... ``nmr_shielding_iso_ppm``) are float64
        ``[n_atoms]``; the vector roles ``magnetic_moment`` and
        ``force_eh_bohr`` are float64 ``[n_atoms, 3]``. All roles of one
        section must have the same ``n_atoms``. The units and conventions are
        in the table of spec § 5.7; ``force_eh_bohr`` is the energy gradient,
        not the force.

        ``custom`` adds roles outside the list: ``{name: (array, entry)}``
        writes the role ``custom.<name>`` (``name`` is ``[A-Za-z0-9_]+``; a
        leading ``custom.`` is accepted), a scalar for a 1-D array and a vector
        for an ``[n_atoms, 3]`` one. ``entry`` is its ``meta`` description,
        ``{"units": ..., "description": ...}`` and optionally ``"kind"``
        (``"scalar"`` or ``"vector"``, which must agree with the array shape;
        it is derived from the shape when absent) and ``x_*`` keys. Instead of
        a tuple, a name may map to a bare array whose entry comes from
        ``meta={name: entry}``. Every custom role needs exactly one entry and
        ``meta`` may not name a role that ``custom`` lacks. The writer emits
        the ``meta`` member (``$defs/AtomPropertiesMetaPayload``) itself.
        """
        self._check_new_id(section_id)
        given = {
            "mulliken_charge": mulliken_charge,
            "loewdin_charge": loewdin_charge,
            "hirshfeld_charge": hirshfeld_charge,
            "iao_charge": iao_charge,
            "spin_population": spin_population,
            "cm5_charge": cm5_charge,
            "chelpg_charge": chelpg_charge,
            "resp_charge": resp_charge,
            "esp_charge": esp_charge,
            "nbo_charge": nbo_charge,
            "mbis_charge": mbis_charge,
            "bader_charge": bader_charge,
            "ddec_charge": ddec_charge,
            "valence": valence,
            "bonded_valence": bonded_valence,
            "fukui_plus": fukui_plus,
            "fukui_minus": fukui_minus,
            "fukui_zero": fukui_zero,
            "dual_descriptor": dual_descriptor,
            "atomic_energy_eh": atomic_energy_eh,
            "atomic_volume_bohr3": atomic_volume_bohr3,
            "nmr_shielding_iso_ppm": nmr_shielding_iso_ppm,
            "magnetic_moment": magnetic_moment,
            "force_eh_bohr": force_eh_bohr,
        }
        assert tuple(given) == _ATOM_PROPERTY_ROLES
        arrays: dict[str, np.ndarray] = {}
        for role in _ATOM_PROPERTY_ROLES:
            if given[role] is None:
                continue
            arrays[role] = _atom_property_array(
                role, given[role], vector=role in _ATOM_VECTOR_ROLES)
        meta_doc: dict[str, dict[str, Any]] = {}
        for role, (data, entry) in _custom_atom_roles(custom, meta).items():
            arr = _real_array(data, f"atom_properties {role}")
            arrays[role] = _atom_property_array(role, arr, vector=arr.ndim == 2)
            meta_doc[role] = _custom_meta_entry(role, entry, arr.ndim == 2)
        if not arrays:
            raise ValueError("atom_properties needs at least one array")
        lengths = {role: arr.shape[0] for role, arr in arrays.items()}
        if len(set(lengths.values())) != 1:
            raise ValueError(f"atom_properties roles differ in length: {lengths}")
        members = {role: self._binary_member(f"atom_properties/{role}.bin",
                                             arr, "float64")
                   for role, arr in arrays.items()}
        if meta_doc:
            members["meta"] = self._json_member("atom_properties/meta.json",
                                                meta_doc)
        return self._push({"id": self._reserve_id(section_id),
                          "kind": "atom_properties", "members": members}, label=label)

    # -- trajectories / reactions / scans --------------------------------

    def add_trajectory(self, metadata: Mapping[str, Any], coords: np.ndarray, *,
                      lattice: Any = None,
                      velocities: np.ndarray | None = None,
                      gradient: np.ndarray | None = None,
                      times: Sequence[float] | None = None,
                      energies: Sequence[float] | None = None,
                      cell_volumes: Sequence[float] | None = None,
                      temperatures: Sequence[float] | None = None,
                      precision: str | None = None,
                      units: Mapping[str, str] | None = None,
                      label: str | None = None, section_id: str = "trajectory") -> str:
        """Add a ``trajectory`` section (spec § 5.5). ``coords`` is
        [n_frames, n_atoms, 3] in Å.

        ``metadata`` names the atoms (``symbols``, ``atomic_numbers`` or
        ``atoms``) and may carry ``energies``, ``kind`` (``"md"`` for
        molecular dynamics), ``pbc``, ``ensemble`` (``"nve"``, ``"nvt"``,
        ``"npt"``, ``"nph"``, ``"custom"``) and ``timestep_fs``.

        The optional per-frame members, each in frame order:

        * ``lattice`` — the cell in the **structure convention**: rows in Å,
          as ``add_structure(lattice_vectors=...)`` writes them, and *not* the
          columns in bohr of ``add_reaction_path`` (spec § 4.3, "Lattice
          conventions"). A :class:`Lattice` writes ``[3, 3]`` (a fixed cell),
          a list/tuple of them ``[n_frames, 3, 3]`` (one cell per frame);
          an array is written as given and must already be rows in Å.
          ``metadata["pbc"]`` is required with it, and a periodic ``pbc``
          entry requires it.
        * ``velocities`` — [n_frames, n_atoms, 3] in Å/fs.
        * ``gradient`` — [n_frames, n_atoms, 3], the energy gradient dE/dR in
          hartree/bohr (the forces are its negative), as ``add_forces``.
        * ``times`` [n_frames] in fs, ``energies`` [n_frames] in hartree (the
          binary form of ``metadata["energies"]``, which it must match within
          1e-10 hartree when both are given), ``cell_volumes`` [n_frames] in
          Å³ and ``temperatures`` [n_frames] in kelvin.

        ``precision="float32"`` writes ``coords`` as float32 and sets
        ``metadata["precision"]``; ``"float64"`` (the default) writes float64.
        A ``precision`` in ``metadata`` is honoured the same way and must not
        contradict the argument.

        ``units`` (optional) is written as the metadata's ``units``
        declaration, e.g. ``{"energies": "hartree", "coords": "angstrom"}``:
        unit names of spec § 4.3, keyed by a metadata field or a member role.
        It declares and never converts, so a key whose unit the spec fixes
        (``energies``, ``coords``, ``lattice``, ``times``, ``timestep_fs``,
        ``cell_volumes``, ``temperatures``) can only restate it; an unknown
        unit name raises ``ValueError``, and so does a ``units`` given both
        here and in ``metadata``. Every check runs before anything is
        written.
        """
        self._check_new_id(section_id)
        arr = np.asarray(coords, dtype=np.float64)
        if arr.ndim != 3 or arr.shape[2] != 3:
            raise ValueError("trajectory coords must be [n_frames, n_atoms, 3]")
        meta = dict(metadata)
        declared = meta.get("precision")
        for value, what in ((precision, "precision"),
                            (declared, "metadata precision")):
            if value is not None and value not in ("float64", "float32"):
                raise ValueError(f"trajectory {what} must be 'float64' or "
                                 f"'float32', got {value!r}")
        if precision is not None and declared is not None and declared != precision:
            raise ValueError(f"trajectory precision {precision!r} contradicts "
                             f"metadata precision {declared!r}")
        if precision is not None:
            meta["precision"] = precision
        dtype = meta.get("precision", "float64")
        lat = None
        if lattice is not None:
            lat = _trajectory_lattice(lattice)
            _check_lattice_frames(lat, arr.shape[0], "trajectory")
        pbc = meta.get("pbc")
        if pbc is not None and (not isinstance(pbc, (list, tuple)) or len(pbc) != 3
                                or not all(isinstance(x, bool) for x in pbc)):
            raise ValueError(f"trajectory metadata pbc must be three booleans, "
                             f"got {pbc!r}")
        if lat is not None and pbc is None:
            raise ValueError("trajectory lattice needs metadata['pbc']: pbc is "
                             "REQUIRED when lattice is present (spec § 5.5)")
        if lat is None and pbc is not None and any(pbc):
            raise ValueError(f"trajectory metadata pbc {list(pbc)} marks a "
                             f"periodic axis, so lattice is required (spec § 5.5)")
        frames = _trajectory_frame_arrays(
            "trajectory", arr,
            {"velocities": velocities, "gradient": gradient, "times": times,
             "energies": energies, "cell_volumes": cell_volumes,
             "temperatures": temperatures}, meta)
        meta = _with_units(meta, units, "trajectory", "TrajectoryMetadataPayload")
        _json_bytes(meta)   # refuse a payload that is not JSON before writing
        sid = self._reserve_id(section_id)
        members = {
            "metadata": self._json_member(f"{sid}/metadata.json", meta),
            "coords": self._binary_member(f"{sid}/coords.bin", arr, dtype),
        }
        if lat is not None:
            members["lattice"] = self._binary_member(f"{sid}/lattice.bin", lat,
                                                     "float64")
        for role, data in frames:
            members[role] = self._binary_member(f"{sid}/{role}.bin", data,
                                                "float64")
        return self._push({"id": sid, "kind": "trajectory", "members": members},
                          label=label)

    def add_reaction_path(self, metadata: Mapping[str, Any], coords: np.ndarray, *,
                         lattice: np.ndarray | Lattice | Sequence[Lattice] | None = None,
                         velocities: np.ndarray | None = None,
                         gradient: np.ndarray | None = None,
                         times: Sequence[float] | None = None,
                         energies: Sequence[float] | None = None,
                         label: str | None = None,
                         section_id: str = "reaction_path") -> str:
        """Add a ``reaction.path`` section (trajectory layout + waypoint
        annotations carried in ``metadata``).

        ``lattice`` (optional) marks the path **periodic**: float64 with shape
        ``[3, 3]`` (fixed cell) or ``[n_frames, 3, 3]`` (variable cell), columns
        = a, b, c, in **bohr**. Its presence — not the manifest version — is how
        a consumer detects a periodic reaction path (spec § 5.5). A path without
        ``lattice`` is unambiguously molecular.

        This is *not* the convention of ``add_structure(lattice_vectors=...)``
        or of ``add_trajectory(lattice=...)`` (rows, Å); see spec § 4.3,
        "Lattice conventions". Pass a :class:`Lattice` (fixed cell) or a
        list/tuple of them (one per frame) and it is converted to columns in
        bohr for you; a plain array is written as given and must already be
        columns in bohr.

        ``velocities`` (Å/fs) and ``gradient`` (the energy gradient in
        hartree/bohr) are [n_frames, n_atoms, 3]; ``times`` (fs) and
        ``energies`` (hartree, matching ``metadata["energies"]`` within 1e-10
        when both are given) are [n_frames], as for :meth:`add_trajectory`.
        The coordinates are written as float64.
        """
        self._check_new_id(section_id)
        arr = np.asarray(coords, dtype=np.float64)
        if arr.ndim != 3 or arr.shape[2] != 3:
            raise ValueError("reaction.path coords must be [n_frames, n_atoms, 3]")
        meta = dict(metadata)
        if meta.get("precision", "float64") != "float64":
            raise ValueError("add_reaction_path writes float64 coords, so "
                             "metadata precision must be absent or 'float64'")
        # Convert and check everything before any member is registered, so a
        # rejected call leaves the writer untouched.
        lat = None
        if lattice is not None:
            lat = _reaction_path_lattice(lattice)
            _check_lattice_frames(lat, arr.shape[0], "reaction.path")
        frames = _trajectory_frame_arrays(
            "reaction.path", arr,
            {"velocities": velocities, "gradient": gradient, "times": times,
             "energies": energies}, meta)
        meta = _with_units(meta, None, "reaction.path",
                           "ReactionPathMetadataPayload")
        _json_bytes(meta)
        sid = self._reserve_id(section_id)
        members = {
            "metadata": self._json_member(f"{sid}/metadata.json", meta),
            "coords": self._binary_member(f"{sid}/coords.bin", arr, "float64"),
        }
        if lat is not None:
            members["lattice"] = self._binary_member(
                f"{sid}/lattice.bin", lat, "float64")
        for role, data in frames:
            members[role] = self._binary_member(f"{sid}/{role}.bin", data,
                                                "float64")
        return self._push({"id": sid, "kind": "reaction.path",
                           "members": members}, label=label)

    def add_reaction_waypoints(self, trajectory_ref: str,
                              waypoints: Mapping[str, Any], *,
                              label: str | None = None,
                              section_id: str = "reaction_waypoints") -> str:
        """Add a ``reaction.waypoints`` annotation over an existing trajectory."""
        self._check_new_id(section_id)
        if trajectory_ref not in self._ids:
            raise ValueError(f"trajectory_ref {trajectory_ref!r} not yet added")
        member = self._json_member(f"{section_id}/waypoints.json", dict(waypoints))
        return self._push({"id": self._reserve_id(section_id),
                          "kind": "reaction.waypoints",
                          "trajectory_ref": trajectory_ref,
                          "members": {"waypoints": member}}, label=label)

    def add_reaction_scheme(self, steps: Sequence[Mapping[str, Any]], *,
                            label: str | None = None,
                            section_id: str = "reaction_scheme") -> str:
        """Add a ``reaction.scheme`` section (spec § 5.14).

        ``steps`` lists the steps in reading order. Each is a mapping with
        ``id`` (unique among the steps), ``reactants`` and ``products`` (at
        least one each) and optionally ``label``, ``reagents``, ``arrow``
        (``"forward"``, the default, ``"equilibrium"``, ``"retrosynthetic"``
        or ``"resonance"``), ``conditions`` (the text lines at the arrow),
        ``yield_percent`` (0 to 100), ``yield_text`` (as printed) and
        ``energies`` (the step's reaction energies with their evidence,
        ``$defs/ReactionStepEnergies``, spec § 5.14), plus keys beginning
        ``x_``. A participant is the id of a ``structure``
        section already added, or ``{"structure_ref": id, "coefficient":
        n}``; within one role of a step a structure appears once.

        Roles are chemical: ``reactants`` are consumed and ``products``
        formed in the forward direction whatever the arrow, so the target
        of a retrosynthetic step is among its ``products``, and an
        ``energies`` term of a reactant has a negative coefficient. The
        energies are written as given and checked as a validator checks
        them: the terms match the participants, the deltas their sums, and
        a ``run_ref`` names an added ``run.record``. ``label`` titles the
        scheme. A wrong step raises ``ValueError`` before anything is added.
        """
        self._check_new_id(section_id)
        if isinstance(steps, (str, Mapping)) or not isinstance(steps, Iterable):
            raise ValueError("reaction.scheme steps must be a list of steps")
        out: list[dict[str, Any]] = []
        seen: set[str] = set()
        for k, step in enumerate(steps):
            what = f"reaction.scheme step {k}"
            _check_keys(step, ("id", "label", "reactants", "products",
                               "reagents", "arrow", "conditions",
                               "yield_percent", "yield_text", "energies"), what)
            sid = step.get("id")
            if not isinstance(sid, str) or not _ID_RE.match(sid):
                raise ValueError(f"{what}: id must match [A-Za-z0-9_.-]+, "
                                 f"got {sid!r}")
            if sid in seen:
                raise ValueError(f"{what}: id {sid!r} is used by another step")
            seen.add(sid)
            entry: dict[str, Any] = {"id": sid}
            if step.get("label") is not None:
                entry["label"] = _checked_text(step["label"], f"{what} label")
            for role in ("reactants", "products", "reagents"):
                if role == "reagents" and step.get(role) is None:
                    continue
                entry[role] = _participants(step.get(role), f"{what} {role}",
                                            required=role != "reagents")
                for participant in entry[role]:
                    self._check_ref(f"{what} {role} structure_ref",
                                    participant["structure_ref"], "structure")
            if step.get("arrow") is not None:
                entry["arrow"] = _choice(step["arrow"], _ARROWS, f"{what} arrow")
            if step.get("conditions") is not None:
                lines = step["conditions"]
                if isinstance(lines, str) or not all(isinstance(t, str)
                                                     for t in lines):
                    raise ValueError(f"{what}: conditions must be a list of "
                                     f"strings, one per line")
                entry["conditions"] = list(lines)
            if step.get("yield_percent") is not None:
                y = step["yield_percent"]
                _checked_number(y, f"{what} yield_percent")
                if not 0 <= y <= 100:
                    raise ValueError(f"{what}: yield_percent must be from 0 to "
                                     f"100, got {y!r}")
                entry["yield_percent"] = int(y) if _is_int(y) else float(y)
            if step.get("yield_text") is not None:
                entry["yield_text"] = _checked_text(step["yield_text"],
                                                    f"{what} yield_text")
            if step.get("energies") is not None:
                entry["energies"] = self._checked_step_energies(
                    entry, step["energies"], what)
            entry.update({k: v for k, v in step.items() if str(k).startswith("x_")})
            out.append(entry)
        if not out:
            raise ValueError("reaction.scheme needs at least one step")
        sid = self._reserve_id(section_id)
        member = self._json_member(f"{sid}/scheme.json", {"steps": out})
        return self._push({"id": sid, "kind": "reaction.scheme",
                           "members": {"scheme": member}}, label=label)

    def _checked_step_energies(self, entry: dict[str, Any], energies: Any,
                               what: str) -> dict[str, Any]:
        """A step's ``energies`` (spec § 5.14), checked against the validator's
        rules and, when ``jsonschema`` is installed, against
        ``$defs/ReactionStepEnergies``."""
        from .validate import _load_payload_validator, _reaction_energies_errors
        if not isinstance(energies, Mapping):
            raise ValueError(f"{what}: energies must be a mapping")
        payload = json.loads(json.dumps(energies, allow_nan=False))
        missing = [k for k in ("version", "deltas", "terms") if k not in payload]
        if missing:
            raise ValueError(f"{what}: energies needs {', '.join(missing)} "
                             f"(spec § 5.14)")
        validator = _load_payload_validator("ReactionStepEnergies")
        if validator is not None:
            error = next(iter(validator.iter_errors(payload)), None)
            if error is not None:
                where = "/".join(str(p) for p in error.absolute_path)
                raise ValueError(f"{what} energies{'/' + where if where else ''}: "
                                 f"{error.message} (spec § 5.14)")
        kinds = {s["id"]: s["kind"] for s in self._sections}
        problems = _reaction_energies_errors(what, {**entry, "energies": payload}, kinds)
        if problems:
            raise ValueError(problems[0][1])
        return payload

    def add_scan_surface(self, metadata: Mapping[str, Any], axis_a: Sequence[float],
                        axis_b: Sequence[float], energies: np.ndarray, *,
                        geometries: np.ndarray | None = None,
                        units: Mapping[str, str] | None = None,
                        label: str | None = None,
                        section_id: str = "scan_surface") -> str:
        """Add a ``scan.surface`` 2-D relaxed-scan energy surface.

        ``units`` (optional) is written as the metadata's ``units``
        declaration, e.g. ``{"axis_a": "bohr", "axis_b": "radian",
        "energies": "hartree"}``: unit names of spec § 4.3, keyed by a member
        role or a metadata field. It declares the units the values are stored
        in and converts nothing; an unknown unit name raises ``ValueError``,
        and so does a ``units`` given both here and in ``metadata``."""
        self._check_new_id(section_id)
        meta = _with_units(dict(metadata), units, "scan.surface",
                           "ScanSurfaceMetadataPayload")
        members = {
            "metadata": self._json_member(f"{section_id}/metadata.json", meta),
            "axis_a": self._binary_member(f"{section_id}/axis_a.bin",
                                          np.asarray(axis_a), "float64"),
            "axis_b": self._binary_member(f"{section_id}/axis_b.bin",
                                          np.asarray(axis_b), "float64"),
            "energies": self._binary_member(f"{section_id}/energies.bin",
                                            np.asarray(energies), "float64"),
        }
        if geometries is not None:
            members["geometries"] = self._binary_member(
                f"{section_id}/geometries.bin", np.asarray(geometries), "float64")
        return self._push({"id": self._reserve_id(section_id), "kind": "scan.surface",
                          "members": members}, label=label)

    # -- vibrations / spectra --------------------------------------------

    def add_vibrations(self, metadata: Mapping[str, Any], displacements: np.ndarray,
                      *, hessian_ref: str | None = None,
                      label: str | None = None, section_id: str = "vibrations") -> str:
        """Add a ``vibrations`` section. ``displacements`` is float64
        [n_modes, n_atoms, 3]; ``metadata.frequencies`` in cm⁻¹.

        ``hessian_ref`` (optional) is the id of the ``hessian`` section the
        modes were derived from (spec § 5.5); it must already be added."""
        self._check_new_id(section_id)
        arr = np.asarray(displacements, dtype=np.float64)
        if arr.ndim != 3 or arr.shape[2] != 3:
            raise ValueError("displacements must be [n_modes, n_atoms, 3]")
        if hessian_ref is not None and not any(
                s["id"] == hessian_ref and s["kind"] == "hessian"
                for s in self._sections):
            raise ValueError(f"vibrations hessian_ref {hessian_ref!r} names no "
                             f"hessian section added so far")
        members = {
            "metadata": self._json_member(f"{section_id}/metadata.json",
                                          dict(metadata)),
            "displacements": self._binary_member(
                f"{section_id}/displacements.bin", arr, "float64"),
        }
        section: dict[str, Any] = {"id": self._reserve_id(section_id),
                                   "kind": "vibrations"}
        if hessian_ref is not None:
            section["hessian_ref"] = hessian_ref
        section["members"] = members
        return self._push(section, label=label)

    # -- energy derivatives ----------------------------------------------

    def add_forces(self, gradient: np.ndarray, *,
                   stress: np.ndarray | None = None,
                   meta: Mapping[str, Any] | None = None,
                   label: str | None = None, section_id: str = "forces") -> str:
        """Add a ``forces`` section (spec § 5.10).

        ``gradient`` is the energy gradient dE/dR, float64 ``[n_atoms, 3]``
        in hartree/bohr, with the rows in the atom order of the structure;
        the forces are its negative, so pass the gradient, not the forces.
        ``stress`` (optional, periodic systems) is the ``[3, 3]`` stress
        tensor (1/V) dE/d(strain) in hartree/bohr³. ``meta`` (optional) is the
        ``$defs/ForcesMetaPayload`` object: ``structure_ref``, ``method``,
        ``convention`` (``"gradient"``) and ``units``
        (``{"energy": "hartree", "length": "bohr"}``). A wrong shape or a
        non-real array raises ``ValueError`` before anything is written.
        """
        self._check_new_id(section_id)
        grad = _real_array(gradient, "forces gradient")
        if grad.ndim != 2 or grad.shape[1] != 3:
            raise ValueError("forces gradient must be [n_atoms, 3], got shape "
                             f"{list(grad.shape)}")
        st = None
        if stress is not None:
            st = _real_array(stress, "forces stress")
            if st.shape != (3, 3):
                raise ValueError("forces stress must be [3, 3], got shape "
                                 f"{list(st.shape)}")
        if meta is not None and not isinstance(meta, Mapping):
            raise ValueError("forces meta must be a mapping")
        sid = self._reserve_id(section_id)
        members: dict[str, Any] = {
            "gradient": self._binary_member(f"{sid}/gradient.bin", grad, "float64")}
        if st is not None:
            members["stress"] = self._binary_member(f"{sid}/stress.bin", st,
                                                    "float64")
        if meta is not None:
            members["meta"] = self._json_member(f"{sid}/meta.json", dict(meta))
        return self._push({"id": sid, "kind": "forces", "members": members},
                          label=label)

    def add_hessian(self, matrix: np.ndarray, *,
                    meta: Mapping[str, Any] | None = None,
                    label: str | None = None, section_id: str = "hessian") -> str:
        """Add a ``hessian`` section (spec § 5.11).

        ``matrix`` is the Cartesian Hessian d2E/(dx_i dx_j), float64
        ``[3N, 3N]`` in hartree/bohr², index ``3 * atom + xyz``; write a
        symmetric matrix ((H + H.T) / 2 for a finite-difference one). ``meta``
        (optional) is the ``$defs/HessianMetaPayload`` object:
        ``structure_ref``, ``mass_weighted``, ``projected``,
        ``atom_masses_amu`` (N entries; required when ``mass_weighted``) and
        ``units``. A non-square matrix, an extent not divisible by 3, or masses
        that do not match raise ``ValueError`` before anything is written.
        """
        self._check_new_id(section_id)
        mat = _real_array(matrix, "hessian matrix")
        if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
            raise ValueError("hessian matrix must be square [3N, 3N], got "
                             f"shape {list(mat.shape)}")
        if mat.shape[0] % 3:
            raise ValueError(f"hessian matrix extent {mat.shape[0]} is not "
                             f"divisible by 3")
        if meta is not None:
            if not isinstance(meta, Mapping):
                raise ValueError("hessian meta must be a mapping")
            masses = meta.get("atom_masses_amu")
            if meta.get("mass_weighted") is True and masses is None:
                raise ValueError("hessian meta: a mass-weighted matrix needs "
                                 "atom_masses_amu")
            if masses is not None and len(masses) != mat.shape[0] // 3:
                raise ValueError(
                    f"hessian meta atom_masses_amu has {len(masses)} entries "
                    f"but the matrix is for {mat.shape[0] // 3} atoms")
        sid = self._reserve_id(section_id)
        members: dict[str, Any] = {
            "matrix": self._binary_member(f"{sid}/matrix.bin", mat, "float64")}
        if meta is not None:
            members["meta"] = self._json_member(f"{sid}/meta.json", dict(meta))
        return self._push({"id": sid, "kind": "hessian", "members": members},
                          label=label)

    # -- excited states -----------------------------------------------------

    def add_excited_states(self, states: Sequence[Mapping[str, Any]], *,
                           meta: Mapping[str, Any] | None = None,
                           wavefunction_ref: str | None = None,
                           label: str | None = None,
                           section_id: str = "excited_states") -> str:
        """Add an ``excited_states`` section (spec § 5.12).

        ``states`` lists the states in ascending order of ``index`` (an int
        >= 1, unique), each a mapping with the fields of
        ``$defs/ExcitedState``: ``index`` and ``excitation_energy`` are
        required; ``total_energy``, ``multiplicity``, ``symmetry``,
        ``oscillator_strength``, ``rotatory_strength_cgs`` (10^-40 esu² cm²),
        ``transition_dipole`` and ``transition_magnetic_dipole`` (three
        numbers each, atomic units), ``contributions`` (``{from_mo, to_mo,
        spin?, weight}``), ``nto_ref``, ``density_ref`` and ``x_`` keys are
        optional. ``excitation_energy`` and ``total_energy`` are ``{value,
        units}`` objects with a unit name, or numbers, which are written in
        ``ev`` and ``hartree`` respectively. ``nto_ref`` must name a
        ``wavefunction.gto`` section and ``density_ref`` a
        ``volume.transition_density`` or ``volume.difference`` section, both
        already added.

        ``meta`` (optional) is the ``$defs/ExcitedStatesMeta`` object
        (``method``, ``reference_state``, ``n_states``, ``tda``,
        ``spin_flip``, ``units``); ``n_states`` must equal the number of
        states. ``wavefunction_ref`` (optional) names the ``wavefunction.gto``
        section of the reference state, already added; every contribution's
        ``from_mo`` and ``to_mo`` must then be one of its orbital rows.
        Anything else raises ``ValueError`` before anything is written.
        """
        self._check_new_id(section_id)
        by_id = {s["id"]: s for s in self._sections}

        def check_ref(what: str, ref: Any, kinds: tuple[str, ...]) -> None:
            target = by_id.get(ref) if isinstance(ref, str) else None
            if target is None or target["kind"] not in kinds:
                found = ("no such section" if target is None
                         else f"a {target['kind']!r} section")
                raise ValueError(f"{what} {ref!r} must name an added "
                                 f"{' or '.join(kinds)} section, found {found} "
                                 f"(spec § 5.12)")

        rows: dict[str, int] | None = None
        if wavefunction_ref is not None:
            check_ref("excited_states: wavefunction_ref", wavefunction_ref,
                      ("wavefunction.gto",))
            wf_members = by_id[wavefunction_ref]["members"]
            rows = {key: int(wf_members[role]["shape"][0])
                    for role, key in (("mo_coefficients", ""),
                                      ("mo_coefficients_alpha", "alpha"),
                                      ("mo_coefficients_beta", "beta"))
                    if role in wf_members}
        if isinstance(states, (str, bytes, Mapping)):
            raise ValueError("excited_states: states must be a sequence of "
                             "state mappings")
        out_states = [_excited_state(state, k, rows, check_ref)
                      for k, state in enumerate(states)]
        for k in range(1, len(out_states)):
            if out_states[k]["index"] <= out_states[k - 1]["index"]:
                raise ValueError(
                    f"excited_states: states[{k}].index "
                    f"{out_states[k]['index']} follows index "
                    f"{out_states[k - 1]['index']}; indices are unique and "
                    f"ascending (spec § 5.12)")
        payload: dict[str, Any] = {"states": out_states}
        if meta is not None:
            if not isinstance(meta, Mapping):
                raise ValueError("excited_states: meta must be a mapping")
            body = _with_units(dict(meta), None, "excited_states meta",
                               "ExcitedStatesMeta")
            n_states = body.get("n_states")
            if n_states is not None and n_states != len(out_states):
                raise ValueError(f"excited_states: meta.n_states {n_states!r} "
                                 f"but {len(out_states)} states (spec § 5.12)")
            payload["meta"] = body
        sid = self._reserve_id(section_id)
        section: dict[str, Any] = {"id": sid, "kind": "excited_states"}
        if wavefunction_ref is not None:
            section["wavefunction_ref"] = wavefunction_ref
        section["members"] = {"states": self._json_member(f"{sid}/states.json",
                                                          payload)}
        return self._push(section, label=label)

    # -- property kinds (spec § 5.13) ------------------------------------

    def _add_properties(self, kind: str, section_id: str,
                        payload: dict[str, Any], *, label: str | None = None) -> str:
        sid = self._reserve_id(section_id)
        member = self._json_member(f"{sid}/properties.json", payload)
        return self._push({"id": sid, "kind": kind,
                           "members": {"properties": member}}, label=label)

    def add_electric_properties(self, properties: Mapping[str, Any],
                                label: str | None = None,
                                section_id: str = "electric") -> str:
        """Add a ``properties.electric`` section (spec § 5.13,
        ``$defs/ElectricPropertiesPayload``).

        ``properties`` is the payload; at least one of its blocks is required.
        Arrays may be NumPy arrays; they are written as nested lists. Blocks
        and units (nothing is converted):

        * ``origin``: a conventional name or a 3-vector in bohr;
        * ``dipole``: ``{"vector": [x, y, z], "total"?, "units"?}`` in debye;
        * ``quadrupole``: ``{"matrix": 3x3 symmetric, "traceless": bool,
          "convention": "buckingham" | "cartesian"}`` in debye·Å;
        * ``octupole``: ``{"tensor": [3, 3, 3] or 27 numbers, "convention"?}``
          in debye·Å²;
        * ``polarizability``: ``{"matrix": 3x3, "frequency_ev"?, "isotropic"?}``
          in atomic units;
        * ``hyperpolarizability``: ``{"tensor": [3, 3, 3] or 27 numbers,
          "frequencies_ev"?: [w1, w2]}`` in atomic units;
        * ``units``: restates ``dipole`` (debye) and ``origin`` (bohr).

        A dipole that the root ``dipole_moment`` also states must agree with
        it within 1e-6 debye (the validator checks, so the order of the two
        calls does not matter). Raises ``ValueError`` before anything is
        written on an unknown key, a wrong shape or a non-finite number.
        """
        self._check_new_id(section_id)
        _check_keys(properties, _ELECTRIC_KEYS, "properties.electric")
        out: dict[str, Any] = {}
        what = "properties.electric"
        if "origin" in properties:
            origin = properties["origin"]
            if isinstance(origin, str):
                if not origin:
                    raise ValueError(f"{what} origin must not be empty")
                out["origin"] = origin
            else:
                out["origin"] = _real_vector(origin, 3, f"{what} origin")
        if "dipole" in properties:
            d = properties["dipole"]
            _check_keys(d, ("vector", "total", "units"), f"{what} dipole")
            if "vector" not in d:
                raise ValueError(f"{what} dipole requires a vector")
            block: dict[str, Any] = {
                "vector": _real_vector(d["vector"], 3, f"{what} dipole vector")}
            if "total" in d:
                block["total"] = _finite_number(d["total"], f"{what} dipole total")
                if block["total"] < 0:
                    raise ValueError(f"{what} dipole total must be >= 0")
            if "units" in d:
                if d["units"] != "debye":
                    raise ValueError(f"{what} dipole units must be 'debye'; "
                                     f"spec § 5.13 fixes it (got {d['units']!r})")
                block["units"] = "debye"
            out["dipole"] = {**block, **{k: v for k, v in d.items()
                                         if str(k).startswith("x_")}}
        if "quadrupole" in properties:
            q = properties["quadrupole"]
            _check_keys(q, ("matrix", "traceless", "convention"),
                        f"{what} quadrupole")
            for key in ("matrix", "traceless", "convention"):
                if key not in q:
                    raise ValueError(f"{what} quadrupole requires {key}")
            matrix = _real_matrix(q["matrix"], 3, 3, f"{what} quadrupole matrix")
            _require_symmetric(matrix, f"{what} quadrupole matrix")
            if not isinstance(q["traceless"], (bool, np.bool_)):
                raise ValueError(f"{what} quadrupole traceless must be a boolean")
            convention = _choice(q["convention"], _MULTIPOLE_CONVENTIONS,
                                 f"{what} quadrupole convention")
            if convention == "buckingham" and not q["traceless"]:
                raise ValueError(f"{what} quadrupole: the buckingham "
                                 f"convention is traceless")
            out["quadrupole"] = {
                "matrix": matrix, "traceless": bool(q["traceless"]),
                "convention": convention,
                **{k: v for k, v in q.items() if str(k).startswith("x_")}}
        if "octupole" in properties:
            o = properties["octupole"]
            _check_keys(o, ("tensor", "convention"), f"{what} octupole")
            if "tensor" not in o:
                raise ValueError(f"{what} octupole requires a tensor")
            block = {"tensor": _real_tensor3(o["tensor"], f"{what} octupole tensor")}
            if "convention" in o:
                block["convention"] = _choice(o["convention"],
                                              _MULTIPOLE_CONVENTIONS,
                                              f"{what} octupole convention")
            out["octupole"] = {**block, **{k: v for k, v in o.items()
                                           if str(k).startswith("x_")}}
        if "polarizability" in properties:
            p = properties["polarizability"]
            _check_keys(p, ("matrix", "frequency_ev", "isotropic"),
                        f"{what} polarizability")
            if "matrix" not in p:
                raise ValueError(f"{what} polarizability requires a matrix")
            block = {"matrix": _real_matrix(p["matrix"], 3, 3,
                                            f"{what} polarizability matrix")}
            if "frequency_ev" in p:
                block["frequency_ev"] = _finite_number(
                    p["frequency_ev"], f"{what} polarizability frequency_ev")
                if block["frequency_ev"] < 0:
                    raise ValueError(f"{what} polarizability frequency_ev must "
                                     f"be >= 0")
            if "isotropic" in p:
                block["isotropic"] = _finite_number(
                    p["isotropic"], f"{what} polarizability isotropic")
            out["polarizability"] = {**block, **{k: v for k, v in p.items()
                                                 if str(k).startswith("x_")}}
        if "hyperpolarizability" in properties:
            b = properties["hyperpolarizability"]
            _check_keys(b, ("tensor", "frequencies_ev"),
                        f"{what} hyperpolarizability")
            if "tensor" not in b:
                raise ValueError(f"{what} hyperpolarizability requires a tensor")
            block = {"tensor": _real_tensor3(b["tensor"],
                                             f"{what} hyperpolarizability tensor")}
            if "frequencies_ev" in b:
                block["frequencies_ev"] = _real_vector(
                    b["frequencies_ev"], 2,
                    f"{what} hyperpolarizability frequencies_ev")
            out["hyperpolarizability"] = {**block, **{
                k: v for k, v in b.items() if str(k).startswith("x_")}}
        if not any(name in out for name in _ELECTRIC_BLOCKS):
            raise ValueError(f"{what} needs at least one of "
                             f"{', '.join(_ELECTRIC_BLOCKS)}")
        if "units" in properties:
            out["units"] = properties["units"]
        out = _with_units(out, None, what, "ElectricPropertiesPayload")
        out.update({k: v for k, v in properties.items() if str(k).startswith("x_")})
        return self._add_properties("properties.electric", section_id, out,
                                    label=label)

    def add_elastic_properties(self, stiffness: Any, *, compliance: Any = None,
                               moduli: Mapping[str, float] | None = None,
                               average: str | None = None,
                               units: Mapping[str, str] | None = None,
                               label: str | None = None,
                               section_id: str = "elastic") -> str:
        """Add a ``properties.elastic`` section (spec § 5.13,
        ``$defs/ElasticPropertiesPayload``).

        ``stiffness`` is the symmetric 6 x 6 stiffness matrix in Voigt notation
        (order xx, yy, zz, yz, xz, xy) in GPa. ``compliance`` (optional) is the
        6 x 6 compliance matrix in 1/GPa. ``moduli`` (optional) maps any of
        ``bulk_modulus_gpa``, ``shear_modulus_gpa``, ``youngs_modulus_gpa``
        (GPa) and ``poisson_ratio`` to a number, and ``average`` names the
        scheme they come from: ``"voigt"``, ``"reuss"`` or ``"hill"`` (a
        producer SHOULD give it with any modulus). ``units`` declares units
        (``stiffness`` and the moduli ``"gpa"``, ``poisson_ratio``
        ``"dimensionless"``); nothing is converted. A matrix that is not
        6 x 6 and symmetric to within 1e-6 of its largest entry raises
        ``ValueError`` before anything is written.
        """
        self._check_new_id(section_id)
        what = "properties.elastic"
        stiff = _real_matrix(stiffness, 6, 6, f"{what} stiffness")
        _require_symmetric(stiff, f"{what} stiffness")
        out: dict[str, Any] = {"stiffness": stiff}
        if compliance is not None:
            comp = _real_matrix(compliance, 6, 6, f"{what} compliance")
            _require_symmetric(comp, f"{what} compliance")
            out["compliance"] = comp
        if moduli is not None:
            if not isinstance(moduli, Mapping):
                raise ValueError(f"{what} moduli must be a mapping")
            unknown = sorted(str(k) for k in moduli if k not in _ELASTIC_MODULI)
            if unknown:
                raise ValueError(f"{what} moduli: unknown keys {unknown}; use "
                                 f"{', '.join(_ELASTIC_MODULI)}")
            for key in _ELASTIC_MODULI:
                if key in moduli:
                    out[key] = _finite_number(moduli[key], f"{what} {key}")
        if average is not None:
            out["average"] = _choice(average, _ELASTIC_AVERAGES,
                                     f"{what} average")
        out = _with_units(out, units, what, "ElasticPropertiesPayload")
        return self._add_properties("properties.elastic", section_id, out,
                                    label=label)

    def add_dielectric_properties(self, electronic: Any, *, ionic: Any = None,
                                  static: Any = None, born_charges: Any = None,
                                  piezoelectric: Any = None,
                                  units: Mapping[str, str] | None = None,
                                  label: str | None = None,
                                  section_id: str = "dielectric") -> str:
        """Add a ``properties.dielectric`` section (spec § 5.13,
        ``$defs/DielectricPropertiesPayload``).

        ``electronic`` is the 3 x 3 electronic (high-frequency) dielectric
        tensor; ``ionic`` and ``static`` (optional, 3 x 3) are the ionic
        contribution and the static tensor, dimensionless; when all three are
        given, ``static`` must equal ``electronic + ionic`` within 1e-6.
        ``born_charges`` (optional) is an ``[n_atoms, 3, 3]`` array, one Born
        tensor per atom in units of e, ``piezoelectric`` (optional) the 3 x 6
        piezoelectric stress tensor in C/m² (Voigt columns xx, yy, zz, yz, xz,
        xy). ``units`` declares units (the three tensors ``"dimensionless"``,
        ``born_charges`` ``"e"``); nothing is converted. That there is one
        Born tensor per atom of the section's structure is checked by the
        validator. Raises ``ValueError`` before anything is written on a wrong
        shape or a ``static`` that is not ``electronic + ionic``.
        """
        self._check_new_id(section_id)
        what = "properties.dielectric"
        out: dict[str, Any] = {
            "electronic": _real_matrix(electronic, 3, 3, f"{what} electronic")}
        if ionic is not None:
            out["ionic"] = _real_matrix(ionic, 3, 3, f"{what} ionic")
        if static is not None:
            out["static"] = _real_matrix(static, 3, 3, f"{what} static")
        if "ionic" in out and "static" in out:
            for i in range(3):
                for j in range(3):
                    total = out["electronic"][i][j] + out["ionic"][i][j]
                    if abs(out["static"][i][j] - total) > _PROPERTY_TOLERANCE:
                        raise ValueError(
                            f"{what} static[{i}][{j}] is {out['static'][i][j]!r} "
                            f"but electronic + ionic is {total!r}; static must "
                            f"equal their sum within 1e-6 (spec § 5.13)")
        if born_charges is not None:
            arr = _real_array(born_charges, f"{what} born_charges")
            if arr.ndim != 3 or arr.shape[1:] != (3, 3) or arr.shape[0] < 1:
                raise ValueError(f"{what} born_charges must be [n_atoms, 3, 3], "
                                 f"got shape {list(arr.shape)}")
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"{what} born_charges must be finite")
            out["born_charges"] = arr.astype(np.float64).tolist()
        if piezoelectric is not None:
            out["piezoelectric"] = _real_matrix(piezoelectric, 3, 6,
                                                f"{what} piezoelectric")
        out = _with_units(out, units, what, "DielectricPropertiesPayload")
        return self._add_properties("properties.dielectric", section_id, out,
                                    label=label)

    _SPECTRA_KINDS = frozenset({
        "spectra.ir", "spectra.raman", "spectra.uvvis", "spectra.ecd",
        "spectra.vcd", "spectra.nmr", "spectra.epr", "spectra.generic",
    })

    def add_spectrum(self, kind: str, *,
                    frequencies: Sequence[float] | None = None,
                    intensities: Sequence[float] | None = None,
                    payload: Mapping[str, Any] | None = None,
                    x_units: str | None = None,
                    y_units: str | None = None,
                    peaks: Sequence[Mapping[str, Any]] | None = None,
                    broadening: Mapping[str, Any] | None = None,
                    label: str | None = None,
                    section_id: str | None = None) -> str:
        """Add a ``spectra.*`` section (spec § 5.4).

        For the 1-D kinds (``spectra.ir`` / ``.raman`` / ``.uvvis`` / ``.ecd`` /
        ``.vcd`` / ``.generic``, ``$defs/SpectrumPayload``) pass
        ``frequencies`` + ``intensities``: the x values, in the kind's § 4.3
        unit unless ``x_units`` declares another, and the y values. ``x_units``
        / ``y_units`` are unit names of spec § 4.3 (an unknown one raises
        ``ValueError``), declared, never converted; ``peaks`` (``{x, y?,
        width?, label?, assignment?}``) and ``broadening`` (``{type, fwhm?,
        units?}``) are written as given. Alternatively pass the whole
        ``payload`` dict, which is how the other optional blocks (``x``,
        ``temperature_k``, ``raman``, ``uvvis``) are written and the only form
        for ``spectra.nmr`` and ``spectra.epr``. A 1-D payload's
        ``frequencies``, ``intensities`` and ``x`` must have one length, and
        its ``states_ref``, which makes the peaks' ``assignment.state_index``
        values state indices (spec § 5.4), must name an ``excited_states``
        section already added; anything else raises ``ValueError``.
        """
        if section_id is not None:
            self._check_new_id(section_id)
        if kind not in self._SPECTRA_KINDS:
            raise ValueError(f"{kind!r} is not a spectra kind")
        one_d = kind not in ("spectra.nmr", "spectra.epr")
        if payload is None:
            if frequencies is None or intensities is None:
                raise ValueError(f"{kind}: pass frequencies+intensities or payload")
            payload = {"frequencies": [float(x) for x in frequencies],
                       "intensities": [float(x) for x in intensities]}
        body = dict(payload)
        extras = {"x_units": x_units, "y_units": y_units, "peaks": peaks,
                  "broadening": broadening}
        for key, value in extras.items():
            if value is None:
                continue
            if not one_d:
                raise ValueError(f"{kind}: {key} applies to the 1-D spectra "
                                 f"kinds only")
            if key in body:
                raise ValueError(f"{kind}: {key} given both as an argument and "
                                 f"in payload")
            if key == "peaks":
                value = [dict(p) for p in value]
            elif key == "broadening":
                value = dict(value)
            body[key] = value
        if one_d:
            for key in ("x_units", "y_units"):
                if key in body:
                    _checked_unit_name(body[key], f"{kind}: {key}")
            broad = body.get("broadening")
            if isinstance(broad, Mapping) and "units" in broad:
                _checked_unit_name(broad["units"], f"{kind}: broadening units")
            freqs = body.get("frequencies")
            if (not isinstance(freqs, (list, tuple))
                    or not isinstance(body.get("intensities"), (list, tuple))):
                raise ValueError(f"{kind}: the payload needs frequencies and "
                                 f"intensities arrays (spec § 5.4)")
            for key in ("intensities", "x"):
                other = body.get(key)
                if other is None:
                    continue
                if not isinstance(other, (list, tuple)):
                    raise ValueError(f"{kind}: {key} must be an array")
                if len(other) != len(freqs):
                    raise ValueError(
                        f"{kind}: {key} has {len(other)} entries but frequencies "
                        f"has {len(freqs)}; they must have the same length "
                        f"(spec § 5.4)")
            if "states_ref" in body:
                ref = body["states_ref"]
                target = next((s for s in self._sections if s["id"] == ref), None)
                if target is None or target["kind"] != "excited_states":
                    raise ValueError(f"{kind}: states_ref {ref!r} must name an "
                                     f"added excited_states section "
                                     f"(spec § 5.4)")
        sid = section_id or kind.replace(".", "_")
        self._check_new_id(sid)
        member = self._json_member(f"spectra/{sid}.json", body)
        section: dict[str, Any] = {"id": self._reserve_id(sid), "kind": kind,
                                   "members": {"spectrum": member}}
        if label is not None:
            section["label"] = label
        return self._push(section)

    # -- bands / dos ------------------------------------------------------

    def add_bands(self, kpath: Mapping[str, Any], eigenvalues: np.ndarray, *,
                 projections: np.ndarray | None = None,
                 occupations: np.ndarray | None = None,
                 label: str | None = None, section_id: str = "bands") -> str:
        """Add a ``bands`` section (spec § 5.6).

        ``eigenvalues`` is float64 [n_spin, n_k, n_bands] in eV. Optional
        ``projections`` (fat bands) is float64 [n_spin, n_k, n_bands,
        n_channels], or the legacy [n_k, n_bands, n_channels], which means
        n_spin = 1. Optional ``occupations`` is float64 [n_spin, n_k, n_bands]
        (electrons per state). Shapes must agree with ``eigenvalues``.
        """
        self._check_new_id(section_id)
        eig = np.asarray(eigenvalues, dtype=np.float64)
        if eig.ndim != 3:
            raise ValueError("eigenvalues must be [n_spin, n_kpoints, n_bands]")
        proj = occ = None
        if projections is not None:
            proj = np.asarray(projections, dtype=np.float64)
            if proj.ndim == 4:
                if proj.shape[:3] != eig.shape:
                    raise ValueError(
                        f"projections {proj.shape} must be [n_spin, n_k, n_bands, "
                        f"n_channels] with [n_spin, n_k, n_bands] = {eig.shape}")
            elif proj.ndim == 3:
                if proj.shape[:2] != eig.shape[1:]:
                    raise ValueError(
                        f"projections {proj.shape} must be [n_k, n_bands, "
                        f"n_channels] with [n_k, n_bands] = {eig.shape[1:]}")
            else:
                raise ValueError("projections must be [n_spin, n_k, n_bands, "
                                 "n_channels] or [n_k, n_bands, n_channels]")
        if occupations is not None:
            occ = np.asarray(occupations, dtype=np.float64)
            if occ.shape != eig.shape:
                raise ValueError(f"occupations {occ.shape} must have the shape "
                                 f"of eigenvalues {eig.shape}")
        members = {
            "kpath": self._json_member(f"{section_id}/kpath.json", dict(kpath)),
            "eigenvalues": self._binary_member(f"{section_id}/eigenvalues.bin",
                                               eig, "float64"),
        }
        if proj is not None:
            members["projections"] = self._binary_member(
                f"{section_id}/projections.bin", proj, "float64")
        if occ is not None:
            members["occupations"] = self._binary_member(
                f"{section_id}/occupations.bin", occ, "float64")
        return self._push({"id": self._reserve_id(section_id), "kind": "bands",
                          "members": members}, label=label)

    def add_dos_total(self, energies: Sequence[float], dos: np.ndarray, *,
                     smearing: float | None = None, smearing_type: str | None = None,
                     fermi_energy_ev: float | None = None,
                     n_electrons: float | None = None, n_spin: int = 1,
                     label: str | None = None, section_id: str = "dos_total") -> str:
        """Add a ``dos.total`` section. Its metadata rides as optional fields on
        the section object, not as members (spec § 5.6): ``n_spin`` (1 or 2),
        ``smearing`` (broadening width in eV, ≥ 0), ``smearing_type`` (e.g.
        ``"gaussian"``), ``fermi_energy_ev`` and ``n_electrons``."""
        self._check_new_id(section_id)
        members = {
            "energies": self._binary_member(f"dos/{section_id}_energies.bin",
                                            np.asarray(energies), "float64"),
            "dos": self._binary_member(f"dos/{section_id}.bin",
                                       np.asarray(dos), "float64"),
        }
        section: dict[str, Any] = {"id": self._reserve_id(section_id),
                                   "kind": "dos.total", "members": members,
                                   "n_spin": int(n_spin)}
        for k, v in (("smearing", smearing), ("smearing_type", smearing_type),
                     ("fermi_energy_ev", fermi_energy_ev),
                     ("n_electrons", n_electrons)):
            if v is not None:
                section[k] = v
        return self._push(section, label=label)

    def add_dos_projected(self, energies: Sequence[float], projections: np.ndarray,
                         channels: Sequence[Mapping[str, Any]], *, n_spin: int = 1,
                         fermi_energy_ev: float | None = None,
                         label: str | None = None,
                         section_id: str = "dos_projected") -> str:
        """Add a ``dos.projected`` section. ``channels`` labels the rows of
        ``projections``, one object per row with the optional keys ``atom``,
        ``l``, ``m``, ``spin`` and ``label``; like ``n_spin`` and
        ``fermi_energy_ev`` it rides as a field on the section object
        (spec § 5.6). ``projections`` is [n_spin, n_channels, n_points], or
        the legacy [n_channels, n_points], which means n_spin = 1; a rank-3
        array's spin axis must equal ``n_spin``."""
        self._check_new_id(section_id)
        proj = np.asarray(projections)
        if proj.ndim not in (2, 3):
            raise ValueError("projections must be [n_spin, n_channels, n_points] "
                             "or [n_channels, n_points]")
        if proj.ndim == 3 and proj.shape[0] != int(n_spin):
            raise ValueError(f"projections {proj.shape} has a spin axis of "
                             f"{proj.shape[0]}, but n_spin is {n_spin}")
        members = {
            "energies": self._binary_member(f"dos/{section_id}_energies.bin",
                                            np.asarray(energies), "float64"),
            "projections": self._binary_member(f"dos/{section_id}_projections.bin",
                                               proj, "float64"),
        }
        section: dict[str, Any] = {"id": self._reserve_id(section_id),
                                   "kind": "dos.projected", "members": members,
                                   "n_spin": int(n_spin),
                                   "channels": [dict(c) for c in channels]}
        if fermi_energy_ev is not None:
            section["fermi_energy_ev"] = fermi_energy_ev
        return self._push(section, label=label)

    def add_dos_coop_cohp(self, kind: str, energies: Sequence[float],
                         projections: np.ndarray, integrated: np.ndarray,
                         meta: Mapping[str, Any], *, label: str | None = None,
                                                     section_id: str | None = None
                         ) -> str:
        """Add a ``dos.coop`` or ``dos.cohp`` section (bonding analysis). Here
        ``meta`` IS a member (JSON), unlike dos.total/projected."""
        if section_id is not None:
            self._check_new_id(section_id)
        if kind not in ("dos.coop", "dos.cohp"):
            raise ValueError("kind must be dos.coop or dos.cohp")
        sid = section_id or kind.replace(".", "_")
        self._check_new_id(sid)
        members = {
            "energies": self._binary_member(f"dos/{sid}_energies.bin",
                                            np.asarray(energies), "float64"),
            "projections": self._binary_member(f"dos/{sid}_projections.bin",
                                               np.asarray(projections), "float64"),
            "integrated": self._binary_member(f"dos/{sid}_integrated.bin",
                                              np.asarray(integrated), "float64"),
            "meta": self._json_member(f"dos/{sid}_meta.json", dict(meta)),
        }
        return self._push({"id": self._reserve_id(sid), "kind": kind,
                          "members": members}, label=label)

    # -- symmetry / scf history / citations ------------------------------

    def add_structure_symmetry(self, data: Mapping[str, Any], *,
                              label: str | None = None,
                              section_id: str = "symmetry") -> str:
        """Add a ``structure.symmetry`` section; ``data`` is the payload
        (``$defs/SymmetryPayload``, spec § 5.7), written as given."""
        self._check_new_id(section_id)
        member = self._json_member(f"structure/symmetry.json", dict(data))
        return self._push({"id": self._reserve_id(section_id),
                          "kind": "structure.symmetry",
                          "members": {"data": member}}, label=label)

    def add_scf_history(self, iterations: Sequence[Mapping[str, Any]], *,
                       label: str | None = None,
                       section_id: str = "scf_history") -> str:
        """Add an ``scf_history`` section (``$defs/ScfHistoryPayload``, spec
        § 5.7). Each iteration is ``{iter, energy_eh, delta_e?, diis_error?,
        rms_density_change?}``."""
        self._check_new_id(section_id)
        member = self._json_member(f"{section_id}/iterations.json",
                                   {"iterations": [dict(i) for i in iterations]})
        return self._push({"id": self._reserve_id(section_id), "kind": "scf_history",
                          "members": {"iterations": member}}, label=label)

    def add_citations(self, bibtex: str, *, label: str | None = None,
                                            section_id: str = "citations") -> str:
        """Add a ``citations`` section from a BibTeX string."""
        self._check_new_id(section_id)
        data = bibtex.encode("utf-8") if isinstance(bibtex, str) else bytes(bibtex)
        member = self._opaque_member(f"{section_id}/references.bib", data)
        return self._push({"id": self._reserve_id(section_id), "kind": "citations",
                          "members": {"references": member}}, label=label)

    def add_run_record(self, *, program: str,
                       input_text: str | bytes | None = None,
                       log_text: str | bytes | None = None,
                       program_version: str | None = None,
                       command: str | None = None,
                       exit_status: int | None = None,
                       started_utc: str | None = None,
                       finished_utc: str | None = None,
                       sequence: int | None = None,
                       files: Mapping[str, Any] | None = None,
                       attachments: Mapping[str, bytes] | None = None,
                       job_ref: str | None = None,
                       run_status: str | None = None,
                       label: str | None = None, section_id: str = "run_record") -> str:
        """Add a ``run.record`` section — the self-contained record of one
        program invocation.

        ``program`` names the code the input/log belong to (e.g. ``"orca"``),
        which may differ from the archive-level ``source.program`` when the
        QVF is written by a converter. At least one of ``input_text`` /
        ``log_text`` is required; both are stored as opaque UTF-8 bytes.
        ``attachments`` maps a role suffix to raw bytes, stored as members
        named ``attachment.<suffix>``; describe original filenames in the
        optional ``files`` index (role -> {filename, description?, media_type?}).

        ``job_ref`` names the ``job.spec`` section this run executed (it
        must already be added) and ``run_status`` gives the run's outcome,
        ``"running"``, ``"converged"`` or ``"failed"``: in a job set they say
        which job ran and how it ended (spec § 5.8, § 5.9).
        """
        self._check_new_id(section_id)
        if not program:
            raise ValueError("run.record requires a non-empty program")
        if job_ref is not None:
            self._check_ref("run.record job_ref", job_ref, "job.spec")
        if run_status is not None:
            _choice(run_status, ("running", "converged", "failed"),
                    "run.record run_status")
        if input_text is None and log_text is None:
            raise ValueError("run.record requires input_text and/or log_text")
        sid = self._reserve_id(section_id)
        members: dict[str, Any] = {}

        def _opaque(role: str, path: str, payload: str | bytes) -> None:
            data = (payload.encode("utf-8") if isinstance(payload, str)
                    else bytes(payload))
            members[role] = self._opaque_member(path, data)

        if input_text is not None:
            _opaque("input", f"{section_id}/input.txt", input_text)
        if log_text is not None:
            _opaque("log", f"{section_id}/log.txt", log_text)
        for suffix, blob in (attachments or {}).items():
            _opaque(f"attachment.{suffix}",
                    f"{section_id}/attachments/{suffix}", blob)
        if files is not None:
            members["files"] = self._json_member(f"{section_id}/files.json",
                                                 dict(files))
        section: dict[str, Any] = {"id": sid, "kind": "run.record",
                                   "program": program, "members": members}
        if program_version is not None:
            section["program_version"] = program_version
        if command is not None:
            section["command"] = command
        if exit_status is not None:
            section["exit_status"] = int(exit_status)
        if started_utc is not None:
            section["started_utc"] = started_utc
        if finished_utc is not None:
            section["finished_utc"] = finished_utc
        if sequence is not None:
            section["sequence"] = int(sequence)
        if job_ref is not None:
            section["job_ref"] = job_ref
        if run_status is not None:
            section["run_status"] = run_status
        return self._push(section, label=label)

    def add_job_spec(self, *, job_type: str,
                     method: str | None = None,
                     basis: str | None | _Unset = _UNSET,
                     functional: str | None | _Unset = _UNSET,
                     charge: int | None = None,
                     multiplicity: int | None = None,
                     kpoints: Sequence[int] | None = None,
                     tasks: Sequence[str] | None = None,
                     structure_ref: str | None = None,
                     solvation: Mapping[str, Any] | None = None,
                     dispersion: str | None = None,
                     relativistic: str | None = None,
                     scf: Mapping[str, Any] | None = None,
                     geometry: Mapping[str, Any] | None = None,
                     resources: Mapping[str, Any] | None = None,
                     options: Mapping[str, Any] | None = None,
                     label: str | None = None, section_id: str = "job_spec") -> str:
        """Add a ``job.spec`` section — the declarative specification of the
        calculation this archive *requests* (QVF spec § 5.9).

        Omitted ``basis`` and ``functional`` leave the choice to the engine.
        Pass ``None`` to write JSON null (the method takes none); strings,
        including the deprecated empty spelling, are written unchanged.

        ``job_type`` is ``"molecular"`` or ``"periodic"`` and selects which
        class of engine runs the job; the geometry itself (including
        periodicity) lives in the archive's ``structure`` section, which
        ``structure_ref`` names (a section id; it may be added later, but one
        that is already added must be a ``structure``). The typed blocks are
        mappings checked against the schema: ``solvation`` ``{model,
        solvent?, epsilon?}`` (``epsilon`` > 0), ``scf`` ``{convergence?,
        max_iterations?, guess?, damping?}`` (the convergence threshold in
        hartree, > 0), ``geometry`` ``{convergence?, max_steps?,
        coordinates?}`` (the threshold in hartree/bohr, > 0; ``coordinates``
        ``"cartesian"``, ``"internal"`` or ``"redundant"``) and
        ``resources`` ``{n_cores?, memory_mb?, walltime_s?}``, which are
        hints, not requirements. ``dispersion`` (e.g. ``"d3bj"``) and
        ``relativistic`` (e.g. ``"x2c"``) are non-empty strings. A wrong
        value raises ``ValueError`` before anything is added. Engine keywords
        beyond the typed fields go under ``options``. Pair with
        ``set_provenance(run_status="pending")`` for an archive describing a
        job that has not yet run.

        ``structure_ref`` is written on the section object, where a reader
        finds it in the manifest, and in the payload, as runners have read
        it so far (spec § 5.9). An archive may carry several ``job.spec``
        sections, a job set: give each its own ``section_id`` and, when the
        archive holds more than one structure, its ``structure_ref``.
        """
        self._check_new_id(section_id)
        if job_type not in ("molecular", "periodic"):
            raise ValueError(
                f"job.spec job_type must be 'molecular' or 'periodic', "
                f"got {job_type!r}")
        payload: dict[str, Any] = {"job_type": job_type}
        if method is not None:
            payload["method"] = str(method)
        for name, value in (("basis", basis), ("functional", functional)):
            if value is _UNSET:
                continue
            if value is not None and not isinstance(value, str):
                raise ValueError(f"job.spec {name} must be a string or None")
            payload[name] = value
        if charge is not None:
            payload["charge"] = int(charge)
        if multiplicity is not None:
            payload["multiplicity"] = int(multiplicity)
        if kpoints is not None:
            mesh = [int(n) for n in kpoints]
            if len(mesh) != 3 or any(n < 1 for n in mesh):
                raise ValueError(
                    f"job.spec kpoints must be three integers >= 1, "
                    f"got {list(kpoints)!r}")
            payload["kpoints"] = mesh
        if tasks is not None:
            payload["tasks"] = [str(t) for t in tasks]
        if structure_ref is not None:
            ref = _checked_text(structure_ref, "job.spec structure_ref")
            if not _ID_RE.match(ref):
                raise ValueError(f"job.spec structure_ref {ref!r} is not a "
                                 f"section id ([A-Za-z0-9_.-]+)")
            for section in self._sections:
                if section["id"] == ref and section["kind"] != "structure":
                    raise ValueError(
                        f"job.spec structure_ref {ref!r} names a "
                        f"{section['kind']!r} section, not a 'structure' "
                        f"section (spec § 5.9)")
            payload["structure_ref"] = ref
        if solvation is not None:
            payload["solvation"] = _job_spec_block("solvation", solvation)
        if dispersion is not None:
            payload["dispersion"] = _checked_text(dispersion,
                                                  "job.spec dispersion")
        if relativistic is not None:
            payload["relativistic"] = _checked_text(relativistic,
                                                    "job.spec relativistic")
        for name, block in (("scf", scf), ("geometry", geometry),
                            ("resources", resources)):
            if block is not None:
                payload[name] = _job_spec_block(name, block)
        if options is not None:
            payload["options"] = dict(options)
        sid = self._reserve_id(section_id)
        member = self._json_member(f"{sid}/spec.json", payload)
        section: dict[str, Any] = {"id": sid, "kind": "job.spec"}
        if structure_ref is not None:
            section["structure_ref"] = payload["structure_ref"]
        section["members"] = {"spec": member}
        return self._push(section, label=label)

    # -- periodic-only kinds ---------------------------------------------

    def add_fermi_surface(self, mesh: Mapping[str, Any], energies: np.ndarray, *,
                         label: str | None = None,
                         section_id: str = "fermi_surface") -> str:
        """Add a ``fermi_surface`` section. ``energies`` is
        [nk1, nk2, nk3, n_bands] float32/64 in eV (signed E(k)-E_F)."""
        self._check_new_id(section_id)
        arr = np.asarray(energies)
        if arr.ndim != 4:
            raise ValueError("fermi_surface energies must be 4-D")
        members = {
            "mesh": self._json_member(f"fermi/{section_id}_mesh.json", dict(mesh)),
            "energies": self._binary_member(f"fermi/{section_id}.bin", arr,
                                            "float64" if arr.dtype == np.float64
                                            else "float32"),
        }
        return self._push({"id": self._reserve_id(section_id), "kind": "fermi_surface",
                          "members": members}, label=label)

    def add_phonon_bands(self, qpath: Mapping[str, Any], frequencies: np.ndarray, *,
                        eigenvectors: np.ndarray | None = None,
                        label: str | None = None,
                        section_id: str = "phonon_bands") -> str:
        """Add a ``phonon_bands`` section (spec § 5.6). ``frequencies`` is
        [n_q, n_modes] in cm⁻¹.

        ``eigenvectors`` (optional) is ``[n_q, n_modes, n_atoms, 3]`` when real,
        or ``[n_q, n_modes, n_atoms, 3, 2]`` with the last axis ``[re, im]``
        when complex; a complex NumPy array of the real shape is split into
        that form for you. Its leading ``[n_q, n_modes]`` must match
        ``frequencies``, and it is written as float64. A wrong rank or trailing
        axis, a non-numeric dtype or a mismatch raises ``ValueError`` before
        anything is written.
        """
        self._check_new_id(section_id)
        eig = (_phonon_eigenvectors(eigenvectors, frequencies)
               if eigenvectors is not None else None)
        members = {
            "qpath": self._json_member(f"phonons/{section_id}_qpath.json",
                                       dict(qpath)),
            "frequencies": self._binary_member(f"phonons/{section_id}_freq.bin",
                                               np.asarray(frequencies), "float64"),
        }
        if eig is not None:
            members["eigenvectors"] = self._binary_member(
                f"phonons/{section_id}_eig.bin", eig, "float64")
        return self._push({"id": self._reserve_id(section_id), "kind": "phonon_bands",
                          "members": members}, label=label)

    def add_phonon_dos(self, meta: Mapping[str, Any], frequencies: Sequence[float],
                      dos: np.ndarray, *, projected: np.ndarray | None = None,
                      label: str | None = None, section_id: str = "phonon_dos") -> str:
        """Add a ``phonon_dos`` section (spec § 5.6). ``frequencies`` is the
        [n_points] grid in cm⁻¹ and ``dos`` the total phonon DOS on it.

        ``projected`` (optional) is the real float64 ``[n_channels, n_points]``
        DOS per channel, on the same grid. Label its rows with a ``channels``
        array in ``meta``, one object per row (e.g. ``{"atom": 0, "label":
        "Mg"}``); when ``meta`` carries ``channels``, its length must equal
        ``n_channels``. A wrong rank, a non-real dtype or a mismatch raises
        ``ValueError`` before anything is written.
        """
        self._check_new_id(section_id)
        proj = None
        if projected is not None:
            proj = _real_array(projected, "phonon_dos projected")
            if proj.ndim != 2:
                raise ValueError("phonon_dos projected must be "
                                 "[n_channels, n_points], got shape "
                                 f"{list(proj.shape)}")
            freq_shape = tuple(np.shape(frequencies))
            if len(freq_shape) != 1 or proj.shape[1] != freq_shape[0]:
                raise ValueError(
                    f"phonon_dos projected has {proj.shape[1]} points per "
                    f"channel but frequencies has shape {list(freq_shape)}")
            channels = meta.get("channels")
            if isinstance(channels, (list, tuple)) and len(channels) != proj.shape[0]:
                raise ValueError(
                    f"phonon_dos meta channels has {len(channels)} entries but "
                    f"projected has {proj.shape[0]} rows")
        members = {
            "meta": self._json_member(f"phonons/{section_id}_meta.json", dict(meta)),
            "frequencies": self._binary_member(f"phonons/{section_id}_freq.bin",
                                               np.asarray(frequencies), "float64"),
            "dos": self._binary_member(f"phonons/{section_id}_dos.bin",
                                       np.asarray(dos), "float64"),
        }
        if proj is not None:
            members["projected"] = self._binary_member(
                f"phonons/{section_id}_proj.bin", proj, "float64")
        return self._push({"id": self._reserve_id(section_id), "kind": "phonon_dos",
                          "members": members}, label=label)

    def add_equation_of_state(self, volumes: Sequence[float],
                             energies: Sequence[float], fit: Mapping[str, Any], *,
                             label: str | None = None, section_id: str = "eos") -> str:
        """Add an ``equation_of_state`` section. ``volumes`` Å³, ``energies`` eV,
        ``fit`` carries model + V0/E0/B0/B0_prime."""
        self._check_new_id(section_id)
        members = {
            "volumes": self._binary_member(f"eos/{section_id}_volumes.bin",
                                           np.asarray(volumes), "float64"),
            "energies": self._binary_member(f"eos/{section_id}_energies.bin",
                                            np.asarray(energies), "float64"),
            "fit": self._json_member(f"eos/{section_id}_fit.json", dict(fit)),
        }
        return self._push({"id": self._reserve_id(section_id),
                          "kind": "equation_of_state", "members": members}, label=label)

    def add_topology_qtaim(self, critical_points: Mapping[str, Any], *,
                          units: Mapping[str, str] | None = None,
                          label: str | None = None, section_id: str = "qtaim") -> str:
        """Add a ``topology.qtaim`` section (spec § 5.7, ``$defs/QtaimPayload``).

        ``critical_points`` is the payload: ``{"points": [{type, position, rho?,
        laplacian?, ellipticity?, atom_pair?, …}], "bond_paths"?: [{atoms,
        path}]}``. Positions are Cartesian in Å unless ``units`` declares
        another unit, e.g. ``units={"position": "bohr"}``; ``units`` is written
        as the payload's ``units`` object and nothing is converted. Its values
        are unit names of spec § 4.3; an unknown one raises ``ValueError``.
        """
        self._check_new_id(section_id)
        payload = _with_units(dict(critical_points), units, "topology.qtaim",
                              "QtaimPayload")
        member = self._json_member(f"topology/{section_id}.json", payload)
        return self._push({"id": self._reserve_id(section_id), "kind": "topology.qtaim",
                          "members": {"critical_points": member}}, label=label)

    # -- vendor escape hatch ---------------------------------------------

    def add_preview(self, thumbnail: bytes, *,
                    thumbnail_media_type: str | None = None,
                    image: bytes | str | None = None,
                    image_media_type: str | None = None,
                    label: str | None = None,
                    section_id: str = "preview") -> str:
        """Add a ``preview`` section, a picture of the archive (spec § 5.15).

        ``thumbnail`` is PNG or JPEG bytes, the picture a file manager
        shows; its longest side should be 256 to 1024 pixels. ``image`` is
        an optional larger or vector rendering: SVG (bytes, or a ``str``,
        which is encoded as UTF-8), PNG or JPEG. Media types are sniffed
        from the bytes unless given, and checked against them. The picture
        of the archive is the thumbnail of the first ``preview`` section.
        An SVG image must be self-contained and inert: no external
        resources, no scripts. Both members are stored uncompressed.
        """
        self._check_new_id(section_id)
        thumb = bytes(thumbnail)
        thumb_type = _image_media_type(thumb, ("image/png", "image/jpeg"),
                                       "thumbnail", thumbnail_media_type)
        picture: bytes | None = None
        picture_type = None
        if image is not None:
            picture = image.encode("utf-8") if isinstance(image, str) else bytes(image)
            picture_type = _image_media_type(
                picture, ("image/svg+xml", "image/png", "image/jpeg"), "image",
                image_media_type)
        sid = self._reserve_id(section_id)
        extension = {"image/png": "png", "image/jpeg": "jpg",
                     "image/svg+xml": "svg"}
        members = {"thumbnail": self._opaque_member(
            f"{sid}/thumbnail.{extension[thumb_type]}", thumb, thumb_type)}
        if picture is not None and picture_type is not None:
            members["image"] = self._opaque_member(
                f"{sid}/image.{extension[picture_type]}", picture, picture_type)
        for member in members.values():
            self._stored.add(member["path"])
        return self._push({"id": sid, "kind": "preview", "members": members},
                          label=label)

    def add_vendor_section(self, kind: str, *,
                          json_members: Mapping[str, Any] | None = None,
                          binary_members: Mapping[str, tuple] | None = None,
                          opaque_members: Mapping[str, Any] | None = None,
                          critical: bool = False, schema_uri: str | None = None,
                          label: str | None = None,
                          section_id: str | None = None) -> str:
        """Add an ``x_<vendor>.*`` section. ``json_members`` maps role -> object;
        ``binary_members`` maps role -> ``(array, dtype)``; ``opaque_members``
        maps role -> bytes, or ``(bytes, media_type)``, stored without dtype
        and shape (spec § 4.2, § 7.3), e.g. ``{"drawing": (svg_bytes,
        "image/svg+xml")}``. If ``critical`` is set, remember to declare the
        namespace via :meth:`set_extensions`."""
        if section_id is not None:
            self._check_new_id(section_id)
        if not _VENDOR_RE.match(kind):
            raise ValueError(f"vendor kind {kind!r} must match x_<vendor>.*")
        roles = [*(json_members or {}), *(binary_members or {}),
                 *(opaque_members or {})]
        repeated = sorted({r for r in roles if roles.count(r) > 1})
        if repeated:
            raise ValueError(f"vendor section: member role(s) {repeated} given "
                             f"more than once")
        for role, blob in (opaque_members or {}).items():
            media_type = blob[1] if isinstance(blob, tuple) else None
            if media_type is not None and not (isinstance(media_type, str)
                                               and _MEDIA_TYPE_RE.fullmatch(media_type)):
                raise ValueError(f"vendor section member {role!r}: media type "
                                 f"{media_type!r} is not type/subtype (spec § 4.2)")
        sid = section_id or kind.replace(".", "_").replace("x_", "x_", 1)
        self._check_new_id(sid)
        prefix = sid
        members: dict[str, Any] = {}
        for role, obj in (json_members or {}).items():
            members[role] = self._json_member(f"{prefix}/{role}.json", obj)
        for role, spec in (binary_members or {}).items():
            arr, dtype = spec
            members[role] = self._binary_member(f"{prefix}/{role}.bin", arr, dtype)
        for role, blob in (opaque_members or {}).items():
            data, media_type = (blob if isinstance(blob, tuple) else (blob, None))
            members[role] = self._opaque_member(f"{prefix}/{role}.bin",
                                                bytes(data), media_type)
        if not members:
            raise ValueError("vendor section needs at least one member")
        section: dict[str, Any] = {"id": self._reserve_id(sid), "kind": kind,
                                   "members": members}
        if critical:
            section["critical"] = True
        if schema_uri is not None:
            section["schema_uri"] = schema_uri
        if label is not None:
            section["label"] = label
        return self._push(section)

    # -- output -----------------------------------------------------------

    def build_manifest(self) -> dict[str, Any]:
        """Return the manifest object (dict) without writing anything."""
        manifest: dict[str, Any] = {
            "qvf_version": QVF_VERSION,
            "schema_uri": _SCHEMA_URI,
            "byte_order": _BYTE_ORDER,
            "source": dict(self._source),
        }
        metadata = dict(self._root.get("metadata", {}))
        if self.auto_metadata:
            if self._auto_values is None:   # made once, then reused
                self._auto_values = {
                    "created_utc": _datetime.now(_timezone.utc).strftime(
                        "%Y-%m-%dT%H:%M:%SZ"),
                    "uuid": str(_uuid.uuid4()),
                }
            for key, value in self._auto_values.items():
                metadata.setdefault(key, value)
        if metadata:
            manifest["metadata"] = metadata
        for key in ("provenance", "thermochemistry", "dipole_moment",
                    "constraints", "energies", "extensions", "viewer_defaults"):
            if key in self._root:
                manifest[key] = self._root[key]
        manifest["sections"] = self._sections
        return manifest

    # -- output and in-place update (spec § 2.1) ---------------------------

    #: The archive a writer from :meth:`open_for_update` appends to: a file,
    #: or the bytes of an archive opened from memory.
    _update_path: str | None = None
    _update_bytes: bytes | None = None
    #: Root keys of that archive the writer does not model, kept verbatim.
    _update_root: dict[str, Any] | None = None
    #: ``metadata.modified_utc`` as the archive has it, and the time stamp the
    #: update writes instead (made once, so ``to_bytes`` previews ``write``).
    _update_loaded_modified: Any = None
    _update_stamp: str | None = None

    @classmethod
    def open_for_update(cls, path: str | bytes | bytearray | memoryview
                        ) -> "QvfWriter":
        """Open an existing archive for an in-place update (spec § 2.3).

        ``path`` is the archive's file, or its bytes. The manifest is loaded
        (``source``, root blocks, ``sections``); the ``add_*`` methods then
        add sections and :meth:`remove_section` removes one (remove and add
        under the same id to replace it). :meth:`write` appends the new
        members and a new ``manifest.json`` (STORED, last) to the same file;
        :meth:`to_bytes` returns the updated archive and leaves the file
        alone. An archive opened from bytes is written with
        ``write(path)``, to a file of the caller's choosing.

        The previous manifest, and the members of removed or replaced
        sections, stay in the archive as superseded entries (a
        ``W_ORPHAN_MEMBER`` warning each); readers use the last manifest.
        Pass ``repack=True`` to ``write`` or ``to_bytes`` to drop them (see
        :func:`repack`), which a program that saves the same document many
        times should do. Existing members keep their bytes and checksums,
        and a new member never reuses an existing entry name: one whose
        default path is taken gets a fresh one. Root keys the writer does not
        model are kept as they are; ``auto_metadata`` is off, so
        ``created_utc`` and ``uuid`` stay those of the archive, and
        ``metadata.modified_utc`` is set to the time of the update unless
        :meth:`set_metadata` gives one. An archive stamped with the withdrawn
        ``qvf_version`` 2 is read as version 1, and the update writes 1 and
        the v1 ``schema_uri`` (spec § 3.1).

        Raises ``ValueError`` for an archive that repeats an entry name other
        than ``manifest.json``, has no ``manifest.json``, a ``qvf_version``
        other than 1 or 2, or a ``byte_order`` other than ``"little"``.
        """
        import io
        import os
        data: bytes | None = None
        if isinstance(path, (bytes, bytearray, memoryview)):
            data = bytes(path)
            source_file: Any = io.BytesIO(data)
        else:
            source_file = path
        with zipfile.ZipFile(source_file) as zf:
            infos = zf.infolist()
            dups = _repeated_names(infos)
            if dups:
                raise ValueError(f"E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry "
                                 f"name(s) {dups!r}; refusing to update")
            try:
                info = zf.getinfo("manifest.json")      # the last one (spec § 2.3)
            except KeyError:
                raise ValueError(f"not a QVF archive (no manifest.json): "
                                 f"{path}") from None
            if info.file_size > _MAX_MANIFEST_BYTES:
                raise ValueError("manifest.json is over the 256 MiB JSON cap")
            manifest = json.loads(zf.read(info))
            names = [i.filename for i in infos]
        version = manifest.get("qvf_version") if isinstance(manifest, dict) else None
        if isinstance(version, bool) or version not in (1, 2):
            raise ValueError("only qvf_version 1 archives (or the withdrawn 2, "
                             "read as 1) can be updated")
        source = manifest.get("source")
        if (not isinstance(source, dict) or not isinstance(source.get("program"), str)
                or not isinstance(source.get("version"), str)):
            raise ValueError("manifest.json has no usable source block")
        sections = manifest.get("sections")
        if not isinstance(sections, list) or not all(
                isinstance(s, dict) and isinstance(s.get("id"), str)
                and isinstance(s.get("members"), dict) for s in sections):
            raise ValueError("manifest.json has no well-formed sections array")
        if manifest.get("byte_order", "little") != "little":
            raise ValueError("byte_order is not 'little' (spec § 2.2); refusing to update")
        w = cls(source["program"], source["version"],
                str(source.get("calculation", "")), auto_metadata=False)
        w._source = dict(source)            # kept verbatim, extra fields included
        w._sections = list(sections)
        w._ids = {s["id"] for s in sections}
        kept: dict[str, Any] = {}
        for key, value in manifest.items():
            if key in ("qvf_version", "source", "sections"):
                continue
            if key == "schema_uri" and version == 2:
                continue                  # the update writes the v1 URI
            if key in _MODELLED_ROOT_BLOCKS:
                w._root[key] = value
            else:
                kept[key] = value
        # Every existing entry name and member path is taken; None marks a
        # member that is already in the archive.
        taken = set(names)
        for s in sections:
            for member in s["members"].values():
                if isinstance(member, dict) and isinstance(member.get("path"), str):
                    taken.add(member["path"])
        w._files = {name: None for name in sorted(taken)}  # type: ignore[misc]
        if data is not None:
            w._update_bytes = data
        else:
            w._update_path = os.fspath(path)
        w._update_root = kept
        metadata = manifest.get("metadata")
        if isinstance(metadata, dict):
            w._update_loaded_modified = metadata.get("modified_utc")
        return w

    def _final_manifest(self) -> dict[str, Any]:
        """The manifest to write; an update stamps ``metadata.modified_utc``
        (spec § 2.3) unless the caller set a value of its own."""
        manifest = self.build_manifest()
        if self._update_path is None and self._update_bytes is None:
            return manifest
        if self._update_root:
            manifest.update(self._update_root)
        metadata = dict(manifest.get("metadata", {}))
        current = metadata.get("modified_utc")
        if current is None or current == self._update_loaded_modified:
            if self._update_stamp is None:
                self._update_stamp = _datetime.now(_timezone.utc).strftime(
                    "%Y-%m-%dT%H:%M:%SZ")
            metadata["modified_utc"] = self._update_stamp
            manifest["metadata"] = metadata
        return manifest

    def _write_entries(self, zf: zipfile.ZipFile, compression: int,
                       force_zip64: bool) -> None:
        """The pending members, in the order added, then manifest.json."""
        manifest_bytes = _json_bytes(self._final_manifest())
        for path, data in self._files.items():
            if data is None:              # already in the archive being updated
                continue
            method = zipfile.ZIP_STORED if path in self._stored else compression
            _zip_write(zf, path, data, method, force_zip64)
        # manifest.json STOREd, for cheap table-of-contents reads, and last.
        _zip_write(zf, zipfile.ZipInfo("manifest.json"), manifest_bytes,
                   zipfile.ZIP_STORED, force_zip64)

    def _append_into(self, fileobj: Any, compression: int, force_zip64: bool) -> None:
        """Append the pending members and a new manifest to the archive in
        ``fileobj`` (zipfile mode ``"a"``, spec § 2.3). The previous manifest
        stays listed as a superseded entry; readers use the last one."""
        import warnings
        with zipfile.ZipFile(fileobj, "a", compression=compression) as zf:
            dups = _repeated_names(zf.infolist())
            if dups:
                raise ValueError(f"E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry "
                                 f"name(s) {dups!r}; refusing to update")
            with warnings.catch_warnings():
                # zipfile warns about the repeated manifest.json, which § 2.3
                # prescribes; no other name can repeat (_add_file refuses).
                warnings.filterwarnings("ignore", message="Duplicate name: 'manifest.json'")
                self._write_entries(zf, compression, force_zip64)

    def to_bytes(self, *, compression: int = zipfile.ZIP_DEFLATED,
                 force_zip64: bool = False, repack: bool = False) -> bytes:
        """Serialize the archive to an in-memory ZIP and return its bytes.

        Members come first, in the order they were added, then
        ``manifest.json``, STOREd and last (spec § 2.1). ``force_zip64``
        writes a ZIP64 extra field into every local header (``zipfile`` adds
        ZIP64 central-directory fields and end records only when a size, an
        offset or the entry count needs them). For a writer from
        :meth:`open_for_update` this is the updated archive; the file itself
        is left untouched. ``repack=True`` returns it without the superseded
        entries (spec § 2.3), as :func:`repack` would.
        """
        import io
        if self._update_path is not None or self._update_bytes is not None:
            if self._update_bytes is not None:
                buf = io.BytesIO(self._update_bytes)
            else:
                with open(self._update_path, "rb") as fh:   # type: ignore[arg-type]
                    buf = io.BytesIO(fh.read())
            self._append_into(buf, compression, force_zip64)
            if repack:
                buf.seek(0)
                out = io.BytesIO()
                _repack_zip(buf, out)
                return out.getvalue()
            return buf.getvalue()
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", compression=compression) as zf:
            self._write_entries(zf, compression, force_zip64)
        return buf.getvalue()

    def write(self, path: str | None = None, *,
              compression: int = zipfile.ZIP_DEFLATED,
              force_zip64: bool = False, repack: bool = False) -> str:
        """Write ``{path}`` (``.qvf`` appended if missing). Returns the path.

        A writer from :meth:`open_for_update` on a file appends to its own
        archive; ``path`` may then be omitted (or name that archive). The
        writer stays usable: members added afterwards are appended by the
        next ``write``. ``repack=True`` then rewrites the file without the
        superseded entries (spec § 2.3), as :func:`repack` does. A writer
        opened from bytes writes the updated archive to ``path``.
        """
        import os
        if self._update_bytes is not None:
            if path is None:
                raise TypeError("write() needs a path for an archive opened "
                                "from bytes")
            if not str(path).endswith(".qvf"):
                path = f"{path}.qvf"
            data = self.to_bytes(compression=compression,
                                 force_zip64=force_zip64, repack=repack)
            with open(path, "wb") as fh:
                fh.write(data)
            # The written archive is the new base of further updates.
            self._update_bytes = data
            self._files = {p: None for p in self._files}  # type: ignore[misc]
            self._mark_written(self._final_manifest())
            return path
        if self._update_path is not None:
            if path is not None and not _same_file(os.fspath(path), self._update_path):
                raise ValueError(f"this writer updates {self._update_path!r}; "
                                 f"call write() without another path")
            written = self._final_manifest()
            with open(self._update_path, "r+b") as fh:
                self._append_into(fh, compression, force_zip64)
            for p, d in list(self._files.items()):
                if d is not None:
                    self._files[p] = None  # type: ignore[assignment]
            self._mark_written(written)
            if repack:
                _repack_file(self._update_path, self._update_path)
            return self._update_path
        if path is None:
            raise TypeError("write() needs a path unless the writer came from "
                            "open_for_update()")
        if not str(path).endswith(".qvf"):
            path = f"{path}.qvf"
        with open(path, "wb") as fh:
            with zipfile.ZipFile(fh, "w", compression=compression) as zf:
                self._write_entries(zf, compression, force_zip64)
        return path


    def _mark_written(self, written: dict[str, Any]) -> None:
        """After an update is written: a later write() is another update,
        which stamps a new ``modified_utc``."""
        if "metadata" in written:
            self._root["metadata"] = dict(written["metadata"])
            self._update_loaded_modified = written["metadata"].get("modified_utc")
        self._update_stamp = None


#: Root keys build_manifest() writes from the writer's own state.
_MODELLED_ROOT_BLOCKS = frozenset({
    "metadata", "provenance", "thermochemistry", "dipole_moment", "constraints",
    "energies", "extensions", "viewer_defaults"})
_MAX_MANIFEST_BYTES = 256 * 1024 ** 2


def _repeated_names(infos: Iterable[zipfile.ZipInfo]) -> list[str]:
    """Entry names that occur more than once, except manifest.json, whose
    superseded copies an in-place update leaves (spec § 2.3)."""
    from collections import Counter
    counts = Counter(i.filename for i in infos)
    return sorted(n for n, k in counts.items() if k > 1 and n != "manifest.json")


def _same_file(a: str, b: str) -> bool:
    import os
    try:
        return os.path.samefile(a, b)
    except OSError:
        return os.path.normpath(a) == os.path.normpath(b)


def _zip_write(zf: zipfile.ZipFile, name: Any, data: bytes, compression: int,
               force_zip64: bool) -> None:
    """One whole entry, as ``ZipFile.writestr`` writes it; with
    ``force_zip64`` the local header carries a ZIP64 extra field."""
    import time
    if isinstance(name, zipfile.ZipInfo):
        info = name
    else:
        info = zipfile.ZipInfo(name, date_time=time.localtime(time.time())[:6])
        info.external_attr = 0o600 << 16
    info.compress_type = compression
    if not force_zip64:
        zf.writestr(info, data)
        return
    info.file_size = len(data)
    with zf.open(info, "w", force_zip64=True) as fh:
        fh.write(data)


def repack(src: str, dst: str) -> str:
    """Rewrite the archive ``src`` as ``dst`` without orphans (spec § 2.1).

    Keeps the entries the current (last) manifest references, in archive
    order and with their compression method, then that ``manifest.json``
    (STORED, last). Orphans and superseded manifests are dropped (spec
    § 2.3); member bytes and the manifest are copied unchanged, so every
    sha256 and size stays valid. ``src`` and ``dst`` may name the same file. Returns
    ``dst``.
    """
    import os
    in_place = os.path.exists(dst) and _same_file(os.fspath(src), os.fspath(dst))
    target = f"{dst}.repack-tmp" if in_place else dst
    _repack_zip(src, target)
    if in_place:
        os.replace(target, dst)
    return dst


#: :func:`repack` under a name that the ``repack`` keyword of
#: :meth:`QvfWriter.write` does not shadow.
_repack_file = repack


def _repack_zip(src: Any, dst: Any) -> None:
    """Write the archive ``src`` (a path or a file object) to ``dst`` from its
    current manifest: the referenced entries, then ``manifest.json``."""
    import shutil
    with zipfile.ZipFile(src) as zin:
        dups = _repeated_names(zin.infolist())
        if dups:
            raise ValueError(f"E_ZIP_DUPLICATE_ENTRY: duplicate ZIP entry name(s) "
                             f"{dups!r}; refusing to repack")
        manifest_info = zin.getinfo("manifest.json")    # the last one (spec § 2.3)
        if manifest_info.file_size > _MAX_MANIFEST_BYTES:
            raise ValueError("manifest.json is over the 256 MiB JSON cap")
        raw = zin.read(manifest_info)
        manifest = json.loads(raw)
        referenced = set()
        for s in manifest.get("sections", []) if isinstance(manifest, dict) else []:
            members = s.get("members") if isinstance(s, dict) else None
            for m in (members.values() if isinstance(members, dict) else []):
                if isinstance(m, dict) and isinstance(m.get("path"), str):
                    referenced.add(m["path"])
        with zipfile.ZipFile(dst, "w") as zout:
            for info in zin.infolist():
                if info.filename == "manifest.json" or info.filename not in referenced:
                    continue
                out = zipfile.ZipInfo(info.filename, date_time=info.date_time)
                out.compress_type = info.compress_type
                out.external_attr = info.external_attr
                out.file_size = info.file_size    # lets zipfile pick ZIP64 when needed
                with zin.open(info) as fsrc, zout.open(out, "w") as fdst:
                    shutil.copyfileobj(fsrc, fdst, 1 << 20)
            m_info = zipfile.ZipInfo("manifest.json", date_time=manifest_info.date_time)
            zout.writestr(m_info, raw, compress_type=zipfile.ZIP_STORED)
