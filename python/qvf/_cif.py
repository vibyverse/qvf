"""Core crystallographic CIF 1.1 interchange; see docs/converters.md.

Import reads one data block's cell, fractional atom sites and space-group
operators and writes every atom of the unit cell; export writes a periodic
structure as a P 1 CIF. mmCIF (PDBx) dictionaries are not interpreted.
"""
from __future__ import annotations

from fractions import Fraction
import math
import re

from . import __version__
from .convert import (_atom, _export_report, _lattice_inverse, _new_output, _number,
                      _pdb_cell_vectors, _pdb_scale, _select, _text, _validated)
from ._legacy_writer import _ELEMENT_SYMBOLS
from .reader import Reader
from .writer import Writer

_KEYWORDS = ("data_", "loop_", "save_", "global_", "stop_")
_OPERATIONS = ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
_HALL = ("_space_group_name_hall", "_symmetry_space_group_name_hall")
_HM = ("_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m")
_NUMBER = ("_space_group_it_number", "_symmetry_int_tables_number")
_CELL = ("_cell_length_a", "_cell_length_b", "_cell_length_c",
         "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma")
_SITE = ("_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x",
         "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy",
         "_atom_site_u_iso_or_equiv", "_atom_site_b_iso_or_equiv",
         "_atom_site_adp_type", "_atom_site_thermal_displace_type",
         "_atom_site_disorder_assembly", "_atom_site_disorder_group")
# The anisotropic displacement loop (CIF core `_atom_site_aniso.U_11` ...,
# `B_11` ...), components in the CIF order 11, 22, 33, 12, 13, 23.
_ANISO_ORDER = ("11", "22", "33", "12", "13", "23")
_ANISO_LABEL = "_atom_site_aniso_label"
_LABEL = re.compile(r"[A-Za-z][A-Za-z0-9_+-]*")


def _tokens(text):
    """(quoted?, value) tokens of CIF 1.1 syntax: comments, quoted values and
    semicolon text fields; a quote closes only before whitespace."""
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            body = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                body.append(lines[i]); i += 1
            if i == len(lines):
                raise ValueError("CIF semicolon text field is not closed")
            yield True, "\n".join(body).strip("\n")
            line = lines[i][1:]
        pos, n = 0, len(line)
        while pos < n:
            c = line[pos]
            if c.isspace():
                pos += 1
            elif c == "#":
                break
            elif c in "'\"":
                end = pos + 1
                while True:
                    end = line.find(c, end)
                    if end < 0:
                        raise ValueError("CIF quoted value is not closed on its line")
                    if end + 1 == n or line[end + 1].isspace():
                        break
                    end += 1
                yield True, line[pos + 1:end]
                pos = end + 1
            else:
                end = pos
                while end < n and not line[end].isspace():
                    end += 1
                yield False, line[pos:end]
                pos = end
        i += 1


def _blocks(text):
    """[(name, items {tag: token}, loops [(tags, rows)])]; tags lowercase."""
    if any(ord(c) < 32 and c not in "\t\r\n" for c in text):
        raise ValueError("CIF requires text without control characters")
    tokens = list(_tokens(text))
    blocks, k = [], 0

    def keyword(token):
        return not token[0] and (token[1].startswith("_")
                                 or token[1].lower().startswith(_KEYWORDS))

    while k < len(tokens):
        quoted, value = tokens[k]
        low = value.lower()
        if not quoted and low.startswith("data_"):
            if len(value) == 5 or any(b[0].lower() == low[5:] for b in blocks):
                raise ValueError("CIF data blocks need distinct nonempty names")
            blocks.append((value[5:], {}, [])); k += 1
            continue
        if not quoted and low.startswith(_KEYWORDS[2:]):
            raise ValueError(f"CIF {value} is not supported")
        if not blocks:
            raise ValueError("CIF has content before its first data_ block")
        _, items, loops = blocks[-1]
        seen = set(items) | {t for tags, _ in loops for t in tags}
        if not quoted and low == "loop_":
            k += 1; tags = []
            while k < len(tokens) and not tokens[k][0] and tokens[k][1].startswith("_"):
                tags.append(tokens[k][1].lower()); k += 1
            values = []
            while k < len(tokens) and not keyword(tokens[k]):
                values.append(tokens[k]); k += 1
            if not tags or len(values) % len(tags):
                raise ValueError("CIF loop values must fill whole rows of its tags")
            if seen & set(tags) or len(set(tags)) != len(tags):
                raise ValueError("CIF repeats a tag within one data block")
            loops.append((tags, [values[r:r + len(tags)]
                                 for r in range(0, len(values), len(tags))]))
        elif not quoted and value.startswith("_"):
            if k + 1 == len(tokens) or keyword(tokens[k + 1]):
                raise ValueError(f"CIF tag {value} has no value")
            if low in seen:
                raise ValueError("CIF repeats a tag within one data block")
            items[low] = tokens[k + 1]; k += 2
        else:
            raise ValueError(f"CIF value {value!r} has no tag")
    return blocks


def _reciprocal_lengths(inverse):
    """a*_j, the reciprocal-cell lengths, from the inverse of the lattice (rows
    are cell vectors): the reciprocal vector b_j is column j of the inverse."""
    return [math.sqrt(sum(inverse[k][j] ** 2 for k in range(3))) for j in range(3)]


def _sym6(u):
    """The symmetric 3 x 3 matrix of the components [11, 22, 33, 12, 13, 23]."""
    return [[u[0], u[3], u[4]], [u[3], u[1], u[5]], [u[4], u[5], u[2]]]


def _six(m):
    return [m[0][0], m[1][1], m[2][2], m[0][1], m[0][2], m[1][2]]


def _congruent(t, m):
    """t m tᵀ for 3 x 3 matrices."""
    tm = [[sum(t[i][k] * m[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    return [[sum(tm[i][k] * t[j][k] for k in range(3)) for j in range(3)] for i in range(3)]


def _rotate_crystal_u(u, rotation, nstar):
    """The crystal-frame tensor (CIF U^ij, spec § 5.1, Crystallographic atom
    sites) of the copy that the operation with fractional rotation R makes:
    N⁻¹ R N U N Rᵀ N⁻¹, N = diag(a*). The fractional-basis tensor N U N
    transforms as R (N U N) Rᵀ, like x' = R x (Grosse-Kunstleve & Adams,
    J. Appl. Cryst. 35, 477 (2002), doi:10.1107/S0021889802008580)."""
    scaled = [[nstar[i] * v * nstar[j] for j, v in enumerate(row)]
              for i, row in enumerate(_sym6(u))]
    turned = _congruent(rotation, scaled)
    return [round(v, 10) + 0.0 for v in
            _six([[turned[i][j] / (nstar[i] * nstar[j]) for j in range(3)] for i in range(3)])]


def _positive_definite(u):
    m = _sym6(u)
    minor = m[0][0] * m[1][1] - m[0][1] ** 2
    det = (m[0][0] * (m[1][1] * m[2][2] - m[1][2] ** 2) - m[0][1] * (m[0][1] * m[2][2] - m[1][2] * m[0][2])
           + m[0][2] * (m[0][1] * m[1][2] - m[1][1] * m[0][2]))
    return m[0][0] > 0 and minor > 0 and det > 0


def _code(token):
    """A CIF code value (a disorder assembly or group), else None."""
    return None if _missing(token) or not token[1].strip() else token[1].strip()


def _aniso_rows(block):
    """{site label: (type, [six components])} of the anisotropic loop, U
    components preferred over B (CIF core `_atom_site_aniso.U_ij`, `B_ij`)."""
    found = _column(block, _ANISO_LABEL)
    if not found:
        return {}
    rows, at = found
    tags = next(t for t, r in block[2] if r is rows)
    out = {}
    for letter in ("u", "b"):
        names = [f"_atom_site_aniso_{letter}_{ij}" for ij in _ANISO_ORDER]
        if not all(n in tags for n in names):
            continue
        for row in rows:
            label = row[at][1]
            values = [_cif_number(row[tags.index(n)], n) for n in names]
            if label not in out and None not in values:
                out[label] = ("Uani" if letter == "u" else "Bani", values)
    return out


def _missing(token):
    return token is None or (not token[0] and token[1] in ("?", "."))


def _cif_number(token, label):
    """A CIF number, its standard uncertainty in parentheses dropped."""
    if _missing(token):
        return None
    match = re.fullmatch(r"([+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)"
                         r"(?:\([0-9]+\))?", token[1])
    if not match:
        raise ValueError(f"CIF {label} requires a number")
    return _number(match.group(1))


def _item(block, tags):
    """The first single-valued item of ``tags``, else None."""
    return next((block[1][t] for t in tags if t in block[1]), None)


def _column(block, tag):
    """(loop rows, column index) of ``tag`` in a loop, else None."""
    for tags, rows in block[2]:
        if tag in tags:
            return rows, tags.index(tag)
    return None


def _values(block, tag):
    if tag in block[1]:
        return [block[1][tag]]
    found = _column(block, tag)
    return [row[found[1]] for row in found[0]] if found else []


def _operation(text):
    """(R, t) of a symmetry operation such as -x+1/2,y,-z+1/2."""
    parts = text.replace(" ", "").lower().split(",")
    if len(parts) != 3:
        raise ValueError(f"CIF symmetry operation {text!r} needs three components")
    r = [[0, 0, 0] for _ in range(3)]
    t = [Fraction(0)] * 3
    for i, part in enumerate(parts):
        terms = re.findall(r"[+-]?[^+-]+", part)
        if not terms or "".join(terms) != part:
            raise ValueError(f"CIF symmetry operation {text!r} is malformed")
        for term in terms:
            sign = -1 if term.startswith("-") else 1
            body = term.lstrip("+-")
            if body in ("x", "y", "z"):
                r[i]["xyz".index(body)] += sign
            elif re.fullmatch(r"[0-9]+/[1-9][0-9]*|[0-9]*\.[0-9]+|[0-9]+", body):
                t[i] += sign * Fraction(body)
            else:
                raise ValueError(f"CIF symmetry operation {text!r} is malformed")
    det = (r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
           - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
           + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]))
    if abs(det) != 1:
        raise ValueError(f"CIF symmetry operation {text!r} is not a lattice symmetry")
    return tuple(tuple(row) for row in r), tuple(v % 1 for v in t)


def _compose(a, b):
    (ra, ta), (rb, tb) = a, b
    r = tuple(tuple(sum(ra[i][k] * rb[k][j] for k in range(3)) for j in range(3))
              for i in range(3))
    return r, tuple((sum(ra[i][k] * tb[k] for k in range(3)) + ta[i]) % 1 for i in range(3))


def _operators(block):
    """(operators, source): the listed operations, else those of the Hall
    symbol, else the identity when the group is P 1."""
    listed = next((_values(block, tag) for tag in _OPERATIONS if _values(block, tag)), [])
    identity = (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (Fraction(0),) * 3)
    if listed:
        if any(_missing(v) for v in listed):
            raise ValueError("CIF symmetry operations must not be unknown")
        ops = [_operation(v[1]) for v in listed]
        if len(set(ops)) != len(ops) or identity not in ops:
            raise ValueError("CIF symmetry operations must be distinct and include x,y,z")
        known = set(ops)
        if any(_compose(a, b) not in known for a in ops for b in ops):
            raise ValueError("CIF symmetry operations do not form a group")
        return ops, "listed"
    hall = _item(block, _HALL)
    if not _missing(hall):
        from ._spacegroups import operators
        try:
            return operators(hall[1]), "hall"
        except (ValueError, KeyError, IndexError) as error:
            raise ValueError(f"CIF Hall symbol {hall[1]!r} is not supported: {error}") from None
    hm, number = _item(block, _HM), _item(block, _NUMBER)
    if (not _missing(hm) and hm[1].replace(" ", "") == "P1") or \
            (not _missing(number) and number[1] == "1"):
        return [identity], "p1"
    raise ValueError("CIF lists no symmetry operations or Hall symbol, and the space "
                     "group is not P 1")


def _element(token):
    """(element, charge given?) of an _atom_site_type_symbol such as Fe3+."""
    match = None if _missing(token) else re.fullmatch(r"([A-Za-z]{1,2})([0-9]*[+-][0-9]*)?",
                                                      token[1])
    symbol = match.group(1).capitalize() if match else None
    if symbol not in _ELEMENT_SYMBOLS[1:]:
        raise ValueError("CIF atom sites require an element _atom_site_type_symbol "
                         "(no element is guessed from labels)")
    return symbol, bool(match.group(2))


def _select_block(blocks, block):
    with_sites = [b for b in blocks if _column(b, "_atom_site_label")]
    if block is not None:
        chosen = [b for b in blocks if b[0].lower() == block.lower()]
        if not chosen:
            raise ValueError(f"CIF has no data block {block!r}")
        if chosen[0] not in with_sites:
            raise ValueError(f"CIF data block {block!r} has no atom sites")
        return chosen[0]
    if len(with_sites) != 1:
        raise ValueError("select a CIF data block (block=) when zero or several "
                         "blocks have atom sites")
    return with_sites[0]


def import_cif(source, destination, *, block=None):
    blocks = _blocks(_text(source))
    if not blocks:
        raise ValueError("CIF has no data block")
    chosen = _select_block(blocks, block)
    cell = [_cif_number(_item(chosen, (tag,)), tag) for tag in _CELL]
    if any(v is None for v in cell) or not all(v > 0 for v in cell[:3]) \
            or not all(0 < v < 180 for v in cell[3:]):
        raise ValueError("CIF requires positive cell lengths and angles in (0, 180)")
    lattice = _pdb_cell_vectors(*cell)
    ops, source_kind = _operators(chosen)
    rows, _ = _column(chosen, "_atom_site_label")
    tags = next(t for t, r in chosen[2] if r is rows)
    column = {tag: tags.index(tag) if tag in tags else None for tag in _SITE}
    if any(column[t] is None for t in _SITE[2:5]):
        raise ValueError("CIF atom sites require fractional coordinates "
                         "(Cartesian sites are not supported)")
    if column["_atom_site_type_symbol"] is None:
        raise ValueError("CIF atom sites require an element _atom_site_type_symbol "
                         "(no element is guessed from labels)")
    if not rows:
        raise ValueError("CIF has no atom sites")
    scale = _pdb_scale(lattice)
    nstar = _reciprocal_lengths(_lattice_inverse(lattice))
    aniso = _aniso_rows(chosen)
    adp_dropped, anisotropic, disordered = 0, 0, 0

    def cartesian(f):
        return [sum(f[i] * lattice[i][k] for i in range(3)) for k in range(3)]

    def apart(f, g):
        return math.dist(cartesian([a - b - round(a - b) for a, b in zip(f, g)]), [0, 0, 0])

    atoms, merged, charges, partial = [], 0, 0, 0
    float_ops = [(r, [float(v) for v in t]) for r, t in ops]
    for row in rows:
        def get(tag):
            return None if column[tag] is None else row[column[tag]]
        label = get("_atom_site_label")
        if _missing(label):
            raise ValueError("CIF atom sites require labels")
        symbol, charged = _element(get("_atom_site_type_symbol"))
        charges += charged
        f = [_cif_number(get(t), t) for t in _SITE[2:5]]
        if None in f:
            raise ValueError(f"CIF atom site {label[1]} has no fractional coordinates")
        occupancy = _cif_number(get("_atom_site_occupancy"), "occupancy")
        if occupancy is not None and not 0 < occupancy <= 1:
            raise ValueError(f"CIF atom site {label[1]} occupancy must lie in (0, 1]")
        partial += occupancy is not None and occupancy < 1
        u = _cif_number(get("_atom_site_u_iso_or_equiv"), "U_iso")
        b = _cif_number(get("_atom_site_b_iso_or_equiv"), "B_iso")
        if b is None and u is not None:
            b = 8 * math.pi ** 2 * u
        # Displacement parameters and disorder grouping (spec § 5.1,
        # Crystallographic atom sites; CIF core _atom_site.ADP_type,
        # _atom_site_aniso.U_ij / B_ij, _atom_site.disorder_assembly and
        # _atom_site.disorder_group).
        declared = _code(get("_atom_site_adp_type") or get("_atom_site_thermal_displace_type"))
        adp = None
        if label[1] in aniso and (declared or "").lower() not in ("uiso", "biso"):
            kind, components = aniso[label[1]]
            if _positive_definite(components):
                adp = {"type": kind, "frame": "crystal", "aniso": components}
                if kind == "Bani":
                    equiv = b                      # B(equiv), or 8π² U(equiv)
                else:
                    equiv = u if u is not None or b is None else b / (8 * math.pi ** 2)
                if equiv is not None and equiv > 0:
                    adp["iso"] = equiv
            else:
                adp_dropped += 1                   # keep the isotropic value below
        if adp is None and (u is not None or b is not None):
            recorded_b = _cif_number(get("_atom_site_b_iso_or_equiv"), "B_iso")
            biso = recorded_b is not None and ((declared or "").lower() == "biso" or u is None)
            value = recorded_b if biso else u
            if value is not None and value > 0:
                adp = {"type": "Biso" if biso else "Uiso", "iso": value}
            else:
                adp_dropped += 1
        group = _code(get("_atom_site_disorder_group"))
        assembly = _code(get("_atom_site_disorder_assembly")) if group else None
        anisotropic += adp is not None and "aniso" in adp
        disordered += group is not None
        images = []
        rotations = []
        for r, t in float_ops:
            g = [round(sum(r[i][k] * f[k] for k in range(3)) + t[i], 10) % 1.0
                 for i in range(3)]
            if any(apart(g, h) < 0.1 for h in images):
                merged += 1
                continue
            images.append(g)
            rotations.append(r)
        for g, r in zip(images, rotations):
            atom = _atom(symbol, [round(v, 6) + 0.0 for v in cartesian(g)])
            atom["atom_name"] = label[1]
            if occupancy is not None:
                atom["occupancy"] = occupancy
            if b is not None:
                atom["b_factor"] = b
            if adp is not None:
                atom["adp"] = dict(adp)
                if "aniso" in adp:
                    atom["adp"]["aniso"] = _rotate_crystal_u(adp["aniso"], r, nstar)
            if group is not None:
                if assembly is not None:
                    atom["disorder_assembly"] = assembly
                atom["disorder_group"] = group
            atoms.append(atom)
    symmetry = {}
    hm, number, hall = _item(chosen, _HM), _item(chosen, _NUMBER), _item(chosen, _HALL)
    if not _missing(hm):
        symmetry["space_group_symbol"] = " ".join(hm[1].split())
    if not _missing(number):
        if not re.fullmatch(r"[0-9]+", number[1]) or not 1 <= int(number[1]) <= 230:
            raise ValueError("CIF space-group number must lie in 1-230")
        symmetry["space_group_number"] = int(number[1])
    if not _missing(hall):
        symmetry["hall"] = " ".join(hall[1].split())
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "CIF crystal structure import") as writer:
            writer.add_section({"id": "structure", "kind": "structure",
                                "label": f"data_{chosen[0]}"})
            writer.add_json("structure", "structure", {"atoms": atoms, "pbc": [True] * 3,
                                                       "lattice_vectors": lattice})
            if symmetry:
                writer.add_section({"id": "symmetry", "kind": "structure.symmetry",
                                    "structure_ref": "structure"})
                writer.add_json("symmetry", "data", symmetry)
        findings = _validated(temporary)
    used = (set(_CELL) | set(_SITE) | set(_OPERATIONS) | set(_HALL) | set(_HM) | set(_NUMBER)
            | {_ANISO_LABEL} | {f"_atom_site_aniso_{x}_{ij}" for x in "ub" for ij in _ANISO_ORDER})
    present = set(chosen[1]) | {t for tags, _ in chosen[2] for t in tags}
    return {"format": "cif", "output": str(destination), "block": chosen[0],
            "blocks": [b[0] for b in blocks], "sections": ["structure"]
            + (["symmetry"] if symmetry else []), "findings": findings,
            "sites": len(rows), "operators": len(ops), "operator_source": source_kind,
            "atoms": len(atoms), "special_position_atoms_merged": merged,
            "partial_occupancy_sites": partial, "charges_omitted": charges,
            "anisotropic_sites": anisotropic, "disordered_sites": disordered,
            "displacement_parameters_dropped": adp_dropped,
            "space_group": symmetry.get("space_group_symbol"),
            "unused_tags": len(present - used),
            "scope": "Cell, fractional atom sites (label, element, occupancy, isotropic "
                     "and anisotropic displacement parameters, disorder assembly and "
                     "group) and the unit cell they generate under the space group, in "
                     "the standard orientation (a along x, b in the xy-plane); each copy "
                     "carries the displacement tensor of its own orientation. A tensor "
                     "that is not positive definite falls back to the isotropic value, "
                     "a value that is not positive is dropped; both are counted. No "
                     "bonds, disorder resolution, charge or job is inferred. Keep the "
                     "original CIF."}


def _export_adp(adp, transpose, nstar):
    """(CIF ADP_type, U_iso_or_equiv or None, crystal-frame U^ij or None) of an
    atom's ``adp``, written as U as the IUCr recommends (B = 8π² U); a
    Cartesian tensor is turned into the CIF frame, U = N⁻¹ A⁻¹ U_cart A⁻ᵀ N⁻¹
    with A the lattice vectors as columns (spec § 5.1, Crystallographic atom
    sites). None without ``adp``."""
    if not isinstance(adp, dict):
        return None
    scale = 1 / (8 * math.pi ** 2) if adp.get("type") in ("Biso", "Bani") else 1.0
    iso = adp.get("iso")
    iso = None if iso is None else iso * scale
    u = adp.get("aniso")
    if u is None:
        return ("Uiso", iso, None)
    u = [v * scale for v in u]
    if adp.get("frame") == "cartesian":
        frac = _congruent(transpose, _sym6(u))
        u = _six([[frac[i][j] / (nstar[i] * nstar[j]) for j in range(3)] for i in range(3)])
    return ("Uani", iso, u)


def export_cif(source, destination, *, section=None):
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"structure"})
        geometry = reader.json(selected.id, "structure")
        lattice = geometry.get("lattice_vectors")
        if not all(geometry.get("pbc", [False] * 3)) or lattice is None:
            raise ValueError("CIF export requires a structure periodic on all three axes")
        rows = [[_number(v) for v in row] for row in lattice]
        inverse = _lattice_inverse(rows)
        lengths = [math.sqrt(sum(v * v for v in row)) for row in rows]

        def angle(u, v):
            dot = sum(x * y for x, y in zip(rows[u], rows[v])) / (lengths[u] * lengths[v])
            return math.degrees(math.acos(max(-1.0, min(1.0, dot))))

        cell = lengths + [angle(1, 2), angle(0, 2), angle(0, 1)]
        atoms = geometry["atoms"]
        if not atoms:
            raise ValueError("CIF export requires at least one atom")
        names = [a.get("atom_name") for a in atoms]
        named = all(isinstance(n, str) and _LABEL.fullmatch(n) for n in names) \
            and len(set(names)) == len(names)
        counts, labels = {}, []
        for a in atoms:
            symbol = a["symbol"]
            counts[symbol] = counts.get(symbol, 0) + 1
            labels.append(f"{symbol}{counts[symbol]}")
        if named:
            labels = names
        has_occupancy = any(a.get("occupancy") is not None for a in atoms)
        has_b = any(a.get("b_factor") is not None for a in atoms)
        nstar = _reciprocal_lengths(inverse)
        transpose = [[inverse[k][i] for k in range(3)] for i in range(3)]
        adps = [_export_adp(a.get("adp"), transpose, nstar) for a in atoms]
        has_adp = any(x is not None for x in adps)
        has_group = any(a.get("disorder_group") is not None for a in atoms)
        has_assembly = any(a.get("disorder_assembly") is not None for a in atoms)
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", selected.id)
        lines = [f"# CIF written by qvf.convert {__version__}: P 1 cell, no symmetry applied",
                 f"data_{name}"]
        lines += [f"{tag:<20s}{value:.6f}" for tag, value in
                  zip(("_cell_length_a", "_cell_length_b", "_cell_length_c",
                       "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"), cell)]
        lines += ["_space_group_name_H-M_alt 'P 1'", "_space_group_IT_number 1", "loop_",
                  "_space_group_symop_operation_xyz", "x,y,z", "loop_", "_atom_site_label",
                  "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y",
                  "_atom_site_fract_z"]
        lines += ["_atom_site_occupancy"] * has_occupancy + ["_atom_site_B_iso_or_equiv"] * has_b
        lines += ["_atom_site_U_iso_or_equiv", "_atom_site_adp_type"] * has_adp
        lines += ["_atom_site_disorder_assembly"] * has_assembly
        lines += ["_atom_site_disorder_group"] * has_group
        error = 0.0
        for index, (label, atom) in enumerate(zip(labels, atoms)):
            x = [_number(v) for v in atom["position"]]
            f = [sum(x[k] * inverse[k][i] for k in range(3)) for i in range(3)]
            text = [f"{v:.8f}" for v in f]
            back = [sum(float(text[i]) * rows[i][k] for i in range(3)) for k in range(3)]
            error = max(error, math.dist(back, x))
            fields = [label, atom["symbol"]] + text
            if has_occupancy:
                fields.append("?" if atom.get("occupancy") is None else f"{atom['occupancy']:.4f}")
            if has_b:
                fields.append("?" if atom.get("b_factor") is None else f"{atom['b_factor']:.4f}")
            if has_adp:
                adp = adps[index]
                fields += (["?", "?"] if adp is None else
                           ["?" if adp[1] is None else f"{adp[1]:.6f}", adp[0]])
            if has_assembly:
                fields.append(atom.get("disorder_assembly") or ".")
            if has_group:
                fields.append(atom.get("disorder_group") or ".")
            lines.append(" ".join(fields))
        tensors = [(label, adp[2]) for label, adp in zip(labels, adps)
                   if adp is not None and adp[2] is not None]
        if tensors:
            lines += ["loop_", "_atom_site_aniso_label"]
            lines += [f"_atom_site_aniso_U_{ij}" for ij in _ANISO_ORDER]
            lines += [" ".join([label] + [f"{v:.6f}" for v in u]) for label, u in tensors]
        with _new_output(destination) as temporary:
            temporary.write_text("\n".join(lines) + "\n", encoding="ascii")
        report = _export_report(reader, selected, {selected.id}, destination, "cif")
        report.update(atoms=len(atoms), labels="atom_name" if named else "generated",
                      max_rounding_error_angstrom=error,
                      generated="P 1 with the identity operator only; a reader places the "
                                "cell in its own standard orientation, so Cartesian "
                                "coordinates may be rotated relative to the archive.")
        return report
