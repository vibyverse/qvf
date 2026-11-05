"""Bounded molecular Molden interchange; see docs/converters.md for conventions.

The format documentation defines the component orders. The normalized-AO
interpretation is independently exercised against PySCF's Molden reader/writer.
No reference-program source is incorporated here.
"""
from __future__ import annotations

import math
import re

from . import __version__
from .convert import (_MAX_VALUES, _atom, _binary, _new_output, _number, _select,
                      _structure, _text, _validated)
from ._legacy_writer import BOHR_PER_ANGSTROM, primitive_norm
from .reader import Reader
from .writer import Writer

_SPINS = {"restricted", "restricted_open", "unrestricted"}
_FLAGS = {"5d", "6d", "5d7f", "5d10f", "7f", "10f", "9g", "15g"}

# Producer dialects: each maps what a program writes onto the standard Molden
# meaning (unit-normalized contractions and AOs, Molden phases), which then
# maps onto QVF Appendix A. The table in docs/converters.md lists the sources.
# The producer-specific rules follow the documented corrections of IOData
# 1.0.1, iodata/formats/molden.py (_fix_obasis_orca, _fix_obasis_turbomole,
# _fix_mo_coeffs_cfour, _fix_mo_coeffs_psi4); no IOData code is used here.
# Unlike IOData, nothing is guessed from orbital norms: the caller names the
# dialect, and a contraction that is still not normalized is refused.
DIALECTS = ("standard", "renormalize", "orca", "turbomole", "cfour", "psi4-1.3", "xtb")


def _shell_norm(shell):
    """Axial contracted-AO norm, Appendix A.1; no per-component scaling."""
    es = shell["exponents"]; cs = shell["coefficients"]
    if not 1 <= len(es) <= 512 or len(es) != len(cs):
        raise ValueError("Molden requires 1-512 primitives per contraction")
    if any(_number(e) <= 0 for e in es):
        raise ValueError("Molden exponents must be positive")
    # Ratio form avoids overflowing alpha_i * alpha_j or alpha_i + alpha_j.
    norm2 = math.fsum(_number(a) * _number(b) *
                     (2 * math.sqrt(min(x, y) / max(x, y)) /
                      (1 + min(x, y) / max(x, y))) ** (shell["l"] + 1.5)
                     for x, a in zip(es, cs) for y, b in zip(es, cs))
    if not math.isfinite(norm2) or norm2 <= 1e-14:
        raise ValueError("Molden requires a finite, nonzero contraction norm")
    return math.sqrt(norm2)


# Molden's Cartesian component orders (exponents of x, y, z).
_CARTESIAN = {
    2: ["xx", "yy", "zz", "xy", "xz", "yz"],
    3: ["xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"],
    4: ["xxxx", "yyyy", "zzzz", "xxxy", "xxxz", "yyyx", "yyyz", "zzzx", "zzzy", "xxyy",
        "xxzz", "yyzz", "xxyz", "yyxz", "zzxy"],
}


def _double_factorial(n):
    return math.prod(range(n, 0, -2)) if n > 0 else 1


def _layout(shell):
    """QVF indices in Molden order; QVF AO / unit-normalized Molden AO."""
    angular = shell["l"]; pure = shell.get("pure", True)
    if angular == 0:
        return [0], [1.0]
    if angular == 1:
        return ([2, 0, 1] if pure else [0, 1, 2]), [1.0] * 3
    if angular > 4:
        raise ValueError("Molden conversion supports s, p, sp, d, f and g shells only")
    if pure:
        # Molden m = 0, +1, -1, ..., +l, -l; QVF m = -l, ..., +l.
        ms = [0] + [s * k for k in range(1, angular + 1) for s in (1, -1)]
        return [m + angular for m in ms], [1.0] * len(ms)
    qvf = [(lx, ly, angular - lx - ly) for lx in range(angular, -1, -1)
           for ly in range(angular - lx, -1, -1)]
    indices, factors = [], []
    for name in _CARTESIAN[angular]:
        powers = tuple(name.count(a) for a in "xyz")
        indices.append(qvf.index(powers))
        # QVF Cartesian AOs carry the axial normalization (Appendix A.1).
        factors.append(math.sqrt(math.prod(_double_factorial(2 * n - 1) for n in powers)
                                 / _double_factorial(2 * angular - 1)))
    return indices, factors


def _sections(text):
    sections = {}; name = None; omitted_comments = 0
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith('#'):
            omitted_comments += 1
            continue
        header = re.fullmatch(r"\[([^]]+)\](.*)", line)
        if header:
            name = header[1].lower().strip()
            if (not sections and name != "molden format") or name in sections:
                raise ValueError("Molden needs its format header first and unique sections")
            if name in {"core", "pseudo"}:
                raise ValueError(f"Molden [{header[1]}] marks ECP cores; the ECP parameters are not in "
                                 "the file, so QVF basis.ecp cannot be built and the orbitals would "
                                 "be wrong without it")
            if name not in {"molden format", "atoms", "gto", "mo", "title"} | _FLAGS:
                raise ValueError(f"unsupported Molden section {name!r}")
            suffix = header[2].strip().lower().strip('()')
            if (name == "atoms" and suffix not in {"au", "angs"}) or (name != "atoms" and suffix):
                raise ValueError("Molden needs explicit Atoms AU/Angs units and plain section headers")
            sections[name] = [suffix, []]
        elif name is None:
            raise ValueError("Molden needs its format header first")
        else:
            sections[name][1].append(line)
    if not {"molden format", "atoms", "gto", "mo"} <= sections.keys():
        raise ValueError("Molden requires Atoms, GTO and MO sections")
    if any(sections[flag][1] for flag in _FLAGS & sections.keys()):
        raise ValueError("Molden angular flags must have no payload")
    flags = _FLAGS & sections.keys()
    pure_d = bool({"5d", "5d7f", "5d10f"} & flags)
    # [5D] also means spherical f; [5D10F] keeps f Cartesian.
    pure_f = bool({"5d", "5d7f", "7f"} & flags)
    pure_g = "9g" in flags
    if (pure_d and "6d" in flags) or (pure_f and {"10f", "5d10f"} & flags) \
            or (pure_g and "15g" in flags):
        raise ValueError("conflicting Molden shell conventions")
    return sections, {2: pure_d, 3: pure_f, 4: pure_g}, omitted_comments


def _primitive_scale(dialect, angular, pure, exponent):
    """Written / standard contraction coefficient of one primitive."""
    if dialect == "orca":
        # IOData _fix_obasis_orca: ORCA multiplies by the norm of the Cartesian
        # primitive (0,0,0), (1,0,0), (1,1,0), (1,1,1), (2,1,1) for s to g,
        # which is N_i * sqrt((2l-1)!! / prod (2l_k-1)!!) with Appendix A.1 N_i.
        product = {0: 1, 1: 1, 2: 1, 3: 1, 4: 3}[angular]
        return primitive_norm(exponent, angular) * math.sqrt(_double_factorial(2 * angular - 1) / product)
    if dialect == "xtb":
        # QVF Appendix A.1. xTB printmold writes basis%cont with the axial
        # primitive normalization already multiplied into it.
        return primitive_norm(exponent, angular)
    if dialect == "turbomole" and not pure and angular >= 2:
        # IOData _fix_obasis_turbomole: tm2molden Cartesian d/f/g contraction
        # coefficients are 1/sqrt(3), 1/sqrt(15), 1/sqrt(105) of the
        # normalized ones.
        return 1 / math.sqrt(_double_factorial(2 * angular - 1))
    return 1.0


def _component_factors(dialect, shell, factors):
    """Adjust the standard per-component factors (file AO / QVF AO)."""
    angular, pure = shell["l"], shell["pure"]
    if dialect == "xtb":
        # xTB cube.f90/primval uses R_atom-r; QVF uses r-R_atom.
        # printmold emits those MO coefficients unchanged. Reverse odd-l
        # phases while preserving Molden's Cartesian component scaling.
        return [factor * (-1) ** angular for factor in factors]
    if dialect == "orca" and pure and angular >= 3:
        # IOData orca_conventions: in Molden order c0 c1 s1 c2 s2 c3 s3 c4 s4,
        # ORCA's c3, s3 (f, g) and c4, s4 (g) are the negated Molden functions.
        return [-factor if k >= 5 else factor for k, factor in enumerate(factors)]
    if not pure and angular >= 2 and dialect == "cfour":
        # IOData _fix_mo_coeffs_cfour (CFOUR up to 2.1): written Cartesian MO
        # coefficients are the Molden ones times 1/sqrt(prod (2l_k-1)!!), so
        # the file AO is a common 1/sqrt((2l-1)!!) multiple of the QVF AO.
        return [1 / math.sqrt(_double_factorial(2 * angular - 1))] * len(factors)
    if not pure and angular >= 2 and dialect == "psi4-1.3":
        # IOData _fix_mo_coeffs_psi4 (Psi4 <= 1.3.2): written Cartesian MO
        # coefficients are the Molden ones times sqrt((2l-1)!! / prod
        # (2l_k-1)!!), which are coefficients on QVF's axially normalized AOs.
        return [1.0] * len(factors)
    return factors


def _atom_line(fields, index, factor):
    """One [Atoms] row: label, index, atomic number, x, y, z."""
    if len(fields) != 6 or not re.fullmatch(r"[0-9]+", fields[1]) or int(fields[1]) != index:
        raise ValueError("Molden atoms require sequential 1-based indices and six fields")
    charge = _number(fields[2])
    # Molpro and others write labels such as O1 or H12; the element is the
    # leading letters, which must agree with the atomic number.
    letters = re.match(r"[A-Za-z]*", fields[0])[0].title()
    atom = _atom(letters, [_number(v) * factor for v in fields[3:]])
    if charge != atom["atomic_number"]:
        raise ValueError(f"Molden atom {index} ({fields[0]}) has nuclear charge {fields[2]}, not "
                         f"{atom['atomic_number']}: an ECP effective core charge, ghost or point "
                         "charge, which this import cannot represent")
    return atom


def _parse(text, spin, dialect="standard"):
    if spin not in _SPINS:
        raise ValueError("Molden import requires explicit restricted, restricted_open or unrestricted spin")
    if dialect not in DIALECTS:
        raise ValueError(f"unknown Molden dialect {dialect!r}; choose one of {', '.join(DIALECTS)}")
    sections, pure, comments = _sections(text)
    title = " ".join(sections.get("title", [None, []])[1] + sections["molden format"][1]).lower()
    atoms = []
    factor = 1 / BOHR_PER_ANGSTROM if sections["atoms"][0] == "au" else 1
    for line in sections["atoms"][1]:
        atoms.append(_atom_line(line.split(), len(atoms) + 1, factor))
    if not atoms:
        raise ValueError("Molden has no atoms")
    rows = sections["gto"][1]; index = 0; center = None; centers = set()
    shells = []; norms = []; values = 0
    while index < len(rows):
        fields = rows[index].split(); index += 1
        if re.fullmatch(r"[0-9]+", fields[0]):
            center = int(fields[0]) - 1
            if (len(fields) not in (1, 2) or (len(fields) == 2 and fields[1] != "0")
                    or not 0 <= center < len(atoms) or center in centers):
                raise ValueError("Molden GTO requires unique valid atom blocks")
            centers.add(center)
            continue
        label = fields[0].lower()
        if center is None or label not in {"s", "p", "sp", "d", "f", "g"} or len(fields) not in (2, 3):
            raise ValueError("Molden conversion supports s, p, sp, d, f and g atom-centered "
                             "contractions only")
        if len(fields) == 3 and _number(fields[2]) != 1:
            raise ValueError("non-unit Molden shell scale is ambiguous and unsupported")
        count = int(fields[1]); values += count
        if not 1 <= count <= 512 or index + count > len(rows) or values > _MAX_VALUES:
            raise ValueError("invalid, truncated or oversized Molden contraction")
        primitives = [[_number(v) for v in row.split()] for row in rows[index:index + count]]
        index += count
        if any(len(row) != (3 if label == "sp" else 2) for row in primitives):
            raise ValueError("Molden primitive rows must match their contraction type")
        for column, letter in enumerate(label if label == "sp" else [label], 1):
            angular = "spdfg".index(letter)
            shell = {"center": center, "l": angular, "pure": pure.get(angular, False),
                     "exponents": [row[0] for row in primitives],
                     "coefficients": [row[column] for row in primitives]}
            if dialect == "orca" and angular >= 2 and not shell["pure"]:
                raise ValueError("ORCA dialect: ORCA writes spherical d, f and g only; a Cartesian "
                                 "shell has no documented ORCA convention")
            scales = [_primitive_scale(dialect, angular, shell["pure"], e) for e in shell["exponents"]]
            shell["coefficients"] = [c / s for c, s in zip(shell["coefficients"], scales)]
            norm = _shell_norm(shell)
            if dialect != "renormalize" and abs(norm * norm - 1) > 1e-5:
                hint = (" The title names ORCA, whose contractions need dialect 'orca'."
                        if "orca" in title and dialect != "orca" else "")
                raise ValueError(f"Molden contraction {len(shells) + 1} ({letter}) has squared norm "
                                 f"{norm * norm:.8g} under dialect {dialect!r}; input contractions must "
                                 "be normalized within 1e-5. Name the producer dialect, or use "
                                 "'renormalize' to rescale like the Molden program." + hint)
            # Remove text-rounding drift exactly as a normalized-AO reader does.
            shell["coefficients"] = [c / norm for c in shell["coefficients"]]
            norms.append(norm); shells.append(shell)
    if not shells or centers != set(range(len(atoms))):
        raise ValueError("Molden requires a basis block for every atom")
    order = []; scales = []
    for shell in shells:
        indices, factors = _layout(shell); offset = len(order)
        factors = _component_factors(dialect, shell, factors)
        order.extend(offset + i for i in indices); scales.extend(factors)
    nao = len(order); orbitals = []; current = None; coefficients = {}

    def finish():
        if current is None:
            return
        if current.keys() != {"sym", "ene", "spin", "occup"} or set(coefficients) != set(range(1, nao + 1)):
            raise ValueError("Molden MO requires Sym/Ene/Spin/Occup and each AO coefficient exactly once")
        row = [0.0] * nao
        for source, (target, scale) in enumerate(zip(order, scales), 1):
            row[target] = coefficients[source] / scale
        orbitals.append((current, row))

    for line in sections["mo"][1]:
        match = re.fullmatch(r"([A-Za-z]+)\s*=\s*(.+)", line)
        if match:
            key = match[1].lower(); value = match[2].strip()
            if key not in {"sym", "ene", "spin", "occup"}:
                raise ValueError("unsupported or repeated Molden MO header")
            # A header after coefficients starts the next orbital, so the
            # four keys may come in any order (the format fixes none).
            if current is None or coefficients:
                finish(); current = {}; coefficients = {}
            if key in current:
                raise ValueError("unsupported or repeated Molden MO header")
            current[key] = _number(value) if key in {"ene", "occup"} else value
        else:
            fields = line.split()
            if current is None or len(fields) != 2 or not re.fullmatch(r"[0-9]+", fields[0]):
                raise ValueError("Molden MO coefficient requires an index and one real value")
            ao = int(fields[0])
            if not 1 <= ao <= nao or ao in coefficients:
                raise ValueError("Molden AO index is out of range or repeated")
            coefficients[ao] = _number(fields[1])
        if (len(orbitals) + 1) * nao > _MAX_VALUES:
            raise ValueError("Molden coefficient limit exceeded")
    finish()
    if not orbitals:
        raise ValueError("Molden has no orbitals")
    channels = {"alpha": [], "beta": []}
    limit = 1 if spin == "unrestricted" else 2
    for number, (metadata, row) in enumerate(orbitals, 1):
        channel = metadata["spin"].lower()
        if channel not in channels:
            raise ValueError("Molden spin label or occupation conflicts with the declared spin")
        if not 0 <= metadata["occup"] <= limit:
            raise ValueError(f"Molden orbital {number} ({metadata['spin']}) has occupation "
                             f"{metadata['occup']:g}, outside 0 to {limit} for {spin} spin: QVF "
                             "electron occupations cannot represent negative or excess (for "
                             "example natural-orbital) occupations, and none is clipped")
        channels[channel].append((metadata, row))
    if not channels["alpha"] or (bool(channels["beta"]) != (spin == "unrestricted")):
        raise ValueError("Molden single sets must be Alpha; unrestricted requires both Alpha and Beta")
    return atoms, shells, channels, nao, norms, comments, sections


def import_molden(source, destination, *, spin, dialect="standard"):
    atoms, shells, channels, nao, norms, comments, sections = _parse(_text(source), spin, dialect)
    metadata = {"spin": spin, "orbital_kind": "custom", "n_ao": nao}
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "Molden orbital import") as writer:
            _structure(writer, atoms, [False] * 3, None)
            writer.add_section({"id": "wavefunction", "kind": "wavefunction.gto", "structure_ref": "structure"})
            writer.add_json("wavefunction", "basis", {"structure_ref": "structure", "n_ao": nao, "shells": shells})
            for channel, orbitals in channels.items():
                if not orbitals:
                    continue
                payload = {"n_mo": len(orbitals), "energies": [m["ene"] for m, _ in orbitals],
                           "occupations": [m["occup"] for m, _ in orbitals],
                           "symmetry_labels": [m["sym"] for m, _ in orbitals]}
                if spin == "unrestricted":
                    metadata[channel] = payload
                else:
                    metadata.update(payload)
                role = "mo_coefficients" + ("_" + channel if spin == "unrestricted" else "")
                _binary(writer, "wavefunction", role, (v for _, row in orbitals for v in row), [len(orbitals), nao])
            writer.add_json("wavefunction", "mo_metadata", metadata)
        findings = _validated(temporary)
    return {"format": "molden", "output": str(destination), "spin": spin, "dialect": dialect, "n_ao": nao,
            "findings": findings, "contraction_norms": norms, "omitted_comments": comments,
            "omitted_text_lines": sum(len(sections.get(name, [None, []])[1]) for name in ("title", "molden format")),
            "scope": "Molecular real s, p, sp, d, f and g orbitals only. Orbital kind is custom: Molden does not establish its meaning. No method, charge, multiplicity, electron count, overlap matrix or job is inferred."}


def export_molden(source, destination, *, section=None):
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"wavefunction.gto"})
        basis = reader.json(selected.id, "basis"); metadata = reader.json(selected.id, "mo_metadata")
        sid = selected.to_json().get("structure_ref") or basis["structure_ref"]
        geometry = reader.json(sid, "structure")
        if geometry.get("lattice_vectors") is not None or any(geometry.get("pbc", [False] * 3)):
            raise ValueError("Molden conversion supports molecular geometry only")
        spin = metadata["spin"]
        if (spin not in _SPINS or metadata.get("complex") or metadata.get("coefficient_encoding", "real") != "real"
                or metadata.get("occupation_semantics", "electron_occupation") != "electron_occupation"
                or any(metadata.get("k_point", [0, 0, 0])) or basis.get("ecp")):
            raise ValueError("Molden refuses ECP, complex/spinor, transition-weight and non-Gamma orbitals")
        # Atomic units avoid a reader's different angstrom/bohr constant
        # changing the geometry on which these atomic-unit orbitals are defined.
        lines = ["[Molden Format]", "[Atoms] AU"]
        for i, atom in enumerate(geometry["atoms"], 1):
            checked = _atom(atom["symbol"], atom["position"])
            lines.append(f"{checked['symbol']} {i} {checked['atomic_number']} " + ' '.join(format(_number(v * BOHR_PER_ANGSTROM), '.17g') for v in checked["position"]))
        lines.append("[GTO]"); order = []; factors = []; norms = []; centers = set(); center = None
        types = {2: set(), 3: set(), 4: set()}
        for shell in basis["shells"]:
            shell = dict(shell, pure=shell.get("pure", basis.get("pure", True)))
            indices, scales = _layout(shell); offset = len(order); norm = _shell_norm(shell)
            order.extend(offset + i for i in indices); factors.extend(norm * s for s in scales); norms.append(norm)
            if shell["l"] in types:
                types[shell["l"]].add(shell["pure"])
            if center != shell["center"]:
                center = shell["center"]
                if center in centers:
                    raise ValueError("Molden export requires contiguous basis shells per atom")
                centers.add(center); lines.append(f"{center + 1} 0")
            lines.append(f"{'spdfg'[shell['l']]} {len(shell['exponents'])} 1.00")
            lines.extend(f"{e:.17g} {c / norm:.17g}" for e, c in zip(shell["exponents"], shell["coefficients"]))
        if any(len(t) > 1 for t in types.values()) or not order \
                or centers != set(range(len(geometry["atoms"]))):
            raise ValueError("Molden requires one convention per d, f and g and a basis for "
                             "every atom")
        pure_d, pure_f = True in types[2], True in types[3]
        cart_f = False in types[3]
        flags = ["[5D10F]" if pure_d and cart_f else "[5D]" if pure_d and not cart_f
                 else "[7F]" if pure_f else "[6D]"]
        if types[4]:
            flags.append("[9G]" if True in types[4] else "[15G]")
        lines.extend(flags + ["[MO]"])
        generated_labels = 0; nao = len(order)
        for channel in (["alpha", "beta"] if spin == "unrestricted" else [None]):
            payload = metadata[channel] if channel else metadata
            role = "mo_coefficients" + ("_" + channel if channel else "")
            member = selected.member(role)
            if math.prod(member.shape) > _MAX_VALUES:
                raise ValueError("Molden coefficient limit exceeded")
            matrix = reader.tensor(selected.id, role).tolist()
            if not matrix or any(len(row) != nao for row in matrix):
                raise ValueError("Molden requires nonempty real MO rows matching the basis")
            labels = payload.get("symmetry_labels")
            if labels is None:
                labels = ["A"] * len(matrix); generated_labels += len(matrix)
            for row, energy, occupation, label in zip(matrix, payload["energies"], payload["occupations"], labels):
                if (not label.strip() or any(ord(c) < 32 for c in label)
                        or label.splitlines() != [label]):
                    raise ValueError("Molden symmetry labels must be nonempty single-line text")
                if not 0 <= occupation <= (1 if spin == "unrestricted" else 2):
                    raise ValueError("Molden requires electron occupations in the declared spin range")
                lines.extend([f"Sym= {label}", f"Ene= {_number(energy):.17g}",
                              f"Spin= {(channel or 'alpha').capitalize()}", f"Occup= {_number(occupation):.17g}"])
                lines.extend(f"{i} {_number(row[index] * scale):.17g}" for i, (index, scale) in enumerate(zip(order, factors), 1))
        with _new_output(destination) as temporary:
            temporary.write_text('\n'.join(lines) + '\n', encoding="utf-8")
        return {"format": "molden", "output": str(destination), "section": selected.id,
                "spin": spin, "n_ao": nao, "contraction_norms": norms, "generated_symmetry_labels": generated_labels,
                "omitted_sections": [s.id for s in reader.sections if s.id not in {selected.id, sid}],
                "omitted_members": [role for role in selected.members if role not in
                                    {"basis", "mo_metadata", "mo_coefficients", "mo_coefficients_alpha", "mo_coefficients_beta"}],
                "scope": "Selected geometry, basis, energies, electron occupations, labels and real coefficients only. Orbital kind, method, charge, multiplicity, identities, bonds, jobs, overlap/density matrices and other QVF metadata are not retained. Supply the reported spin on reimport. Keep the original archive."}
