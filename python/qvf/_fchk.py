"""Gaussian formatted checkpoint (fchk) interchange; see docs/converters.md.

The layout follows Gaussian's formchk documentation (https://gaussian.com/formchk/):
a title line, a type/method/basis line, then named fields, each a header
``name (A40), type I/R/C/H/L, value`` or ``N=count`` followed by the values
(integers 6I12 per line, reals 5E16.8 per line). The AO conventions are the
ones documented for fchk by IOData 1.0.1 (iodata/formats/fchk.py, CONVENTIONS,
with L2-normalized Cartesian primitives) and by Schlegel and Frisch,
Int. J. Quantum Chem. 54, 83 (1995) for the real solid harmonics. No Gaussian,
IOData or cclib code is incorporated here.

Gaussian's AO order within a shell:

* p: x, y, z.
* Cartesian d: xx, yy, zz, xy, xz, yz. Cartesian f: xxx, yyy, zzz, xyy,
  xxy, xxz, xzz, yzz, yyz, xyz. Cartesian g and higher: the reverse of the
  QVF (libint) order, zzzz, yzzz, yyzz, ... , xxxx.
* Pure shells: m = 0, +1, -1, +2, -2, ... (c0, c1, s1, c2, s2, ...), with the
  QVF real solid harmonics (no Condon-Shortley phase, no extra signs).

Every Gaussian AO is unit-normalized: a Cartesian component is normalized on
its own, so d_xy is sqrt(3) times the QVF d_xy of Appendix A.1, whose
primitives carry the axial normalization. Contraction coefficients refer to
normalized primitives.
"""
from __future__ import annotations

import math
import re

from . import __version__
from .convert import _MAX_VALUES, _atom, _binary, _new_output, _number, _select, _structure, _text, _validated
from ._legacy_writer import BOHR_PER_ANGSTROM, _ELEMENT_SYMBOLS
from .adapters.pyscf import cartesian_powers, real_solid_harmonics
from .reader import Reader
from .units import HARTREE_IN_EV
from .writer import Writer

_SPINS = {"restricted", "restricted_open", "unrestricted"}
# A header line: name, type and either "N=  count" or a scalar value. Names
# start in column 1; numeric data lines are right-aligned and start blank.
_HEADER = re.compile(r"(?P<name>\S(?:.*?\S)?)\s+(?P<type>[IRCHL])\s+"
                     r"(?:N=\s*(?P<count>[0-9]+)|(?P<value>\S+))\s*")
_CARTESIAN = {
    2: ["xx", "yy", "zz", "xy", "xz", "yz"],
    3: ["xxx", "yyy", "zzz", "xyy", "xxy", "xxz", "xzz", "yzz", "yyz", "xyz"],
}
_USED = {"Number of atoms", "Charge", "Multiplicity", "Number of electrons",
         "Number of alpha electrons", "Number of beta electrons", "Number of basis functions",
         "Number of independent functions", "Atomic numbers", "Nuclear charges",
         "Current cartesian coordinates", "Number of contracted shells",
         "Number of primitive shells", "Pure/Cartesian d shells", "Pure/Cartesian f shells",
         "Highest angular momentum", "Largest degree of contraction", "Shell types",
         "Number of primitives per shell", "Shell to atom map", "Primitive exponents",
         "Contraction coefficients", "P(S=P) Contraction coefficients",
         "Coordinates of each shell", "Total Energy", "Alpha Orbital Energies",
         "Alpha MO coefficients", "Beta Orbital Energies", "Beta MO coefficients",
         "Total SCF Density", "Spin SCF Density", "Mulliken Charges", "SCF Energy"}


def _double_factorial(n):
    return math.prod(range(n, 0, -2)) if n > 0 else 1


def _n_ao(l, pure):
    return 2 * l + 1 if pure else (l + 1) * (l + 2) // 2


def _layout(l, pure):
    """For each fchk AO of a shell in file order: its QVF index in the shell
    and the factor QVF AO / fchk AO (both on the same unit contraction)."""
    if l == 0:
        return [0], [1.0]
    if l == 1:
        # Gaussian p is x, y, z; QVF pure p is y, z, x, Cartesian p x, y, z.
        return ([2, 0, 1] if pure else [0, 1, 2]), [1.0] * 3
    if pure:
        ms = [0] + [s * k for k in range(1, l + 1) for s in (1, -1)]
        return [m + l for m in ms], [1.0] * len(ms)
    qvf = cartesian_powers(l)
    if l in _CARTESIAN:
        powers = [tuple(name.count(a) for a in "xyz") for name in _CARTESIAN[l]]
    else:
        powers = qvf[::-1]
    factors = [math.sqrt(math.prod(_double_factorial(2 * n - 1) for n in p) / _double_factorial(2 * l - 1))
               for p in powers]
    return [qvf.index(p) for p in powers], factors


def _contraction_norm(l, exponents, coefficients):
    """sqrt(S_axial) of Appendix A.1 (ratio form, no overflow)."""
    total = math.fsum(a * b * (2 * math.sqrt(min(x, y) / max(x, y)) / (1 + min(x, y) / max(x, y))) ** (l + 1.5)
                      for x, a in zip(exponents, coefficients) for y, b in zip(exponents, coefficients))
    if not math.isfinite(total) or total <= 1e-14:
        raise ValueError("fchk requires finite, nonzero contraction norms")
    return math.sqrt(total)


def _real(token):
    # Fortran E16.8 drops the E of a three-digit exponent: 1.23456789-100.
    match = re.fullmatch(r"([+-]?[0-9.]+)([+-][0-9]{3})", token)
    return _number(match[1] + "E" + match[2] if match else token)


def _integer(token):
    if not re.fullmatch(r"[+-]?[0-9]+", token):
        raise ValueError(f"fchk integer field holds {token!r}")
    return int(token)


def _fields(text):
    lines = text.splitlines()
    if len(lines) < 3:
        raise ValueError("fchk needs a title line, a type/method/basis line and fields")
    second = lines[1]
    header = {"title": lines[0].strip(), "job_type": second[:10].strip(),
              "method": second[10:40].strip(), "basis": second[40:].strip()}
    fields = {}; skipped = []; index = 2; total = 0
    while index < len(lines):
        line = lines[index].rstrip(); index += 1
        if not line:
            continue
        match = _HEADER.fullmatch(line)
        if not match or line[0].isspace():
            raise ValueError(f"fchk line {index}: expected a field header, found {line[:60]!r}")
        name, kind = match["name"], match["type"]
        if name in fields or name in skipped:
            raise ValueError(f"fchk field {name!r} is repeated")
        if kind not in "IR":
            # Text (C, H) and logical (L) fields: nothing this import reads.
            skipped.append(name)
            if match["count"] is not None:
                while index < len(lines) and not (_HEADER.fullmatch(lines[index].rstrip())
                                                  and not lines[index][:1].isspace()):
                    index += 1
            continue
        convert = _integer if kind == "I" else _real
        if match["count"] is None:
            fields[name] = convert(match["value"])
            continue
        count = int(match["count"]); total += count
        if count > _MAX_VALUES or total > 4 * _MAX_VALUES:
            raise ValueError("fchk value limit exceeded")
        values = []
        while len(values) < count:
            if index >= len(lines):
                raise ValueError(f"fchk field {name!r} is truncated")
            values.extend(convert(t) for t in lines[index].split()); index += 1
        if len(values) != count:
            raise ValueError(f"fchk field {name!r} has {len(values)} values, not {count}")
        fields[name] = values
    return header, fields, skipped


def _get(fields, name, kind=None, length=None):
    if name not in fields:
        raise ValueError(f"fchk import requires the field {name!r}")
    value = fields[name]
    if (kind == "array") != isinstance(value, list) and kind is not None:
        raise ValueError(f"fchk field {name!r} must be {'an array' if kind == 'array' else 'a scalar'}")
    if length is not None and len(value) != length:
        raise ValueError(f"fchk field {name!r} has {len(value)} values, expected {length}")
    return value


def _unpack(packed, n):
    """Gaussian's lower-triangle row-major packing: (0,0), (1,0), (1,1), (2,0), ..."""
    matrix = [[0.0] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1):
            matrix[i][j] = matrix[j][i] = packed[k]; k += 1
    return matrix


def _overlap_1d(la, lb, pa, pb, p):
    """Obara-Saika table of int (x-A)^i (x-B)^j exp(-p (x-P)^2) dx / sqrt(pi/p)."""
    table = [[0.0] * (lb + 1) for _ in range(la + 1)]; table[0][0] = 1.0; h = 0.5 / p
    for i in range(la + 1):
        for j in range(lb + 1):
            if i == j == 0:
                continue
            if i:
                value = pa * table[i - 1][j] + h * ((i - 1) * (table[i - 2][j] if i > 1 else 0)
                                                   + j * (table[i - 1][j - 1] if j else 0))
            else:
                value = pb * table[i][j - 1] + h * ((j - 1) * (table[i][j - 2] if j > 1 else 0))
            table[i][j] = value
    return table


def _overlap(shells, centres):
    """AO overlap of a QVF basis (Appendix A.1 primitives, A.2 order)."""
    from ._legacy_writer import primitive_norm
    info = []; offset = 0
    for shell in shells:
        l, pure = shell["l"], shell["pure"]
        prims = [(a, c * primitive_norm(a, l)) for a, c in zip(shell["exponents"], shell["coefficients"])]
        info.append((l, pure, centres[shell["center"]], prims, cartesian_powers(l),
                     real_solid_harmonics(l) if pure and l > 1 else None, offset))
        offset += _n_ao(l, pure)
    size = offset; matrix = [[0.0] * size for _ in range(size)]
    for ia, (la, pure_a, ra, pa, powers_a, ta, oa) in enumerate(info):
        for ib in range(ia + 1):
            lb, pure_b, rb, pb, powers_b, tb, ob = info[ib]
            ab2 = sum((x - y) ** 2 for x, y in zip(ra, rb))
            block = [[0.0] * len(powers_b) for _ in powers_a]
            for a, ca in pa:
                for b, cb in pb:
                    p = a + b; scale = ca * cb * math.exp(-a * b / p * ab2) * (math.pi / p) ** 1.5
                    if scale == 0.0:
                        continue
                    centre = [(a * x + b * y) / p for x, y in zip(ra, rb)]
                    tables = [_overlap_1d(la, lb, centre[k] - ra[k], centre[k] - rb[k], p) for k in range(3)]
                    for i, (ax, ay, az) in enumerate(powers_a):
                        row = block[i]
                        for j, (bx, by, bz) in enumerate(powers_b):
                            row[j] += scale * tables[0][ax][bx] * tables[1][ay][by] * tables[2][az][bz]
            if ta is not None:
                block = [[sum(ta[k][m] * block[k][j] for k in range(len(powers_a))) for j in range(len(powers_b))]
                         for m in range(2 * la + 1)]
            elif pure_a and la == 1:
                block = [block[1], block[2], block[0]]
            if tb is not None:
                block = [[sum(row[k] * tb[k][m] for k in range(len(powers_b))) for m in range(2 * lb + 1)]
                         for row in block]
            elif pure_b and lb == 1:
                block = [[row[1], row[2], row[0]] for row in block]
            for i, row in enumerate(block):
                for j, value in enumerate(row):
                    matrix[oa + i][ob + j] = matrix[ob + j][oa + i] = value
    return matrix


def _occupations(n_mo, n_alpha, n_beta, spin):
    if spin == "unrestricted":
        return ([1.0] * n_alpha + [0.0] * (n_mo - n_alpha), [1.0] * n_beta + [0.0] * (n_mo - n_beta))
    return [2.0] * n_beta + [1.0] * (n_alpha - n_beta) + [0.0] * (n_mo - n_alpha), None


def _parse(text):
    header, fields, skipped = _fields(text)
    numbers = _get(fields, "Atomic numbers", "array")
    natoms = len(numbers)
    if not natoms or fields.get("Number of atoms", natoms) != natoms:
        raise ValueError("fchk atom count disagrees with its Atomic numbers")
    if any(not 1 <= z < len(_ELEMENT_SYMBOLS) for z in numbers):
        raise ValueError("fchk atomic numbers must name elements; ghost and dummy atoms are not supported")
    if any(name.startswith("ECP-") for name in fields) or \
            [float(z) for z in numbers] != [float(q) for q in fields.get("Nuclear charges", numbers)]:
        raise ValueError("fchk has effective core potentials (nuclear charges below the atomic numbers "
                         "or ECP- fields); QVF basis.ecp cannot be built from them, and the orbitals "
                         "would be wrong without it")
    coordinates = _get(fields, "Current cartesian coordinates", "array", 3 * natoms)
    centres = [coordinates[3 * i:3 * i + 3] for i in range(natoms)]
    atoms = [_atom(_ELEMENT_SYMBOLS[z], [v / BOHR_PER_ANGSTROM for v in centre])
             for z, centre in zip(numbers, centres)]
    types = _get(fields, "Shell types", "array")
    counts = _get(fields, "Number of primitives per shell", "array", len(types))
    owners = _get(fields, "Shell to atom map", "array", len(types))
    exponents = _get(fields, "Primitive exponents", "array")
    first = _get(fields, "Contraction coefficients", "array", len(exponents))
    second = fields.get("P(S=P) Contraction coefficients")
    if not types or sum(counts) != len(exponents) or any(n < 1 for n in counts):
        raise ValueError("fchk shell primitive counts disagree with the exponents")
    if second is not None and (not isinstance(second, list) or len(second) != len(exponents)):
        raise ValueError("fchk P(S=P) contraction coefficients must match the exponents")
    if -1 in types and second is None:
        raise ValueError("fchk SP shells need P(S=P) Contraction coefficients")
    if any(not 1 <= a <= natoms for a in owners):
        raise ValueError("fchk shell to atom map names no atom")
    if any(abs(t) > 7 for t in types):
        raise ValueError("QVF shells have l <= 7")
    # Gaussian's declared conventions must agree with the shell-type signs:
    # 0 pure, 1 Cartesian, for d and for f and higher.
    for flag, ls in (("Pure/Cartesian d shells", {2}), ("Pure/Cartesian f shells", set(range(3, 8)))):
        if flag in fields and fields[flag] in (0, 1):
            if any(abs(t) in ls and (t < 0) != (fields[flag] == 0) for t in types):
                raise ValueError(f"fchk shell types disagree with {flag!r}")
    shell_coordinates = fields.get("Coordinates of each shell")
    if shell_coordinates is not None:
        if not isinstance(shell_coordinates, list) or len(shell_coordinates) != 3 * len(types):
            raise ValueError("fchk Coordinates of each shell must have three values per shell")
        for k, owner in enumerate(owners):
            if any(abs(a - b) > 1e-6 for a, b in zip(shell_coordinates[3 * k:3 * k + 3], centres[owner - 1])):
                raise ValueError(f"fchk shell {k + 1} is not centred on atom {owner}")
    has_pure = any(t < -1 for t in types)
    shells = []; norms = []; order = []; scales = []; start = 0
    for code, count, owner in zip(types, counts, owners):
        es = exponents[start:start + count]
        if any(e <= 0 for e in es):
            raise ValueError("fchk primitive exponents must be positive")
        # s and p shells follow the file's d/f convention; they are the same functions.
        parts = [(0, has_pure, first), (1, has_pure, second)] if code == -1 else \
                [(abs(code), code < -1 or (code in (0, 1) and has_pure), first)]
        for l, pure, source in parts:
            cs = source[start:start + count]
            norm = _contraction_norm(l, es, cs)
            if abs(norm * norm - 1) > 1e-5:
                raise ValueError(f"fchk contraction {len(shells) + 1} (l = {l}) has squared norm "
                                 f"{norm * norm:.8g}; Gaussian writes normalized contractions, and an "
                                 "unnormalized one is ambiguous")
            indices, factors = _layout(l, pure); offset = sum(_n_ao(s["l"], s["pure"]) for s in shells)
            order.extend(offset + i for i in indices); scales.extend(factors)
            # Remove text rounding exactly as a normalized-AO reader does.
            shells.append({"center": owner - 1, "l": l, "pure": pure, "exponents": list(es),
                           "coefficients": [c / norm for c in cs]})
            norms.append(norm)
        start += count
    nao = len(order)
    if _get(fields, "Number of basis functions") != nao:
        raise ValueError(f"fchk Number of basis functions disagrees with the {nao} AOs of its shells")
    n_electrons = _get(fields, "Number of electrons"); n_alpha = _get(fields, "Number of alpha electrons")
    n_beta = _get(fields, "Number of beta electrons")
    charge = _get(fields, "Charge"); multiplicity = _get(fields, "Multiplicity")
    if (n_alpha + n_beta != n_electrons or n_beta < 0 or n_alpha < n_beta
            or sum(numbers) - n_electrons != charge or n_alpha - n_beta + 1 != multiplicity):
        raise ValueError("fchk charge, multiplicity and electron counts disagree (alpha >= beta expected)")
    energies = _get(fields, "Alpha Orbital Energies", "array"); n_mo = len(energies)
    if fields.get("Number of independent functions", n_mo) != n_mo or not 0 < n_mo <= nao or n_alpha > n_mo:
        raise ValueError("fchk orbital count disagrees with Number of independent functions or the electrons")
    unrestricted = "Beta MO coefficients" in fields
    spin = "unrestricted" if unrestricted else ("restricted" if n_alpha == n_beta else "restricted_open")
    channels = {}
    for channel in (["alpha", "beta"] if unrestricted else ["alpha"]):
        prefix = channel.capitalize()
        values = _get(fields, prefix + " MO coefficients", "array", n_mo * nao)
        channel_energies = _get(fields, prefix + " Orbital Energies", "array", n_mo)
        rows = []
        for r in range(n_mo):
            source = values[r * nao:(r + 1) * nao]; row = [0.0] * nao
            for value, target, scale in zip(source, order, scales):
                row[target] = value / scale
            rows.append(row)
        channels[channel] = (channel_energies, rows)
    if not unrestricted and ("Beta Orbital Energies" in fields):
        raise ValueError("fchk has Beta Orbital Energies without Beta MO coefficients")
    alpha_occ, beta_occ = _occupations(n_mo, n_alpha, n_beta, spin)
    occupations = {"alpha": alpha_occ, "beta": beta_occ}
    # Densities in the QVF basis: P_qvf[I,J] = P_file[k,l] / (f_k f_l).
    densities = {}; packed_size = nao * (nao + 1) // 2

    def transform(packed, name):
        if not isinstance(packed, list) or len(packed) != packed_size:
            raise ValueError(f"fchk {name} must hold the {packed_size} values of a packed lower triangle")
        source = _unpack(packed, nao); matrix = [[0.0] * nao for _ in range(nao)]
        for k, (i, fi) in enumerate(zip(order, scales)):
            row = source[k]; target = matrix[i]
            for l, (j, fj) in enumerate(zip(order, scales)):
                target[j] = row[l] / (fi * fj)
        return matrix

    omitted = sorted(set(fields) - _USED)
    if "Total SCF Density" in fields:
        total = transform(fields["Total SCF Density"], "Total SCF Density")
        if spin == "unrestricted":
            if "Spin SCF Density" in fields:
                spin_density = transform(fields["Spin SCF Density"], "Spin SCF Density")
                densities["density_matrix_alpha"] = [[(a + b) / 2 for a, b in zip(r, s)] for r, s in zip(total, spin_density)]
                densities["density_matrix_beta"] = [[(a - b) / 2 for a, b in zip(r, s)] for r, s in zip(total, spin_density)]
            else:
                omitted.append("Total SCF Density")
        else:
            densities["density_matrix"] = total
            if "Spin SCF Density" in fields:
                omitted.append("Spin SCF Density")
    elif "Spin SCF Density" in fields:
        omitted.append("Spin SCF Density")
    deviation = None
    if densities:
        # The derived occupations must reproduce the file's SCF density:
        # otherwise the occupied orbitals are not the lowest n (non-aufbau).
        deviation = 0.0
        for role, matrix in densities.items():
            channel = "beta" if role.endswith("beta") else "alpha"
            occ = occupations[channel]; rows = channels[channel][1]
            occupied = [(n, rows[i]) for i, n in enumerate(occ) if n]
            for i in range(nao):
                for j in range(i + 1):
                    value = sum(n * row[i] * row[j] for n, row in occupied)
                    deviation = max(deviation, abs(value - matrix[i][j]))
        largest = max(abs(v) for matrix in densities.values() for row in matrix for v in row)
        if deviation > 1e-5 * max(1.0, largest):
            raise ValueError(f"the electron counts do not reproduce the fchk Total SCF Density (maximum "
                             f"deviation {deviation:.3g}): the occupied orbitals are not the lowest ones, "
                             "and occupations cannot be derived")
    charges = fields.get("Mulliken Charges")
    if charges is not None and (not isinstance(charges, list) or len(charges) != natoms):
        raise ValueError("fchk Mulliken Charges must have one value per atom")
    return dict(header=header, atoms=atoms, centres=centres, shells=shells, norms=norms, nao=nao,
                spin=spin, channels=channels, occupations=occupations, densities=densities,
                deviation=deviation, charge=charge, multiplicity=multiplicity, n_electrons=n_electrons,
                energy=fields.get("Total Energy"), mulliken=charges, omitted=omitted, skipped=skipped,
                scf_energy=fields.get("SCF Energy"))


def import_fchk(source, destination, *, overlap=True):
    data = _parse(_text(source)); nao = data["nao"]; spin = data["spin"]
    header = data["header"]
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "Gaussian fchk import") as writer:
            _structure(writer, data["atoms"], [False] * 3, None)
            writer.add_section({"id": "wavefunction", "kind": "wavefunction.gto", "structure_ref": "structure"})
            writer.add_json("wavefunction", "basis", {"structure_ref": "structure", "n_ao": nao,
                                                     "n_electrons": data["n_electrons"],
                                                     "shells": data["shells"]})
            metadata = {"spin": spin, "orbital_kind": "custom", "n_ao": nao}
            for channel, (energies, rows) in data["channels"].items():
                payload = {"n_mo": len(rows), "energies": energies,
                           "occupations": data["occupations"][channel]}
                if spin == "unrestricted":
                    metadata[channel] = payload
                else:
                    metadata.update(payload)
                role = "mo_coefficients" + ("_" + channel if spin == "unrestricted" else "")
                _binary(writer, "wavefunction", role, (v for row in rows for v in row), [len(rows), nao])
            writer.add_json("wavefunction", "mo_metadata", metadata)
            if overlap:
                matrix = _overlap(data["shells"], data["centres"])
                _binary(writer, "wavefunction", "overlap", (v for row in matrix for v in row), [nao, nao])
            for role, matrix in data["densities"].items():
                _binary(writer, "wavefunction", role, (v for row in matrix for v in row), [nao, nao])
            if data["mulliken"] is not None:
                writer.add_section({"id": "atom_properties", "kind": "atom_properties",
                                    "structure_ref": "structure"})
                _binary(writer, "atom_properties", "mulliken_charge", data["mulliken"], [len(data["atoms"])])
            provenance = {"charge": data["charge"], "multiplicity": data["multiplicity"],
                          "n_electrons": data["n_electrons"]}
            if header["basis"]:
                provenance["basis"] = header["basis"]
            writer.set_provenance(**provenance)
            if data["energy"] is not None:
                writer.set_energies(data["energy"], method=header["method"].lower() or None)
        findings = _validated(temporary)
    omitted = list(data["omitted"])
    # An SCF energy that differs from Total Energy (a post-SCF run) has no
    # method name in the file to store it under.
    if data["scf_energy"] is not None and (data["energy"] is None or abs(data["scf_energy"] - data["energy"]) > 1e-8):
        omitted.append("SCF Energy")
    return {"format": "fchk", "output": str(destination), "spin": spin, "n_ao": nao,
            "findings": findings, "contraction_norms": data["norms"],
            "density_members": sorted(data["densities"]), "density_deviation": data["deviation"],
            "overlap": bool(overlap), "omitted_fields": sorted(set(omitted)),
            "omitted_text_fields": data["skipped"],
            "omitted_header": {k: v for k, v in header.items() if v and k != "basis"},
            "scope": "Molecular real orbitals, s to l = 7. Occupations are derived from the electron "
                     "counts in the file's orbital order (aufbau) and are checked against the SCF "
                     "density when the file has one. Orbital kind is custom: the file does not say "
                     "whether the orbitals are canonical. The overlap is computed from the basis; "
                     "ECPs, post-SCF densities, gradients, Hessians and other fields are not imported."}


def _aufbau(occupations, values, name):
    """Count of each occupation value in a leading aufbau block (values descending)."""
    counts = []; index = 0
    for value in values:
        start = index
        while index < len(occupations) and abs(occupations[index] - value) <= 1e-10:
            index += 1
        counts.append(index - start)
    if any(abs(n) > 1e-10 for n in occupations[index:]):
        raise ValueError(f"fchk has no occupations: the {name} occupations must be integral and in "
                         f"aufbau order ({', '.join(f'{v:g}' for v in values)}, then 0) in the row order")
    return counts


def _real_array(name, values):
    values = list(values)
    lines = [f"{name:<40}   R   N={len(values):12d}"]
    for k in range(0, len(values), 5):
        lines.append("".join(f"{_number(v):16.8E}" for v in values[k:k + 5]))
    return lines


def _integer_array(name, values):
    values = list(values)
    lines = [f"{name:<40}   I   N={len(values):12d}"]
    for k in range(0, len(values), 6):
        lines.append("".join(f"{int(v):12d}" for v in values[k:k + 6]))
    return lines


def _scalar(name, value):
    if isinstance(value, int):
        return [f"{name:<40}   I     {value:12d}"]
    return [f"{name:<40}   R     {_number(value):22.15E}"]


def _hartree(entry):
    if not isinstance(entry, dict) or not isinstance(entry.get("value"), (int, float)):
        return None
    units = str(entry.get("units", "hartree")).lower()
    factor = {"hartree": 1.0, "eh": 1.0, "ha": 1.0, "ev": 1 / HARTREE_IN_EV}.get(units)
    return None if factor is None else entry["value"] * factor


def export_fchk(source, destination, *, section=None):
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"wavefunction.gto"})
        basis = reader.json(selected.id, "basis"); metadata = reader.json(selected.id, "mo_metadata")
        sid = selected.to_json().get("structure_ref") or basis["structure_ref"]
        geometry = reader.json(sid, "structure")
        if geometry.get("lattice_vectors") is not None or any(geometry.get("pbc", [False] * 3)):
            raise ValueError("fchk export supports molecular geometry only")
        spin = metadata["spin"]
        if (spin not in _SPINS or metadata.get("complex") or metadata.get("coefficient_encoding", "real") != "real"
                or metadata.get("occupation_semantics", "electron_occupation") != "electron_occupation"
                or any(metadata.get("k_point", [0, 0, 0])) or basis.get("ecp")):
            raise ValueError("fchk export refuses ECP, complex/spinor, transition-weight and non-Gamma orbitals")
        atoms = [_atom(a["symbol"], a["position"]) for a in geometry["atoms"]]
        numbers = [a["atomic_number"] for a in atoms]
        centres = [[v * BOHR_PER_ANGSTROM for v in a["position"]] for a in atoms]
        default_pure = basis.get("pure", True)
        shells = [dict(s, pure=s.get("pure", default_pure)) for s in basis["shells"]]
        for name, ls in (("d", {2}), ("f and higher", set(range(3, 8)))):
            if len({s["pure"] for s in shells if s["l"] in ls}) > 1:
                raise ValueError(f"fchk has one pure/Cartesian convention for {name} shells; the basis mixes them")
        # Gaussian groups shells by atom; a stable sort keeps the in-atom order.
        offsets = []; offset = 0
        for shell in shells:
            offsets.append(offset); offset += _n_ao(shell["l"], shell["pure"])
        nao = offset
        ranked = sorted(range(len(shells)), key=lambda k: shells[k]["center"])
        order = []; factors = []; norms = []
        for k in ranked:
            shell = shells[k]
            norm = _contraction_norm(shell["l"], shell["exponents"], shell["coefficients"]); norms.append(norm)
            indices, scales = _layout(shell["l"], shell["pure"])
            order.extend(offsets[k] + i for i in indices); factors.extend(norm * s for s in scales)
        channels = ["alpha", "beta"] if spin == "unrestricted" else [None]
        matrices = {}; payloads = {}
        for channel in channels:
            payload = metadata[channel] if channel else metadata
            role = "mo_coefficients" + ("_" + channel if channel else "")
            member = selected.member(role)
            if math.prod(member.shape) > _MAX_VALUES:
                raise ValueError("fchk coefficient limit exceeded")
            matrix = reader.tensor(selected.id, role).tolist()
            if not matrix or any(len(row) != nao for row in matrix):
                raise ValueError("fchk export requires nonempty real MO rows matching the basis")
            matrices[channel] = matrix; payloads[channel] = payload
        n_mo = len(matrices[channels[0]])
        if any(len(m) != n_mo for m in matrices.values()):
            raise ValueError("fchk stores the same number of alpha and beta orbitals; the channels differ")
        if spin == "unrestricted":
            n_alpha = _aufbau(payloads["alpha"]["occupations"], [1], "alpha")[0]
            n_beta = _aufbau(payloads["beta"]["occupations"], [1], "beta")[0]
        else:
            doubly, singly = _aufbau(metadata["occupations"], [2, 1], spin)
            if spin == "restricted" and singly:
                raise ValueError("restricted orbitals with singly occupied rows: declare restricted_open")
            n_alpha, n_beta = doubly + singly, doubly
        if n_alpha < n_beta:
            raise ValueError("fchk requires at least as many alpha as beta electrons")
        n_electrons = n_alpha + n_beta; charge = sum(numbers) - n_electrons
        multiplicity = n_alpha - n_beta + 1
        provenance = reader.raw_manifest.get("provenance") or {}
        stated = {"n_electrons": n_electrons, "charge": charge, "multiplicity": multiplicity}
        for key, value in stated.items():
            if key in provenance and provenance[key] != value:
                raise ValueError(f"provenance {key} {provenance[key]} disagrees with the occupations ({value})")
        if basis.get("n_electrons", n_electrons) != n_electrons:
            raise ValueError("basis n_electrons disagrees with the occupations")
        energies_block = reader.raw_manifest.get("energies") or {}
        total = _hartree(energies_block.get("total"))
        scf = _hartree(provenance.get("scf_energy"))
        method = energies_block.get("method") or provenance.get("method") or "unknown"
        basis_label = provenance.get("basis") or "unknown"
        for text in (method, basis_label):
            if any(ord(c) < 32 for c in str(text)) or len(str(text)) > 30:
                raise ValueError("fchk method and basis names must be single-line text of at most 30 characters")
        title = (reader.raw_manifest.get("source") or {}).get("calculation") or "QVF export"
        title = "".join(c if 32 <= ord(c) < 127 else " " for c in str(title))[:72]
        lines = [title or "QVF export", f"{'SP':<10}{str(method).upper():<30}{basis_label:<30}".rstrip()]
        lines += _scalar("Number of atoms", len(atoms)) + _scalar("Charge", charge)
        lines += _scalar("Multiplicity", multiplicity) + _scalar("Number of electrons", n_electrons)
        lines += _scalar("Number of alpha electrons", n_alpha) + _scalar("Number of beta electrons", n_beta)
        lines += _scalar("Number of basis functions", nao) + _scalar("Number of independent functions", n_mo)
        lines += _integer_array("Atomic numbers", numbers)
        lines += _real_array("Nuclear charges", [float(z) for z in numbers])
        lines += _real_array("Current cartesian coordinates", [v for c in centres for v in c])
        ordered = [shells[k] for k in ranked]
        types = [(-1 if s["pure"] and s["l"] > 1 else 1) * s["l"] for s in ordered]
        d_flag = next((0 if s["pure"] else 1 for s in ordered if s["l"] == 2), 0)
        f_flag = next((0 if s["pure"] else 1 for s in ordered if s["l"] >= 3), 0)
        lines += _scalar("Number of contracted shells", len(ordered))
        lines += _scalar("Number of primitive shells", sum(len(s["exponents"]) for s in ordered))
        lines += _scalar("Pure/Cartesian d shells", d_flag) + _scalar("Pure/Cartesian f shells", f_flag)
        lines += _scalar("Highest angular momentum", max(s["l"] for s in ordered))
        lines += _scalar("Largest degree of contraction", max(len(s["exponents"]) for s in ordered))
        lines += _integer_array("Shell types", types)
        lines += _integer_array("Number of primitives per shell", [len(s["exponents"]) for s in ordered])
        lines += _integer_array("Shell to atom map", [s["center"] + 1 for s in ordered])
        lines += _real_array("Primitive exponents", [e for s in ordered for e in s["exponents"]])
        lines += _real_array("Contraction coefficients",
                             [c / n for s, n in zip(ordered, norms) for c in s["coefficients"]])
        lines += _real_array("Coordinates of each shell", [v for s in ordered for v in centres[s["center"]]])
        if total is not None:
            lines += _scalar("Total Energy", float(total))
        if scf is not None:
            lines += _scalar("SCF Energy", float(scf))

        def file_rows(matrix):
            return [[row[i] * f for i, f in zip(order, factors)] for row in matrix]

        rows = {channel: file_rows(matrix) for channel, matrix in matrices.items()}
        for channel in channels:
            prefix = (channel or "alpha").capitalize()
            lines += _real_array(prefix + " Orbital Energies", payloads[channel]["energies"])
            lines += _real_array(prefix + " MO coefficients", [v for row in rows[channel] for v in row])

        def packed(matrix):
            return [matrix[i][j] * factors[a] * factors[b]
                    for a, i in enumerate(order) for b, j in enumerate(order[:a + 1])]

        def from_orbitals(pairs):
            matrix = [[0.0] * nao for _ in range(nao)]
            for n, row in pairs:
                for a in range(nao):
                    if row[a]:
                        target = matrix[a]; scaled = n * row[a]
                        for b in range(a + 1):
                            target[b] += scaled * row[b]
            return [[matrix[max(a, b)][min(a, b)] for b in range(nao)] for a in range(nao)]

        members = set(selected.members)
        if spin == "unrestricted":
            if {"density_matrix_alpha", "density_matrix_beta"} <= members:
                pa = reader.tensor(selected.id, "density_matrix_alpha").tolist()
                pb = reader.tensor(selected.id, "density_matrix_beta").tolist()
                density_source = "archive"
            else:
                pa = from_orbitals([(n, r) for n, r in zip(payloads["alpha"]["occupations"], matrices["alpha"]) if n])
                pb = from_orbitals([(n, r) for n, r in zip(payloads["beta"]["occupations"], matrices["beta"]) if n])
                density_source = "orbitals"
            total_density = [[a + b for a, b in zip(x, y)] for x, y in zip(pa, pb)]
            spin_density = [[a - b for a, b in zip(x, y)] for x, y in zip(pa, pb)]
        else:
            pairs = list(zip(metadata["occupations"], matrices[None]))
            if "density_matrix" in members:
                total_density = reader.tensor(selected.id, "density_matrix").tolist(); density_source = "archive"
            else:
                total_density = from_orbitals([(n, r) for n, r in pairs if n]); density_source = "orbitals"
            spin_density = from_orbitals([(1.0, r) for n, r in pairs if abs(n - 1) <= 1e-10]) \
                if spin == "restricted_open" else None
        lines += _real_array("Total SCF Density", packed(total_density))
        if spin_density is not None:
            lines += _real_array("Spin SCF Density", packed(spin_density))
        charges = None; used = {selected.id, sid}
        candidates = [s for s in reader.sections if s.kind == "atom_properties" and "mulliken_charge" in s.members
                      and (s.to_json().get("structure_ref") or sid) == sid]
        if len(candidates) == 1:
            charges = reader.tensor(candidates[0].id, "mulliken_charge").tolist(); used.add(candidates[0].id)
            lines += _real_array("Mulliken Charges", charges)
        with _new_output(destination) as temporary:
            temporary.write_text("\n".join(lines) + "\n", encoding="ascii")
        return {"format": "fchk", "output": str(destination), "section": selected.id, "spin": spin,
                "n_ao": nao, "n_mo": n_mo, "charge": charge, "multiplicity": multiplicity,
                "contraction_norms": norms, "density_source": density_source,
                "mulliken_charges": charges is not None, "total_energy": total is not None,
                "omitted_sections": [s.id for s in reader.sections if s.id not in used],
                "omitted_members": [role for role in selected.members if role not in
                                    {"basis", "mo_metadata", "mo_coefficients", "mo_coefficients_alpha",
                                     "mo_coefficients_beta", "density_matrix", "density_matrix_alpha",
                                     "density_matrix_beta"}],
                "scope": "Geometry, basis, orbital energies, real coefficients, SCF densities, the total "
                         "and SCF energies and Mulliken charges only, with 9 significant digits in "
                         "array fields (Gaussian E16.8). Charge and multiplicity are derived from the "
                         "occupations, which fchk does not store. Symmetry labels, orbital kind, "
                         "identities, bonds, jobs, the overlap and other QVF data are not retained. "
                         "Keep the original archive."}
