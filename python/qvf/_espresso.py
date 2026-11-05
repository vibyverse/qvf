"""Quantum ESPRESSO pw.x structure interchange; see docs/converters.md.

Import reads the cell and atoms of a pw.x input (any ibrav of INPUT_PW, or a
space_group input with crystal_sg sites, expanded to the conventional cell)
or the final structure, total energy, forces and stress of a pw.x XML data
file (data-file-schema.xml). Export writes a structure periodic on all three
axes as a pw.x input fragment with ``ibrav = 0``. pp.x cube files import
through the cube reader with the periodic, endpoint-free FFT-grid
convention, and pp.x filplot files with the reader below.

The filplot file is the intermediate file of pp.x (``filplot`` in INPUT_PP,
PP/Doc/INPUT_PP.def), written and read by ``plot_io`` in Modules/plot_io.f90
with these Fortran records:

    title                                         (a)
    nr1x nr2x nr3x nr1 nr2 nr3 nat ntyp           (8i8)
    ibrav celldm(1:6)                             (i6,6f12.8)
    at(:,1); at(:,2); at(:,3)                     list-directed, ibrav = 0 only
    gcutm dual ecutwfc plot_num                   (3f20.10,i6)
    nt atm(nt) zv(nt)        one line per species (i4,3x,a2,3x,f5.2)
    na tau(1:3,na) ityp(na)  one line per atom    (i4,3x,3f15.9,3x,i2)
    plot(1 : nr1x*nr2x*nr3)                       (5(1pe17.9))

celldm(1) (alat) is in bohr, at and tau in units of alat. For ibrav /= 0 the
cell is latgen's lattice of ibrav and celldm. The values run x fastest over
the padded FFT dimensions nr1x and nr2x, and the grid points are
r = (i/nr1) a1 + (j/nr2) a2 + (k/nr3) a3, without the endpoint plane.
pp.x writes densities in e/bohr^3 and potentials in Ry as the potential
energy of an electron (INPUT_PP, plot_num). The file has no version marker,
so a header that does not fit these records is refused.
"""
from __future__ import annotations

import math
import re
import struct
import xml.etree.ElementTree as ET

from . import __version__
from ._legacy_writer import BOHR_PER_ANGSTROM
from ._textgeom import cartesian, element_of_label, export_geometry, fixed, volume, \
    write_structure
from .convert import (_MAX_VALUES, _binary, _export_report, _new_output, _number, _select,
                      _text, _validated)
from .reader import Reader
from .writer import Writer

_CARDS = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS", "OCCUPATIONS",
          "CONSTRAINTS", "ATOMIC_VELOCITIES", "ATOMIC_FORCES", "ADDITIONAL_K_POINTS",
          "SOLVENTS", "HUBBARD")
_REAL = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?")
_KEY = re.compile(r"([A-Za-z][A-Za-z0-9_%]*)\s*(?:\(\s*([0-9]+(?:\s*,\s*[0-9]+)*)\s*\))?\s*=")


# -- pw.x input ------------------------------------------------------------

def _value(token):
    if token[:1] in "'\"":
        return token[1:-1].replace(token[0] * 2, token[0])
    lowered = token.lower()
    if lowered in (".true.", ".t.", "t", "true"):
        return True
    if lowered in (".false.", ".f.", "f", "false"):
        return False
    if re.fullmatch(r"[+-]?[0-9]+", token):
        return int(token)
    if _REAL.fullmatch(token):
        return _number(token)
    raise ValueError(f"pw.x namelist value {token!r} is not a number, logical or string")


def _namelist_body(text, index):
    """({key: value}, index after the closing slash) of the namelist body
    starting at ``index``; ``!`` comments run to the end of the line."""
    values, i = {}, index
    while True:
        while i < len(text) and (text[i].isspace() or text[i] == ","):
            i += 1
        if i == len(text):
            raise ValueError("pw.x namelist is not closed by '/'")
        if text[i] == "!":
            i = text.find("\n", i)
            i = len(text) if i < 0 else i
            continue
        if text[i] == "/":
            return values, i + 1
        if text.startswith("&end", i) or text.startswith("&END", i):
            return values, i + 4
        match = _KEY.match(text, i)
        if not match:
            raise ValueError("pw.x namelist requires key = value assignments and a closing '/'")
        name = match.group(1).lower()
        if match.group(2):
            name += "(" + ",".join(v.strip() for v in match.group(2).split(",")) + ")"
        i = match.end()
        while i < len(text) and text[i] in " \t":
            i += 1
        if i < len(text) and text[i] in "'\"":
            quote, j = text[i], i + 1
            while True:
                j = text.find(quote, j)
                if j < 0:
                    raise ValueError("pw.x namelist string is not closed")
                if text[j + 1:j + 2] == quote:
                    j += 2
                    continue
                break
            token, i = text[i:j + 1], j + 1
        else:
            match = re.compile(r"[^\s,/!]+").match(text, i)
            if not match:
                raise ValueError(f"pw.x namelist key {name!r} has no value")
            token, i = match.group(), match.end()
        if name in values:
            raise ValueError(f"pw.x namelist key {name!r} is given twice")
        values[name] = _value(token)


def _split(text):
    """({namelist: {key: value}}, [(card, option, [lines])])."""
    namelists, i = {}, 0
    while True:
        while i < len(text) and text[i].isspace():
            i += 1
        if i < len(text) and text[i] in "!#":
            i = text.find("\n", i)
            i = len(text) if i < 0 else i
            continue
        if i == len(text) or text[i] != "&":
            break
        match = re.compile(r"&([A-Za-z]+)").match(text, i)
        if not match:
            raise ValueError("pw.x namelist requires a name after '&'")
        name = match.group(1).lower()
        if name in namelists:
            raise ValueError(f"pw.x namelist &{name} is given twice")
        namelists[name], i = _namelist_body(text, match.end())
    cards = []
    for line in text[i:].splitlines():
        stripped = line.strip()
        if not stripped or stripped[0] in "!#":
            continue
        words = stripped.split(None, 1)
        if words[0].upper() in _CARDS:
            option = words[1].strip() if len(words) > 1 else ""
            option = option.strip("{}() \t").lower()
            cards.append((words[0].upper(), option, []))
        elif not cards:
            raise ValueError(f"pw.x input has text that is neither a namelist nor a card: "
                             f"{stripped[:40]!r}")
        else:
            cards[-1][2].append(stripped.split("!")[0].split("#")[0])
    return namelists, cards


def _expression(token):
    """A real number or a QE card expression with + - * / ^ and parentheses."""
    tokens = re.findall(r"(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eEdD][+-]?[0-9]+)?|[-+*/^()]",
                        token)
    if "".join(tokens) != token:
        raise ValueError(f"pw.x position {token!r} is not a number or arithmetic expression")
    position = 0

    def take():
        nonlocal position
        position += 1
        return tokens[position - 1]

    def peek():
        return tokens[position] if position < len(tokens) else None

    def atom():
        t = take() if peek() is not None else None
        if t == "(":
            value = add()
            if take() != ")":
                raise ValueError
            return value
        if t in ("-", "+"):
            return -atom() if t == "-" else atom()
        if t is None or not t[0].isdigit() and t[0] != ".":
            raise ValueError
        return _number(t)

    def power():
        base = atom()
        if peek() == "^":
            take()
            return base ** power()
        return base

    def mul():
        value = power()
        while peek() in ("*", "/"):
            value = value * power() if take() == "*" else value / power()
        return value

    def add():
        value = mul()
        while peek() in ("+", "-"):
            value = value + mul() if take() == "+" else value - mul()
        return value

    try:
        value = add()
        if position != len(tokens):
            raise ValueError
        return _number(value)
    except (ValueError, ZeroDivisionError, IndexError, OverflowError):
        raise ValueError(f"pw.x position {token!r} is not a valid arithmetic expression") \
            from None


def _celldm(system, ibrav):
    """[None, alat in angstrom, b/a, c/a, celldm(4), celldm(5), celldm(6)]
    from celldm(1..6) (bohr) or A, B, C, cosAB, cosAC, cosBC (angstrom)."""
    dm = {k: system[f"celldm({k})"] for k in range(1, 7) if f"celldm({k})" in system}
    abc = {k: system[k] for k in ("a", "b", "c", "cosab", "cosac", "cosbc") if k in system}
    if any(not isinstance(v, (int, float)) or isinstance(v, bool)
           for v in list(dm.values()) + list(abc.values())):
        raise ValueError("pw.x celldm and A, B, C, cosAB, cosAC, cosBC must be numbers")
    if dm and abc:
        raise ValueError("pw.x input gives both celldm and A, B, C: use one")
    if any(f"celldm({k})" in system for k in range(7, 100)) or "celldm(0)" in system:
        raise ValueError("pw.x celldm has indices 1 to 6")
    if dm:
        c = [None, dm.get(1, 0) / BOHR_PER_ANGSTROM] + [dm.get(k, 0.0) for k in range(2, 7)]
        source = "celldm"
    elif abc:
        a = abc.get("a", 0.0)
        c = [None, a, abc.get("b", 0.0) / a if a else 0.0, abc.get("c", 0.0) / a if a else 0.0,
             0.0, 0.0, 0.0]
        if ibrav in (14, 0):
            c[4:7] = [abc.get("cosbc", 0.0), abc.get("cosac", 0.0), abc.get("cosab", 0.0)]
        elif ibrav in (-12, -13):
            c[5] = abc.get("cosac", 0.0)
        else:
            c[4] = abc.get("cosab", 0.0)
        source = "A"
    else:
        return None, None
    return c, source


def _bravais(ibrav, c):
    """The lattice rows of ``ibrav`` in units of alat (INPUT_PW, latgen)."""
    ba, ca, c4, c5, c6 = c[2:7]
    needs_b = ibrav in (8, 9, -9, 91, 10, 11, 12, -12, 13, -13, 14)
    needs_c = ibrav in (4, 6, 7) or needs_b
    if needs_b and ba <= 0 or needs_c and ca <= 0:
        raise ValueError(f"pw.x ibrav = {ibrav} requires positive b/a and c/a")
    r3 = math.sqrt(3.0)
    if ibrav == 1:
        return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    if ibrav == 2:
        return [[-0.5, 0, 0.5], [0, 0.5, 0.5], [-0.5, 0.5, 0]]
    if ibrav == 3:
        return [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]]
    if ibrav == -3:
        return [[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]
    if ibrav == 4:
        return [[1, 0, 0], [-0.5, r3 / 2, 0], [0, 0, ca]]
    if ibrav in (5, -5):
        if not -0.5 < c4 < 1:
            raise ValueError("pw.x ibrav = 5 or -5 requires -1/2 < cos(alpha) < 1")
        tx, ty, tz = math.sqrt((1 - c4) / 2), math.sqrt((1 - c4) / 6), math.sqrt((1 + 2 * c4) / 3)
        if ibrav == 5:
            return [[tx, -ty, tz], [0, 2 * ty, tz], [-tx, -ty, tz]]
        u, v = (tz - 2 * math.sqrt(2) * ty) / r3, (tz + math.sqrt(2) * ty) / r3
        return [[u, v, v], [v, u, v], [v, v, u]]
    if ibrav == 6:
        return [[1, 0, 0], [0, 1, 0], [0, 0, ca]]
    if ibrav == 7:
        return [[0.5, -0.5, ca / 2], [0.5, 0.5, ca / 2], [-0.5, -0.5, ca / 2]]
    if ibrav == 8:
        return [[1, 0, 0], [0, ba, 0], [0, 0, ca]]
    if ibrav == 9:
        return [[0.5, ba / 2, 0], [-0.5, ba / 2, 0], [0, 0, ca]]
    if ibrav == -9:
        return [[0.5, -ba / 2, 0], [0.5, ba / 2, 0], [0, 0, ca]]
    if ibrav == 91:
        return [[1, 0, 0], [0, ba / 2, -ca / 2], [0, ba / 2, ca / 2]]
    if ibrav == 10:
        return [[0.5, 0, ca / 2], [0.5, ba / 2, 0], [0, ba / 2, ca / 2]]
    if ibrav == 11:
        return [[0.5, ba / 2, ca / 2], [-0.5, ba / 2, ca / 2], [-0.5, -ba / 2, ca / 2]]
    if ibrav in (12, 13):
        if not -1 < c4 < 1:
            raise ValueError(f"pw.x ibrav = {ibrav} requires -1 < cos(gamma) < 1")
        b = [ba * c4, ba * math.sqrt(1 - c4 * c4), 0]
        if ibrav == 12:
            return [[1, 0, 0], b, [0, 0, ca]]
        return [[0.5, 0, -ca / 2], b, [0.5, 0, ca / 2]]
    if ibrav in (-12, -13):
        if not -1 < c5 < 1:
            raise ValueError(f"pw.x ibrav = {ibrav} requires -1 < cos(beta) < 1")
        v3 = [ca * c5, 0, ca * math.sqrt(1 - c5 * c5)]
        if ibrav == -12:
            return [[1, 0, 0], [0, ba, 0], v3]
        return [[0.5, ba / 2, 0], [-0.5, ba / 2, 0], v3]
    if ibrav == 14:
        cos_a, cos_b, cos_g = c4, c5, c6
        if not all(-1 < v < 1 for v in (cos_a, cos_b, cos_g)):
            raise ValueError("pw.x ibrav = 14 requires cosines between -1 and 1")
        sin_g = math.sqrt(1 - cos_g * cos_g)
        root = 1 + 2 * cos_a * cos_b * cos_g - cos_a ** 2 - cos_b ** 2 - cos_g ** 2
        if root <= 0:
            raise ValueError("pw.x ibrav = 14 angles do not describe a cell")
        return [[1, 0, 0], [ba * cos_g, ba * sin_g, 0],
                [ca * cos_b, ca * (cos_a - cos_b * cos_g) / sin_g, ca * math.sqrt(root) / sin_g]]
    raise ValueError(f"pw.x ibrav = {ibrav} is not a Bravais lattice index of INPUT_PW")


def _card(cards, name):
    found = [c for c in cards if c[0] == name]
    if len(found) > 1:
        raise ValueError(f"pw.x card {name} is given twice")
    return found[0] if found else None


def _input_cell(system, cards):
    """(lattice rows in angstrom, alat in angstrom, cell report)."""
    ibrav = system.get("ibrav")
    if type(ibrav) is not int:
        raise ValueError("pw.x &SYSTEM requires an integer ibrav")
    c, source = _celldm(system, ibrav)
    card = _card(cards, "CELL_PARAMETERS")
    if ibrav != 0:
        if card is not None:
            raise ValueError("pw.x CELL_PARAMETERS requires ibrav = 0")
        if c is None or c[1] <= 0:
            raise ValueError(f"pw.x ibrav = {ibrav} requires celldm(1) or A")
        alat = c[1]
        lattice = [[v * alat for v in row] for row in _bravais(ibrav, c)]
        return lattice, alat, {"ibrav": ibrav, "cell_source": source}
    if card is None:
        raise ValueError("pw.x ibrav = 0 requires a CELL_PARAMETERS card")
    option = card[1]
    if option not in ("alat", "bohr", "angstrom", ""):
        raise ValueError(f"pw.x CELL_PARAMETERS unit {option!r} is not alat, bohr or angstrom")
    if c is not None and c[1] <= 0:
        raise ValueError("pw.x celldm(1) or A must be positive")
    if option == "":
        # Deprecated: alat when celldm(1) or A is given, else bohr (QE).
        option = "alat" if c is not None else "bohr"
        deprecated = True
    else:
        deprecated = False
    if option == "alat" and c is None:
        raise ValueError("pw.x CELL_PARAMETERS alat requires celldm(1) or A")
    if option != "alat" and c is not None:
        raise ValueError(f"pw.x CELL_PARAMETERS {option} must not be combined with celldm(1) "
                         "or A (the lattice parameter would be given twice)")
    if len(card[2]) < 3:
        raise ValueError("pw.x CELL_PARAMETERS requires three vectors")
    rows = []
    for line in card[2][:3]:
        tokens = line.split()
        if len(tokens) != 3 or not all(_REAL.fullmatch(t) for t in tokens):
            raise ValueError("pw.x CELL_PARAMETERS vectors require three numbers")
        rows.append([_number(t) for t in tokens])
    if len(card[2]) > 3:
        raise ValueError("pw.x CELL_PARAMETERS has more than three vectors")
    factor = {"alat": c[1] if c else None, "bohr": 1 / BOHR_PER_ANGSTROM, "angstrom": 1.0}[option]
    lattice = [[v * factor for v in row] for row in rows]
    # Without celldm(1) or A, alat is the length of the first vector (QE).
    alat = c[1] if c is not None else math.sqrt(sum(v * v for v in lattice[0]))
    report = {"ibrav": 0, "cell_units": option, "cell_units_defaulted": deprecated}
    if source:
        report["cell_source"] = source
    return lattice, alat, report


_R_GROUPS = (146, 148, 155, 160, 161, 166, 167)


# Conventional lattice of a space-group input, as the pw.x Bravais index whose
# celldm (or A, B, C, cosAB, cosAC, cosBC) convention gives its parameters.
def _space_group_ibrav(number, uniqueb, rhombohedral):
    if number <= 2:
        return 14
    if number <= 15:
        return -12 if uniqueb else 12
    if number <= 74:
        return 8
    if number <= 142:
        return 6
    if number <= 194:
        return 5 if number in _R_GROUPS and rhombohedral else 4
    return 1


def _space_group(system, cards):
    """The conventional cell and operators of a pw.x space_group input."""
    from ._spacegroups import ITA_HALL_ORIGIN_2, ita_hall, operators
    number = system["space_group"]
    if type(number) is not int or not 1 <= number <= 230:
        raise ValueError("pw.x space_group must be an integer from 1 to 230")
    uniqueb = system.get("uniqueb", False)
    rhombohedral = system.get("rhombohedral", True)
    origin = system.get("origin_choice", 1)
    if type(uniqueb) is not bool or type(rhombohedral) is not bool:
        raise ValueError("pw.x uniqueb and rhombohedral must be logicals")
    if origin not in (1, 2) or type(origin) is not int:
        raise ValueError("pw.x origin_choice must be 1 or 2")
    if _card(cards, "CELL_PARAMETERS") is not None:
        raise ValueError("pw.x space_group input takes the cell from celldm or A, B, C, "
                         "not from CELL_PARAMETERS")
    ibrav = _space_group_ibrav(number, uniqueb, rhombohedral)
    c, source = _celldm(system, ibrav)
    if c is None or c[1] <= 0:
        raise ValueError("pw.x space_group input requires celldm(1) or A")
    lattice = [[v * c[1] for v in row] for row in _bravais(ibrav, c)]
    hall = ita_hall(number, unique_axis="b" if uniqueb else "c", origin_choice=origin,
                    rhombohedral=rhombohedral)
    report = {"space_group": number, "hall": hall, "conventional_ibrav": ibrav,
              "cell_source": source, "uniqueb": uniqueb, "rhombohedral": rhombohedral,
              "origin_choice": origin,
              "settings_ignored": [k for k, applies in (
                  ("ibrav", False), ("origin_choice", number in ITA_HALL_ORIGIN_2),
                  ("rhombohedral", number in _R_GROUPS), ("uniqueb", 3 <= number <= 15))
                  if k in system and not applies]}
    return {"lattice": lattice, "alat": c[1], "ops": operators(hall), "report": report}


def _expand_sites(sites, group, lattice):
    """(atoms, frozen, report) of the unit cell the space group generates from
    the crystal_sg sites. Images closer than 0.01 angstrom are one atom."""
    from ._spacegroups import orbit
    atoms, frozen, merged = [], [], 0
    for symbol, f, fixed_site in sites:
        images, dropped = orbit(f, group["ops"], lattice, 0.01)
        merged += dropped
        for g in images:
            x = cartesian(g, lattice)
            if any(math.dist(x, y) < 0.01 for _, y in _wrapped(atoms, g, lattice)):
                raise ValueError("pw.x crystal_sg sites generate overlapping atoms")
            if fixed_site:
                frozen.append(len(atoms))
            atoms.append((symbol, x, g))
    return [(s, x) for s, x, _ in atoms], frozen, {
        "sites": len(sites), "operators": len(group["ops"]),
        "special_position_images_merged": merged}


def _wrapped(atoms, g, lattice):
    """(symbol, the image of each earlier atom nearest to fractional g)."""
    for symbol, _, h in atoms:
        d = [b + round(a - b) for a, b in zip(g, h)]
        yield symbol, cartesian(d, lattice)


def _input_structure(text):
    namelists, cards = _split(text)
    system = namelists.get("system")
    if system is None:
        raise ValueError("pw.x input requires a &SYSTEM namelist")
    nat, ntyp = system.get("nat"), system.get("ntyp")
    if type(nat) is not int or type(ntyp) is not int or nat < 1 or ntyp < 1:
        raise ValueError("pw.x &SYSTEM requires positive integer nat and ntyp")
    group = _space_group(system, cards) if "space_group" in system else None
    if group is None:
        lattice, alat, cell = _input_cell(system, cards)
    else:
        lattice, alat, cell = group["lattice"], group["alat"], group["report"]
    if volume(lattice) < 1e-12:
        raise ValueError("pw.x lattice vectors are linearly dependent")
    species_card = _card(cards, "ATOMIC_SPECIES")
    if species_card is None or len(species_card[2]) != ntyp:
        raise ValueError("pw.x ATOMIC_SPECIES requires ntyp species lines")
    species = {}
    for line in species_card[2]:
        tokens = line.split()
        if len(tokens) < 3:
            raise ValueError("pw.x ATOMIC_SPECIES lines need label, mass and pseudopotential")
        if tokens[0] in species:
            raise ValueError(f"pw.x species {tokens[0]!r} is given twice")
        species[tokens[0]] = element_of_label(tokens[0], "pw.x")
    positions = _card(cards, "ATOMIC_POSITIONS")
    if positions is None:
        raise ValueError("pw.x input requires an ATOMIC_POSITIONS card")
    option = positions[1] or "alat"
    if (option == "crystal_sg") != (group is not None):
        raise ValueError("pw.x ATOMIC_POSITIONS crystal_sg and the &SYSTEM space_group "
                         "require each other")
    if option not in ("alat", "bohr", "angstrom", "crystal", "crystal_sg"):
        raise ValueError(f"pw.x ATOMIC_POSITIONS unit {option!r} is not alat, bohr, angstrom "
                         "or crystal")
    if len(positions[2]) != nat:
        raise ValueError("pw.x ATOMIC_POSITIONS requires nat position lines")
    atoms, frozen, partial, sites = [], [], 0, []
    for k, line in enumerate(positions[2]):
        tokens = line.split()
        if option == "crystal_sg" and len(tokens) > 1 and re.fullmatch(r"[0-9]+[A-Za-z]",
                                                                        tokens[1]):
            raise ValueError(f"pw.x crystal_sg Wyckoff label {tokens[1]!r} is not supported: "
                             "give the site as x y z in crystal coordinates")
        if len(tokens) not in (4, 7):
            raise ValueError("pw.x position lines need a label, three coordinates and "
                             "optionally three if_pos flags")
        if tokens[0] not in species:
            raise ValueError(f"pw.x position label {tokens[0]!r} is not in ATOMIC_SPECIES")
        x = [_expression(t) for t in tokens[1:4]]
        if len(tokens) == 7:
            if any(t not in ("0", "1") for t in tokens[4:]):
                raise ValueError("pw.x if_pos flags must be 0 or 1")
            if tokens[4:] == ["0"] * 3:
                frozen.append(k)
            elif "0" in tokens[4:]:
                partial += 1
        if option == "crystal_sg":
            sites.append((species[tokens[0]], x, k in frozen))
            continue
        if option == "crystal":
            x = cartesian(x, lattice)
        else:
            factor = {"alat": alat, "bohr": 1 / BOHR_PER_ANGSTROM, "angstrom": 1.0}[option]
            x = [v * factor for v in x]
        atoms.append((species[tokens[0]], x))
    if group is not None:
        atoms, frozen, expansion = _expand_sites(sites, group, lattice)
        cell.update(expansion)
    labels = {k: v for k, v in species.items()}
    omitted = sorted({c[0] for c in cards} - {"ATOMIC_SPECIES", "ATOMIC_POSITIONS",
                                               "CELL_PARAMETERS"})
    title = namelists.get("control", {}).get("title")
    return {"atoms": atoms, "lattice": lattice, "frozen": frozen, "partial": partial,
            "label": title if isinstance(title, str) and title.strip() else None,
            "report": {**cell, "alat_angstrom": alat, "positions_units": option,
                       "positions_units_defaulted": not positions[1],
                       "species_labels": labels, "namelists": sorted(namelists),
                       "cards_omitted": omitted}}


def import_espresso_in(source, destination):
    head = _input_structure(_text(source))
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "pw.x input structure import") \
                as writer:
            write_structure(writer, head["atoms"], head["lattice"], head["label"],
                            head["frozen"])
        findings = _validated(temporary)
    return {"format": "espresso-in", "output": str(destination), "atoms": len(head["atoms"]),
            **head["report"], "cell_volume_angstrom3": volume(head["lattice"]),
            "frozen_atoms": len(head["frozen"]), "partially_fixed_atoms_omitted": head["partial"],
            "findings": findings,
            "scope": "Cell, species elements, positions and atoms fixed on all three axes (as "
                     "constraints.frozen_atoms) only; namelist settings, pseudopotentials, "
                     "masses, k-points, magnetization, charge and the job are not stored. Keep "
                     "the original pw.x input."}


def export_espresso_in(source, destination, *, section=None):
    _validated(source)
    with Reader(source) as reader:
        selected = _select(reader, section, {"structure"})
        atoms, lattice, frozen, omitted, label = export_geometry(reader, selected, "pw.x",
                                                                 periodic_only=True)
        elements = list(dict.fromkeys(s for s, _ in atoms))
        title = label.replace("'", "").replace('"', "")
        lines = ["&CONTROL", f"  title = '{title}'", "/", "&SYSTEM", "  ibrav = 0",
                 f"  nat = {len(atoms)}", f"  ntyp = {len(elements)}",
                 "  ! Add ecutwfc and the other settings of the calculation here.", "/",
                 "&ELECTRONS", "/", "ATOMIC_SPECIES"]
        # A mass of 0 makes pw.x use its own table of atomic masses.
        lines += [f"  {s:<3s} 0.0 {s}.UPF" for s in elements]
        lines.append("CELL_PARAMETERS angstrom")
        lines += ["  " + " ".join(f"{fixed(v):>22s}" for v in row) for row in lattice]
        lines.append("ATOMIC_POSITIONS angstrom")
        for k, (symbol, x) in enumerate(atoms):
            flags = (" 0 0 0" if k in frozen else " 1 1 1") if frozen else ""
            lines.append(f"  {symbol:<3s}" + " ".join(f"{fixed(v):>22s}" for v in x) + flags)
        with _new_output(destination) as temporary:
            temporary.write_text("\n".join(lines) + "\n", encoding="ascii")
        report = _export_report(reader, selected, {selected.id}, destination, "espresso-in")
        report.update(atoms=len(atoms), species=elements, frozen_atoms=len(frozen),
                      constraints_not_written=omitted,
                      generated="&CONTROL title, &SYSTEM ibrav = 0, nat and ntyp, one species "
                                "per element with mass 0.0 (pw.x then uses its own masses) and "
                                "the placeholder pseudopotential <element>.UPF, CELL_PARAMETERS "
                                "and ATOMIC_POSITIONS in angstrom with 15 decimals; add "
                                "ecutwfc, pseudopotentials and k-points before running pw.x.")
        return report


# -- pw.x XML data file ----------------------------------------------------

def _local(element):
    return element.tag.rsplit("}", 1)[-1]


def _child(element, name, required=True):
    found = [c for c in element if _local(c) == name]
    if len(found) > 1:
        raise ValueError(f"pw.x XML has more than one <{name}> in <{_local(element)}>")
    if not found and required:
        raise ValueError(f"pw.x XML requires <{name}> in <{_local(element)}>")
    return found[0] if found else None


def _numbers(element, count, label):
    tokens = (element.text or "").split()
    if len(tokens) != count or not all(_REAL.fullmatch(t) for t in tokens):
        raise ValueError(f"pw.x XML <{label}> requires {count} numbers")
    return [_number(t) for t in tokens]


def _flag(element):
    if element is None:
        return None
    text = (element.text or "").strip().lower()
    if text not in ("true", "false"):
        raise ValueError(f"pw.x XML <{_local(element)}> must be true or false")
    return text == "true"


def import_espresso_xml(source, destination):
    text = _text(source)
    if re.search(r"<!\s*(DOCTYPE|ENTITY)", text, re.IGNORECASE):
        raise ValueError("pw.x XML must not contain a DOCTYPE or entity declarations")
    try:
        root = ET.fromstring(text)
    except ET.ParseError as error:
        raise ValueError(f"pw.x XML is not well-formed: {error}") from None
    if _local(root) != "espresso":
        raise ValueError("pw.x XML requires the qes:espresso root element")
    output = _child(root, "output")
    structure = _child(output, "atomic_structure")
    nat = structure.get("nat", "")
    if not re.fullmatch(r"[0-9]+", nat) or int(nat) < 1:
        raise ValueError("pw.x XML atomic_structure requires a positive nat")
    nat = int(nat)
    cell = _child(structure, "cell")
    lattice = [[v / BOHR_PER_ANGSTROM for v in _numbers(_child(cell, f"a{i}"), 3, f"a{i}")]
               for i in (1, 2, 3)]
    if volume(lattice) < 1e-12:
        raise ValueError("pw.x XML lattice vectors are linearly dependent")
    cartesian_block = _child(structure, "atomic_positions", required=False)
    crystal_block = _child(structure, "crystal_positions", required=False)
    if (cartesian_block is None) == (crystal_block is None):
        raise ValueError("pw.x XML atomic_structure requires atomic_positions or "
                         "crystal_positions (Wyckoff positions are not supported)")
    block = cartesian_block if cartesian_block is not None else crystal_block
    rows = [a for a in block if _local(a) == "atom"]
    if len(rows) != nat or len(rows) != len(list(block)):
        raise ValueError("pw.x XML requires one <atom> per atom of nat")
    atoms, labels = [], {}
    for row in rows:
        name = row.get("name", "")
        labels[name] = element_of_label(name, "pw.x XML")
        x = _numbers(row, 3, "atom")
        x = cartesian(x, lattice) if block is crystal_block else [v / BOHR_PER_ANGSTROM
                                                                 for v in x]
        atoms.append((labels[name], x))
    frozen, partial = [], 0
    inputs = _child(root, "input", required=False)
    free = _child(inputs, "free_positions", required=False) if inputs is not None else None
    if free is not None:
        flags = (free.text or "").split()
        if len(flags) != 3 * nat or any(f not in ("0", "1") for f in flags):
            raise ValueError("pw.x XML free_positions requires 0 or 1 for each coordinate")
        for k in range(nat):
            mask = flags[3 * k:3 * k + 3]
            if mask == ["0"] * 3:
                frozen.append(k)
            elif "0" in mask:
                partial += 1
    convergence = _child(output, "convergence_info", required=False)
    scf = opt = None
    if convergence is not None:
        for name in ("scf_conv", "opt_conv"):
            block_ = _child(convergence, name, required=False)
            value = _flag(_child(block_, "convergence_achieved", required=False)) \
                if block_ is not None else None
            scf, opt = (value, opt) if name == "scf_conv" else (scf, value)
    results = scf is not False
    energy_block = _child(output, "total_energy", required=False)
    energy = None
    if energy_block is not None and _child(energy_block, "etot", required=False) is not None:
        energy = _numbers(_child(energy_block, "etot"), 1, "etot")[0]
    forces = _child(output, "forces", required=False)
    gradient = [-v for v in _numbers(forces, 3 * nat, "forces")] if forces is not None else None
    stress_block = _child(output, "stress", required=False)
    # pw.x writes sigma = -(1/V) dE/d(strain) (positive pressure under
    # compression); QVF stores (1/V) dE/d(strain).
    stress = [-v for v in _numbers(stress_block, 9, "stress")] if stress_block is not None \
        else None
    if stress is not None and any(abs(stress[3 * i + j] - stress[3 * j + i]) > 1e-9 * max(
            1e-12, max(abs(v) for v in stress)) for i in range(3) for j in range(i)):
        raise ValueError("pw.x XML stress is not symmetric")
    exit_status = _child(root, "exit_status", required=False)
    used = {"atomic_structure", "total_energy", "forces", "stress", "convergence_info"}
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "pw.x XML result import") as writer:
            write_structure(writer, atoms, lattice, None, frozen)
            if results and energy is not None:
                writer.set_energies(energy)
            if results and gradient is not None:
                writer.add_section({"id": "forces", "kind": "forces",
                                    "structure_ref": "structure"})
                _binary(writer, "forces", "gradient", gradient, [nat, 3])
                if stress is not None:
                    _binary(writer, "forces", "stress", stress, [3, 3])
        findings = _validated(temporary)
    dft = _child(output, "dft", required=False)
    functional = _child(dft, "functional", required=False) if dft is not None else None
    return {"format": "espresso-xml", "output": str(destination), "atoms": nat,
            "species_labels": labels, "bravais_index": structure.get("bravais_index"),
            "cell_volume_angstrom3": volume(lattice),
            "scf_converged": scf, "optimization_converged": opt,
            "exit_status": (exit_status.text or "").strip() if exit_status is not None else None,
            "functional": (functional.text or "").strip() if functional is not None else None,
            "energy_hartree": energy if results else None,
            "forces": results and gradient is not None,
            "stress": results and gradient is not None and stress is not None,
            "stress_omitted_without_forces": results and gradient is None and stress is not None,
            "results_omitted_unconverged": not results,
            "frozen_atoms": len(frozen), "partially_fixed_atoms_omitted": partial,
            "output_elements_omitted": sorted({_local(c) for c in output} - used),
            "findings": findings,
            "scope": "Final structure, total energy, forces and stress only; band structure, "
                     "magnetization, symmetry and input settings are not stored, and results "
                     "of an SCF that did not converge are omitted. Keep the original data "
                     "directory."}


# -- pp.x cube -------------------------------------------------------------

def import_espresso_cube(source, destination, *, kind="volume.density"):
    from .convert import import_cube
    report = import_cube(source, destination, kind=kind, periodic=(True, True, True),
                         includes_endpoint=False)
    with Reader(destination) as reader:
        grid = reader.json("volume", "grid")
        lattice = reader.json("structure", "structure")["lattice_vectors"]
        total = sum(v for (v,) in struct.iter_unpack("<d", reader.tensor("volume", "data").data))
    cell = volume(lattice)
    report.update(format="espresso-cube", shape=grid["shape"], includes_endpoint=False,
                  origin_bohr=grid["origin"], cell_volume_angstrom3=cell,
                  integral=total / math.prod(grid["shape"]) * cell * BOHR_PER_ANGSTROM ** 3,
                  scope="pp.x writes the FFT grid without the endpoint plane, so the grid is "
                        "periodic with includes_endpoint false and the cell is shape times the "
                        "voxel vectors. Values are stored as written (e/bohr^3 for plot_num "
                        "0, where integral is the electron count); the kind is asserted by the "
                        "caller, not read from the file.")
    return report


# -- pp.x filplot ----------------------------------------------------------

# plot_num -> (kind, label, factor, declared unit). Values are multiplied by
# factor. pp.x writes densities in e/bohr^3 and potentials as the potential
# energy of an electron in Ry: as an electrostatic potential (hartree/e) that
# is -V/2, as an electron potential energy in hartree V/2.
_PLOT_NUM = {
    0: ("volume.density", "Charge density (pp.x plot_num 0)", 1.0, "e_per_bohr3"),
    1: ("volume.generic", "Total local potential V_bare + V_H + V_xc, electron potential "
        "energy (pp.x plot_num 1)", 0.5, "hartree"),
    2: ("volume.generic", "Local ionic potential V_bare, electron potential energy "
        "(pp.x plot_num 2)", 0.5, "hartree"),
    6: ("volume.spin", "Spin polarization rho_up - rho_down (pp.x plot_num 6)", 1.0,
        "e_per_bohr3"),
    8: ("volume.elf", "Electron localization function (pp.x plot_num 8)", 1.0,
        "dimensionless"),
    9: ("volume.difference", "Charge density minus superposition of atomic densities "
        "(pp.x plot_num 9)", 1.0, "e_per_bohr3"),
    11: ("volume.potential", "Electrostatic potential from V_bare + V_H (pp.x plot_num 11)",
         -0.5, "hartree_per_e"),
    12: ("volume.potential", "Sawtooth electric-field potential (pp.x plot_num 12)", -0.5,
         "hartree_per_e"),
    17: ("volume.density", "All-electron valence charge density (pp.x plot_num 17)", 1.0,
         "e_per_bohr3"),
    19: ("volume.rdg", "Reduced density gradient (pp.x plot_num 19)", 1.0, "dimensionless"),
    21: ("volume.density", "All-electron charge density (pp.x plot_num 21)", 1.0,
         "e_per_bohr3"),
}
_FORTRAN_REAL = re.compile(r"([+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+))(?:[eEdD]?([+-][0-9]+))?")


def _fortran_real(token):
    """A Fortran E-format real. Ew.d drops the E of a three-digit exponent."""
    match = _FORTRAN_REAL.fullmatch(token)
    if not match:
        raise ValueError(f"pp.x filplot value {token!r} is not a number")
    mantissa, exponent = match.groups()
    return _number(float(mantissa + ("e" + exponent if exponent else "")))


def _fields(line, count, widths, label):
    """``count`` whitespace-separated fields, else the fixed-width fields of
    the Fortran format (adjacent fields may touch)."""
    tokens = line.split()
    if len(tokens) == count:
        return tokens
    tokens, position = [], 0
    for width in widths:
        tokens.append(line[position:position + width].strip()); position += width
    if any(not t for t in tokens) or line[position:].strip():
        raise ValueError(f"pp.x filplot {label} line needs {count} fields")
    return tokens


def import_espresso_filplot(source, destination):
    lines = _text(source).splitlines()
    if len(lines) < 5:
        raise ValueError("pp.x filplot header is truncated")
    try:
        dims = [int(v) for v in _fields(lines[1], 8, [8] * 8, "dimension")]
        head = _fields(lines[2], 7, [6] + [12] * 6, "ibrav/celldm")
        ibrav, celldm = int(head[0]), [_number(v) for v in head[1:]]
    except ValueError as error:
        raise ValueError(f"pp.x filplot header is not the plot_io layout: {error}") from None
    nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp = dims
    if min(dims) < 1 or nr1 > nr1x or nr2 > nr2x or nr3 > nr3x:
        raise ValueError("pp.x filplot needs positive counts and FFT dimensions nr <= nrx")
    if celldm[0] <= 0:
        raise ValueError("pp.x filplot celldm(1) must be positive")
    alat = celldm[0] / BOHR_PER_ANGSTROM
    index = 3
    if ibrav == 0:
        rows = []
        for line in lines[3:6]:
            tokens = line.split()
            if len(tokens) != 3 or not all(_REAL.fullmatch(t) for t in tokens):
                raise ValueError("pp.x filplot ibrav = 0 requires three lattice vectors in alat")
            rows.append([_number(t) for t in tokens])
        index = 6
    else:
        rows = _bravais(ibrav, [None, alat] + celldm[1:])
    lattice = [[v * alat for v in row] for row in rows]
    if volume(lattice) < 1e-12:
        raise ValueError("pp.x filplot lattice vectors are linearly dependent")
    tokens = lines[index].split() if index < len(lines) else []
    if len(tokens) != 4 or not re.fullmatch(r"-?[0-9]+", tokens[3]):
        raise ValueError("pp.x filplot requires the gcutm, dual, ecutwfc, plot_num line")
    gcutm, dual, ecutwfc = (_number(t) for t in tokens[:3])
    plot_num = int(tokens[3])
    index += 1
    if len(lines) < index + ntyp + nat:
        raise ValueError("pp.x filplot species or atom lines are truncated")
    species, valence = {}, {}
    for line in lines[index:index + ntyp]:
        tokens = line.split()
        if len(tokens) != 3 or not tokens[0].isdigit() or int(tokens[0]) in species:
            raise ValueError("pp.x filplot species lines need index, label and valence")
        species[int(tokens[0])] = (tokens[1], element_of_label(tokens[1], "pp.x filplot"))
        valence[tokens[1]] = _number(tokens[2])
    if sorted(species) != list(range(1, ntyp + 1)):
        raise ValueError("pp.x filplot species indices must be 1 to ntyp")
    index += ntyp
    atoms = []
    for line in lines[index:index + nat]:
        tokens = line.split()
        if len(tokens) != 5 or not tokens[4].isdigit() or int(tokens[4]) not in species:
            raise ValueError("pp.x filplot atom lines need index, x, y, z (alat) and species")
        atoms.append((species[int(tokens[4])][1], [_number(t) * alat for t in tokens[1:4]]))
    index += nat
    count = nr1x * nr2x * nr3
    if count > _MAX_VALUES:
        raise ValueError("pp.x filplot grid exceeds the value limit")
    values = [_fortran_real(t) for line in lines[index:] for t in line.split()]
    if len(values) != count:
        raise ValueError(f"pp.x filplot requires nr1x * nr2x * nr3 = {count} values, "
                         f"found {len(values)}")
    kind, label, factor, unit = _PLOT_NUM.get(
        plot_num, ("volume.generic", f"pp.x plot_num {plot_num}, values as written", 1.0, None))
    shape = [nr1, nr2, nr3]
    grid = {"origin": [0.0, 0.0, 0.0],
            "voxel_vectors": [[v * BOHR_PER_ANGSTROM / shape[a] for v in row]
                              for a, row in enumerate(lattice)],
            "shape": shape, "periodic": [True] * 3, "includes_endpoint": False,
            "structure_ref": "structure"}
    if unit:
        grid["units"] = {"data": unit}
    # Fortran order x + nr1x * (y + nr2x * z), QVF order z + nr3 * (y + nr2 * x).
    order = (x + nr1x * (y + nr2x * z) for x in range(nr1) for y in range(nr2)
             for z in range(nr3))
    with _new_output(destination) as temporary:
        with Writer(temporary, "qvf.convert", __version__, "pp.x filplot import") as writer:
            write_structure(writer, atoms, lattice, lines[0].strip() or None, [])
            writer.add_section({"id": "volume", "kind": kind, "label": label,
                                "structure_ref": "structure"})
            writer.add_json("volume", "grid", grid)
            _binary(writer, "volume", "data", (values[i] * factor + 0.0 for i in order), shape)
        findings = _validated(temporary)
    cell = volume(lattice)
    used = [values[x + nr1x * (y + nr2x * z)] for z in range(nr3) for y in range(nr2)
            for x in range(nr1)]
    mean = sum(used) / len(used)
    report = {"format": "espresso-filplot", "output": str(destination), "plot_num": plot_num,
              "kind": kind, "units": unit, "factor": factor, "shape": shape,
              "fft_dimensions": [nr1x, nr2x, nr3x], "ibrav": ibrav, "alat_angstrom": alat,
              "atoms": nat, "species_labels": {v[0]: v[1] for v in species.values()},
              "valence_charges_omitted": valence, "gcutm": gcutm, "dual": dual,
              "ecutwfc_ry": ecutwfc, "includes_endpoint": False,
              "cell_volume_angstrom3": cell, "findings": findings,
              "scope": "Structure and one scalar field on the FFT grid, periodic without the "
                       "endpoint plane. The kind and unit follow from plot_num. A plot_num "
                       "outside the table is volume.generic with the values as written. "
                       "Keep the original filplot."}
    if unit == "e_per_bohr3":
        report["integral"] = mean * cell * BOHR_PER_ANGSTROM ** 3
    return report
