"""Crystallographic symmetry operators from Hall symbols, standard library only.

A small Hall-symbol parser (S. R. Hall, Acta Cryst. A37, 517-525, 1981), the
Hall symbols of the 65 chiral (Sohncke) space groups under their wwPDB
CRYST1 names, with the settings macromolecular PDB files commonly use, and
the Hall symbols of all 230 space groups by ITA number (``ita_hall``). An
operator is (R, t): R a 3x3 integer matrix on fractional coordinates and t a
fractional translation in [0, 1) as Fractions.
"""

from fractions import Fraction

# wwPDB CRYST1 space-group symbol -> Hall symbol. "H" names the hexagonal
# setting of a rhombohedral group and "R" the rhombohedral axes.
PDB_HALL = {
    "P 1": "P 1", "P 1 2 1": "P 2y", "P 1 21 1": "P 2yb", "C 1 2 1": "C 2y",
    "I 1 2 1": "I 2y", "A 1 2 1": "A 2y", "P 2 2 2": "P 2 2", "P 2 2 21": "P 2c 2",
    "P 2 21 2": "P 2 2b", "P 21 2 2": "P 2a 2a", "P 21 21 2": "P 2 2ab",
    "P 21 2 21": "P 2ac 2ac", "P 2 21 21": "P 2bc 2", "P 21 21 21": "P 2ac 2ab",
    "C 2 2 21": "C 2c 2", "C 2 2 2": "C 2 2", "F 2 2 2": "F 2 2", "I 2 2 2": "I 2 2",
    "I 21 21 21": "I 2b 2c", "P 4": "P 4", "P 41": "P 4w", "P 42": "P 4c",
    "P 43": "P 4cw", "I 4": "I 4", "I 41": "I 4bw", "P 4 2 2": "P 4 2",
    "P 4 21 2": "P 4ab 2ab", "P 41 2 2": "P 4w 2c", "P 41 21 2": "P 4abw 2nw",
    "P 42 2 2": "P 4c 2", "P 42 21 2": "P 4n 2n", "P 43 2 2": "P 4cw 2c",
    "P 43 21 2": "P 4nw 2abw", "I 4 2 2": "I 4 2", "I 41 2 2": "I 4bw 2bw",
    "P 3": "P 3", "P 31": "P 31", "P 32": "P 32", "H 3": "R 3", "R 3": "P 3*",
    "P 3 1 2": "P 3 2", "P 3 2 1": 'P 3 2"', "P 31 1 2": "P 31 2 (0 0 4)",
    "P 31 2 1": 'P 31 2"', "P 32 1 2": "P 32 2 (0 0 2)", "P 32 2 1": 'P 32 2"',
    "H 3 2": 'R 3 2"', "R 3 2": "P 3* 2", "P 6": "P 6", "P 61": "P 61", "P 65": "P 65",
    "P 62": "P 62", "P 64": "P 64", "P 63": "P 6c", "P 6 2 2": "P 6 2",
    "P 61 2 2": "P 61 2 (0 0 5)", "P 65 2 2": "P 65 2 (0 0 1)",
    "P 62 2 2": "P 62 2 (0 0 4)", "P 64 2 2": "P 64 2 (0 0 2)", "P 63 2 2": "P 6c 2c",
    "P 2 3": "P 2 2 3", "F 2 3": "F 2 2 3", "I 2 3": "I 2 2 3", "P 21 3": "P 2ac 2ab 3",
    "I 21 3": "I 2b 2c 3", "P 4 3 2": "P 4 2 3", "P 42 3 2": "P 4n 2 3",
    "F 4 3 2": "F 4 2 3", "F 41 3 2": "F 4d 2 3", "I 4 3 2": "I 4 2 3",
    "P 43 3 2": "P 4acd 2ab 3", "P 41 3 2": "P 4bd 2ab 3", "I 41 3 2": "I 4bd 2c 3",
    # Short monoclinic symbols some programs write.
    "P 2": "P 2y", "P 21": "P 2yb", "C 2": "C 2y", "I 2": "I 2y",
}

_HALF, _THIRD, _QUARTER = Fraction(1, 2), Fraction(1, 3), Fraction(1, 4)
_CENTRING = {
    "P": [], "A": [(0, _HALF, _HALF)], "B": [(_HALF, 0, _HALF)], "C": [(_HALF, _HALF, 0)],
    "I": [(_HALF, _HALF, _HALF)], "R": [(2 * _THIRD, _THIRD, _THIRD), (_THIRD, 2 * _THIRD, 2 * _THIRD)],
    "F": [(0, _HALF, _HALF), (_HALF, 0, _HALF), (_HALF, _HALF, 0)],
}
_SHIFT = {"a": (_HALF, 0, 0), "b": (0, _HALF, 0), "c": (0, 0, _HALF), "n": (_HALF, _HALF, _HALF),
          "u": (_QUARTER, 0, 0), "v": (0, _QUARTER, 0), "w": (0, 0, _QUARTER),
          "d": (_QUARTER, _QUARTER, _QUARTER)}
_AXIS = {
    "x": {1: ((1, 0, 0), (0, 1, 0), (0, 0, 1)), 2: ((1, 0, 0), (0, -1, 0), (0, 0, -1)),
          3: ((1, 0, 0), (0, 0, -1), (0, 1, -1)), 4: ((1, 0, 0), (0, 0, -1), (0, 1, 0)),
          6: ((1, 0, 0), (0, 1, -1), (0, 1, 0))},
    "y": {1: ((1, 0, 0), (0, 1, 0), (0, 0, 1)), 2: ((-1, 0, 0), (0, 1, 0), (0, 0, -1)),
          3: ((-1, 0, 1), (0, 1, 0), (-1, 0, 0)), 4: ((0, 0, 1), (0, 1, 0), (-1, 0, 0)),
          6: ((0, 0, 1), (0, 1, 0), (-1, 0, 1))},
    "z": {1: ((1, 0, 0), (0, 1, 0), (0, 0, 1)), 2: ((-1, 0, 0), (0, -1, 0), (0, 0, 1)),
          3: ((0, -1, 0), (1, -1, 0), (0, 0, 1)), 4: ((0, -1, 0), (1, 0, 0), (0, 0, 1)),
          6: ((1, -1, 0), (1, 0, 0), (0, 0, 1))},
}
# Two-fold axes along a face diagonal, relative to the preceding axis.
_DIAGONAL = {
    ("z", "'"): ((0, -1, 0), (-1, 0, 0), (0, 0, -1)), ("z", '"'): ((0, 1, 0), (1, 0, 0), (0, 0, -1)),
    ("x", "'"): ((-1, 0, 0), (0, 0, -1), (0, -1, 0)), ("x", '"'): ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
    ("y", "'"): ((0, 0, -1), (0, -1, 0), (-1, 0, 0)), ("y", '"'): ((0, 0, 1), (0, -1, 0), (1, 0, 0)),
}
_BODY_DIAGONAL = ((0, 0, 1), (1, 0, 0), (0, 1, 0))
_UNIT = {"x": 0, "y": 1, "z": 2}


def _mod(t):
    return tuple(v % 1 for v in t)


def _apply(r, t):
    return tuple(sum(r[i][k] * t[k] for k in range(3)) for i in range(3))


def _compose(a, b):
    """The operator a after b."""
    (ra, ta), (rb, tb) = a, b
    r = tuple(tuple(sum(ra[i][k] * rb[k][j] for k in range(3)) for j in range(3))
              for i in range(3))
    return r, _mod(tuple(x + y for x, y in zip(_apply(ra, tb), ta)))


def _matrix(token, index, previous):
    """(operator, axis) of one Hall matrix symbol such as 2ac, -4bw or 3*."""
    improper = token.startswith("-")
    body = token.lstrip("-")
    order = int(body[0])
    if order not in (1, 2, 3, 4, 6):
        raise ValueError(f"unsupported Hall rotation {token!r}")
    rest = body[1:]
    axis = rest[0] if rest[:1] in ("x", "y", "z", "'", '"', "*") else None
    if axis:
        rest = rest[1:]
    if axis is None:
        if index == 0:
            axis = "z"
        elif index == 1 and order == 2:
            axis = "x" if previous[0] in (2, 4) else "'"
        elif index == 2 and order == 3:
            axis = "*"
        elif order != 1:
            raise ValueError(f"Hall rotation {token!r} needs an explicit axis")
    if axis in ("'", '"'):
        r = _DIAGONAL[(previous[1], axis)]
    elif axis == "*":
        r = _BODY_DIAGONAL
    else:
        r = _AXIS[axis or "z"][order]
    t = [Fraction(0)] * 3
    for symbol in rest:
        if symbol.isdigit():
            if axis not in _UNIT:
                raise ValueError(f"Hall screw {token!r} needs an x, y or z axis")
            t[_UNIT[axis]] += Fraction(int(symbol), order)
        elif symbol in _SHIFT:
            t = [a + b for a, b in zip(t, _SHIFT[symbol])]
        else:
            raise ValueError(f"unsupported Hall translation in {token!r}")
    if improper:
        r = tuple(tuple(-v for v in row) for row in r)
    return (r, _mod(t)), (order, axis if axis in _UNIT else previous[1])


def operators(hall):
    """The space group's operators in a fixed order: the identity first, then
    the point operators in generation order, each block repeated for every
    lattice centring vector."""
    text, _, basis = hall.partition("(")
    tokens = text.split()
    lattice = tokens[0]
    centric = lattice.startswith("-")
    if lattice.lstrip("-") not in _CENTRING:
        raise ValueError(f"unsupported Hall lattice {lattice!r}")
    centring = [(Fraction(0),) * 3] + [tuple(Fraction(v) for v in c)
                                       for c in _CENTRING[lattice.lstrip("-")]]
    generators, previous = [], (None, "z")
    for index, token in enumerate(tokens[1:]):
        op, previous = _matrix(token, index, previous)
        generators.append(op)
    identity = (((1, 0, 0), (0, 1, 0), (0, 0, 1)), (Fraction(0),) * 3)
    if centric:
        generators.append((((-1, 0, 0), (0, -1, 0), (0, 0, -1)), (Fraction(0),) * 3))
    if basis:
        # Change of basis by an origin shift v (in twelfths): (R, t + v - R v).
        v = tuple(Fraction(int(x), 12) for x in basis.rstrip(") ").split())
        generators = [(r, _mod(tuple(a + b - c for a, b, c in zip(t, v, _apply(r, v)))))
                      for r, t in generators]
    generators += [(identity[0], c) for c in centring[1:]]
    group, seen = [identity], {identity}
    for op in group:
        for g in generators:
            new = _compose(g, op)
            if new not in seen:
                seen.add(new); group.append(new)
                if len(group) > 192:
                    raise ValueError(f"Hall symbol {hall!r} does not close")
    points, covered = [], set()
    for r, t in group:
        if (r, t) not in covered:
            points.append((r, t))
            covered |= {(r, _mod(tuple(a + b for a, b in zip(t, c)))) for c in centring}
    return [(r, _mod(tuple(a + b for a, b in zip(t, c)))) for c in centring for r, t in points]


# Hall symbols of the 230 space groups in the settings of International
# Tables A (Hall 1981, and the 530-setting table of R. W. Grosse-Kunstleve,
# Acta Cryst. A55, 383-395, 1999). ITA_HALL[n - 1] is the standard setting of
# group n: unique axis b and cell choice 1 for monoclinic groups, origin
# choice 1 where ITA gives two origins, hexagonal axes for R groups. The
# dictionaries hold the alternatives that pw.x space_group input selects.
ITA_HALL = (
    'P 1', '-P 1', 'P 2y', 'P 2yb', 'C 2y', 'P -2y', 'P -2yc', 'C -2y', 'C -2yc', '-P 2y',
    '-P 2yb', '-C 2y', '-P 2yc', '-P 2ybc', '-C 2yc', 'P 2 2', 'P 2c 2', 'P 2 2ab',
    'P 2ac 2ab', 'C 2c 2', 'C 2 2', 'F 2 2', 'I 2 2', 'I 2b 2c', 'P 2 -2', 'P 2c -2',
    'P 2 -2c', 'P 2 -2a', 'P 2c -2ac', 'P 2 -2bc', 'P 2ac -2', 'P 2 -2ab', 'P 2c -2n',
    'P 2 -2n', 'C 2 -2', 'C 2c -2', 'C 2 -2c', 'A 2 -2', 'A 2 -2b', 'A 2 -2a', 'A 2 -2ab',
    'F 2 -2', 'F 2 -2d', 'I 2 -2', 'I 2 -2c', 'I 2 -2a', '-P 2 2', 'P 2 2 -1n', '-P 2 2c',
    'P 2 2 -1ab', '-P 2a 2a', '-P 2a 2bc', '-P 2ac 2', '-P 2a 2ac', '-P 2 2ab', '-P 2ab 2ac',
    '-P 2c 2b', '-P 2 2n', 'P 2 2ab -1ab', '-P 2n 2ab', '-P 2ac 2ab', '-P 2ac 2n', '-C 2c 2',
    '-C 2ac 2', '-C 2 2', '-C 2 2c', '-C 2a 2', 'C 2 2 -1ac', '-F 2 2', 'F 2 2 -1d', '-I 2 2',
    '-I 2 2c', '-I 2b 2c', '-I 2b 2', 'P 4', 'P 4w', 'P 4c', 'P 4cw', 'I 4', 'I 4bw', 'P -4',
    'I -4', '-P 4', '-P 4c', 'P 4ab -1ab', 'P 4n -1n', '-I 4', 'I 4bw -1bw', 'P 4 2',
    'P 4ab 2ab', 'P 4w 2c', 'P 4abw 2nw', 'P 4c 2', 'P 4n 2n', 'P 4cw 2c', 'P 4nw 2abw',
    'I 4 2', 'I 4bw 2bw', 'P 4 -2', 'P 4 -2ab', 'P 4c -2c', 'P 4n -2n', 'P 4 -2c', 'P 4 -2n',
    'P 4c -2', 'P 4c -2ab', 'I 4 -2', 'I 4 -2c', 'I 4bw -2', 'I 4bw -2c', 'P -4 2', 'P -4 2c',
    'P -4 2ab', 'P -4 2n', 'P -4 -2', 'P -4 -2c', 'P -4 -2ab', 'P -4 -2n', 'I -4 -2',
    'I -4 -2c', 'I -4 2', 'I -4 2bw', '-P 4 2', '-P 4 2c', 'P 4 2 -1ab', 'P 4 2 -1n',
    '-P 4 2ab', '-P 4 2n', 'P 4ab 2ab -1ab', 'P 4ab 2n -1ab', '-P 4c 2', '-P 4c 2c',
    'P 4n 2c -1n', 'P 4n 2 -1n', '-P 4c 2ab', '-P 4n 2n', 'P 4n 2n -1n', 'P 4n 2ab -1n',
    '-I 4 2', '-I 4 2c', 'I 4bw 2bw -1bw', 'I 4bw 2aw -1bw', 'P 3', 'P 31', 'P 32', 'R 3',
    '-P 3', '-R 3', 'P 3 2', 'P 3 2"', 'P 31 2 (0 0 4)', 'P 31 2"', 'P 32 2 (0 0 2)',
    'P 32 2"', 'R 3 2"', 'P 3 -2"', 'P 3 -2', 'P 3 -2"c', 'P 3 -2c', 'R 3 -2"', 'R 3 -2"c',
    '-P 3 2', '-P 3 2c', '-P 3 2"', '-P 3 2"c', '-R 3 2"', '-R 3 2"c', 'P 6', 'P 61', 'P 65',
    'P 62', 'P 64', 'P 6c', 'P -6', '-P 6', '-P 6c', 'P 6 2', 'P 61 2 (0 0 5)',
    'P 65 2 (0 0 1)', 'P 62 2 (0 0 4)', 'P 64 2 (0 0 2)', 'P 6c 2c', 'P 6 -2', 'P 6 -2c',
    'P 6c -2', 'P 6c -2c', 'P -6 2', 'P -6c 2', 'P -6 -2', 'P -6c -2c', '-P 6 2', '-P 6 2c',
    '-P 6c 2', '-P 6c 2c', 'P 2 2 3', 'F 2 2 3', 'I 2 2 3', 'P 2ac 2ab 3', 'I 2b 2c 3',
    '-P 2 2 3', 'P 2 2 3 -1n', '-F 2 2 3', 'F 2 2 3 -1d', '-I 2 2 3', '-P 2ac 2ab 3',
    '-I 2b 2c 3', 'P 4 2 3', 'P 4n 2 3', 'F 4 2 3', 'F 4d 2 3', 'I 4 2 3', 'P 4acd 2ab 3',
    'P 4bd 2ab 3', 'I 4bd 2c 3', 'P -4 2 3', 'F -4 2 3', 'I -4 2 3', 'P -4n 2 3', 'F -4a 2 3',
    'I -4bd 2c 3', '-P 4 2 3', 'P 4 2 3 -1n', '-P 4n 2 3', 'P 4n 2 3 -1n', '-F 4 2 3',
    '-F 4a 2 3', 'F 4d 2 3 -1d', 'F 4d 2 3 -1ad', '-I 4 2 3', '-I 4bd 2c 3',
)
ITA_HALL_ORIGIN_2 = {
    48: '-P 2ab 2bc', 50: '-P 2ab 2b', 59: '-P 2ab 2a', 68: '-C 2a 2ac', 70: '-F 2uv 2vw',
    85: '-P 4a', 86: '-P 4bc', 88: '-I 4ad', 125: '-P 4a 2b', 126: '-P 4a 2bc',
    129: '-P 4a 2a', 130: '-P 4a 2ac', 133: '-P 4ac 2b', 134: '-P 4ac 2bc', 137: '-P 4ac 2a',
    138: '-P 4ac 2ac', 141: '-I 4bd 2', 142: '-I 4bd 2c', 201: '-P 2ab 2bc 3',
    203: '-F 2uv 2vw 3', 222: '-P 4a 2bc 3', 224: '-P 4bc 2bc 3', 227: '-F 4vw 2vw 3',
    228: '-F 4ud 2vw 3',
}
ITA_HALL_UNIQUE_C = {
    3: 'P 2', 4: 'P 2c', 5: 'B 2', 6: 'P -2', 7: 'P -2a', 8: 'B -2', 9: 'B -2ab', 10: '-P 2',
    11: '-P 2c', 12: '-B 2', 13: '-P 2a', 14: '-P 2ac', 15: '-B 2ab',
}
ITA_HALL_RHOMBOHEDRAL = {
    146: 'P 3*', 148: '-P 3*', 155: 'P 3* 2', 160: 'P 3* -2', 161: 'P 3* -2n', 166: '-P 3* 2',
    167: '-P 3* 2n',
}


def ita_hall(number, *, unique_axis="b", origin_choice=1, rhombohedral=False):
    """The Hall symbol of space group ``number`` (1-230) in the setting that
    ``unique_axis`` ("b" or "c", monoclinic groups only), ``origin_choice``
    (1 or 2, groups with two ITA origins only) and ``rhombohedral`` (R groups
    only) select. For unique axis c, centred groups use the B-centred cell
    (ITA cell choice 2) and primitive groups cell choice 1."""
    if type(number) is not int or not 1 <= number <= 230:
        raise ValueError("space-group number must be an integer from 1 to 230")
    if unique_axis == "c" and number in ITA_HALL_UNIQUE_C:
        return ITA_HALL_UNIQUE_C[number]
    if origin_choice == 2 and number in ITA_HALL_ORIGIN_2:
        return ITA_HALL_ORIGIN_2[number]
    if rhombohedral and number in ITA_HALL_RHOMBOHEDRAL:
        return ITA_HALL_RHOMBOHEDRAL[number]
    return ITA_HALL[number - 1]


def orbit(fractional, ops, lattice, tolerance):
    """(images, merged): the distinct images in [0, 1) of a fractional
    position under ``ops`` with lattice rows ``lattice``, and how many images
    fell within ``tolerance`` (Cartesian, the unit of ``lattice``) of an
    earlier one, wrapping across cell boundaries."""
    def apart(f, g):
        d = [a - b - round(a - b) for a, b in zip(f, g)]
        x = [sum(d[i] * lattice[i][k] for i in range(3)) for k in range(3)]
        return sum(v * v for v in x) ** 0.5

    images, merged = [], 0
    for r, t in ops:
        g = [round(sum(r[i][k] * fractional[k] for k in range(3)) + float(t[i]), 10) % 1.0
             for i in range(3)]
        if any(apart(g, h) < tolerance for h in images):
            merged += 1
        else:
            images.append(g)
    return images, merged
