"""Generate the synthetic Gaussian fchk fixtures of test_convert_fchk.py.

Run from the repository root:

    python python/tests/fchk_fixtures/generate.py

The files are built from the documented fchk layout (Gaussian formchk
documentation, https://gaussian.com/formchk/): field headers as
``name (A40), 3 blanks, type, 3 blanks, N=count`` or a scalar, integer arrays
6I12 and real arrays 5E16.8 per line. No Gaussian output is used. Geometry,
basis and coefficients are synthetic. Contractions are normalized, and the SCF
densities are built from the stated occupations so that the files are
self-consistent. Deliberately independent of qvf._fchk: this writer and its
AO counts are written from the format description alone.
"""
from __future__ import annotations

import math
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent


def n_functions(shell_type):
    l = abs(shell_type)
    if shell_type == -1:
        return 4                      # SP: s, px, py, pz
    return 2 * l + 1 if shell_type < -1 else (l + 1) * (l + 2) // 2


def normalized(l, exponents, raw):
    """Scale contraction coefficients of normalized primitives to a unit AO."""
    total = sum(a * b * (2 * math.sqrt(x * y) / (x + y)) ** (l + 1.5)
                for x, a in zip(exponents, raw) for y, b in zip(exponents, raw))
    return [c / math.sqrt(total) for c in raw]


# Each fixture: atoms (Z, bohr), shells (type, atom 1-based, exponents, s or
# l coefficients, P(S=P) coefficients for SP), electrons and spin layout.
FIXTURES = {
    "rhf_pure": dict(
        title="synthetic RHF hydroxide, pure d/f/g and SP shells",
        method="RHF", basis="Gen",
        atoms=[(8, (0.0, 0.1, 0.2)), (1, (0.3, 1.4, -1.1))],
        shells=[(0, 1, [130.7, 23.8, 6.44], [0.154, 0.535, 0.444], None),
                (-1, 1, [5.03, 1.17], [-0.0999, 0.3995], [0.1559, 0.6077]),
                (-2, 1, [1.1], [1.0], None),
                (-3, 1, [0.9], [1.0], None),
                (-4, 1, [0.7], [1.0], None),
                (0, 2, [3.42, 0.62], [0.15, 0.53], None),
                (1, 2, [0.8], [1.0], None)],
        charge=-1, multiplicity=1, n_alpha=5, n_beta=5, unrestricted=False,
        energy=-75.123456789012, mulliken=[-1.1234, 0.1234], flags=(0, 0)),
    "uhf_cart": dict(
        title="synthetic UHF CH doublet, Cartesian d/f/g and SP shells",
        method="UHF", basis="Gen",
        atoms=[(6, (0.0, 0.0, 0.0)), (1, (0.4, -0.7, 2.0))],
        shells=[(0, 1, [71.6, 13.0, 3.53], [0.154, 0.535, 0.444], None),
                (-1, 1, [2.94, 0.68], [-0.1, 0.4], [0.156, 0.608]),
                (2, 1, [0.8, 0.25], [0.6, 0.5], None),
                (3, 1, [0.6], [1.0], None),
                (4, 1, [0.5], [1.0], None),
                (0, 2, [1.3], [1.0], None),
                (1, 2, [0.75], [1.0], None)],
        charge=0, multiplicity=2, n_alpha=4, n_beta=3, unrestricted=True,
        energy=-38.2701234567, mulliken=[-0.2, 0.2], flags=(1, 1)),
    "rohf_5d10f": dict(
        title="synthetic ROHF NH triplet, pure d with Cartesian f",
        method="ROHF", basis="Gen",
        atoms=[(7, (0.2, -0.1, 0.0)), (1, (-0.5, 0.3, 1.9))],
        shells=[(0, 1, [99.1, 18.1, 4.9], [0.154, 0.535, 0.444], None),
                (-1, 1, [3.78, 0.88], [-0.1, 0.4], [0.156, 0.608]),
                (-2, 1, [0.9], [1.0], None),
                (3, 1, [0.7], [1.0], None),
                (0, 2, [1.1], [1.0], None)],
        charge=0, multiplicity=3, n_alpha=5, n_beta=3, unrestricted=False,
        energy=None, mulliken=None, flags=(0, 1)),
}


def coefficients(seed, n_mo, n_ao):
    """Diagonally dominant synthetic MO coefficients, MO-major (row = MO)."""
    rng = random.Random(seed)
    return [[(1.0 if i == j else 0.0) + rng.uniform(-0.3, 0.3) for j in range(n_ao)] for i in range(n_mo)]


def packed_density(rows, occupations):
    n = len(rows[0])
    return [sum(o * row[i] * row[j] for o, row in zip(occupations, rows)) for i in range(n) for j in range(i + 1)]


def header(name, kind, count=None, value=None):
    if count is not None:
        return f"{name:<40}   {kind}   N={count:12d}"
    return f"{name:<40}   {kind}     {value:12d}" if kind == "I" else f"{name:<40}   {kind}     {value:22.15E}"


def integers(name, values):
    out = [header(name, "I", len(values))]
    out += ["".join(f"{v:12d}" for v in values[k:k + 6]) for k in range(0, len(values), 6)]
    return out


def reals(name, values):
    out = [header(name, "R", len(values))]
    out += ["".join(f"{v:16.8E}" for v in values[k:k + 5]) for k in range(0, len(values), 5)]
    return out


def resolved(spec):
    """Shells with normalized coefficients and the AO count."""
    shells = []
    for shell_type, atom, exponents, first, second in spec["shells"]:
        l = 0 if shell_type == -1 else abs(shell_type)
        first = normalized(l, exponents, first)
        second = normalized(1, exponents, second) if second is not None else None
        shells.append((shell_type, atom, exponents, first, second))
    return shells, sum(n_functions(s[0]) for s in shells)


def orbitals(name):
    """The fixture's MO coefficient rows (file AO order) per spin channel."""
    spec = FIXTURES[name]; shells, n_ao = resolved(spec)
    channels = {"alpha": coefficients(len(name), n_ao, n_ao)}
    if spec["unrestricted"]:
        channels["beta"] = coefficients(len(name) + 100, n_ao, n_ao)
    return shells, channels


def build(name):
    spec = FIXTURES[name]; shells, channels = orbitals(name)
    n_ao = len(channels["alpha"][0]); n_mo = n_ao
    atoms = spec["atoms"]; n_e = spec["n_alpha"] + spec["n_beta"]
    lines = [spec["title"], f"{'SP':<10}{spec['method']:<30}{spec['basis']:<30}".rstrip()]
    lines.append(header("Number of atoms", "I", value=len(atoms)))
    lines += integers("Info1-9", [9, 9, 0, 0, 0, 110, 1, 18, -502])
    lines.append(header("Route", "C", 2))
    lines.append("#P Gen GFIn" + " " * 13)
    for key, value in (("Charge", spec["charge"]), ("Multiplicity", spec["multiplicity"]),
                       ("Number of electrons", n_e), ("Number of alpha electrons", spec["n_alpha"]),
                       ("Number of beta electrons", spec["n_beta"]),
                       ("Number of basis functions", n_ao), ("Number of independent functions", n_mo)):
        lines.append(header(key, "I", value=value))
    lines += integers("Atomic numbers", [z for z, _ in atoms])
    lines += reals("Nuclear charges", [float(z) for z, _ in atoms])
    lines += reals("Current cartesian coordinates", [v for _, xyz in atoms for v in xyz])
    lines.append(header("Number of contracted shells", "I", value=len(shells)))
    lines.append(header("Number of primitive shells", "I", value=sum(len(s[2]) for s in shells)))
    lines.append(header("Pure/Cartesian d shells", "I", value=spec["flags"][0]))
    lines.append(header("Pure/Cartesian f shells", "I", value=spec["flags"][1]))
    lines.append(header("Highest angular momentum", "I", value=max(abs(s[0]) for s in shells)))
    lines.append(header("Largest degree of contraction", "I", value=max(len(s[2]) for s in shells)))
    lines += integers("Shell types", [s[0] for s in shells])
    lines += integers("Number of primitives per shell", [len(s[2]) for s in shells])
    lines += integers("Shell to atom map", [s[1] for s in shells])
    lines += reals("Primitive exponents", [e for s in shells for e in s[2]])
    lines += reals("Contraction coefficients", [c for s in shells for c in s[3]])
    if any(s[4] for s in shells):
        lines += reals("P(S=P) Contraction coefficients",
                       [c for s in shells for c in (s[4] or [0.0] * len(s[2]))])
    lines += reals("Coordinates of each shell", [v for s in shells for v in atoms[s[1] - 1][1]])
    if spec["energy"] is not None:
        lines.append(header("Total Energy", "R", value=spec["energy"]))
    occupied = {}
    if spec["unrestricted"]:
        occupied["alpha"] = [1.0] * spec["n_alpha"]; occupied["beta"] = [1.0] * spec["n_beta"]
    else:
        occupied["alpha"] = [2.0] * spec["n_beta"] + [1.0] * (spec["n_alpha"] - spec["n_beta"])
    for channel, rows in channels.items():
        prefix = channel.capitalize()
        lines += reals(prefix + " Orbital Energies", [-1.0 + 0.1 * k for k in range(n_mo)])
        lines += reals(prefix + " MO coefficients", [v for row in rows for v in row])
    if spec["unrestricted"]:
        pa = packed_density(channels["alpha"], occupied["alpha"])
        pb = packed_density(channels["beta"], occupied["beta"])
        total = [a + b for a, b in zip(pa, pb)]; spin = [a - b for a, b in zip(pa, pb)]
    else:
        total = packed_density(channels["alpha"], occupied["alpha"])
        singly = spec["n_alpha"] - spec["n_beta"]
        spin = packed_density(channels["alpha"][spec["n_beta"]:spec["n_alpha"]], [1.0] * singly) if singly else None
    lines += reals("Total SCF Density", total)
    if spin is not None:
        lines += reals("Spin SCF Density", spin)
    if spec["mulliken"] is not None:
        lines += reals("Mulliken Charges", spec["mulliken"])
    return "\n".join(lines) + "\n"


def main(directory=HERE):
    for name in FIXTURES:
        (Path(directory) / f"{name}.fchk").write_text(build(name), encoding="ascii")


if __name__ == "__main__":
    main()
