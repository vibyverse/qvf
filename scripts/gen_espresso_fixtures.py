#!/usr/bin/env python3
"""Write the synthetic pp.x fixtures of python/tests/fixtures/espresso.

Each case is one analytic periodic field on an FFT grid, written twice: as
the pp.x filplot intermediate file (the records of plot_io in QE
Modules/plot_io.f90) and as the cube file pp.x writes for the same data with
output_format = 6 (no endpoint plane, bohr, z fastest, e13.5 values). The
script uses the standard library only and does not import qvf, so the
fixtures are independent of the reader they test.

    python3 scripts/gen_espresso_fixtures.py          # rewrite the fixtures
    python3 scripts/gen_espresso_fixtures.py --check  # fail if they are stale
"""
import argparse
import math
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "python" / "tests" / "fixtures" / "espresso"


def fortran_e(value, width=17, digits=9):
    """Fortran 1pEw.d: the E is dropped from a three-digit exponent."""
    text = f"{value:.{digits}E}"
    mantissa, exponent = text.split("E")
    if abs(int(exponent)) > 99:
        text = mantissa + exponent
    return text.rjust(width)


def lattice_alat(ibrav, at):
    """Lattice rows in units of alat (INPUT_PW: ibrav 2 fcc, 3 bcc, 0 given)."""
    if ibrav == 2:
        return [[-0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [-0.5, 0.5, 0.0]]
    if ibrav == 3:
        return [[0.5, 0.5, 0.5], [-0.5, 0.5, 0.5], [-0.5, -0.5, 0.5]]
    return at


CASES = {
    # Si, fcc, charge density: 8 electrons in the cell.
    "si_rho": dict(title="synthetic Si charge density", ibrav=2, alat=10.26, at=None,
                   nr=(6, 6, 6), nrx=(6, 6, 6), plot_num=0,
                   species=[("Si", 4.0, 14)], atoms=[(1, (0, 0, 0)), (1, (0.25, 0.25, 0.25))],
                   mean=lambda volume: 8.0 / volume, amplitude=0.02, tiny=None),
    # GaN, hexagonal ibrav = 0, V_bare + V_H in Ry, padded nr1x and nr2x and
    # one value with a three-digit exponent.
    "gan_v": dict(title="synthetic GaN potential", ibrav=0, alat=6.026,
                  at=[[1.0, 0.0, 0.0], [-0.5, math.sqrt(3) / 2, 0.0], [0.0, 0.0, 1.626]],
                  nr=(5, 5, 8), nrx=(6, 7, 8), plot_num=11,
                  species=[("Ga", 13.0, 31), ("N", 5.0, 7)],
                  atoms=[(1, (1 / 3, 2 / 3, 0.0)), (1, (2 / 3, 1 / 3, 0.5)),
                         (2, (1 / 3, 2 / 3, 0.377)), (2, (2 / 3, 1 / 3, 0.877))],
                  mean=lambda volume: -0.35, amplitude=0.8, tiny=(2, 3, 5)),
    # Fe, bcc, spin polarization with zero mean.
    "fe_spin": dict(title="", ibrav=3, alat=5.42, at=None, nr=(4, 4, 4), nrx=(4, 4, 4),
                    plot_num=6, species=[("Fe1", 8.0, 26)], atoms=[(1, (0, 0, 0))],
                    mean=lambda volume: 0.0, amplitude=0.05, tiny=None),
}


def field(case, volume):
    nr = case["nr"]
    values = {}
    for x in range(nr[0]):
        for y in range(nr[1]):
            for z in range(nr[2]):
                u, v, w = x / nr[0], y / nr[1], z / nr[2]
                values[x, y, z] = case["mean"](volume) + case["amplitude"] * (
                    math.cos(2 * math.pi * u) * math.cos(2 * math.pi * v)
                    * math.cos(2 * math.pi * w) + 0.5 * math.sin(2 * math.pi * (u + 2 * w)))
    if case["tiny"]:
        values[case["tiny"]] = 1.25e-120
    return values


def render(case):
    rows = lattice_alat(case["ibrav"], case["at"])
    alat = case["alat"]
    a = [[c * alat for c in row] for row in rows]
    volume = abs(a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                 - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                 + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    values = field(case, volume)
    nr, nrx = case["nr"], case["nrx"]
    tau = [(t, [sum(f[i] * rows[i][k] for i in range(3)) for k in range(3)])
           for t, f in case["atoms"]]
    celldm = [alat, 0, 0, 0, 0, 0]
    out = [case["title"],
           "".join(f"{n:8d}" for n in (*nrx, *nr, len(tau), len(case["species"]))),
           f"{case['ibrav']:6d}" + "".join(f"{c:12.8f}" for c in celldm)]
    if case["ibrav"] == 0:
        out += ["  " + "  ".join(f"{c:.15f}" for c in row) for row in rows]
    out.append(f"{60.0:20.10f}{4.0:20.10f}{30.0:20.10f}{case['plot_num']:6d}")
    out += [f"{n:4d}   {label[:2]:<2s}   {zv:5.2f}"
            for n, (label, zv, _) in enumerate(case["species"], 1)]
    out += [f"{n:4d}   " + "".join(f"{c:15.9f}" for c in x) + f"   {t:2d}"
            for n, (t, x) in enumerate(tau, 1)]
    flat = [values.get((x, y, z), 0.0) for z in range(nr[2]) for y in range(nrx[1])
            for x in range(nrx[0])]
    for i in range(0, len(flat), 5):
        out.append("".join(fortran_e(v) for v in flat[i:i + 5]))
    filplot = "\n".join(out) + "\n"
    cube = [f" Cubefile created from PWScf calculation: {case['title']}", " Cube file",
            f"{len(tau):5d}{0.0:12.6f}{0.0:12.6f}{0.0:12.6f}"]
    cube += [f"{n:5d}" + "".join(f"{c / n:12.6f}" for c in row) for n, row in zip(nr, a)]
    cube += [f"{case['species'][t - 1][2]:5d}{case['species'][t - 1][1]:12.6f}"
             + "".join(f"{c * alat:12.6f}" for c in x) for t, x in tau]
    for x in range(nr[0]):
        for y in range(nr[1]):
            row = [values[x, y, z] for z in range(nr[2])]
            for i in range(0, len(row), 6):
                cube.append("".join(f"{v:13.5E}" for v in row[i:i + 6]))
    return filplot, "\n".join(cube) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="fail if a fixture is stale")
    args = parser.parse_args()
    stale = []
    for name, case in CASES.items():
        for suffix, text in zip(("filplot", "cube"), render(case)):
            path = OUT / f"{name}.{suffix}"
            if args.check:
                if not path.exists() or path.read_text() != text:
                    stale.append(path.name)
            else:
                OUT.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
    if stale:
        raise SystemExit("stale fixtures: " + ", ".join(stale))


if __name__ == "__main__":
    main()
