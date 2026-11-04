"""QVF archive inspection, validation and explicit-output maintenance."""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import sys
import zipfile

from . import __version__, validate
from .reader import Reader
from . import operations


def _json(value):
    print(json.dumps(value, ensure_ascii=False, indent=2))


def _parser():
    parser = argparse.ArgumentParser(prog="qvf", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("validate", help="validate cumulatively through a selected stage")
    command.add_argument("archive")
    command.add_argument("--level", choices=[name.lower() for name in validate.Level.__members__], default="all")
    command.add_argument("--json", action="store_true")
    for name, help_text in (("info", "show sections, kinds and member sizes"),
                            ("cat", "print the manifest or one JSON member")):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("archive")
        if name == "info":
            command.add_argument("--json", action="store_true")
        else:
            command.add_argument("section", nargs="?")
            command.add_argument("member", nargs="?")
    command = commands.add_parser("extract", help="verify and extract one declared member")
    command.add_argument("archive")
    command.add_argument("section")
    command.add_argument("member")
    command.add_argument("-o", "--output", required=True, help="new output file (never overwritten)")
    command = commands.add_parser("diff", help="compare decoded content, independent of ZIP layout")
    command.add_argument("left")
    command.add_argument("right")
    command.add_argument("--atol", type=float, default=0.0)
    command.add_argument("--rtol", type=float, default=0.0)
    command.add_argument("--ignore-metadata", action="store_true", help="ignore root metadata (UUID/timestamps/authors)")
    command.add_argument("--max-differences", type=int, default=20)
    command.add_argument("--json", action="store_true")
    command = commands.add_parser("repack", help="copy current members to a new archive, dropping orphans")
    command.add_argument("archive")
    command.add_argument("-o", "--output", required=True)
    command.add_argument("--compression", choices=("preserve", "store", "deflate"), default="preserve")
    command = commands.add_parser("scrub", help="report and remove selected personal content into a new archive")
    command.add_argument("archive")
    command.add_argument("-o", "--output")
    command.add_argument("--dry-run", action="store_true", help="report changes without writing")
    command.add_argument("--path-pattern", action="append", default=[], metavar="REGEX")
    command.add_argument("--user-pattern", action="append", default=[], metavar="REGEX")
    command.add_argument("--replacement", default="[redacted]", help="literal replacement for every pattern match")
    command = commands.add_parser("conformance", help="run an explicit positive corpus and optional siblings")
    command.add_argument("corpus")
    command.add_argument("--negative", metavar="DIR")
    command.add_argument("--no-negative", action="store_true")
    command.add_argument("--legacy", metavar="DIR")
    command.add_argument("--no-legacy", action="store_true")
    command = commands.add_parser("convert", help="convert XYZ/extXYZ, PDB, XSF, scalar cube or Molden through QVF")
    command.add_argument("input")
    command.add_argument("-o", "--output", required=True)
    command.add_argument("--from", dest="input_format", choices=("qvf", "xyz", "extxyz", "cube", "xsf", "pdb", "cif", "mmcif", "poscar", "chgcar", "molden"), required=True)
    command.add_argument("--to", dest="output_format", choices=("qvf", "xyz", "extxyz", "cube", "xsf", "pdb", "cif", "mmcif", "poscar", "chgcar", "molden"), required=True)
    command.add_argument("--section", help="QVF section to export")
    command.add_argument("--block", help="CIF/mmCIF import: data block name (without data_)")
    command.add_argument("--species", nargs="+", metavar="ELEMENT",
                         help="POSCAR/CHGCAR import: elements of a VASP 4 file, one per count")
    command.add_argument("--expand-ncs", action="store_true",
                         help="PDB import: also generate MTRIX copies the file does not give")
    command.add_argument("--expand-symmetry", action="store_true",
                         help="PDB import: also write the unit cell with every space-group copy")
    command.add_argument("--kind", help="scalar volume kind asserted on cube/XSF import (default volume.generic)")
    command.add_argument("--periodic", nargs=3, choices=("true", "false"), help="cube/XSF scalar grid periodic axes")
    command.add_argument("--no-endpoint", action="store_true", help="cube import spans n voxels rather than n-1")
    command.add_argument("--spin", choices=("restricted", "restricted_open", "unrestricted"),
                         help="required spin interpretation for Molden import")
    command.add_argument("--molden-dialect", default=None,
                         choices=("standard", "renormalize", "orca", "turbomole", "cfour", "psi4-1.3", "xtb"),
                         help="Molden import: producer convention (default standard)")
    command = commands.add_parser("prepare-job", help="preserve one geometry/drawing in an explicit pending job")
    command.add_argument("archive")
    command.add_argument("-o", "--output", required=True, help="new archive (never overwritten)")
    command.add_argument("--job-type", choices=("molecular", "periodic"), required=True)
    command.add_argument("--method", required=True)
    command.add_argument("--task", choices=("single_point", "optimize", "hessian"), required=True)
    command.add_argument("--charge", type=int, default=0)
    command.add_argument("--multiplicity", type=int, default=1)
    command.add_argument("--kpoints", nargs=3, type=int, metavar=("N1", "N2", "N3"))
    for name in ("basis", "functional"):
        group = command.add_mutually_exclusive_group()
        group.add_argument("--" + name, default=argparse.SUPPRESS)
        group.add_argument("--no-" + name, dest=name, action="store_const", const=None,
                           default=argparse.SUPPRESS, help=f"write explicit null for {name}")
    command = commands.add_parser("export-xtb", help="package molecular xTB outputs as a new QVF")
    command.add_argument("results", help="xtbout.json from the selected calculation")
    command.add_argument("--geometry", required=True, help="matching final XYZ or molecular coord")
    command.add_argument("-o", "--output", required=True, help="new archive (never overwritten)")
    for option in ("trajectory", "wbo", "molden", "modes"):
        command.add_argument("--" + option, help="explicit companion file from the same calculation")
    command.add_argument("--spin", choices=("restricted", "restricted_open", "unrestricted"),
                         help="required with --molden")
    command.add_argument("--cube", nargs=2, action="append", default=[], metavar=("FILE", "KIND"),
                         help="scalar cube and explicit QVF quantity, e.g. density.cub volume.density")
    command.add_argument("--charge", type=int)
    command.add_argument("--multiplicity", type=int)
    command.add_argument("--run-status", choices=("converged", "failed"),
                         help="status established independently; omitted by default")
    command.add_argument("--canonical-orbitals", action="store_true",
                         help="assert that the Molden orbitals are untransformed canonical xTB orbitals")
    command = commands.add_parser("run-xtb", help="run an xTB bundle once and export a converged QVF")
    command.add_argument("request", help="version-1 xtb-job.json beside coord and xcontrol")
    command.add_argument("--binary", required=True, help="absolute registered xTB executable")
    commands.add_parser("schema", help="print the packaged normative JSON Schema")
    return parser


def _execute(args):
    if args.command == "schema":
        print(Path(validate._SCHEMA_PATH).read_text(encoding="utf-8"), end="")
    elif args.command == "validate":
        report = validate.run(args.archive, level=args.level)
        if args.json:
            _json(report)
        else:
            print(f"{args.level}: {'OK' if report['ok'] else 'FAILED'} ({report['n_sections']} sections)")
            print(validate._context_summary(report))
            for finding in report["findings"]:
                print(f"{finding['severity']} {finding['code']} {finding['where']}: {finding['message']}")
        return 0 if report["ok"] else 1
    elif args.command in ("info", "cat"):
        with Reader(args.archive) as reader:
            if args.command == "cat":
                if (args.section is None) != (args.member is None):
                    raise ValueError("cat requires both SECTION and MEMBER, or neither for the manifest")
                _json(reader.raw_manifest if args.section is None else reader.json(args.section, args.member))
            else:
                sections = []
                for section in reader.sections:
                    members = []
                    for role, member in (section.members or {}).items():
                        info = reader._zip.getinfo(member.path)
                        members.append(dict(member.to_json(), role=role, size=info.file_size,
                                            compressed_size=info.compress_size))
                    sections.append({"id": section.id, "kind": section.kind, "members": members})
                value = {"qvf_version": reader.raw_manifest.get("qvf_version"),
                         "source": reader.raw_manifest.get("source"), "sections": sections}
                if args.json:
                    _json(value)
                else:
                    print(f"QVF {value['qvf_version']} | {len(sections)} sections | {value['source']}")
                    for section in sections:
                        print(f"  {section['id']} [{section['kind']}]")
                        for member in section["members"]:
                            tensor = f" {member['dtype']}{member.get('shape', [])}" if "dtype" in member else ""
                            print(f"    {member['role']}: {member['path']} ({member['format']}{tensor}, {member['size']} bytes)")
    elif args.command == "extract":
        operations.extract(args.archive, args.section, args.member, args.output)
        print(f"Extracted {args.section}/{args.member} to {args.output}")
    elif args.command == "diff":
        result = operations.diff(args.left, args.right, atol=args.atol, rtol=args.rtol,
                                 ignore_metadata=args.ignore_metadata, max_differences=args.max_differences)
        if args.json:
            _json(result)
        else:
            print("Equal" if result["equal"] else f"Different ({result['difference_count']} differences)")
            for item in result["differences"]:
                print(f"  {item['where']}: {item['message']}")
            if result["truncated"]:
                print("  additional differences omitted; raise --max-differences to display more")
        return 0 if result["equal"] else 1
    elif args.command == "repack":
        methods = {"preserve": None, "store": zipfile.ZIP_STORED, "deflate": zipfile.ZIP_DEFLATED}
        operations.repack(args.archive, args.output, compression=methods[args.compression])
        print(f"Repacked to {args.output}")
    elif args.command == "scrub":
        if args.dry_run and args.output:
            raise ValueError("--dry-run and --output are mutually exclusive")
        _json(operations.scrub(args.archive, args.output, dry_run=args.dry_run,
                              patterns=args.path_pattern + args.user_pattern, replacement=args.replacement))
    elif args.command == "convert":
        from . import convert
        grid_options = args.kind is not None or args.periodic is not None
        if args.input_format not in ("cube", "xsf") and grid_options:
            raise ValueError("--kind and --periodic apply only to cube/XSF import")
        if args.input_format != "cube" and args.no_endpoint:
            raise ValueError("--no-endpoint applies only to cube import")
        if args.input_format not in ("cif", "mmcif") and args.block is not None:
            raise ValueError("--block applies only to CIF/mmCIF import")
        if args.input_format not in ("poscar", "chgcar") and args.species is not None:
            raise ValueError("--species applies only to POSCAR/CHGCAR import")
        if args.input_format != "pdb" and (args.expand_symmetry or args.expand_ncs):
            raise ValueError("--expand-symmetry and --expand-ncs apply only to PDB import")
        if args.input_format != "molden" and args.spin is not None:
            raise ValueError("--spin applies only to Molden import")
        if args.input_format != "molden" and args.molden_dialect is not None:
            raise ValueError("--molden-dialect applies only to Molden import")
        if args.input_format == "qvf" and args.output_format in ("xyz", "extxyz", "cube", "xsf", "pdb", "cif", "mmcif", "poscar", "chgcar", "molden"):
            if args.output_format == "cube":
                result = convert.export_cube(args.input, args.output, section=args.section)
            elif args.output_format == "xsf":
                result = convert.export_xsf(args.input, args.output, section=args.section)
            elif args.output_format == "pdb":
                result = convert.export_pdb(args.input, args.output, section=args.section)
            elif args.output_format == "mmcif":
                result = convert.export_mmcif(args.input, args.output, section=args.section)
            elif args.output_format == "chgcar":
                result = convert.export_chgcar(args.input, args.output, section=args.section)
            elif args.output_format == "poscar":
                result = convert.export_poscar(args.input, args.output, section=args.section)
            elif args.output_format == "cif":
                result = convert.export_cif(args.input, args.output, section=args.section)
            elif args.output_format == "molden":
                result = convert.export_molden(args.input, args.output, section=args.section)
            else:
                result = convert.export_xyz(args.input, args.output, section=args.section,
                                            extended=args.output_format == "extxyz")
        elif args.output_format == "qvf" and args.input_format != "qvf":
            if args.section is not None:
                raise ValueError("--section applies only to QVF export")
            if args.input_format == "cube":
                result = convert.import_cube(args.input, args.output, kind=args.kind or "volume.generic",
                    periodic=tuple(v == "true" for v in args.periodic) if args.periodic else (False,) * 3,
                    includes_endpoint=not args.no_endpoint)
            elif args.input_format == "xsf":
                result = convert.import_xsf(args.input, args.output, kind=args.kind,
                    periodic=tuple(v == "true" for v in args.periodic) if args.periodic else None)
            elif args.input_format == "pdb":
                result = convert.import_pdb(args.input, args.output,
                                            expand_symmetry=args.expand_symmetry,
                                            expand_ncs=args.expand_ncs)
            elif args.input_format == "mmcif":
                result = convert.import_mmcif(args.input, args.output, block=args.block)
            elif args.input_format == "chgcar":
                result = convert.import_chgcar(args.input, args.output, species=args.species)
            elif args.input_format == "poscar":
                result = convert.import_poscar(args.input, args.output, species=args.species)
            elif args.input_format == "cif":
                result = convert.import_cif(args.input, args.output, block=args.block)
            elif args.input_format == "molden":
                result = convert.import_molden(args.input, args.output, spin=args.spin,
                                               dialect=args.molden_dialect or "standard")
            else:
                result = convert.import_xyz(args.input, args.output)
        else:
            raise ValueError("conversion must have QVF as exactly one endpoint")
        _json(result)
    elif args.command == "export-xtb":
        from .adapters.xtb import write_xtb
        _json(write_xtb(args.results, args.geometry, args.output,
                        trajectory=args.trajectory, wbo=args.wbo, molden=args.molden,
                        spin=args.spin, modes=args.modes, cubes=args.cube,
                        charge=args.charge, multiplicity=args.multiplicity,
                        run_status=args.run_status, canonical_orbitals=args.canonical_orbitals))
    elif args.command == "run-xtb":
        from .adapters.xtb_run import run_xtb
        _json(run_xtb(args.request, binary=args.binary))
    elif args.command == "prepare-job":
        from .jobs import prepare_job
        names = {name: getattr(args, name) for name in ("basis", "functional") if hasattr(args, name)}
        _json(prepare_job(args.archive, args.output, job_type=args.job_type,
                          method=args.method, task=args.task, charge=args.charge,
                          multiplicity=args.multiplicity, kpoints=args.kpoints, **names))
    elif args.command == "conformance":
        from .conformance import run
        return run(args.corpus, args.negative, negative=not args.no_negative,
                   legacy_dir=args.legacy, legacy=False if args.no_legacy else (True if args.legacy else None))
    return 0


def main(argv=None):
    """Console entry point: 0 success/equal, 1 failed/different, 2 usage/I/O."""
    args = _parser().parse_args(argv)
    try:
        return _execute(args)
    except (OSError, ValueError, TypeError, KeyError, zipfile.BadZipFile, RuntimeError, re.error) as exc:
        print(f"qvf: {exc}", file=sys.stderr)
        return 2
