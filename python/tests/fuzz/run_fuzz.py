"""Time-budgeted Hypothesis fuzz campaign over one target.

    python python/tests/fuzz/run_fuzz.py --target mixed --seconds 1200

Runs Hypothesis in rounds until the budget is spent, checks every input
against the invariants of harness.py, and saves the first input of each new
problem signature to OUT/<target>/ (default build/fuzz-py/, ignored by Git).
Prints the number of inputs, the slowest input and the largest peak
allocation, and exits non-zero when any problem was found.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from hypothesis import HealthCheck, given, settings  # noqa: E402

import harness  # noqa: E402
import strategies  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--target", choices=sorted(strategies.TARGETS), default="mixed")
    ap.add_argument("--seconds", type=float, default=60)
    ap.add_argument("--round", type=int, default=500, help="inputs per Hypothesis round")
    ap.add_argument("--out", type=Path, default=harness.HERE.parents[2] / "build" / "fuzz-py")
    ap.add_argument("--no-reader", action="store_true", help="validator only")
    args = ap.parse_args(argv)
    warnings.simplefilter("ignore")
    out = args.out / args.target
    out.mkdir(parents=True, exist_ok=True)
    limits = harness.DEFAULT if args.target == "bomb" else harness.SMALL
    stats = {"inputs": 0, "slowest": 0.0, "peak": 0, "signatures": {}}
    deadline = time.monotonic() + args.seconds

    @settings(max_examples=args.round, deadline=None, database=None,
              suppress_health_check=list(HealthCheck))
    @given(strategies.inputs(args.target))
    def one(case):
        if time.monotonic() > deadline:
            return
        name, plan, byte_seed, raw = case
        problems, elapsed, peak = harness.check(raw, limits=limits,
                                                reader=not args.no_reader)
        stats["inputs"] += 1
        stats["slowest"] = max(stats["slowest"], elapsed)
        stats["peak"] = max(stats["peak"], peak)
        for p in problems:
            key = p.split(":")[0] + ":" + p.split(":", 2)[1][:80] if ":" in p else p
            if key in stats["signatures"]:
                stats["signatures"][key]["count"] += 1
                continue
            digest = hashlib.sha256(raw).hexdigest()[:16]
            (out / f"{digest}.qvf").write_bytes(raw)
            stats["signatures"][key] = {"count": 1, "example": f"{digest}.qvf",
                                        "seed": name, "plan": plan,
                                        "byte_seed": byte_seed, "problem": p}
            print(f"NEW {p}\n    seed={name} plan={plan} byte_seed={byte_seed} "
                  f"-> {digest}.qvf", flush=True)

    t0 = time.monotonic()
    rounds = 0
    harness_errors = 0

    def summary():
        # Written after every round, so an interrupted campaign keeps its numbers.
        s = {"target": args.target, "seconds": round(time.monotonic() - t0, 1),
             "rounds": rounds, "harness_errors": harness_errors,
             "inputs": stats["inputs"], "slowest_s": round(stats["slowest"], 3),
             "peak_bytes": stats["peak"], "problems": stats["signatures"]}
        (out / "summary.json").write_text(json.dumps(s, indent=2))
        return s

    while time.monotonic() < deadline:
        try:
            one()
        except Exception as exc:  # noqa: BLE001 - a harness bug must not end the campaign
            harness_errors += 1
            print(f"harness error in round {rounds}: {exc!r:.300}", flush=True)
        rounds += 1
        summary()
    summary = summary()
    print(json.dumps({k: v for k, v in summary.items() if k != "problems"}))
    print(f"problem signatures: {len(stats['signatures'])}")
    return 1 if stats["signatures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
