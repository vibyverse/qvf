"""Coverage-guided fuzzing of the Python validator and readers with atheris.

    mkdir -p build/fuzz-py/atheris
    python python/tests/fuzz/atheris_fuzz.py -max_total_time=600 -timeout=10 \\
        -rss_limit_mb=2048 build/fuzz-py/atheris conformance/corpus conformance/negative

libFuzzer adds new inputs to the first directory, so make that a scratch one.
Half of the mutations are libFuzzer's byte-level ones, half the
structure-aware operators of mutators.py, chosen from libFuzzer's seed. An
input that breaks an invariant of harness.py raises, which libFuzzer reports
as a crash and saves.
"""
from __future__ import annotations

import random
import sys
import warnings
from pathlib import Path

import atheris

sys.path.insert(0, str(Path(__file__).resolve().parent))

with atheris.instrument_imports(include=["qvf"]):
    import harness

import mutators  # noqa: E402

warnings.simplefilter("ignore")
OPS = list(mutators.OPS)


def custom_mutator(data: bytes, max_size: int, seed: int) -> bytes:
    rng = random.Random(seed)
    if rng.random() < 0.5:
        return atheris.Mutate(data, max_size)
    plan = [(rng.choice(OPS), rng.randrange(2 ** 32)) for _ in range(rng.randint(1, 3))]
    try:
        out = mutators.mutate(data, plan)
    except Exception:  # noqa: BLE001 - a mutator bug must not stop the campaign
        return atheris.Mutate(data, max_size)
    return out[:max_size] if len(out) > max_size else out


def test_one_input(data: bytes) -> None:
    problems, _, _ = harness.check(data, memory=False)
    if problems:
        raise AssertionError("; ".join(problems))


if __name__ == "__main__":
    atheris.Setup(sys.argv, test_one_input, custom_mutator=custom_mutator)
    atheris.Fuzz()
