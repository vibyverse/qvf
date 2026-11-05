"""Opt-in fuzz tests of the Python validator and readers (set QVF_FUZZ=1).

    QVF_FUZZ=1 python -m pytest python/tests/fuzz -q

runs a short Hypothesis pass per target (QVF_FUZZ_EXAMPLES inputs each,
default 200). `run_fuzz.py` runs time-budgeted campaigns, and
`atheris_fuzz.py` a coverage-guided one. See harness.py for the invariants.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(not os.environ.get("QVF_FUZZ"),
                                reason="fuzzing is opt-in: set QVF_FUZZ=1")

sys.path.insert(0, str(Path(__file__).resolve().parent))

if os.environ.get("QVF_FUZZ"):
    from hypothesis import HealthCheck, given, settings

    import harness
    import strategies

    EXAMPLES = int(os.environ.get("QVF_FUZZ_EXAMPLES", "200"))

    def _run(target):
        limits = harness.DEFAULT if target == "bomb" else harness.SMALL

        @settings(max_examples=EXAMPLES, deadline=None, database=None,
                  suppress_health_check=list(HealthCheck))
        @given(strategies.inputs(target))
        def check(case):
            name, plan, byte_seed, raw = case
            problems, _, _ = harness.check(raw, limits=limits)
            assert not problems, (name, plan, byte_seed, problems)

        check()

    @pytest.mark.parametrize("target", sorted(strategies.TARGETS))
    def test_fuzz(target):
        _run(target)

    def test_seeds_are_clean():
        for name, raw in strategies.SEEDS:
            problems, _, _ = harness.check(raw)
            assert not problems, (name, problems)


def test_fuzz_is_opt_in():
    """Placeholder so a plain run of this directory reports the skip."""
