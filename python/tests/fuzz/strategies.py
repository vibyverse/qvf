"""Hypothesis strategies: a seed archive plus a shrinkable list of operators."""
from __future__ import annotations

from hypothesis import strategies as st

import mutators

SEEDS = mutators.seeds()

TARGETS = {
    # ZIP headers, ZIP64 fields, sizes, offsets, methods, descriptors, names.
    "zip": (mutators.ZIP_OPS, 0.0),
    # manifest.json and member payloads: types, numbers, nesting, dtypes,
    # shapes, sha256, sizes, text-level JSON damage.
    "json": (mutators.JSON_OPS, 0.0),
    # Both, plus raw byte noise after serialization.
    "mixed": (mutators.ZIP_OPS + mutators.JSON_OPS, 0.3),
    # Zip bombs and huge declared sizes, run with the default caps.
    "bomb": (["bomb", "declare_huge", "zip64", "usize", "descriptor", "duplicate"], 0.0),
}


@st.composite
def inputs(draw, target: str):
    ops, noise = TARGETS[target]
    name, raw = draw(st.sampled_from(SEEDS))
    plan = draw(st.lists(st.tuples(st.sampled_from(ops), st.integers(0, 2 ** 32 - 1)),
                         min_size=1, max_size=4))
    byte_seed = None
    if noise and draw(st.floats(0, 1)) < noise:
        byte_seed = draw(st.integers(0, 2 ** 32 - 1))
    return name, plan, byte_seed, mutators.mutate(raw, plan, byte_seed)
