// Cross-check: the Python reference validator accepts archives the JS writer
// produced. Skipped when Python with the qvf package and jsonschema is not
// available (pip install -e "python/[numpy,validate]").

import assert from "node:assert/strict";
import { execFileSync, spawnSync } from "node:child_process";
import { mkdtempSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { test } from "node:test";

import { Writer } from "../dist/index.js";
import { FIXED, H2, SOURCE, fullArchive } from "./writerFixtures.mjs";

const PYTHON_DIR = fileURLToPath(new URL("../../python/", import.meta.url));
const PYTHON = process.env.QVF_PYTHON ?? (process.platform === "win32" ? "python" : "python3");
const env = { ...process.env, PYTHONPATH: [PYTHON_DIR, process.env.PYTHONPATH].filter(Boolean).join(process.platform === "win32" ? ";" : ":") };

const probe = spawnSync(PYTHON, ["-c", "import qvf, jsonschema, numpy"], { env, encoding: "utf8" });
const skip = probe.status === 0 ? false : `Python with qvf, jsonschema and numpy is unavailable (${(probe.error?.message ?? probe.stderr ?? "").trim().split("\n").pop()})`;

/**
 * Validate `bytes` with `python -m qvf validate --json`.
 * @param {Uint8Array} bytes
 */
function pythonValidate(bytes) {
  const dir = mkdtempSync(join(tmpdir(), "qvf-js-writer-"));
  try {
    const file = join(dir, "archive.qvf");
    writeFileSync(file, bytes);
    const run = spawnSync(PYTHON, ["-m", "qvf", "validate", "--json", file], { env, encoding: "utf8" });
    assert.ok(run.stdout, `python produced no report: ${run.stderr}`);
    return JSON.parse(run.stdout);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
}

/** @param {any} report */
function assertPythonOk(report) {
  assert.deepEqual(report.findings, []);
  assert.equal(report.ok, true);
  assert.equal(report.validation_context.schema_checked, true);
}

test("the Python validator accepts the full JS-written archive", { skip }, async () => {
  const report = pythonValidate(await (await fullArchive()).toBytes());
  assertPythonOk(report);
  assert.equal(report.n_sections, 19);
});

test("the Python validator accepts forced ZIP64 and stored archives", { skip }, async () => {
  assertPythonOk(pythonValidate(await (await fullArchive({ forceZip64: true })).toBytes()));
  assertPythonOk(pythonValidate(await (await fullArchive({ compression: "store" })).toBytes()));
});

test("the Python reader decodes JS-written tensors bit for bit", { skip }, async () => {
  const w = new Writer(SOURCE, FIXED);
  await w.addStructure(H2);
  w.addSection({ id: "t", kind: "x_test.t" });
  await w.addBinary("t", "f64", [[Math.PI, -0.0], [1e-300, 2 ** 52 + 1]]);
  await w.addBinary("t", "i64", [-(2n ** 63n), 2n ** 63n - 1n], { dtype: "int64" });
  await w.addBinary("t", "u16", new Uint16Array([1, 65535]));
  const dir = mkdtempSync(join(tmpdir(), "qvf-js-writer-"));
  try {
    const file = join(dir, "t.qvf");
    writeFileSync(file, await w.toBytes());
    const script = [
      "import json, sys, qvf",
      "r = qvf.QvfReader(sys.argv[1])",
      "members = next(s for s in r.sections if s['id'] == 't')['members']",
      "out = {}",
      "for role in ('f64', 'i64', 'u16'):",
      "    a = r.read_array(members[role])",
      "    out[role] = [str(a.dtype), list(a.shape), [repr(x) for x in a.ravel().tolist()]]",
      "print(json.dumps(out))",
    ].join("\n");
    const out = JSON.parse(execFileSync(PYTHON, ["-c", script, file], { env, encoding: "utf8" }));
    assert.deepEqual(out.f64, ["float64", [2, 2], ["3.141592653589793", "-0.0", "1e-300", "4503599627370497.0"]]);
    assert.deepEqual(out.i64, ["int64", [2], ["-9223372036854775808", "9223372036854775807"]]);
    assert.deepEqual(out.u16, ["uint16", [2], ["1", "65535"]]);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
