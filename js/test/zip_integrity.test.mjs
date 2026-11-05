import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";

import { QvfReader, validate } from "../dist/index.js";

/**
 * @param {"corpus" | "negative"} kind
 * @param {string} name
 */
const corpus = (kind, name) => new URL(`../../conformance/${kind}/${name}`, import.meta.url);

test("reference-writer zero spin density validates and decodes", async () => {
  const bytes = await readFile(corpus("corpus", "volume_spin_zero.qvf"));
  const report = await validate(bytes);
  assert.equal(report.ok, true, JSON.stringify(report.findings));
  const reader = QvfReader.fromBytes(bytes);
  const data = await reader.readBytes("spin", "data");
  assert.equal(data.length, 64 ** 3 * 8);
  assert.equal(data.every((b) => b === 0), true);
});

for (const name of ["entry_size_understated", "member_too_large", "total_too_large",
                     "crc_mismatch_stored", "crc_mismatch_deflated", "compression_ratio"]) {
  test(`bounded ZIP rejection: ${name}`, async () => {
    /** @type {{
     *   limits?: Partial<import("../dist/index.js").Limits>,
     *   expect: { code: string, severity: "error" | "warning" }[]
     * }} */
    const expected = JSON.parse(await readFile(corpus("negative", `${name}.expected.json`), "utf8"));
    const report = await validate(await readFile(corpus("negative", `${name}.qvf`)), { limits: expected.limits });
    assert.equal(report.ok, false);
    // As the conformance runners read `expect`: every expected error code is
    // reported (the schema checker may add E_SCHEMA for a member without sha256).
    const got = new Set(report.findings.filter((f) => f.severity === "error").map((f) => String(f.code)));
    for (const f of expected.expect.filter((f) => f.severity === "error")) {
      assert.ok(got.has(f.code), `${name}: ${f.code} missing from ${[...got].join(", ")}`);
    }
  });
}

test("a tiny entry declaring gigabytes is refused before allocation", async () => {
  const { ZipArchive } = await import("../dist/index.js");
  const { buildZip } = await import("./helpers.mjs");
  // 4 KiB of zeros deflate to a few bytes; the headers claim 2 GiB and, through
  // ZIP64, 7 GiB, both under the 8 GiB per-member cap.
  for (const fields of [{ usize: 2 ** 31 }, { zip64: true, zip64Size: 7 * 1024 ** 3 }]) {
    const bytes = buildZip([{ name: "x.bin", data: new Uint8Array(4096), ...fields }]);
    const zip = ZipArchive.open(bytes);
    assert.throws(() => zip.read("x.bin"), (err) => {
      assert.equal(/** @type {any} */ (err).code, "E_ZIP_UNREADABLE");
      assert.match(String(err), /compressed bytes can inflate to/);
      return true;
    });
    const report = await validate(bytes);
    assert.equal(report.ok, false);
  }
});
