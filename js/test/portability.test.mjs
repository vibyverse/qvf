// The built package must run in a browser: no Node-only imports or globals.

import assert from "node:assert/strict";
import { readFileSync, readdirSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { test } from "node:test";

const dist = fileURLToPath(new URL("../dist/", import.meta.url));
const src = fileURLToPath(new URL("../src/", import.meta.url));

/** @param {string} dir @param {RegExp} pattern */
function filesIn(dir, pattern) {
  return readdirSync(dir)
    .filter((n) => pattern.test(n))
    .map((n) => ({ name: n, text: readFileSync(dir + n, "utf8") }));
}

/** Specifiers imported by an ES module, from `import ... from "x"` and `export ... from "x"`. */
function specifiers(/** @type {string} */ text) {
  // No quote or '=' before `from`: string data such as "derived_from" in the
  // generated schemaTable.js is not an import.
  return [...text.matchAll(/(?:import|export)\s[^;"'=]*?\bfrom\s*["']([^"']+)["']/g)].map((m) => m[1]);
}

test("dist imports only sibling modules and fflate", () => {
  const files = filesIn(dist, /\.js$/);
  assert.ok(files.length >= 10);
  for (const { name, text } of files) {
    for (const spec of specifiers(text)) {
      assert.ok(spec === "fflate" || /^\.\/[A-Za-z0-9]+\.js$/.test(spec), `${name} imports ${spec}`);
    }
    assert.doesNotMatch(text, /\brequire\s*\(/, `${name} uses require`);
    assert.doesNotMatch(text, /\bnode:/, `${name} mentions node:`);
    assert.doesNotMatch(text, /\bprocess\.|\bBuffer\b|__dirname|__filename/, `${name} uses a Node global`);
  }
});

test("src does not import Node built-ins", () => {
  for (const { name, text } of filesIn(src, /\.ts$/)) {
    assert.doesNotMatch(text, /from\s*["'](?:node:)?(?:fs|path|crypto|zlib|buffer|os|url|stream|util)["']/, name);
  }
});

test("the only runtime dependency is fflate", () => {
  const pkg = JSON.parse(readFileSync(new URL("../package.json", import.meta.url), "utf8"));
  assert.deepEqual(Object.keys(pkg.dependencies), ["fflate"]);
  assert.deepEqual(Object.keys(pkg.devDependencies).sort(), ["@types/node", "typescript"]);
  assert.equal(pkg.type, "module");
  assert.equal(pkg.private, true);
  assert.equal(pkg.name, "@qvf/core");
});

test("a packed tarball ships the build and the repository license files", () => {
  const pkg = JSON.parse(readFileSync(new URL("../package.json", import.meta.url), "utf8"));
  assert.deepEqual(pkg.files, ["dist", "src", "LICENSE", "NOTICE", "README.md"]);
  assert.equal(pkg.scripts.prepack, "npm run clean && npm run build");
  for (const name of ["LICENSE", "NOTICE"]) {
    assert.equal(readFileSync(new URL(`../${name}`, import.meta.url), "utf8"),
                 readFileSync(new URL(`../../${name}`, import.meta.url), "utf8"), name);
  }
});
