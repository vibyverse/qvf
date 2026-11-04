#!/usr/bin/env node
// Seeded, time-budgeted fuzz campaign of the JavaScript validator and reader.
//
//   npm run build
//   node test/fuzz/run_fuzz.mjs --target mixed --seconds 600 [--seed N] [--out DIR]
//   node test/fuzz/run_fuzz.mjs --replay FILE.qvf [--target bomb]
//
// Targets (the mutation classes of python/tests/fuzz/mutators.py):
//   zip    ZIP headers, ZIP64 fields, sizes, offsets, methods, descriptors, names
//   json   manifest.json and member payloads: types, numbers, nesting, dtypes,
//          shapes, sha256, sizes, text-level JSON damage
//   mixed  both, plus raw byte noise after serialization
//   bomb   zip bombs and huge declared sizes, checked with the default caps
//
// Inputs are generated and checked in a worker thread; the main thread is a
// watchdog that terminates a worker stuck on one input for --timeout seconds
// (a hang) and starts a fresh one. The first input of each new problem
// signature is saved to OUT/<target>/ (default ../build/fuzz-js, ignored by
// Git) with its seed archive and plan. Exits 1 when any problem was found.

import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { createHash } from "node:crypto";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { Worker, isMainThread, parentPort, workerData } from "node:worker_threads";

const TARGETS = {
  zip: { ops: "zip", noise: 0 },
  json: { ops: "json", noise: 0 },
  mixed: { ops: "all", noise: 0.3 },
  bomb: { ops: ["bomb", "declareHuge", "zip64", "usize", "descriptor", "duplicate"], noise: 0 },
};

function parseArgs(argv) {
  const args = { target: "mixed", seconds: 60, seed: Date.now() >>> 0, timeout: 10, replay: null,
    out: fileURLToPath(new URL("../../../build/fuzz-js/", import.meta.url)) };
  for (let i = 0; i < argv.length; i++) {
    const k = argv[i].replace(/^--/, "");
    if (!(k in args)) throw new Error(`unknown option ${argv[i]}`);
    args[k] = ["target", "out", "replay"].includes(k) ? argv[++i] : Number(argv[++i]);
  }
  if (!(args.target in TARGETS)) throw new Error(`unknown target ${args.target}`);
  return args;
}

async function workerMain() {
  const harness = await import("./harness.mjs");
  harness.install();
  await harness.load();
  const m = await import("./mutators.mjs");
  const { target, seed, start, replay } = workerData;
  const limits = target === "bomb" ? harness.DEFAULT : harness.SMALL;
  if (replay) {
    const r = await harness.check(readFileSync(replay), { limits });
    parentPort.postMessage({ type: "replay", ...r });
    return;
  }
  const t = TARGETS[target];
  const ops = t.ops === "zip" ? m.ZIP_OPS : t.ops === "json" ? m.JSON_OPS : t.ops === "all" ? Object.keys(m.OPS) : t.ops;
  const all = m.seeds();
  for (let i = start; ; i++) {
    const r = m.rngFrom((seed ^ Math.imul(i, 0x9e3779b1)) >>> 0);
    const [name, raw] = r.choice(all);
    const plan = Array.from({ length: r.randint(1, 4) }, () => [r.choice(ops), r.int(2 ** 32)]);
    const byteSeed = t.noise && r.random() < t.noise ? r.int(2 ** 32) : null;
    let input;
    try {
      input = m.mutate(raw, plan, byteSeed);
    } catch (err) {
      parentPort.postMessage({ type: "harness-error", i, error: String(err?.stack ?? err) });
      continue;
    }
    parentPort.postMessage({ type: "begin", i, name, plan, byteSeed, input });
    const t0 = performance.now();
    const res = await harness.check(input, { limits });
    parentPort.postMessage({ type: "end", i, ms: performance.now() - t0, ...res });
  }
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const worker = (data) => new Worker(fileURLToPath(import.meta.url), { workerData: data,
    resourceLimits: { maxOldGenerationSizeMb: 4096 } });
  if (args.replay) {
    const w = worker({ target: args.target, replay: args.replay });
    w.on("message", (msg) => {
      console.log(JSON.stringify(msg, null, 2));
      process.exitCode = msg.problems.length ? 1 : 0;
    });
    return;
  }
  const out = join(args.out, args.target);
  mkdirSync(out, { recursive: true });
  const stats = { inputs: 0, slowestMs: 0, largestAlloc: 0, hangs: 0, crashes: 0, signatures: {} };
  const deadline = Date.now() + args.seconds * 1000;

  const record = (problem, cur) => {
    const key = problem.replace(/\d{3,}/g, "N").slice(0, 120);
    if (stats.signatures[key]) {
      stats.signatures[key].count++;
      return;
    }
    const digest = createHash("sha256").update(cur.input).digest("hex").slice(0, 16);
    writeFileSync(join(out, `${digest}.qvf`), cur.input);
    stats.signatures[key] = { count: 1, example: `${digest}.qvf`, seed: cur.name, plan: cur.plan,
      byteSeed: cur.byteSeed, problem };
    console.log(`NEW ${problem}\n    seed=${cur.name} plan=${JSON.stringify(cur.plan)} byteSeed=${cur.byteSeed} -> ${digest}.qvf`);
  };

  const started = Date.now();
  // Written every few seconds too, so an interrupted campaign keeps its numbers.
  const writeSummary = () => {
    const summary = { target: args.target, seed: args.seed, seconds: Math.round((Date.now() - started) / 1000),
      inputs: stats.inputs, slowest_ms: Math.round(stats.slowestMs), largest_alloc_bytes: stats.largestAlloc,
      hangs: stats.hangs, crashes: stats.crashes, problems: stats.signatures };
    writeFileSync(join(out, "summary.json"), JSON.stringify(summary, null, 2));
    return summary;
  };
  const ticker = setInterval(writeSummary, 5000);
  let next = 0;
  while (Date.now() < deadline) {
    await new Promise((resolve) => {
      const w = worker({ target: args.target, seed: args.seed, start: next });
      let cur = null;
      let timer = null;
      const stop = () => {
        clearTimeout(timer);
        clearTimeout(budget);
        w.terminate().then(resolve);
      };
      const budget = setTimeout(stop, Math.max(0, deadline - Date.now()));
      w.on("message", (msg) => {
        if (msg.type === "begin") {
          cur = msg;
          timer = setTimeout(() => {
            stats.hangs++;
            record(`timeout: over ${args.timeout} s`, cur);
            next = cur.i + 1;
            stop();
          }, args.timeout * 1000);
        } else if (msg.type === "end") {
          clearTimeout(timer);
          stats.inputs++;
          stats.slowestMs = Math.max(stats.slowestMs, msg.ms);
          stats.largestAlloc = Math.max(stats.largestAlloc, msg.largest);
          for (const p of msg.problems) record(p, cur);
          next = msg.i + 1;
        } else if (msg.type === "harness-error") {
          console.error(`harness error at input ${msg.i}: ${msg.error}`);
        }
      });
      w.on("error", (err) => {
        stats.crashes++;
        if (cur) record(`worker crashed: ${String(err?.message ?? err)}`, cur);
        next = (cur?.i ?? next) + 1;
        stop();
      });
    });
  }
  clearInterval(ticker);
  const summary = writeSummary();
  const { problems, ...head } = summary;
  console.log(JSON.stringify(head));
  console.log(`problem signatures: ${Object.keys(problems).length}`);
  process.exitCode = Object.keys(problems).length ? 1 : 0;
}

if (isMainThread) await main();
else await workerMain();
