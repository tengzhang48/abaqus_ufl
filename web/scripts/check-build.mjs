import assert from "node:assert/strict";
import { readFile, readdir } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const webRoot = resolve(here, "..");
const dist = resolve(webRoot, "dist");

const index = await readFile(resolve(dist, "index.html"), "utf8");
const notices = await readFile(resolve(dist, "THIRD_PARTY_NOTICES.txt"), "utf8");
const files = await readdir(resolve(dist, "assets"));
const base = process.env.VITE_BASE_PATH || "/abaqus_ufl/";

assert.ok(index.includes(`${base}assets/`), `built asset URLs must use the configured base path: ${base}`);
assert.match(notices, /React 19\.2\.8/);
assert.match(notices, /React DOM 19\.2\.8/);
assert.match(notices, /Scheduler 0\.27\.0/);
assert.ok(
  files.every((name) => !/corrosion|corrison/i.test(name)),
  "the provenance-pending corrosion image must not enter the Pages artifact",
);

console.log(`Checked ${files.length} built assets, the Pages base path, and third-party notices.`);
