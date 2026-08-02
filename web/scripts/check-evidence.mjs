import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { access } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const webRoot = resolve(here, "..");
const repoRoot = resolve(webRoot, "..");
const data = JSON.parse(await readFile(resolve(webRoot, "site-data.json"), "utf8"));

const validationPath = resolve(repoRoot, data.validation.sourcePath);
const validation = await readFile(validationPath, "utf8");
const examples = await readFile(resolve(repoRoot, "examples/README.md"), "utf8");
const paperExamples = await readFile(resolve(repoRoot, "paper_examples/README.md"), "utf8");
const citation = await readFile(resolve(repoRoot, "CITATION.cff"), "utf8");
const credits = await readFile(resolve(repoRoot, "CREDITS.md"), "utf8");
const gelRecord = await readFile(resolve(repoRoot, "paper_examples/gel_bilayer/README.md"), "utf8");
const corrosionRecord = await readFile(resolve(repoRoot, "paper_examples/phasefield_corrosion/README.md"), "utf8");
const appSource = await readFile(resolve(webRoot, "src/App.tsx"), "utf8");
const webReadme = await readFile(resolve(webRoot, "README.md"), "utf8");
const notices = await readFile(resolve(webRoot, "public/THIRD_PARTY_NOTICES.txt"), "utf8");

for (const phrase of [
  "137 passed",
  "clone started at revision `0f52533`",
  "committed in `c5d2a59`",
  "all six released examples",
  "All four manuscript generation entry points",
  "The full corrosion, `n=16` Tet4, gel-bilayer, and pasta analyses were not",
  "A successful datacheck verifies input processing",
]) {
  assert.ok(validation.includes(phrase), `validation source is missing: ${phrase}`);
}

assert.equal(data.examples.length, 6, "the public example allowlist must contain six entries");
assert.equal(data.examples.filter((item) => item.target === "UMAT").length, 4);
assert.equal(data.examples.filter((item) => item.target === "UEL").length, 2);
assert.deepEqual(data.paperEvidence.map((item) => item.id), ["tet4", "pasta", "gel", "corrosion"]);
assert.equal(data.paperEvidence.find((item) => item.id === "corrosion").figurePaths.length, 0);
assert.ok(!appSource.includes("phasefield_corrosion/figure"), "corrosion figure must not be bundled while redistribution status is open");
assert.ok(!appSource.includes('caption: "t ='), "pasta exports must not receive inferred time captions");
assert.equal(data.workflow.length, 5);
assert.equal(data.validation.sourceRevision, "0f52533");
assert.equal(data.validation.reportRevision, "c5d2a59");
assert.match(citation, /version: 0\.1\.0/);
assert.match(citation, /family-names: Zhang/);
assert.match(credits, /cd1fb320a90ada8ebb7a9437254549a0d181a0e0/);
assert.match(credits, /neither distribution includes the exact\s+BSD variant, license text, or copyright notice/);
assert.match(corrosionRecord, /text-only record/);
assert.match(gelRecord, /does not retain the full-run solver log/);
assert.doesNotMatch(webReadme, /license-cleared/i);
assert.match(notices, /React 19\.2\.8/);
assert.match(notices, /Copyright \(c\) Meta Platforms/);

const paths = [
  data.validation.sourcePath,
  "CREDITS.md",
  ...data.examples.map((item) => item.path),
  ...data.paperEvidence.flatMap((item) => [item.path, ...item.figurePaths]),
];
for (const path of paths) {
  await access(resolve(repoRoot, path));
}

for (const example of data.examples) {
  const directoryLink = example.path.replace("examples/", "").replace("README.md", "");
  assert.ok(examples.includes(`(${directoryLink})`), `example is not allowlisted: ${example.path}`);
}
for (const paper of data.paperEvidence) {
  assert.ok(paperExamples.includes(`${paper.path.split("/")[1]}/`), `paper package is not indexed: ${paper.path}`);
}

const serialized = JSON.stringify(data).toLowerCase();
for (const forbidden of ["fully validated", "production-ready", "abaqus compatible", "signoff-ready"]) {
  assert.ok(!serialized.includes(forbidden), `unbounded claim in site data: ${forbidden}`);
}

console.log(`Checked ${paths.length} evidence paths and the ${data.validation.date} validation record.`);
