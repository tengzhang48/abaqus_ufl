import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const data = JSON.parse(await readFile(resolve(here, "../site-data.json"), "utf8"));

test("workflow keeps Abaqus analysis user-owned", () => {
  assert.equal(data.workflow.length, 5);
  assert.match(data.workflow.at(-1).detail, /user owns/i);
});

test("scope names UMAT and UEL without merging their boundaries", () => {
  assert.deepEqual(data.scopes.map((item) => item.target), ["UMAT", "UEL"]);
  for (const scope of data.scopes) assert.ok(scope.boundary.length > 40);
});

test("paper evidence order and boundaries are explicit", () => {
  assert.deepEqual(data.paperEvidence.map((item) => item.id), ["tet4", "pasta", "gel", "corrosion"]);
  for (const item of data.paperEvidence) {
    if (item.id !== "corrosion") assert.ok(item.figurePaths.length > 0);
    assert.ok(item.boundary.length > 50);
  }
  const corrosion = data.paperEvidence.find((item) => item.id === "corrosion");
  assert.equal(corrosion.figurePaths.length, 0);
  assert.match(corrosion.boundary, /redistribution status/i);
});

test("fresh validation never conflates datacheck and completed analysis", () => {
  assert.equal(data.validation.sourceRevision, "0f52533");
  assert.equal(data.validation.reportRevision, "c5d2a59");
  assert.equal(data.validation.completeSolves.length, 6);
  assert.equal(data.validation.datachecks.length, 3);
  assert.match(data.validation.boundary, /does not establish nonlinear convergence/i);
});
