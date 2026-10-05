import { readFile } from "node:fs/promises";
import ts from "typescript";

// Import the browser's parser in Node without introducing a second validator
// or requiring a Node release with native TypeScript support.
const source = await readFile(new URL("../src/livebench-report.ts", import.meta.url), "utf8");
const javascript = ts.transpileModule(source, { compilerOptions: { module: ts.ModuleKind.ESNext, target: ts.ScriptTarget.ES2022 } }).outputText;
export const { validateReport } = await import(`data:text/javascript;base64,${Buffer.from(javascript).toString("base64")}`);
