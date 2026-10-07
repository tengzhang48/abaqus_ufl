import { defineConfig, loadEnv } from "vite";
import react from "@vitejs/plugin-react";
import { resolve } from "node:path";

export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), "VITE_");

  return {
    base: env.VITE_BASE_PATH || "/abaqus_ufl/",
    plugins: [react()],
    build: {
      sourcemap: true,
      rolldownOptions: {
        input: {
          home: resolve(process.cwd(), "index.html"),
          livebench: resolve(process.cwd(), "livebench/index.html"),
        },
      },
    },
  };
});
