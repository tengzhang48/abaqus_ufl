import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import LivebenchPage from "./LivebenchPage";
import "./styles.css";

createRoot(document.getElementById("root")!).render(
  <StrictMode>
    <LivebenchPage />
  </StrictMode>,
);
