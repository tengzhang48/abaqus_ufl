import type { MouseEvent } from "react";
import siteData from "../site-data.json";

const home = import.meta.env.BASE_URL;
export const livebenchUrl = `${home}livebench/`;

function closeMenu(event: MouseEvent<HTMLAnchorElement>) {
  event.currentTarget.closest("details")?.removeAttribute("open");
}

export default function SiteHeader({
  livebench = false,
}: {
  livebench?: boolean;
}) {
  const links = [
    ["How it works", `${home}#how-it-works`],
    ["Scope", `${home}#scope`],
    ["Examples", `${home}#examples`],
    ["Livebench", livebenchUrl],
    ["Paper examples", `${home}#paper-evidence`],
    ["Status", `${home}#status`],
  ];
  return (
    <header className="site-header">
      <a className="brand" href={`${home}#top`} aria-label="abaqus ufl home">
        <span className="brand-mark" aria-hidden="true">
          a/u
        </span>
        <span>
          abaqus_<strong>ufl</strong>
        </span>
      </a>
      <nav className="desktop-nav" aria-label="Project website">
        {links.map(([label, href]) => (
          <a
            key={label}
            href={href}
            aria-current={
              livebench && label === "Livebench" ? "page" : undefined
            }
          >
            {label}
          </a>
        ))}
      </nav>
      <a className="header-link" href={siteData.repository.url}>
        GitHub <span aria-hidden="true">↗</span>
      </a>
      <details className="mobile-menu">
        <summary aria-label="Open navigation">
          <span aria-hidden="true">Menu</span>
        </summary>
        <nav aria-label="Mobile project website">
          {links.map(([label, href]) => (
            <a
              key={label}
              href={href}
              onClick={closeMenu}
              aria-current={
                livebench && label === "Livebench" ? "page" : undefined
              }
            >
              {label}
            </a>
          ))}
          <a href={siteData.repository.url} onClick={closeMenu}>
            GitHub
          </a>
        </nav>
      </details>
    </header>
  );
}
