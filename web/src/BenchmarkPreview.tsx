import { useState } from "react";

export default function BenchmarkPreview({
  caseId,
  plot,
  alt,
  eager = false,
}: {
  caseId: string;
  plot: string;
  alt: string;
  eager?: boolean;
}) {
  const [unavailable, setUnavailable] = useState(false);
  return (
    <div className="benchmark-preview">
      {unavailable ? (
        <p className="preview-unavailable">
          Recorded field preview unavailable
        </p>
      ) : (
        <img
          src={`${import.meta.env.BASE_URL}livebench/figures/${caseId}/${plot}-preview.png`}
          width="540"
          height="330"
          alt={alt}
          loading={eager ? "eager" : "lazy"}
          fetchPriority={eager ? "high" : "auto"}
          onError={() => setUnavailable(true)}
        />
      )}
    </div>
  );
}
