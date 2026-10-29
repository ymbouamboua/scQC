# scQC (development version)

* Add optional per-library ambient RNA correction to `run_qc()`, with decontX
  as default when enabled and SoupX as an alternative.
* Default to corrected counts only; optionally retain raw and corrected assays.
* Preserve raw QC metrics, refresh corrected metrics, invalidate stale compact
  analyses, and export per-cell/per-library contamination diagnostics.
