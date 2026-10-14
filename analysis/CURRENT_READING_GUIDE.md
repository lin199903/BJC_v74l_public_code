# Current reading guide for the v74h analysis archive

Read this guide before interpreting historical status headings or cohort notes. All original archive members, scripts and numeric outputs are retained unchanged. The current manuscript and Supplementary Material incorporate the following verified reporting corrections.

1. **Chronology:** `v72_burden/V72_BURDEN_ESTIMAND_CONTRACT.md` is a historical, hashed contract. Its freeze-before-inspection heading refers to the v72 reanalysis; the identical sample-aggregate point estimate was already known in v39. See `v74h_reporting_corrections/CHRONOLOGY_OVERLAY.md`.
2. **Analytical units:** `extension_meta/meta_results.csv` retains historical notes using patient-level language. Those notes are superseded by `v74h_reporting_corrections/METADATA_ERRATA.md` and the companion corrected-notes CSV. Numeric fields and original weights are unchanged. The number of distinct patients underlying GSE205154's 216 core-collapsed analysis units is unresolved.
3. **TREM2 subsampling:** `scripts/run_S22_dilution.py` sets detectable-cell quotas relative to the original compartment, not exact retained-subset prevalence. It removes unselected detectable cells and reaggregates both target and programme. See `v74h_reporting_corrections/METHODS_READING_NOTES.md` and the deterministic denominator table.
4. **Rank pooling/bootstrap:** the current SI §23 supplies ddof=0, sample-count weights, fixed score-construction mapping and degenerate-draw handling. See the methods notes for reproduction details.

The historical contract, archived code and results remain accessible at their original paths for provenance. Current interpretation requires these overlays; isolated old metadata must not be treated as a verified statement of patient independence or prospective outcome ignorance. No new Monte Carlo simulation or core-result re-estimation was performed for v74h.
