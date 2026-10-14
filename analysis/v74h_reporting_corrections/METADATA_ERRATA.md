# Field-level metadata erratum: extension meta-analysis

Date: 2026-10-09. Applies to the notes fields of `extension_meta/meta_results.csv` and inherited patient-level descriptions of the extension.

GSE205154 contributes 216 primary-tumour core-collapsed analysis units; the number of distinct patients is unresolved. GSE21501 contributes 129 eligible arrays/specimens with no established one-to-one specimen-to-patient mapping. The extension retains the archived analysis-unit weights and should be read as a cohort-composition and inferential-unit sensitivity, not confirmed independent-patient evidence.

Historical notes using `patient-level (n=216 patients/specimens)` or calling the complete extension `patient-level cohorts` are superseded by this erratum. `meta_results_notes_corrected.csv` retains every non-notes field exactly and supplies explicit corrected notes. The original CSV remains unchanged at its historical path. This is a reporting correction, not a numerical correction or a patient-cluster reanalysis.
