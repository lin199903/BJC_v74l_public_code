# Provenance timeline

This file records the file-system modification times of the protocol, script and
result files that establish the analysis chronology. Archive packaging rewrites
embedded timestamps, so the values below were read from the working file system
at packaging time and are reproduced here as a durable record.

The chronology relevant to the pre-declaration of the supplementary
direction-consistency rule is the ordering of `run_P22_scorecard.py` (which
carries the rule in its docstring) against the S20 reference demonstrations and
the S20b comparator run.

| Modification time | Archive path |
|---|---|
| 2026-07-22 23:29:46 | `scripts/LOCKED_PROTOCOL_2026-07-22.md` |
| 2026-07-22 23:29:49 | `scripts/analysis_contract.json` |
| 2026-07-22 23:43:04 | `scripts/POSTHOC_SENSITIVITY_ADDENDUM_2026-07-22.md` |
| 2026-08-10 23:51:48 | `scripts/run_P22_scorecard.py` |
| 2026-08-10 23:52:41 | `results/scorecard_summary.json` |
| 2026-08-14 00:43:12 | `results/S20_summary.json` |
| 2026-08-14 19:21:54 | `scripts/run_S20_state_controls.py` |
| 2026-08-14 19:33:11 | `scripts/run_S20b_csf1r.py` |
| 2026-08-14 19:34:10 | `results/S20b_summary.json` |
| 2026-09-10 19:49:15 | `scripts/S24_fine_myeloid_sensitivity/PROTOCOL_LOCKED_BEFORE_RESULTS.md` |

Sequenced reading: `run_P22_scorecard.py` predates the S20 reference
demonstrations and the S20b comparator run, so the supplementary rule was
written before those outcomes were inspected. The rule is applied uniformly
to all pairs, not only to the positive reference pairs.
